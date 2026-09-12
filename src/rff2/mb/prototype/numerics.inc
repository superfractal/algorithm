// Modified by GPT-6 on 2026-09-10, 2026-09-11, 2026-09-12
#include "locate.hpp"
#include "integer_square.hpp"
#include <gmp.h>
#include <omp.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <memory>

// Small owning GMP wrappers: F is a multiprecision float, I an integer,
// and C a complex pair. Precision belongs to each F instance.
struct F {
    mpf_t v;
    F(unsigned b = 64) {
        mpf_init2(v, b);
    }
    F(const F &o) {
        mpf_init2(v, mpf_get_prec(o.v));
        mpf_set(v, o.v);
    }
    F(F &&o) noexcept {
        mpf_init2(v, 64);
        mpf_swap(v, o.v);
    }
    F &operator=(const F &o) {
        if (this != &o) {
            mpf_set_prec(v, mpf_get_prec(o.v));
            mpf_set(v, o.v);
        }
        return *this;
    }
    ~F() {
        mpf_clear(v);
    }
    operator mpf_ptr() {
        return v;
    }
    operator mpf_srcptr() const {
        return v;
    }
};
struct I {
    mpz_t v;
    I() {
        mpz_init(v);
    }
    ~I() {
        mpz_clear(v);
    }
    operator mpz_ptr() {
        return v;
    }
    operator mpz_srcptr() const {
        return v;
    }
};
struct C {
    F r, i;
    C(unsigned b = 64) : r(b), i(b) {};
};
void set(C &o, const C &a) {
    mpf_set(o.r, a.r);
    mpf_set(o.i, a.i);
}
void zero(C &o) {
    mpf_set_ui(o.r, 0);
    mpf_set_ui(o.i, 0);
}
void one(C &o) {
    mpf_set_ui(o.r, 1);
    mpf_set_ui(o.i, 0);
}
void add(C &o, const C &a, const C &b) {
    mpf_add(o.r, a.r, b.r);
    mpf_add(o.i, a.i, b.i);
}
void sub(C &o, const C &a, const C &b) {
    mpf_sub(o.r, a.r, b.r);
    mpf_sub(o.i, a.i, b.i);
}
struct Ops {
    // Reuse scratch storage; publish the result only after reading both inputs.
    // The order of these rounded GMP operations is part of the kernel behavior.
    F a, b, t, sumX, sumY, third;
    const bool threeProducts;
    Ops(unsigned bits) :
        a(bits), b(bits), t(bits), sumX(bits), sumY(bits), third(bits),
        threeProducts(mpf_get_prec(a) >= 512) {};
    void mul(C &o, const C &x, const C &y) {
        if (threeProducts) {
            mpf_add(sumX, x.r, x.i);
            mpf_add(sumY, y.r, y.i);
            mpf_mul(a, x.r, y.r);
            mpf_mul(t, x.i, y.i);
            mpf_mul(third, sumX, sumY);
            mpf_sub(b, third, a);
            mpf_sub(b, b, t);
            mpf_sub(a, a, t);
        } else {
            mpf_mul(a, x.r, y.r);
            mpf_mul(t, x.i, y.i);
            mpf_sub(a, a, t);
            mpf_mul(b, x.r, y.i);
            mpf_mul(t, x.i, y.r);
            mpf_add(b, b, t);
        }
        mpf_set(o.r, a);
        mpf_set(o.i, b);
    }
    void div(C &o, const C &x, const C &y) {
        mpf_mul(a, y.r, y.r);
        mpf_mul(b, y.i, y.i);
        mpf_add(t, a, b);
        if (mpf_sgn((mpf_srcptr) t) == 0)
            throw std::runtime_error("zero Newton derivative");
        mpf_mul(a, x.r, y.r);
        mpf_mul(b, x.i, y.i);
        mpf_add(a, a, b);
        mpf_div(a, a, t);
        F k(mpf_get_prec(a));
        mpf_mul(b, x.i, y.r);
        mpf_mul(k, x.r, y.i);
        mpf_sub(b, b, k);
        mpf_div(b, b, t);
        mpf_set(o.r, a);
        mpf_set(o.i, b);
    }
};
double lg(mpf_srcptr a) {
    long e;
    double x = mpf_get_d_2exp(&e, a);
    return x ? std::log2(std::abs(x)) + e : -1e300;
}
// Exponent-aligned norm: one logarithm rather than two logs, exp2, and log2.
// The tiny component may underflow after alignment, but is then negligible
// for this double-valued precision estimator (not for the orbit itself).
double la(const C &a) {
    long er = 0, ei = 0;
    double r = mpf_get_d_2exp(&er, a.r), i = mpf_get_d_2exp(&ei, a.i);
    if (r == 0 && i == 0)
        return -1e300;
    if (r == 0)
        er = ei;
    if (i == 0)
        ei = er;
    const long e = std::max(er, ei);
    r = std::ldexp(r, int(std::max(-1075L, er - e)));
    i = std::ldexp(i, int(std::max(-1075L, ei - e)));
    return e + 0.5 * std::log2(r * r + i * i);
}
void scale(mpz_ptr a, int oldb, int newb) {
    if (newb < oldb)
        mpz_tdiv_q_2exp(a, a, oldb - newb);
    else if (newb > oldb)
        mpz_mul_2exp(a, a, newb - oldb);
}
// Clamp the estimate before narrowing it to an integer. A zero intermediate
// orbit uses log norm -1e300 and must request full precision, not overflow int.
inline int taperedPrecisionBits(int goalbits, double offset, int fullbits) {
    if (fullbits < 128 || !std::isfinite(offset))
        throw std::runtime_error("invalid taper precision estimate");
    const double required = double(goalbits) - offset + 64;
    if (required <= 128)
        return 128;
    if (required >= fullbits)
        return fullbits;
    return int(std::ceil(required / 64)) * 64;
}
struct Point {
    unsigned long n;
    double logA;
    C z, D;
    Point(unsigned b, unsigned long n_, double l) : n(n_), logA(l), z(b), D(b) {};
};
// A shooting segment starts at z and ends at out. Its tangent transport is
// D_out = fa * D_in + fb. The state product fa excludes the critical z_0 factor;
// fb includes the parameter derivative's +1 term at every step.
// offset and logA store base-2 sensitivity estimates used for precision tapering.
struct Block {
    unsigned long start = 0, end = 0;
    double offset = 0, logA = 0;
    C z, fa, fb, out;
    Block(unsigned bits) : z(bits), fa(bits), fb(bits), out(bits) {};
};
struct KernelResult {
    double logA = 0;
    unsigned long period = 0;
};
// The integer orbit stores z = kr + (xr + i*xi) / 2^b. In the special
// near-minus-two path kr holds the integer real part; tapering rescales the
// residual integers and c together. Derivatives use separate precision budgets.
KernelResult kernel(const C &c,
                    const C &s,
                    unsigned long start,
                    unsigned long length,
                    int goalbits,
                    int abits,
                    int bbits,
                    double offset,
                    bool special,
                    bool pilot,
                    double zoom,
                    C &oz,
                    C &oA,
                    C &oB,
                    std::vector<Point> *points = nullptr) {
    int fullbits = mpf_get_prec(c.r), b = taperedPrecisionBits(goalbits, offset, fullbits);
    F tmp(fullbits);
    C z(std::max(abits, bbits)), A(abits), B(bbits), na(abits), nb(bbits);
    Ops oa(abits), ob(bbits);
    one(A);
    zero(B);
    I xr, xi, rr, ii, ri, t, crfull, cifull, crb, cib;
    mpf_set(tmp, c.r);
    if (special)
        mpf_add_ui(tmp, tmp, 2);
    mpf_mul_2exp(tmp, tmp, fullbits);
    mpz_set_f(crfull, tmp);
    mpf_set(tmp, c.i);
    mpf_mul_2exp(tmp, tmp, fullbits);
    mpz_set_f(cifull, tmp);
    mpz_set(crb, crfull);
    scale(crb, fullbits, b);
    mpz_set(cib, cifull);
    scale(cib, fullbits, b);
    long kr = 0;
    if (special)
        kr = std::lround(mpf_get_d(s.r));
    mpf_set(tmp, s.r);
    if (kr >= 0)
        mpf_sub_ui(tmp, tmp, kr);
    else
        mpf_add_ui(tmp, tmp, -kr);
    mpf_mul_2exp(tmp, tmp, b);
    mpz_set_f(xr, tmp);
    mpf_set(tmp, s.i);
    mpf_mul_2exp(tmp, tmp, b);
    mpz_set_f(xi, tmp);
    double logA = offset;
    auto exportZ = [&](C &dst) {
        mpf_set_z(dst.r, xr);
        mpf_div_2exp(dst.r, dst.r, b);
        if (kr >= 0)
            mpf_add_ui(dst.r, dst.r, kr);
        else
            mpf_sub_ui(dst.r, dst.r, -kr);
        mpf_set_z(dst.i, xi);
        mpf_div_2exp(dst.i, dst.i, b);
    };
    for (unsigned long j = 0; j < length; ++j) {
        if ((j & 255) == 0 && cancelled())
            return {};
        unsigned long n = start + j + 1;
        exportZ(z);
        if (n > 1)
            logA += 1 + la(z);
        ob.mul(nb, z, B);
        mpf_mul_2exp(nb.r, nb.r, 1);
        mpf_mul_2exp(nb.i, nb.i, 1);
        mpf_add_ui(nb.r, nb.r, 1);
        mpf_swap(nb.r, B.r);
        mpf_swap(nb.i, B.i);
        if (!pilot && n > 1) {
            oa.mul(na, z, A);
            mpf_mul_2exp(na.r, na.r, 1);
            mpf_mul_2exp(na.i, na.i, 1);
            mpf_swap(na.r, A.r);
            mpf_swap(na.i, A.i);
        }
        if (n > 1) {
            int nb_ = taperedPrecisionBits(goalbits, logA, fullbits);
            if (nb_ < b - 64 || nb_ > b) {
                scale(xr, b, nb_);
                scale(xi, b, nb_);
                b = nb_;
                mpz_set(crb, crfull);
                scale(crb, fullbits, b);
                mpz_set(cib, cifull);
                scale(cib, fullbits, b);
            }
        }
        integerSquare(rr, ri, xr, xi, ii, t);
        mpz_tdiv_q_2exp(rr, rr, b);
        mpz_tdiv_q_2exp(ri, ri, b - 1);
        if (special) {
            mpz_mul_si(t, xr, 2 * kr);
            mpz_add(rr, rr, t);
            mpz_mul_si(t, xi, 2 * kr);
            mpz_add(ri, ri, t);
            kr = kr * kr - 2;
        }
        mpz_add(xr, rr, crb);
        mpz_add(xi, ri, cib);
        if (special) {
            long ex;
            double mx = mpz_get_d_2exp(&ex, xr);
            double dx = std::ldexp(mx, int(ex - b));
            if (std::abs(dx) > .25 || std::abs(kr) > 2) {
                long k = std::lround(dx + kr);
                mpz_set_si(t, kr - k);
                mpz_mul_2exp(t, t, b);
                mpz_add(xr, xr, t);
                kr = k;
            }
        }
        if (pilot) {
            exportZ(oz);
            double lz = la(oz), ld = la(B);
            if (n % 256 == 0 && points) {
                points->emplace_back(fullbits, n, logA);
                set(points->back().z, oz);
                set(points->back().D, B);
            }
            if (lz - ld < -zoom * std::log2(10.) + 5) {
                set(oB, B);
                return {logA, n};
            }
            if (lz > 16)
                throw std::runtime_error("escape before period");
        }
    }
    exportZ(oz);
    set(oA, A);
    set(oB, B);
    return {logA, start + length};
}
struct EComplex {
    double r = 0, i = 0;
    long e = 0;
};
EComplex normalize(EComplex z) {
    double m = std::max(std::abs(z.r), std::abs(z.i));
    if (m == 0)
        return {};
    int e;
    std::frexp(m, &e);
    z.r = std::ldexp(z.r, -e);
    z.i = std::ldexp(z.i, -e);
    z.e += e;
    return z;
}
EComplex integer_complex(mpz_srcptr xr, mpz_srcptr xi, int bits, long kr) {
    long er = 0, ei = 0;
    double r = mpz_get_d_2exp(&er, xr), i = mpz_get_d_2exp(&ei, xi);
    if (r == 0)
        er = -1000000000;
    if (i == 0)
        ei = -1000000000;
    long e = std::max(er, ei) - bits;
    if (kr)
        e = std::max(e, 2L);
    if (e < -999999999)
        return {};
    return normalize({std::ldexp(r, int(std::clamp(er - bits - e, -2000L, 2000L))) +
                          std::ldexp(double(kr), int(std::clamp(-e, -2000L, 2000L))),
                      std::ldexp(i, int(std::clamp(ei - bits - e, -2000L, 2000L))),
                      e});
}
EComplex tangent(EComplex z, EComplex d) {
    EComplex q{z.r * d.r - z.i * d.i, z.r * d.i + z.i * d.r, z.e + d.e + 1};
    q = normalize(q);
    if (q.r == 0 && q.i == 0)
        return {1, 0, 0};
    if (q.e < -55)
        return {1, 0, 0};
    if (q.e < 55) {
        q.r += std::ldexp(1., -q.e);
        q = normalize(q);
    }
    return q;
}
double elg(EComplex z) {
    return z.r == 0 && z.i == 0 ? -1e300 : z.e + 0.5 * std::log2(z.r * z.r + z.i * z.i);
}
void to_f(C &dst, EComplex z) {
    mpf_set_d(dst.r, z.r);
    mpf_set_d(dst.i, z.i);
    if (z.e >= 0) {
        mpf_mul_2exp(dst.r, dst.r, z.e);
        mpf_mul_2exp(dst.i, dst.i, z.e);
    } else {
        mpf_div_2exp(dst.r, dst.r, -z.e);
        mpf_div_2exp(dst.i, dst.i, -z.e);
    }
}
KernelResult fastpilot(const C &c,
                       int goalbits,
                       bool special,
                       double zoom,
                       C &oz,
                       C &oD,
                       std::vector<Point> &points,
                       unsigned long cap = 100000000,
                       unsigned long knownPeriod = 0) {
    int fullbits = mpf_get_prec(c.r), b = std::min(fullbits, ((goalbits + 127) / 64) * 64);
    F tmp(fullbits);
    I xr, xi, rr, ii, ri, t, crfull, cifull, crb, cib;
    mpf_set(tmp, c.r);
    if (special)
        mpf_add_ui(tmp, tmp, 2);
    mpf_mul_2exp(tmp, tmp, fullbits);
    mpz_set_f(crfull, tmp);
    mpf_set(tmp, c.i);
    mpf_mul_2exp(tmp, tmp, fullbits);
    mpz_set_f(cifull, tmp);
    mpz_set(crb, crfull);
    scale(crb, fullbits, b);
    mpz_set(cib, cifull);
    scale(cib, fullbits, b);
    long kr = 0;
    double logA = 0;
    EComplex z, D;
    auto exportZ = [&](C &dst) {
        mpf_set_z(dst.r, xr);
        mpf_div_2exp(dst.r, dst.r, b);
        if (kr >= 0)
            mpf_add_ui(dst.r, dst.r, kr);
        else
            mpf_sub_ui(dst.r, dst.r, -kr);
        mpf_set_z(dst.i, xi);
        mpf_div_2exp(dst.i, dst.i, b);
    };
    for (unsigned long n = 1; n <= cap; ++n) {
        if ((n & 255) == 1 && cancelled())
            return {};
        if (n > 1)
            logA += 1 + elg(z);
        D = tangent(z, D);
        if (n > 1) {
            int nb_ = taperedPrecisionBits(goalbits, logA, fullbits);
            if (nb_ < b - 64 || nb_ > b) {
                scale(xr, b, nb_);
                scale(xi, b, nb_);
                b = nb_;
                mpz_set(crb, crfull);
                scale(crb, fullbits, b);
                mpz_set(cib, cifull);
                scale(cib, fullbits, b);
            }
        }
        integerSquare(rr, ri, xr, xi, ii, t);
        mpz_tdiv_q_2exp(rr, rr, b);
        mpz_tdiv_q_2exp(ri, ri, b - 1);
        if (special) {
            mpz_mul_si(t, xr, 2 * kr);
            mpz_add(rr, rr, t);
            mpz_mul_si(t, xi, 2 * kr);
            mpz_add(ri, ri, t);
            kr = kr * kr - 2;
        }
        mpz_add(xr, rr, crb);
        mpz_add(xi, ri, cib);
        if (special) {
            long ex;
            double mx = mpz_get_d_2exp(&ex, xr);
            double dx = std::ldexp(mx, int(ex - b));
            if (std::abs(dx) > .25 || std::abs(kr) > 2) {
                long k = std::lround(dx + kr);
                mpz_set_si(t, kr - k);
                mpz_mul_2exp(t, t, b);
                mpz_add(xr, xr, t);
                kr = k;
            }
        }
        z = integer_complex(xr, xi, b, kr);
        double lz = elg(z), ld = elg(D);
        // A completed orbit endpoint is not a segment start. Test completion before
        // sampling checkpoints, including when an inferred period is a multiple of 256.
        if (knownPeriod ? n == knownPeriod : lz - ld < -zoom * std::log2(10.) + 5) {
            exportZ(oz);
            to_f(oD, D);
            return {logA, n};
        }
        if (n % 256 == 0) {
            points.emplace_back(fullbits, n, logA);
            exportZ(points.back().z);
            to_f(points.back().D, D);
        }
        if (lz > 16)
            throw std::runtime_error("escape before period");
    }
    return {logA, 0};
}

EComplex emul2(EComplex z, EComplex a) {
    return normalize({z.r * a.r - z.i * a.i, z.r * a.i + z.i * a.r, z.e + a.e + 1});
}
KernelResult fastblock(const C &c,
                       const C &s,
                       unsigned long start,
                       unsigned long length,
                       int goalbits,
                       double offset,
                       bool special,
                       C &oz,
                       C &oA,
                       C &oB) {
    int fullbits = mpf_get_prec(c.r), b = taperedPrecisionBits(goalbits, offset, fullbits);
    F tmp(fullbits);
    I xr, xi, rr, ii, ri, t, crfull, cifull, crb, cib;
    mpf_set(tmp, c.r);
    if (special)
        mpf_add_ui(tmp, tmp, 2);
    mpf_mul_2exp(tmp, tmp, fullbits);
    mpz_set_f(crfull, tmp);
    mpf_set(tmp, c.i);
    mpf_mul_2exp(tmp, tmp, fullbits);
    mpz_set_f(cifull, tmp);
    mpz_set(crb, crfull);
    scale(crb, fullbits, b);
    mpz_set(cib, cifull);
    scale(cib, fullbits, b);
    long kr = special ? std::lround(mpf_get_d(s.r)) : 0;
    mpf_set(tmp, s.r);
    if (kr >= 0)
        mpf_sub_ui(tmp, tmp, kr);
    else
        mpf_add_ui(tmp, tmp, -kr);
    mpf_mul_2exp(tmp, tmp, b);
    mpz_set_f(xr, tmp);
    mpf_set(tmp, s.i);
    mpf_mul_2exp(tmp, tmp, b);
    mpz_set_f(xi, tmp);
    double logA = offset;
    EComplex z = integer_complex(xr, xi, b, kr), A{1, 0, 0}, B;
    for (unsigned long j = 0; j < length; ++j) {
        if ((j & 255) == 0 && cancelled())
            return {};
        unsigned long n = start + j + 1;
        if (n > 1) {
            logA += 1 + elg(z);
            A = emul2(z, A);
        }
        B = tangent(z, B);
        if (n > 1) {
            int nb_ = taperedPrecisionBits(goalbits, logA, fullbits);
            if (nb_ < b - 64 || nb_ > b) {
                scale(xr, b, nb_);
                scale(xi, b, nb_);
                b = nb_;
                mpz_set(crb, crfull);
                scale(crb, fullbits, b);
                mpz_set(cib, cifull);
                scale(cib, fullbits, b);
            }
        }
        integerSquare(rr, ri, xr, xi, ii, t);
        mpz_tdiv_q_2exp(rr, rr, b);
        mpz_tdiv_q_2exp(ri, ri, b - 1);
        if (special) {
            mpz_mul_si(t, xr, 2 * kr);
            mpz_add(rr, rr, t);
            mpz_mul_si(t, xi, 2 * kr);
            mpz_add(ri, ri, t);
            kr = kr * kr - 2;
        }
        mpz_add(xr, rr, crb);
        mpz_add(xi, ri, cib);
        if (special) {
            long ex;
            double mx = mpz_get_d_2exp(&ex, xr);
            double dx = std::ldexp(mx, int(ex - b));
            if (std::abs(dx) > .25 || std::abs(kr) > 2) {
                long k = std::lround(dx + kr);
                mpz_set_si(t, kr - k);
                mpz_mul_2exp(t, t, b);
                mpz_add(xr, xr, t);
                kr = k;
            }
        }
        z = integer_complex(xr, xi, b, kr);
    }
    mpf_set_z(oz.r, xr);
    mpf_div_2exp(oz.r, oz.r, b);
    if (kr >= 0)
        mpf_add_ui(oz.r, oz.r, kr);
    else
        mpf_sub_ui(oz.r, oz.r, -kr);
    mpf_set_z(oz.i, xi);
    mpf_div_2exp(oz.i, oz.i, b);
    to_f(oA, A);
    to_f(oB, B);
    return {logA, start + length};
}
