// Created by GPT-6 on 2026-09-10
// Modified by GPT-6 on 2026-09-10, 2026-09-12
#pragma once
#include "inverse_samples.hpp"
#include <bit>
#include <cfenv>
#include <immintrin.h>
#include <limits>
#ifdef __FAST_MATH__
#error Outward bounds require strict floating-point arithmetic
#endif
#pragma clang fp contract(off)
namespace upper {
    using namespace inverse_samples;
    inline void environment() {
        if (!std::numeric_limits<double>::is_iec559 || std::numeric_limits<double>::digits != 53 ||
            std::fegetround() != FE_TONEAREST || (_mm_getcsr() & ((1u << 15) | (1u << 6))))
            throw std::runtime_error("unsupported bound arithmetic environment");
    }
    inline double up(double x) {
        if (!(x >= 0) || !std::isfinite(x))
            throw std::runtime_error("invalid positive bound");
        return std::bit_cast<double>(std::bit_cast<uint64_t>(x) + 1);
    }
    struct Positive {
        double m = 0;
        int64_t e = 0;
        Positive() = default;
        Positive(double value, int64_t exponent) {
            if (value == 0)
                return;
            if (!(value > 0) || !std::isfinite(value))
                throw std::runtime_error("invalid bound mantissa");
            int shift = 0;
            m = std::frexp(value, &shift);
            e = exponent + shift;
            if (e < -(int64_t(1) << 50) || e > (int64_t(1) << 50))
                throw std::runtime_error("bound exponent guard");
        }
        static Positive power(int64_t exponent) {
            return Positive(.5, exponent + 1);
        }
        static Positive abs(mpfr_srcptr x) {
            mpfr_exp_t exponent = 0;
            const double m = std::abs(mpfr_get_d_2exp(&exponent, x, MPFR_RNDA));
            return Positive(m, exponent);
        }
        Real real() const {
            if (e < std::numeric_limits<long>::min() || e > std::numeric_limits<long>::max())
                throw std::runtime_error("MPFR shift conversion guard");
            Real value;
            mpfr_set_d(value.x, m, MPFR_RNDU);
            mpfr_mul_2si(value.x, value.x, static_cast<long>(e), MPFR_RNDU);
            return value;
        }
    };
    // Inputs from add/mul lie in [1/4, 2 + one ULP]. Power-of-two
    // scaling is exact here and reproduces frexp without a library call.
    inline Positive normalizeArithmetic(double m, int64_t e) {
        if (m < .5) {
            m *= 2;
            --e;
        } else if (m >= 2) {
            m *= .25;
            e += 2;
        } else if (m >= 1) {
            m *= .5;
            ++e;
        }
        if (e < -(int64_t(1) << 50) || e > (int64_t(1) << 50))
            throw std::runtime_error("bound exponent guard");
        Positive result;
        result.m = m;
        result.e = e;
        return result;
    }
    inline Positive add(Positive a, Positive b) {
        if (!a.m)
            return b;
        if (!b.m)
            return a;
        if (a.e < b.e)
            std::swap(a, b);
        const auto distance = a.e - b.e;
        // For a discarded tiny summand, one upward ULP of a is conservative.
        if (distance > 1074)
            return normalizeArithmetic(up(a.m), a.e);
        const double small = up(std::scalbn(b.m, -static_cast<int>(distance)));
        return normalizeArithmetic(up(a.m + small), a.e);
    }
    inline Positive mul(Positive a, Positive b) {
        if (!a.m || !b.m)
            return {};
        return normalizeArithmetic(up(a.m * b.m), a.e + b.e);
    }
    inline Positive twice(Positive a) {
        if (a.m)
            ++a.e;
        return a;
    }
    inline Positive norm(const C &z) {
        auto a = Positive::abs(z.re.x), b = Positive::abs(z.im.x);
        if (!a.m)
            return b;
        if (!b.m)
            return a;
        if (a.e < b.e)
            std::swap(a, b);
        const auto distance = a.e - b.e;
        const double small = distance > 1074 ? std::numeric_limits<double>::denorm_min()
                                             : up(std::scalbn(b.m, -static_cast<int>(distance)));
        const double aa = up(a.m * a.m), bb = up(small * small);
        return Positive(up(std::sqrt(up(aa + bb))), a.e);
    }
    inline void selftest() {
        environment();
        v4::precision = 256;
        Real a, b, expected;
        uint64_t random = 0x947251;
        for (int i = 0; i < 10000; ++i) {
            auto next = [&] {
                random ^= random << 13;
                random ^= random >> 7;
                random ^= random << 17;
                return random;
            };
            Positive x(.5 + static_cast<double>(next() >> 12) * 0x1p-53,
                       static_cast<int64_t>(next() % 2000001) - 1000000);
            Positive y(.5 + static_cast<double>(next() >> 12) * 0x1p-53,
                       static_cast<int64_t>(next() % 2000001) - 1000000);
            if (i % 3 == 0)
                y.e = x.e + static_cast<int64_t>(next() % 1200) - 600;
            if (i % 19 == 0)
                x = {};
            if (i % 23 == 0)
                y = {};
            a = x.real();
            b = y.real();
            mpfr_add(expected.x, a.x, b.x, MPFR_RNDU);
            auto sum = add(x, y).real();
            if (mpfr_cmp(sum.x, expected.x) < 0)
                throw std::runtime_error("upper addition undershot MPFR");
            mpfr_mul(expected.x, a.x, b.x, MPFR_RNDU);
            auto product = mul(x, y).real();
            if (mpfr_cmp(product.x, expected.x) < 0)
                throw std::runtime_error("upper product undershot MPFR");
            C z;
            z.re = a;
            z.im = b;
            if (i % 2)
                mpfr_neg(z.im.x, z.im.x, MPFR_RNDN);
            mpfr_hypot(expected.x, z.re.x, z.im.x, MPFR_RNDU);
            auto magnitude = norm(z).real();
            if (mpfr_cmp(magnitude.x, expected.x) < 0)
                throw std::runtime_error("upper norm undershot MPFR");
        }
        for (int distance: {0, 1, 52, 53, 54, 1022, 1073, 1074, 1075, 3000000}) {
            Positive x(std::nextafter(1., 0.), 0), y(.5, -distance);
            a = x.real();
            b = y.real();
            mpfr_add(expected.x, a.x, b.x, MPFR_RNDU);
            auto sum = add(x, y).real();
            if (mpfr_cmp(sum.x, expected.x) < 0)
                throw std::runtime_error("boundary addition undershot MPFR");
            C z;
            z.re = a;
            z.im = b;
            mpfr_hypot(expected.x, a.x, b.x, MPFR_RNDU);
            auto magnitude = norm(z).real();
            if (mpfr_cmp(magnitude.x, expected.x) < 0)
                throw std::runtime_error("boundary norm undershot MPFR");
        }
        bool rejected = false;
        std::fesetround(FE_DOWNWARD);
        try {
            environment();
        } catch (const std::runtime_error &) {
            rejected = true;
        }
        std::fesetround(FE_TONEAREST);
        if (!rejected)
            throw std::runtime_error("wrong rounding mode accepted");
        const auto csr = _mm_getcsr();
        _mm_setcsr(csr | (1u << 15));
        rejected = false;
        try {
            environment();
        } catch (const std::runtime_error &) {
            rejected = true;
        }
        _mm_setcsr(csr);
        if (!rejected)
            throw std::runtime_error("flush-to-zero accepted");
        if (std::numeric_limits<long>::max() < (int64_t(1) << 49)) {
            rejected = false;
            try {
                Positive::power(int64_t(1) << 49).real();
            } catch (const std::runtime_error &) {
                rejected = true;
            }
            if (!rejected)
                throw std::runtime_error("large MPFR shift narrowed silently");
        }
        environment();
        v4::precision = 256;
    }
} // namespace upper
