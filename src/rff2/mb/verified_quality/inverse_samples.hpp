// Created by GPT-6 on 2026-09-09
// Modified by GPT-6 on 2026-09-09, 2026-09-10, 2026-09-12
#pragma once
#include "../STMSLimits.hpp"
#include "interval_solver.hpp"
namespace inverse_samples {
    using v4::Real;
    using v4::Clock;
    struct C {
        Real re, im;
    };
    inline C add(const C &a, const C &b) {
        C r;
        mpfr_add(r.re.x, a.re.x, b.re.x, MPFR_RNDN);
        mpfr_add(r.im.x, a.im.x, b.im.x, MPFR_RNDN);
        return r;
    }
    inline C neg(C a) {
        mpfr_neg(a.re.x, a.re.x, MPFR_RNDN);
        mpfr_neg(a.im.x, a.im.x, MPFR_RNDN);
        return a;
    }
    inline C sub(const C &a, const C &b) {
        return add(a, neg(b));
    }
    inline C mul(const C &a, const C &b) {
        C r;
        Real t, u;
        mpfr_mul(t.x, a.re.x, b.re.x, MPFR_RNDN);
        mpfr_mul(u.x, a.im.x, b.im.x, MPFR_RNDN);
        mpfr_sub(r.re.x, t.x, u.x, MPFR_RNDN);
        mpfr_mul(t.x, a.re.x, b.im.x, MPFR_RNDN);
        mpfr_mul(u.x, a.im.x, b.re.x, MPFR_RNDN);
        mpfr_add(r.im.x, t.x, u.x, MPFR_RNDN);
        return r;
    }
    inline bool zero(const C &a) {
        return mpfr_zero_p(a.re.x) && mpfr_zero_p(a.im.x);
    }
    inline bool divide(const C &a, const C &b, C &out) {
        if (zero(b))
            return false;
        C conjugate = b;
        mpfr_neg(conjugate.im.x, conjugate.im.x, MPFR_RNDN);
        auto num = mul(a, conjugate);
        Real d, t;
        mpfr_sqr(d.x, b.re.x, MPFR_RNDN);
        mpfr_sqr(t.x, b.im.x, MPFR_RNDN);
        mpfr_add(d.x, d.x, t.x, MPFR_RNDN);
        mpfr_div(out.re.x, num.re.x, d.x, MPFR_RNDN);
        mpfr_div(out.im.x, num.im.x, d.x, MPFR_RNDN);
        return true;
    }
    inline Real magnitude(const C &a) {
        Real r;
        mpfr_hypot(r.x, a.re.x, a.im.x, MPFR_RNDN);
        return r;
    }
    inline C scaled(C a, const Real &r) {
        mpfr_mul(a.re.x, a.re.x, r.x, MPFR_RNDN);
        mpfr_mul(a.im.x, a.im.x, r.x, MPFR_RNDN);
        return a;
    }
    struct Evaluation {
        C value;
        uint64_t iterations = 0;
        std::string status = "COMPLETE";
    };
    inline Evaluation orbit(const C &c, uint64_t period, const v4::Control &control) {
        Evaluation r;
        Real rr, ii, ri;
        for (uint64_t k = 0; k < period; ++k) {
            if ((k & 255) == 0 && control.stopped()) {
                r.status = "TIME_BUDGET";
                return r;
            }
            // Same rounded operations as add(mul(z,z),c), without per-step allocations.
            mpfr_sqr(rr.x, r.value.re.x, MPFR_RNDN);
            mpfr_sqr(ii.x, r.value.im.x, MPFR_RNDN);
            mpfr_mul(ri.x, r.value.re.x, r.value.im.x, MPFR_RNDN);
            mpfr_sub(rr.x, rr.x, ii.x, MPFR_RNDN);
            mpfr_mul_2ui(ri.x, ri.x, 1, MPFR_RNDN);
            mpfr_add(r.value.re.x, rr.x, c.re.x, MPFR_RNDN);
            mpfr_add(r.value.im.x, ri.x, c.im.x, MPFR_RNDN);
            ++r.iterations;
            if (!mpfr_number_p(r.value.re.x) || !mpfr_number_p(r.value.im.x)) {
                r.status = "NONFINITE";
                return r;
            }
            if ((!mpfr_zero_p(r.value.re.x) &&
                 mpfr_get_exp(r.value.re.x) > stms_limits::exponentLimit(v4::precision)) ||
                (!mpfr_zero_p(r.value.im.x) &&
                 mpfr_get_exp(r.value.im.x) > stms_limits::exponentLimit(v4::precision))) {
                r.status = "EXPONENT_LIMIT";
                return r;
            }
        }
        return r;
    }
    inline bool interpolate(const std::vector<C> &positions,
                            const std::vector<Evaluation> &values,
                            int stride,
                            C &out) {
        out = C{};
        for (size_t j = 0; j < values.size(); j += stride) {
            if (values[j].status != "COMPLETE")
                return false;
            if (zero(values[j].value)) {
                out = positions[j];
                return true;
            }
        }
        for (size_t j = 0; j < values.size(); j += stride) {
            C weight;
            weight.re = Real(1);
            for (size_t k = 0; k < values.size(); k += stride)
                if (k != j) {
                    C ratio;
                    if (!divide(neg(values[k].value), sub(values[j].value, values[k].value), ratio))
                        return false;
                    weight = mul(weight, ratio);
                }
            out = add(out, mul(positions[j], weight));
        }
        return mpfr_number_p(out.re.x) && mpfr_number_p(out.im.x);
    }
    // Fit u(y)=(a0+a1*y+a2*y^2+a3*y^3)/(1+b1*y+...+b4*y^4).
    // y is F_p normalized by the largest sample magnitude. No Taylor coefficients.
    inline bool
    rational(const std::vector<C> &positions, const std::vector<Evaluation> &values, C &out) {
        if (values.size() != 8 || positions.size() != 8)
            return false;
        Real scale;
        for (const auto &s: values) {
            if (s.status != "COMPLETE")
                return false;
            auto m = magnitude(s.value);
            if (mpfr_cmp(m.x, scale.x) > 0)
                scale = m;
        }
        if (mpfr_zero_p(scale.x))
            return false;
        std::vector<std::vector<C>> a(8, std::vector<C>(9));
        for (size_t j = 0; j < 8; ++j) {
            C y = values[j].value;
            mpfr_div(y.re.x, y.re.x, scale.x, MPFR_RNDN);
            mpfr_div(y.im.x, y.im.x, scale.x, MPFR_RNDN);
            C power;
            power.re = Real(1);
            for (int k = 0; k < 4; ++k) {
                a[j][k] = power;
                power = mul(power, y);
            }
            power = y;
            for (int k = 4; k < 8; ++k) {
                a[j][k] = neg(mul(positions[j], power));
                power = mul(power, y);
            }
            a[j][8] = positions[j];
        }
        for (size_t k = 0; k < 8; ++k) {
            size_t pivot = k;
            Real largest = magnitude(a[k][k]);
            for (size_t j = k + 1; j < 8; ++j) {
                auto m = magnitude(a[j][k]);
                if (mpfr_cmp(m.x, largest.x) > 0) {
                    largest = m;
                    pivot = j;
                }
            }
            if (mpfr_zero_p(largest.x))
                return false;
            std::swap(a[k], a[pivot]);
            C divisor = a[k][k];
            for (size_t col = k; col <= 8; ++col) {
                C q;
                if (!divide(a[k][col], divisor, q))
                    return false;
                a[k][col] = std::move(q);
            }
            for (size_t j = 0; j < 8; ++j)
                if (j != k) {
                    C factor = a[j][k];
                    for (size_t col = k; col <= 8; ++col)
                        a[j][col] = sub(a[j][col], mul(factor, a[k][col]));
                }
        }
        out = a[0][8];
        return mpfr_number_p(out.re.x) && mpfr_number_p(out.im.x);
    }
    struct Result {
        std::string status;
        C center, candidate, coarse;
        Real radius, gap, residual;
        bool hasCandidate = false, rationalUsed = false;
        uint64_t iterations = 0;
        std::vector<C> positions;
        std::vector<Evaluation> samples;
        Evaluation validation;
    };
    inline std::unique_ptr<Result> solve(merutilm::rff2::ParallelRenderState &state,
                                         const v4::Box &box,
                                         uint64_t period,
                                         double budget) {
        auto result = std::make_unique<Result>();
        const auto mid = v4::enclosingDisk(box).center;
        result->center.re = mid.re.lo;
        result->center.im = mid.im.lo;
        result->radius = v4::width(box);
        mpfr_div_2ui(result->radius.x, result->radius.x, 1, MPFR_RNDN);
        const auto start = Clock::now();
        v4::Control control{state,
                            start + std::chrono::duration_cast<Clock::duration>(
                                        std::chrono::duration<double>(budget))};
        Real pi;
        mpfr_const_pi(pi.x, MPFR_RNDN);
        for (int j = 0; j < 8; ++j) {
            Real angle;
            mpfr_mul_ui(angle.x, pi.x, j, MPFR_RNDN);
            mpfr_div_ui(angle.x, angle.x, 4, MPFR_RNDN);
            C u;
            mpfr_sin_cos(u.im.x, u.re.x, angle.x, MPFR_RNDN);
            result->positions.push_back(std::move(u));
        }
        result->samples.resize(8);
        std::vector<std::jthread> workers;
        for (int j = 0; j < 8; ++j)
            workers.emplace_back([&, j] {
                result->samples[j] =
                    orbit(add(result->center, scaled(result->positions[j], result->radius)),
                          period,
                          control);
                mpfr_free_cache();
            });
        for (auto &t: workers)
            t.join();
        if (state.interruptRequested())
            return nullptr;
        bool complete = true;
        for (const auto &s: result->samples) {
            result->iterations += s.iterations;
            complete &= s.status == "COMPLETE";
        }
        if (!complete) {
            result->status = "INCOMPLETE_SAMPLES";
            return result;
        }
        C coarse, fine;
        if (!interpolate(result->positions, result->samples, 1, coarse)) {
            result->status = "SINGULAR_INVERSE_INTERPOLATION";
            return result;
        }
        result->rationalUsed = rational(result->positions, result->samples, fine);
        if (!result->rationalUsed)
            fine = coarse; // Degenerate rational system: explicit polynomial fallback.
        result->coarse = add(result->center, scaled(coarse, result->radius));
        result->candidate = add(result->center, scaled(fine, result->radius));
        result->hasCandidate = true;
        result->gap = magnitude(scaled(sub(fine, coarse), result->radius));
        result->validation = orbit(result->candidate, period, control);
        result->iterations += result->validation.iterations;
        if (state.interruptRequested())
            return nullptr;
        if (result->validation.status != "COMPLETE") {
            result->status = "INCOMPLETE_VALIDATION";
            return result;
        }
        result->residual = magnitude(result->validation.value);
        result->status = "CANDIDATE_UNVALIDATED";
        return result;
    }
} // namespace inverse_samples
