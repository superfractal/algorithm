// Created by GPT-6 on 2026-09-10
// Modified by GPT-6 on 2026-09-10, 2026-09-12
#pragma once
#include "upper.hpp"
#include "../STMSLimits.hpp"
namespace disks {
    using namespace inverse_samples;
    struct Disk {
        C value;
        upper::Positive radius;
    };
    struct Ops {
        Real t, u, re, im, sumA, sumB, third;
        mpfr_exp_t largest;
        void reset() {
            largest = -4 * stms_limits::exponentLimit(v4::precision);
        }
        void track(mpfr_srcptr x) {
            if (!mpfr_number_p(x) ||
                (!mpfr_zero_p(x) &&
                 std::abs(mpfr_get_exp(x)) > stms_limits::exponentLimit(v4::precision)))
                throw std::runtime_error("disk exponent guard");
            if (!mpfr_zero_p(x))
                largest = std::max(largest, mpfr_get_exp(x));
        }
        upper::Positive rounding() const {
            return upper::Positive::power(largest - v4::precision + 4);
        }
        void multiply(Disk &out, const Disk &a, const Disk &b) {
            multiplyKnownNorms(out, a, b, upper::norm(a.value), upper::norm(b.value));
        }
        static upper::Positive l1(const C &z) {
            return upper::add(upper::Positive::abs(z.re.x), upper::Positive::abs(z.im.x));
        }
        void multiplyKnownNorms(
            Disk &out, const Disk &a, const Disk &b, upper::Positive normA, upper::Positive normB) {
            auto error =
                upper::add(upper::add(upper::mul(normA, b.radius), upper::mul(normB, a.radius)),
                           upper::mul(a.radius, b.radius));
            upper::Positive inputSumError;
            if (v4::precision >= 4096) {
                reset();
                mpfr_add(sumA.x, a.value.re.x, a.value.im.x, MPFR_RNDN);
                track(sumA.x);
                auto errorA = mpfr_zero_p(sumA.x)
                                  ? upper::Positive{}
                                  : upper::Positive::power(mpfr_get_exp(sumA.x) - v4::precision);
                mpfr_add(sumB.x, b.value.re.x, b.value.im.x, MPFR_RNDN);
                track(sumB.x);
                auto errorB = mpfr_zero_p(sumB.x)
                                  ? upper::Positive{}
                                  : upper::Positive::power(mpfr_get_exp(sumB.x) - v4::precision);
                inputSumError =
                    upper::add(upper::add(upper::mul(upper::Positive::abs(sumA.x), errorB),
                                          upper::mul(upper::Positive::abs(sumB.x), errorA)),
                               upper::mul(errorA, errorB));
                reset();
                mpfr_mul(t.x, a.value.re.x, b.value.re.x, MPFR_RNDN);
                track(t.x);
                mpfr_mul(u.x, a.value.im.x, b.value.im.x, MPFR_RNDN);
                track(u.x);
                mpfr_mul(third.x, sumA.x, sumB.x, MPFR_RNDN);
                track(third.x);
                mpfr_sub(re.x, t.x, u.x, MPFR_RNDN);
                track(re.x);
                mpfr_sub(im.x, third.x, t.x, MPFR_RNDN);
                track(im.x);
                mpfr_sub(im.x, im.x, u.x, MPFR_RNDN);
                track(im.x);
            } else {
                reset();
                mpfr_mul(t.x, a.value.re.x, b.value.re.x, MPFR_RNDN);
                track(t.x);
                mpfr_mul(u.x, a.value.im.x, b.value.im.x, MPFR_RNDN);
                track(u.x);
                mpfr_sub(re.x, t.x, u.x, MPFR_RNDN);
                track(re.x);
                mpfr_mul(t.x, a.value.re.x, b.value.im.x, MPFR_RNDN);
                track(t.x);
                mpfr_mul(u.x, a.value.im.x, b.value.re.x, MPFR_RNDN);
                track(u.x);
                mpfr_add(im.x, t.x, u.x, MPFR_RNDN);
                track(im.x);
            }
            mpfr_set(out.value.re.x, re.x, MPFR_RNDN);
            mpfr_set(out.value.im.x, im.x, MPFR_RNDN);
            out.radius = upper::add(upper::add(error, inputSumError), rounding());
        }
        void add(Disk &out, const Disk &a, const Disk &b) {
            auto error = upper::add(a.radius, b.radius);
            reset();
            mpfr_add(re.x, a.value.re.x, b.value.re.x, MPFR_RNDN);
            track(re.x);
            mpfr_add(im.x, a.value.im.x, b.value.im.x, MPFR_RNDN);
            track(im.x);
            mpfr_set(out.value.re.x, re.x, MPFR_RNDN);
            mpfr_set(out.value.im.x, im.x, MPFR_RNDN);
            out.radius = upper::add(error, rounding());
        }
        void addOne(Disk &out) {
            reset();
            mpfr_add_ui(out.value.re.x, out.value.re.x, 1, MPFR_RNDN);
            track(out.value.re.x);
            out.radius = upper::add(out.radius, rounding());
        }
    };
} // namespace disks
