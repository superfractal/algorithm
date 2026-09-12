// Created by GPT-6 on 2026-09-10
// Modified by GPT-6 on 2026-09-10, 2026-09-12
#pragma once
#include "upper.hpp"
#include "../STMSLimits.hpp"
namespace enclosure {
    using namespace inverse_samples;
    struct Orbit {
        C value;
        upper::Positive radiusBound, criticalBound = upper::Positive::power(0);
        upper::Positive previousCenterNorm;
        Real rr, ii, ri, sum, difference;
        uint64_t iterations = 0;
        std::string status = "COMPLETE";
        Real radiusValue() const {
            return radiusBound.real();
        }
        bool step(const C &c) {
            const auto bits = mpfr_get_prec(value.re.x);
            if (bits > stms_limits::maxVerificationBits) {
                status = "PRECISION_LIMIT";
                return false;
            }
            const auto exponentLimit = stms_limits::exponentLimit(bits);
            if (mpfr_get_emin() >= -4 * exponentLimit || mpfr_get_emax() <= 4 * exponentLimit) {
                status = "EXPONENT_ENVIRONMENT";
                return false;
            }
            for (auto *x:
                 std::initializer_list<mpfr_srcptr>{value.re.x, value.im.x, c.re.x, c.im.x})
                if (!mpfr_number_p(x) ||
                    (!mpfr_zero_p(x) && std::abs(mpfr_get_exp(x)) > exponentLimit)) {
                    status = "STATE_RANGE";
                    return false;
                }
            const auto norm = upper::norm(value);
            previousCenterNorm = norm;
            if (iterations)
                criticalBound =
                    upper::mul(criticalBound, upper::twice(upper::add(norm, radiusBound)));
            radiusBound = upper::mul(radiusBound, upper::add(upper::twice(norm), radiusBound));
            mpfr_exp_t largest = -4 * exponentLimit;
            auto track = [&](mpfr_srcptr x) {
                if (!mpfr_zero_p(x))
                    largest = std::max(largest, mpfr_get_exp(x));
            };
            upper::Positive inputError;
            if (bits >= 4096) {
                mpfr_add(sum.x, value.re.x, value.im.x, MPFR_RNDN);
                mpfr_sub(difference.x, value.re.x, value.im.x, MPFR_RNDN);
                const auto es = mpfr_zero_p(sum.x)
                                    ? upper::Positive{}
                                    : upper::Positive::power(mpfr_get_exp(sum.x) - bits);
                const auto ed = mpfr_zero_p(difference.x)
                                    ? upper::Positive{}
                                    : upper::Positive::power(mpfr_get_exp(difference.x) - bits);
                inputError =
                    upper::add(upper::add(upper::mul(upper::Positive::abs(sum.x), ed),
                                          upper::mul(upper::Positive::abs(difference.x), es)),
                               upper::mul(es, ed));
                mpfr_mul(rr.x, sum.x, difference.x, MPFR_RNDN);
                track(rr.x);
                mpfr_mul(ri.x, value.re.x, value.im.x, MPFR_RNDN);
                track(ri.x);
            } else {
                mpfr_sqr(rr.x, value.re.x, MPFR_RNDN);
                track(rr.x);
                mpfr_sqr(ii.x, value.im.x, MPFR_RNDN);
                track(ii.x);
                mpfr_mul(ri.x, value.re.x, value.im.x, MPFR_RNDN);
                track(ri.x);
                mpfr_sub(rr.x, rr.x, ii.x, MPFR_RNDN);
                track(rr.x);
            }
            mpfr_mul_2ui(ri.x, ri.x, 1, MPFR_RNDN);
            mpfr_add(value.re.x, rr.x, c.re.x, MPFR_RNDN);
            track(value.re.x);
            mpfr_add(value.im.x, ri.x, c.im.x, MPFR_RNDN);
            track(value.im.x);
            // Eight units cover either kernel; the two-product path additionally
            // propagates rounded input-sum errors.
            radiusBound = upper::add(upper::add(radiusBound, inputError),
                                     upper::Positive::power(largest - bits + 3));
            ++iterations;
            if (radiusBound.m &&
                (radiusBound.e > -1 || (radiusBound.e == -1 && radiusBound.m > .5))) {
                status = "ENCLOSURE_TOO_WIDE";
                return false;
            }
            return true;
        }
        Real scaledUpper() {
            return upper::mul(criticalBound, upper::add(upper::norm(value), radiusBound)).real();
        }
    };
    inline Orbit run(const C &c, uint64_t period, const v4::Control &control) {
        upper::environment();
        Orbit orbit;
        for (uint64_t k = 0; k < period; ++k) {
            if ((k & 255) == 0 && control.stopped()) {
                orbit.status = "TIME_BUDGET";
                break;
            }
            if (!orbit.step(c))
                break;
        }
        return orbit;
    }
    inline void selftest() {
        upper::selftest();
        for (auto bits: {64, 163, 919})
            for (auto pair: {std::pair{-.125, .25},
                             std::pair{.125, .03125},
                             std::pair{-2., 0.},
                             std::pair{0., 0.}}) {
                v4::precision = bits;
                C c;
                c.re = Real(pair.first);
                c.im = Real(pair.second);
                Orbit ball;
                for (int k = 0; k < 16; ++k)
                    if (!ball.step(c))
                        throw std::runtime_error("ball fixture incomplete");
                v4::precision = bits + 128;
                v4::Box z, parameter;
                parameter.re = v4::Interval(pair.first, pair.first);
                parameter.im = v4::Interval(pair.second, pair.second);
                for (int k = 0; k < 16; ++k) {
                    v4::Box next;
                    next.re = v4::add(v4::sub(v4::square(z.re), v4::square(z.im)), parameter.re);
                    next.im =
                        v4::add(v4::mul(v4::Interval(2, 2), v4::mul(z.re, z.im)), parameter.im);
                    z = std::move(next);
                }
                Real dx, dy, delta, temp;
                for (auto *x: {z.re.lo.x, z.re.hi.x}) {
                    mpfr_sub(temp.x, x, ball.value.re.x, MPFR_RNDA);
                    mpfr_abs(temp.x, temp.x, MPFR_RNDU);
                    mpfr_max(dx.x, dx.x, temp.x, MPFR_RNDU);
                }
                for (auto *x: {z.im.lo.x, z.im.hi.x}) {
                    mpfr_sub(temp.x, x, ball.value.im.x, MPFR_RNDA);
                    mpfr_abs(temp.x, temp.x, MPFR_RNDU);
                    mpfr_max(dy.x, dy.x, temp.x, MPFR_RNDU);
                }
                mpfr_hypot(delta.x, dx.x, dy.x, MPFR_RNDU);
                if (mpfr_cmp(delta.x, ball.radiusValue().x) > 0)
                    throw std::runtime_error("independent directed interval outside ball");
            }
        // At c=-2 the intentionally conservative local-error budget grows under
        // multiplier four. This must become inconclusive, not silently accept.
        v4::precision = 64;
        C boundary;
        boundary.re = Real(-2);
        Orbit wide;
        for (int k = 0; k < 64 && wide.status == "COMPLETE"; ++k)
            wide.step(boundary);
        if (wide.status != "ENCLOSURE_TOO_WIDE")
            throw std::runtime_error("wide enclosure guard fixture");
        v4::precision = 256;
    }
} // namespace enclosure
