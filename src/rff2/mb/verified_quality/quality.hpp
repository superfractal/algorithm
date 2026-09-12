// Created by GPT-6 on 2026-09-09
// Modified by GPT-6 on 2026-09-09, 2026-09-10, 2026-09-12
#pragma once
#include "../STMSLimits.hpp"
#include "inverse_samples.hpp"
#include <condition_variable>
#include <deque>
#include <mutex>
namespace quality {
    using namespace inverse_samples;
    // Diagnostic state is set by the bridge under its invocation mutex.
    struct Checkpoint {
        uint64_t index;
        C value;
        double offset;
    };
    struct CheckpointProbe {
        uint64_t index;
        double logDistance, normalized;
    };
    inline std::vector<Checkpoint> checkpoints;
    inline std::vector<CheckpointProbe> probes;
    inline size_t nextCheckpoint = 0;
    inline void observe(uint64_t k, const C &z) {
        if (nextCheckpoint >= checkpoints.size() || checkpoints[nextCheckpoint].index != k)
            return;
        const auto &seed = checkpoints[nextCheckpoint++];
        auto distance = magnitude(sub(seed.value, z));
        double logDistance = mpfr_zero_p(distance.x) ? -INFINITY : v4::log10value(distance);
        probes.push_back({k, logDistance, logDistance - seed.offset / std::log2(10.)});
    }
    struct Evaluation : inverse_samples::Evaluation {
        C derivative, critical;
        Real scaledResidual;
        double logZoom = NAN;
        bool usable = false;
        Evaluation() {
            status = "NOT_EVALUATED";
        }
    };
    inline Evaluation evaluateSerial(const C &c, uint64_t period, const v4::Control &control) {
        Evaluation r;
        r.status = "COMPLETE";
        r.critical.re = Real(1);
        Real two(2);
        C one;
        one.re = Real(1);
        for (uint64_t k = 0; k < period; ++k) {
            if ((k & 255) == 0 && control.stopped()) {
                r.status = "TIME_BUDGET";
                return r;
            }
            observe(k, r.value);
            C twice = scaled(r.value, two);
            if (k)
                r.critical = mul(r.critical, twice);
            r.derivative = add(mul(twice, r.derivative), one);
            r.value = add(mul(r.value, r.value), c);
            ++r.iterations;
            for (auto *x: {r.value.re.x,
                           r.value.im.x,
                           r.derivative.re.x,
                           r.derivative.im.x,
                           r.critical.re.x,
                           r.critical.im.x}) {
                if (!mpfr_number_p(x)) {
                    r.status = "NONFINITE";
                    return r;
                }
                if (!mpfr_zero_p(x) &&
                    mpfr_get_exp(x) > stms_limits::exponentLimit(v4::precision)) {
                    r.status = "EXPONENT_LIMIT";
                    return r;
                }
            }
        }
        r.scaledResidual = magnitude(mul(r.critical, r.value));
        Real ad = magnitude(mul(r.critical, r.derivative));
        if (mpfr_sgn(ad.x) > 0)
            r.logZoom = v4::log10value(ad) + 2;
        r.usable =
            std::isfinite(r.logZoom) && r.logZoom > 0 && mpfr_cmp_d(r.scaledResidual.x, 1e-12) < 0;
        return r;
    }
    // Compose affine derivative maps (a,b): d_out=a*d_in+b.
    struct Summary {
        C a, b;
        bool complete = false;
    };
    struct Tile {
        std::vector<C> coefficients;
        size_t count = 0, index = 0;
        Tile() : coefficients(512) {
        }
    };
    inline void multiplyInto(C &out, const C &a, const C &b, Real &t, Real &u, Real &re, Real &im) {
        mpfr_mul(t.x, a.re.x, b.re.x, MPFR_RNDN);
        mpfr_mul(u.x, a.im.x, b.im.x, MPFR_RNDN);
        mpfr_sub(re.x, t.x, u.x, MPFR_RNDN);
        mpfr_mul(t.x, a.re.x, b.im.x, MPFR_RNDN);
        mpfr_mul(u.x, a.im.x, b.re.x, MPFR_RNDN);
        mpfr_add(im.x, t.x, u.x, MPFR_RNDN);
        mpfr_set(out.re.x, re.x, MPFR_RNDN);
        mpfr_set(out.im.x, im.x, MPFR_RNDN);
    }
    inline Evaluation evaluate(const C &c, uint64_t period, const v4::Control &control) {
        if (period <= 512)
            return evaluateSerial(c, period, control);
        Evaluation result;
        if (control.stopped()) {
            result.status = "TIME_BUDGET";
            return result;
        }
        const size_t tileCount = (period - 1 + 511) / 512;
        std::vector<Summary> summaries(tileCount);
        std::vector<std::unique_ptr<Tile>> pool;
        std::deque<Tile *> free, ready;
        for (size_t j = 0; j < std::min<size_t>(8, tileCount); ++j) {
            pool.push_back(std::make_unique<Tile>());
            free.push_back(pool.back().get());
        }
        std::mutex mutex;
        std::condition_variable cv;
        bool finished = false;
        std::atomic<bool> failed = false;
        std::vector<std::jthread> workers;
        workers.reserve(7);
        auto join = [&] {
            {
                std::lock_guard lock(mutex);
                finished = true;
            }
            cv.notify_all();
            for (auto &t: workers)
                t.join();
        };
        try {
            for (size_t j = 0; j < std::min<size_t>(7, tileCount); ++j)
                workers.emplace_back([&] {
                    try {
                        Real t, u, re, im;
                        for (;;) {
                            Tile *tile = nullptr;
                            {
                                std::unique_lock lock(mutex);
                                while (ready.empty() && !finished && !failed) {
                                    cv.wait_for(lock, std::chrono::milliseconds(10));
                                    if (control.stopped())
                                        failed = true;
                                }
                                if (failed) {
                                    cv.notify_all();
                                    break;
                                }
                                if (ready.empty() && finished)
                                    break;
                                tile = ready.front();
                                ready.pop_front();
                            }
                            auto &s = summaries[tile->index];
                            mpfr_set_ui(s.a.re.x, 1, MPFR_RNDN);
                            for (size_t k = 0; k < tile->count; ++k) {
                                if ((k & 63) == 0 && control.stopped()) {
                                    failed = true;
                                    break;
                                }
                                multiplyInto(s.a, tile->coefficients[k], s.a, t, u, re, im);
                                multiplyInto(s.b, tile->coefficients[k], s.b, t, u, re, im);
                                mpfr_add_ui(s.b.re.x, s.b.re.x, 1, MPFR_RNDN);
                                if (!mpfr_number_p(s.a.re.x) || !mpfr_number_p(s.a.im.x) ||
                                    !mpfr_number_p(s.b.re.x) || !mpfr_number_p(s.b.im.x)) {
                                    failed = true;
                                    break;
                                }
                            }
                            s.complete = !failed;
                            {
                                std::lock_guard lock(mutex);
                                free.push_back(tile);
                            }
                            cv.notify_all();
                            if (failed)
                                break;
                        }
                    } catch (...) {
                        failed = true;
                        cv.notify_all();
                    }
                    mpfr_free_cache();
                });
            Tile *current = nullptr;
            size_t nextTile = 0;
            Real rr, ii, ri;
            for (uint64_t k = 0; k < period; ++k) {
                if ((k & 255) == 0 && (control.stopped() || failed)) {
                    failed = true;
                    break;
                }
                observe(k, result.value);
                if (k) {
                    if (!current) {
                        std::unique_lock lock(mutex);
                        while (free.empty() && !failed) {
                            cv.wait_for(lock, std::chrono::milliseconds(10));
                            if (control.stopped())
                                failed = true;
                        }
                        if (failed)
                            break;
                        current = free.front();
                        free.pop_front();
                        current->count = 0;
                        current->index = nextTile++;
                    }
                    auto &m = current->coefficients[current->count++];
                    mpfr_mul_2ui(m.re.x, result.value.re.x, 1, MPFR_RNDN);
                    mpfr_mul_2ui(m.im.x, result.value.im.x, 1, MPFR_RNDN);
                    if (current->count == 512 || k + 1 == period) {
                        {
                            std::lock_guard lock(mutex);
                            ready.push_back(current);
                        }
                        current = nullptr;
                        cv.notify_all();
                    }
                }
                mpfr_sqr(rr.x, result.value.re.x, MPFR_RNDN);
                mpfr_sqr(ii.x, result.value.im.x, MPFR_RNDN);
                mpfr_mul(ri.x, result.value.re.x, result.value.im.x, MPFR_RNDN);
                mpfr_sub(rr.x, rr.x, ii.x, MPFR_RNDN);
                mpfr_mul_2ui(ri.x, ri.x, 1, MPFR_RNDN);
                mpfr_add(result.value.re.x, rr.x, c.re.x, MPFR_RNDN);
                mpfr_add(result.value.im.x, ri.x, c.im.x, MPFR_RNDN);
                ++result.iterations;
                for (auto *x: {result.value.re.x, result.value.im.x})
                    if (!mpfr_number_p(x) ||
                        (!mpfr_zero_p(x) &&
                         mpfr_get_exp(x) > stms_limits::exponentLimit(v4::precision))) {
                        failed = true;
                        break;
                    }
            }
            if (result.iterations != period)
                failed = true;
            join();
        } catch (...) {
            failed = true;
            join();
            throw;
        }
        if (failed) {
            result.status = control.stopped() ? "TIME_BUDGET" : "INCOMPLETE_REDUCTION";
            return result;
        }
        result.critical.re = Real(1);
        result.derivative.re = Real(1);
        for (size_t j = 0; j < summaries.size(); ++j) {
            if ((j & 255) == 0 && control.stopped()) {
                result.status = "TIME_BUDGET";
                return result;
            }
            const auto &s = summaries[j];
            if (!s.complete) {
                result.status = "INCOMPLETE_REDUCTION";
                return result;
            }
            result.derivative = add(mul(s.a, result.derivative), s.b);
            result.critical = mul(s.a, result.critical);
        }
        result.status = "COMPLETE";
        result.scaledResidual = magnitude(mul(result.critical, result.value));
        Real ad = magnitude(mul(result.critical, result.derivative));
        if (mpfr_sgn(ad.x) > 0)
            result.logZoom = v4::log10value(ad) + 2;
        result.usable = std::isfinite(result.logZoom) && result.logZoom > 0 &&
                        mpfr_cmp_d(result.scaledResidual.x, 1e-12) < 0;
        return result;
    }
} // namespace quality
