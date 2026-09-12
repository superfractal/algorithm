// Created by GPT-6 on 2026-09-08
// Independent interval-constraint feasibility experiment, not a production locator.
// Modified by GPT-6 on 2026-09-12
#pragma once
#include <mpfr.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <memory>
#include <thread>
#include <vector>
#include "ParallelRenderState.h"

namespace v4 {
    inline mpfr_prec_t precision =
        256; // Set before launching workers; immutable during each solve.
    struct Real {
        mpfr_t x;
        Real() {
            mpfr_init2(x, precision);
            mpfr_set_zero(x, 1);
        }
        explicit Real(double v) : Real() {
            mpfr_set_d(x, v, MPFR_RNDN);
        }
        Real(const Real &v) : Real() {
            mpfr_set(x, v.x, MPFR_RNDN);
        }
        Real(Real &&v) noexcept {
            mpfr_init2(x, 2);
            mpfr_swap(x, v.x);
        }
        Real &operator=(Real v) {
            mpfr_swap(x, v.x);
            return *this;
        }
        ~Real() {
            mpfr_clear(x);
        }
    };
    struct Interval {
        Real lo, hi;
        Interval() = default;
        Interval(double l, double h) : lo(l), hi(h) {
        }
        bool empty() const {
            return mpfr_cmp(lo.x, hi.x) > 0;
        }
    };
    struct Box {
        Interval re, im;
        Box() = default;
        Box(double rl, double rh, double il, double ih) : re(rl, rh), im(il, ih) {
        }
        bool empty() const {
            return re.empty() || im.empty();
        }
    };
    inline Interval intersect(const Interval &a, const Interval &b) {
        Interval r;
        mpfr_max(r.lo.x, a.lo.x, b.lo.x, MPFR_RNDD);
        mpfr_min(r.hi.x, a.hi.x, b.hi.x, MPFR_RNDU);
        return r;
    }
    inline Box intersect(const Box &a, const Box &b) {
        Box r;
        r.re = intersect(a.re, b.re);
        r.im = intersect(a.im, b.im);
        return r;
    }
    inline Interval add(const Interval &a, const Interval &b) {
        Interval r;
        mpfr_add(r.lo.x, a.lo.x, b.lo.x, MPFR_RNDD);
        mpfr_add(r.hi.x, a.hi.x, b.hi.x, MPFR_RNDU);
        return r;
    }
    inline Interval sub(const Interval &a, const Interval &b) {
        Interval r;
        mpfr_sub(r.lo.x, a.lo.x, b.hi.x, MPFR_RNDD);
        mpfr_sub(r.hi.x, a.hi.x, b.lo.x, MPFR_RNDU);
        return r;
    }
    inline Interval neg(const Interval &a) {
        Interval r;
        mpfr_neg(r.lo.x, a.hi.x, MPFR_RNDD);
        mpfr_neg(r.hi.x, a.lo.x, MPFR_RNDU);
        return r;
    }
    inline Box add(const Box &a, const Box &b) {
        Box r;
        r.re = add(a.re, b.re);
        r.im = add(a.im, b.im);
        return r;
    }
    inline Box sub(const Box &a, const Box &b) {
        Box r;
        r.re = sub(a.re, b.re);
        r.im = sub(a.im, b.im);
        return r;
    }
    inline Box neg(const Box &a) {
        Box r;
        r.re = neg(a.re);
        r.im = neg(a.im);
        return r;
    }
    inline Interval mul(const Interval &a, const Interval &b) {
        Interval r;
        Real tmp;
        mpfr_set_inf(r.lo.x, 1);
        mpfr_set_inf(r.hi.x, -1);
        for (auto x: {a.lo.x, a.hi.x})
            for (auto y: {b.lo.x, b.hi.x}) {
                mpfr_mul(tmp.x, x, y, MPFR_RNDD);
                mpfr_min(r.lo.x, r.lo.x, tmp.x, MPFR_RNDD);
                mpfr_mul(tmp.x, x, y, MPFR_RNDU);
                mpfr_max(r.hi.x, r.hi.x, tmp.x, MPFR_RNDU);
            }
        return r;
    }
    inline Interval square(const Interval &a) {
        Interval r;
        Real x, y;
        mpfr_sqr(x.x, a.lo.x, MPFR_RNDD);
        mpfr_sqr(y.x, a.hi.x, MPFR_RNDD);
        if (mpfr_sgn(a.lo.x) <= 0 && mpfr_sgn(a.hi.x) >= 0)
            mpfr_set_zero(r.lo.x, 1);
        else
            mpfr_min(r.lo.x, x.x, y.x, MPFR_RNDD);
        mpfr_sqr(x.x, a.lo.x, MPFR_RNDU);
        mpfr_sqr(y.x, a.hi.x, MPFR_RNDU);
        mpfr_max(r.hi.x, x.x, y.x, MPFR_RNDU);
        return r;
    }
    inline Interval half(const Interval &a) {
        Interval r;
        mpfr_div_2ui(r.lo.x, a.lo.x, 1, MPFR_RNDD);
        mpfr_div_2ui(r.hi.x, a.hi.x, 1, MPFR_RNDU);
        return r;
    }
    inline Interval nonnegativeSqrt(Interval a) {
        Interval r; // The true radicand is nonnegative; intersection with [0,+inf) is valid.
        if (mpfr_sgn(a.lo.x) < 0)
            mpfr_set_zero(a.lo.x, 1);
        if (mpfr_sgn(a.hi.x) < 0) {
            r.lo = Real(1);
            r.hi = Real(0);
            return r;
        }
        mpfr_sqrt(r.lo.x, a.lo.x, MPFR_RNDD);
        mpfr_sqrt(r.hi.x, a.hi.x, MPFR_RNDU);
        return r;
    }
    inline Box square(const Box &a) {
        Box r;
        r.re = sub(square(a.re), square(a.im));
        r.im = mul(a.re, a.im);
        mpfr_mul_2ui(r.im.lo.x, r.im.lo.x, 1, MPFR_RNDD);
        mpfr_mul_2ui(r.im.hi.x, r.im.hi.x, 1, MPFR_RNDU);
        return r;
    }
    inline bool same(const Interval &a, const Interval &b) {
        return mpfr_equal_p(a.lo.x, b.lo.x) && mpfr_equal_p(a.hi.x, b.hi.x);
    }
    inline bool same(const Box &a, const Box &b) {
        return same(a.re, b.re) && same(a.im, b.im);
    }
    inline bool contains(const Interval &a, const Real &v) {
        return mpfr_cmp(a.lo.x, v.x) <= 0 && mpfr_cmp(v.x, a.hi.x) <= 0;
    }
    inline bool contains(const Box &a, double re, double im) {
        return contains(a.re, Real(re)) && contains(a.im, Real(im));
    }
    inline Interval decimal(const std::string &s) {
        Interval r;
        mpfr_set_str(r.lo.x, s.c_str(), 10, MPFR_RNDD);
        mpfr_set_str(r.hi.x, s.c_str(), 10, MPFR_RNDU);
        return r;
    }
    inline Box hull(const Box &a, const Box &b) {
        if (a.empty())
            return b;
        if (b.empty())
            return a;
        Box r;
        mpfr_min(r.re.lo.x, a.re.lo.x, b.re.lo.x, MPFR_RNDD);
        mpfr_max(r.re.hi.x, a.re.hi.x, b.re.hi.x, MPFR_RNDU);
        mpfr_min(r.im.lo.x, a.im.lo.x, b.im.lo.x, MPFR_RNDD);
        mpfr_max(r.im.hi.x, a.im.hi.x, b.im.hi.x, MPFR_RNDU);
        return r;
    }
    struct Counts {
        uint64_t forward = 0, backward = 0, skipped = 0, branches = 0, hulls = 0;
    };
    inline Box preimage(const Box &w, const Box &limit, Counts &counts) {
        const auto radius = nonnegativeSqrt(add(square(w.re), square(w.im)));
        const auto u = nonnegativeSqrt(half(add(radius, w.re)));
        const auto v = nonnegativeSqrt(half(sub(radius, w.re)));
        Box result(1, 0, 1, 0);
        auto consider = [&](const Interval &imag) {
            Box principal;
            principal.re = u;
            principal.im = imag;
            for (const auto &branch: {principal, neg(principal)}) {
                auto clipped = intersect(branch, limit);
                if (!clipped.empty()) {
                    ++counts.branches;
                    if (!result.empty())
                        ++counts.hulls;
                    result = hull(result, clipped);
                }
            }
        };
        if (mpfr_sgn(w.im.hi.x) >= 0)
            consider(v);
        if (mpfr_sgn(w.im.lo.x) <= 0)
            consider(neg(v));
        return result; // All surviving branches are enclosed, never selected by proximity.
    }
    inline Real width(const Box &b) {
        Real a, c;
        mpfr_sub(a.x, b.re.hi.x, b.re.lo.x, MPFR_RNDU);
        mpfr_sub(c.x, b.im.hi.x, b.im.lo.x, MPFR_RNDU);
        mpfr_max(a.x, a.x, c.x, MPFR_RNDU);
        return a;
    }
    inline double log10value(const Real &r) {
        if (mpfr_zero_p(r.x))
            return -INFINITY;
        long e = 0;
        const double m = mpfr_get_d_2exp(&e, r.x, MPFR_RNDU);
        return std::log10(m) + e * std::log10(2.0);
    }
    inline std::string text(const Real &r) {
        char *p = nullptr;
        mpfr_asprintf(&p,
                      "%.*Re",
                      static_cast<int>(std::ceil(mpfr_get_prec(r.x) * std::log10(2.0))) + 3,
                      r.x);
        std::string s = p;
        mpfr_free_str(p);
        return s;
    }
    using Clock = std::chrono::steady_clock;
    struct Disk {
        Box center;
        Real radius;
    };
    inline Real normUpper(const Box &b) {
        return nonnegativeSqrt(add(square(b.re), square(b.im))).hi;
    }
    inline Disk enclosingDisk(const Box &b) {
        Disk d;
        for (auto pair: {std::pair{&b.re, &d.center.re}, std::pair{&b.im, &d.center.im}}) {
            mpfr_add(pair.second->lo.x, pair.first->lo.x, pair.first->hi.x, MPFR_RNDN);
            mpfr_div_2ui(pair.second->lo.x, pair.second->lo.x, 1, MPFR_RNDN);
            pair.second->hi = pair.second->lo;
        }
        d.radius = normUpper(sub(b, d.center));
        return d;
    }
    inline Box diskBox(const Disk &d) {
        Box b;
        Interval radius;
        mpfr_neg(radius.lo.x, d.radius.x, MPFR_RNDD);
        radius.hi = d.radius;
        b.re = add(d.center.re, radius);
        b.im = add(d.center.im, radius);
        return b;
    }
    inline Disk diskStep(const Disk &d, const Disk &c) {
        // Exact quadratic identity: retain the entire r^2 term and all rounding error.
        auto next = enclosingDisk(add(square(d.center), c.center));
        Real growth, tmp;
        auto magnitude = normUpper(d.center);
        mpfr_mul(growth.x, magnitude.x, d.radius.x, MPFR_RNDU);
        mpfr_mul_2ui(growth.x, growth.x, 1, MPFR_RNDU);
        mpfr_sqr(tmp.x, d.radius.x, MPFR_RNDU);
        mpfr_add(growth.x, growth.x, tmp.x, MPFR_RNDU);
        mpfr_add(growth.x, growth.x, c.radius.x, MPFR_RNDU);
        mpfr_add(next.radius.x, next.radius.x, growth.x, MPFR_RNDU);
        return next;
    }
    inline Disk tightenDisk(Disk d, const Box &b) {
        auto alternative = enclosingDisk(b);
        if (mpfr_cmp(alternative.radius.x, d.radius.x) < 0)
            return alternative;
        return d;
    }
    inline bool sameDisk(const Disk &a, const Disk &b) {
        return same(a.center, b.center) && mpfr_equal_p(a.radius.x, b.radius.x);
    }
    struct Control {
        merutilm::rff2::ParallelRenderState &state;
        Clock::time_point deadline;
        bool stopped() const {
            return state.interruptRequested() || Clock::now() >= deadline;
        }
    };
    // Preserve distinct inverse branches instead of taking their rectangular hull at every step.
    inline std::vector<Box> inverseBranches(const Box &w, const Box &limit, Counts &counts) {
        auto radius = nonnegativeSqrt(add(square(w.re), square(w.im)));
        auto u = nonnegativeSqrt(half(add(radius, w.re)));
        auto v = nonnegativeSqrt(half(sub(radius, w.re)));
        std::vector<Box> branches;
        for (int sign: {1, -1}) {
            if (sign == 1 && mpfr_sgn(w.im.hi.x) < 0)
                continue;
            if (sign == -1 && mpfr_sgn(w.im.lo.x) > 0)
                continue;
            Box principal;
            principal.re = u;
            principal.im = sign == 1 ? v : neg(v);
            for (const auto &candidate: {principal, neg(principal)}) {
                auto b = intersect(candidate, limit);
                if (b.empty())
                    continue;
                bool duplicate = false;
                for (const auto &existing: branches)
                    duplicate |= same(existing, b);
                if (!duplicate) {
                    branches.push_back(std::move(b));
                    ++counts.branches;
                }
            }
        }
        return branches;
    }
    struct BlockResult {
        Box entrance, exit, param;
        Counts counts;
        bool interrupted = false, branchLimited = false, memoryLimited = false;
    };
    inline BlockResult
    contract(const Box &in, const Box &out, const Box &c, uint64_t length, const Control &control) {
        const Box universe(-2, 2, -2, 2);
        BlockResult r;
        r.param = c;
        Box cur = in;
        std::vector<Box> prefix{in};
        Box repeated = universe;
        auto orbitDisk = enclosingDisk(in);
        const auto parameterDisk = enclosingDisk(c);
        // Account conservatively for vector capacity growth as well as MPFR limbs.
        const uint64_t bytesPerEntry = 2 * (4 * ((precision + 63) / 64) * 8 + sizeof(Box) + 512);
        const uint64_t prefixLimit = std::max<uint64_t>(2, (64ULL * 1024 * 1024) / bytesPerEntry);
        for (uint64_t n = 0; n < length; ++n) {
            if (control.stopped()) {
                r.interrupted = true;
                return r;
            }
            auto nextDisk = diskStep(orbitDisk, parameterDisk);
            auto next = intersect(intersect(add(square(cur), c), diskBox(nextDisk)), universe);
            ++r.counts.forward;
            if (!next.empty())
                nextDisk = tightenDisk(std::move(nextDisk), next);
            if (n + 1 == length) {
                r.param = intersect(r.param, sub(out, square(cur)));
                next = intersect(next, out);
            } else if (same(next, cur) && sameDisk(nextDisk, orbitDisk)) {
                // Exact fixed point of this rounded box map: remaining repeated maps are identical.
                r.counts.skipped += length - n - 1;
                r.param = intersect(r.param, sub(out, square(cur)));
                repeated = cur;
                cur = intersect(cur, out);
                break;
            }
            cur = std::move(next);
            orbitDisk = std::move(nextDisk);
            if (prefix.size() >= prefixLimit) {
                r.memoryLimited = true;
                return r;
            }
            prefix.push_back(cur);
            if (cur.empty())
                break;
        }
        r.exit = cur;
        if (cur.empty()) {
            r.entrance = cur;
            return r;
        }
        std::vector<Box> branches{out};
        for (uint64_t n = 0; n < length; ++n) {
            if (control.stopped()) {
                r.interrupted = true;
                return r;
            }
            const bool last = n + 1 == length;
            const auto index = length - n - 1;
            const Box &limit = index < prefix.size() ? prefix[index] : repeated;
            std::vector<Box> nextBranches;
            Box possibleC(1, 0, 1, 0);
            for (const auto &branch: branches) {
                if (control.stopped()) {
                    r.interrupted = true;
                    return r;
                }
                auto candidates = inverseBranches(sub(branch, c), limit, r.counts);
                ++r.counts.backward;
                if (last && !candidates.empty())
                    possibleC = hull(possibleC, intersect(c, sub(branch, square(in))));
                for (auto &candidate: candidates) {
                    bool duplicate = false;
                    for (const auto &existing: nextBranches)
                        duplicate |= same(existing, candidate);
                    if (!duplicate)
                        nextBranches.push_back(std::move(candidate));
                    // Stop explicitly: never discard a possible branch or claim a located center.
                    if (nextBranches.size() > 256) {
                        r.branchLimited = true;
                        return r;
                    }
                }
            }
            if (last)
                r.param = intersect(r.param, possibleC);
            branches = std::move(nextBranches);
            if (branches.empty())
                break;
        }
        r.entrance = Box(1, 0, 1, 0);
        for (const auto &branch: branches) {
            if (!r.entrance.empty())
                ++r.counts.hulls;
            r.entrance = hull(r.entrance, branch);
        }
        return r;
    }
    struct Round {
        int index;
        double seconds, parameterDigits;
        uint64_t changed;
        Counts counts;
    };
    struct Result {
        Box parameter;
        std::string status;
        std::vector<Round> rounds;
        Counts totals;
        double seconds = 0, seedSeconds = 0;
        int segments = 0;
        uint64_t seedEvaluations = 0, seedSkipped = 0, firstBroad = 0;
        std::vector<uint64_t> cuts;
        std::vector<double> seedWidths;
    };
    inline std::unique_ptr<Result> solveConnections(merutilm::rff2::ParallelRenderState &state,
                                                    Box c,
                                                    uint64_t period,
                                                    const Real &target,
                                                    int segmentLimit,
                                                    int workerLimit,
                                                    int maxRounds,
                                                    double budgetSeconds) {
        if (!period || segmentLimit < 1 || workerLimit < 1)
            throw std::invalid_argument("Invalid solver dimensions");
        auto result = std::make_unique<Result>();
        const auto start = Clock::now();
        Control control{state,
                        start + std::chrono::duration_cast<Clock::duration>(
                                    std::chrono::duration<double>(budgetSeconds))};
        const int segments = static_cast<int>(std::min<uint64_t>(period, segmentLimit));
        result->segments = segments;
        std::vector<uint64_t> cuts(segments + 1);
        std::vector<Box> bounds(segments + 1, Box(-2, 2, -2, 2));
        for (int j = 0; j <= segments; ++j)
            cuts[j] = period * j / segments;
        // One certified forward sweep supplies the boundary domains. Within disjoint
        // windows around uniform cuts, select the narrowest enclosure observed.
        Box seed;
        const Box universe(-2, 2, -2, 2);
        const auto seedStart = Clock::now();
        auto seedDisk = enclosingDisk(seed);
        const auto parameterDisk = enclosingDisk(c);
        std::vector<bool> selected(segments + 1, false);
        for (uint64_t n = 1; n <= period; ++n) {
            if (state.interruptRequested())
                return nullptr;
            if (control.stopped()) {
                result->parameter = c;
                result->status = "SEED_TIME_BUDGET";
                result->seconds = std::chrono::duration<double>(Clock::now() - start).count();
                return result;
            }
            auto nextDisk = diskStep(seedDisk, parameterDisk);
            auto nextSeed = intersect(intersect(add(square(seed), c), diskBox(nextDisk)), universe);
            ++result->seedEvaluations;
            if (nextSeed.empty()) {
                result->parameter = nextSeed;
                result->status = "EMPTY_DOMAIN";
                result->seconds = std::chrono::duration<double>(Clock::now() - start).count();
                return result;
            }
            nextDisk = tightenDisk(std::move(nextDisk), nextSeed);
            if (!result->firstBroad && mpfr_cmp_ui(width(nextSeed).x, 1) >= 0)
                result->firstBroad = n;
            const auto j = static_cast<int>((n * segments + period / 2) / period);
            if (j > 0 && j < segments) {
                const auto a = n * segments, b = static_cast<uint64_t>(j) * period;
                if (4 * (a > b ? a - b : b - a) <= period &&
                    (!selected[j] || mpfr_cmp(width(nextSeed).x, width(bounds[j]).x) < 0)) {
                    selected[j] = true;
                    bounds[j] = nextSeed;
                    cuts[j] = n;
                }
            }
            if (same(seed, nextSeed) && sameDisk(seedDisk, nextDisk)) {
                result->seedSkipped = period - n;
                for (int k = 1; k < segments; ++k)
                    if (!selected[k]) {
                        bounds[k] = nextSeed;
                        selected[k] = true;
                    }
                seed = std::move(nextSeed);
                break;
            }
            seed = std::move(nextSeed);
            seedDisk = std::move(nextDisk);
        }
        result->seedSeconds = std::chrono::duration<double>(Clock::now() - seedStart).count();
        if (!contains(seed, 0, 0)) {
            result->parameter = Box(1, 0, 1, 0);
            result->status = "EMPTY_DOMAIN";
            result->seconds = std::chrono::duration<double>(Clock::now() - start).count();
            return result;
        }
        // Very short periods can have empty integer windows: retain a safe universe.
        bounds.front() = Box();
        bounds.back() = Box();
        result->cuts = cuts;
        for (const auto &b: bounds)
            result->seedWidths.push_back(log10value(width(b)));
        for (int j = 1; j < segments; ++j)
            if (cuts[j] == 1)
                bounds[j] = intersect(bounds[j], c);
        const auto initialWidth = width(c);
        int unchanged = 0;
        for (int round = 1; round <= maxRounds; ++round) {
            if (state.interruptRequested())
                return nullptr;
            if (control.stopped()) {
                result->status = "TIME_BUDGET";
                break;
            }
            const auto roundStart = Clock::now();
            std::vector<BlockResult> updates(segments);
            std::atomic<int> next{0};
            std::atomic<bool> interrupted{false};
            std::vector<std::jthread> workers;
            for (int worker = 0; worker < std::min(workerLimit, segments); ++worker)
                workers.emplace_back([&] {
                    for (int j; (j = next.fetch_add(1)) < segments;) {
                        updates[j] =
                            contract(bounds[j], bounds[j + 1], c, cuts[j + 1] - cuts[j], control);
                        if (updates[j].interrupted)
                            interrupted = true;
                    }
                    mpfr_free_cache();
                });
            for (auto &worker: workers)
                worker.join();
            if (state.interruptRequested())
                return nullptr;
            bool memoryLimited = false;
            for (const auto &u: updates)
                memoryLimited |= u.memoryLimited;
            if (memoryLimited) {
                result->status = "PREFIX_MEMORY_LIMIT";
                break;
            }
            bool branchLimited = false;
            for (const auto &u: updates)
                branchLimited |= u.branchLimited;
            if (branchLimited) {
                result->status = "BRANCH_LIMIT";
                break;
            }
            if (interrupted) {
                result->status = "TIME_BUDGET";
                break;
            }
            std::vector<Box> nextBounds = bounds;
            std::vector<Box> reduction;
            reduction.reserve(segments + 1);
            reduction.push_back(c);
            Counts counts;
            for (int j = 0; j < segments; ++j) {
                nextBounds[j] = intersect(nextBounds[j], updates[j].entrance);
                nextBounds[j + 1] = intersect(nextBounds[j + 1], updates[j].exit);
                reduction.push_back(updates[j].param);
                counts.forward += updates[j].counts.forward;
                counts.backward += updates[j].counts.backward;
                counts.skipped += updates[j].counts.skipped;
                counts.branches += updates[j].counts.branches;
                counts.hulls += updates[j].counts.hulls;
            }
            // Pairwise global intersection; local workers only read the previous sweep.
            while (reduction.size() > 1) {
                std::vector<Box> upper;
                for (size_t j = 0; j < reduction.size(); j += 2)
                    upper.push_back(j + 1 < reduction.size()
                                        ? intersect(reduction[j], reduction[j + 1])
                                        : reduction[j]);
                reduction = std::move(upper);
            }
            auto nextC = std::move(reduction.front());
            for (int j = 0; j <= segments; ++j)
                if (cuts[j] == 1) {
                    nextC = intersect(nextC, nextBounds[j]);
                    nextBounds[j] = intersect(nextBounds[j], nextC);
                }
            bool empty = nextC.empty();
            uint64_t changed = 0;
            for (int j = 0; j <= segments; ++j) {
                empty |= nextBounds[j].empty();
                changed += !same(nextBounds[j], bounds[j]);
            }
            const bool cChanged = !same(c, nextC);
            c = std::move(nextC);
            bounds = std::move(nextBounds);
            result->totals.forward += counts.forward;
            result->totals.backward += counts.backward;
            result->totals.skipped += counts.skipped;
            result->totals.branches += counts.branches;
            result->totals.hulls += counts.hulls;
            result->rounds.push_back(
                {round,
                 std::chrono::duration<double>(Clock::now() - roundStart).count(),
                 empty ? 0 : log10value(initialWidth) - log10value(width(c)),
                 changed,
                 counts});
            if (empty) {
                result->status = "EMPTY_DOMAIN";
                break;
            }
            if (mpfr_cmp(width(c).x, target.x) <= 0) {
                result->status = "NARROW_CANDIDATE_UNVALIDATED";
                break;
            }
            unchanged = (!changed && !cChanged) ? unchanged + 1 : 0;
            if (unchanged >= 2) {
                result->status = "STAGNATED";
                break;
            }
        }
        if (result->status.empty())
            result->status = "ROUND_LIMIT";
        result->parameter = std::move(c);
        result->seconds = std::chrono::duration<double>(Clock::now() - start).count();
        return result;
    }
} // namespace v4
