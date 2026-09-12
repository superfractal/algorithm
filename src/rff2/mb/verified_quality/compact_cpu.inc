// Created by GPT-6 on 2026-09-12
// Included within namespace segments. Full-P orbit enclosures are unchanged.
// Derivative centers use 192 bits, with outward conversion/operation errors.
// Modified by GPT-6 on 2026-09-12
// Local RFF_Astra implementation; project license: GNU GPL version 3.
// See ../../../../LICENSE and ../../../../SOURCES_AND_REFERENCES.md.
// Source lineage: local v201 segments.hpp, disks.hpp, upper.hpp, and the
// bounded compact derivative design in gpu_segments.inc; CPU path added v202.
// This application code calls MPFR/GMP APIs; it is not their library source.
inline thread_local bool compactRetryDisabled = false;
namespace compact {
    struct R {
        mpfr_t x;
        R() {
            mpfr_init2(x, 192);
            mpfr_set_zero(x, 1);
        }
        R(const R &) = delete;
        R &operator=(const R &) = delete;
        ~R() {
            mpfr_clear(x);
        }
    };
    struct D {
        R re, im;
        upper::Positive radius;
    };
    inline void copy(D &out, const D &in) {
        mpfr_set(out.re.x, in.re.x, MPFR_RNDN);
        mpfr_set(out.im.x, in.im.x, MPFR_RNDN);
        out.radius = in.radius;
    }
    inline upper::Positive norm(const D &z) {
        auto a = upper::Positive::abs(z.re.x), b = upper::Positive::abs(z.im.x);
        if (!a.m)
            return b;
        if (!b.m)
            return a;
        if (a.e < b.e)
            std::swap(a, b);
        const auto distance = a.e - b.e;
        const double small = distance > 1074
                                 ? std::numeric_limits<double>::denorm_min()
                                 : upper::up(std::scalbn(b.m, -static_cast<int>(distance)));
        const double aa = upper::up(a.m * a.m), bb = upper::up(small * small);
        return upper::Positive(upper::up(std::sqrt(upper::up(aa + bb))), a.e);
    }
    struct Ops {
        R t, u, re, im;
        mpfr_exp_t largest;
        void reset() {
            largest = std::numeric_limits<mpfr_exp_t>::min();
        }
        void track(mpfr_srcptr x) {
            if (!mpfr_number_p(x) ||
                (!mpfr_zero_p(x) &&
                 std::abs(mpfr_get_exp(x)) > stms_limits::exponentLimit(v4::precision)))
                throw std::runtime_error("COMPACT_DERIVATIVE_RANGE");
            if (!mpfr_zero_p(x))
                largest = std::max(largest, mpfr_get_exp(x));
        }
        upper::Positive rounding() const {
            return largest == std::numeric_limits<mpfr_exp_t>::min()
                       ? upper::Positive{}
                       : upper::Positive::power(largest - 192 + 4);
        }
        void coefficient(D &out, const C &value, upper::Positive radius) {
            reset();
            mpfr_set(out.re.x, value.re.x, MPFR_RNDN);
            mpfr_set(out.im.x, value.im.x, MPFR_RNDN);
            mpfr_mul_2ui(out.re.x, out.re.x, 1, MPFR_RNDN);
            track(out.re.x);
            mpfr_mul_2ui(out.im.x, out.im.x, 1, MPFR_RNDN);
            track(out.im.x);
            out.radius = upper::add(upper::twice(radius), rounding());
        }
        void mul(D &out, const D &a, const D &b, upper::Positive na, upper::Positive nb) {
            auto error = upper::add(upper::add(upper::mul(na, b.radius), upper::mul(nb, a.radius)),
                                    upper::mul(a.radius, b.radius));
            reset();
            mpfr_mul(t.x, a.re.x, b.re.x, MPFR_RNDN);
            track(t.x);
            mpfr_mul(u.x, a.im.x, b.im.x, MPFR_RNDN);
            track(u.x);
            mpfr_sub(re.x, t.x, u.x, MPFR_RNDN);
            track(re.x);
            mpfr_mul(t.x, a.re.x, b.im.x, MPFR_RNDN);
            track(t.x);
            mpfr_mul(u.x, a.im.x, b.re.x, MPFR_RNDN);
            track(u.x);
            mpfr_add(im.x, t.x, u.x, MPFR_RNDN);
            track(im.x);
            mpfr_set(out.re.x, re.x, MPFR_RNDN);
            mpfr_set(out.im.x, im.x, MPFR_RNDN);
            out.radius = upper::add(error, rounding());
        }
        void mul(D &out, const D &a, const D &b) {
            mul(out, a, b, norm(a), norm(b));
        }
        void add(D &out, const D &a, const D &b) {
            auto error = upper::add(a.radius, b.radius);
            reset();
            mpfr_add(re.x, a.re.x, b.re.x, MPFR_RNDN);
            track(re.x);
            mpfr_add(im.x, a.im.x, b.im.x, MPFR_RNDN);
            track(im.x);
            mpfr_set(out.re.x, re.x, MPFR_RNDN);
            mpfr_set(out.im.x, im.x, MPFR_RNDN);
            out.radius = upper::add(error, rounding());
        }
        void addOne(D &out) {
            reset();
            mpfr_add_ui(out.re.x, out.re.x, 1, MPFR_RNDN);
            track(out.re.x);
            out.radius = upper::add(out.radius, rounding());
        }
    };
    inline void exportDisk(disks::Disk &out, const D &in) {
        if (mpfr_get_prec(out.value.re.x) < 192 || mpfr_get_prec(out.value.im.x) < 192)
            throw std::runtime_error("COMPACT_EXPORT_PRECISION");
        mpfr_set(out.value.re.x, in.re.x, MPFR_RNDN);
        mpfr_set(out.value.im.x, in.im.x, MPFR_RNDN);
        out.radius = in.radius;
    }
} // namespace compact
inline bool compactCpuParts(const C &c,
                            uint64_t period,
                            const v4::Control &control,
                            const std::vector<upper::Positive> &radii,
                            std::vector<Part> &parts,
                            const std::vector<size_t> &order) {
    const char *mode = std::getenv("STMS_COMPACT_CPU");
    if (compactRetryDisabled || v4::precision <= 192 || (mode && std::string(mode) == "off"))
        return false;
    const auto &seeds = quality::checkpoints;
    std::atomic<size_t> next = 0;
    std::vector<std::jthread> workers;
    for (size_t worker = 0; worker < std::min<size_t>(15, seeds.size()); ++worker)
        workers.emplace_back([&] {
            compact::Ops ops;
            auto group = std::make_unique<compact::D[]>(64);
            std::vector<upper::Positive> norms(64);
            compact::D initial, suffix, constant, groupA;
            for (;;) {
                const auto slot = next.fetch_add(1);
                if (slot >= order.size())
                    break;
                const auto j = order[slot];
                auto &part = parts[j];
                try {
                    upper::environment();
                    part = Part{};
                    enclosure::Orbit orbit;
                    orbit.value = seeds[j].value;
                    orbit.radiusBound = radii[j];
                    orbit.iterations = seeds[j].index;
                    compact::D A, B;
                    mpfr_set_ui(A.re.x, 1, MPFR_RNDN);
                    size_t used = 0;
                    // Compose the same affine recurrence as the full-P route:
                    // groupA is the product, constant the suffix-product sum.
                    // Each rounded compact operation contributes an outward radius.
                    auto flush = [&] {
                        if (!used)
                            return;
                        if (used == 1) {
                            compact::copy(groupA, group[0]);
                            mpfr_set_ui(constant.re.x, 1, MPFR_RNDN);
                            mpfr_set_zero(constant.im.x, 1);
                            constant.radius = {};
                        } else {
                            compact::copy(suffix, group[used - 1]);
                            compact::copy(constant, suffix);
                            ops.addOne(constant);
                            for (size_t k = used - 2; k > 0; --k) {
                                ops.mul(suffix, suffix, group[k], compact::norm(suffix), norms[k]);
                                ops.add(constant, constant, suffix);
                            }
                            ops.mul(groupA, suffix, group[0], compact::norm(suffix), norms[0]);
                        }
                        const auto norm = compact::norm(groupA);
                        ops.mul(A, groupA, A, norm, compact::norm(A));
                        ops.mul(B, groupA, B, norm, compact::norm(B));
                        ops.add(B, B, constant);
                        used = 0;
                    };
                    const uint64_t end = j + 1 < seeds.size() ? seeds[j + 1].index : period;
                    while (orbit.iterations < end) {
                        if ((part.iterations & 255) == 0 && control.stopped()) {
                            orbit.status = "TIME_BUDGET";
                            break;
                        }
                        const bool grouped = orbit.iterations != 0;
                        auto &coefficient = grouped ? group[used] : initial;
                        ops.coefficient(coefficient, orbit.value, orbit.radiusBound);
                        if (grouped) {
                            norms[used] = compact::norm(coefficient);
                            ++used;
                        } else {
                            ops.mul(B, coefficient, B);
                            ops.addOne(B);
                        }
                        if (!orbit.step(c))
                            break;
                        ++part.iterations;
                        if (grouped && used == 64)
                            flush();
                    }
                    if (orbit.iterations == end && orbit.status == "COMPLETE")
                        flush();
                    compact::exportDisk(part.a, A);
                    compact::exportDisk(part.b, B);
                    part.endpoint = std::move(orbit.value);
                    part.radius = orbit.radiusBound;
                    part.critical = orbit.criticalBound;
                    part.status = orbit.status;
                    if (orbit.iterations != end && part.status == "COMPLETE")
                        part.status = "INCOMPLETE";
                } catch (const std::exception &e) {
                    part.status = e.what();
                }
            }
            mpfr_free_cache();
        });
    for (auto &w: workers)
        w.join();
    return true;
}
