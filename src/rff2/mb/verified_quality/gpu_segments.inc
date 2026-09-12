// Created by GPT-6 on 2026-09-11
// Modified by GPT-6 on 2026-09-12
// Included inside namespace segments after Part. All reduced-precision derivative
// operations carry outward error bounds; endpoints retain full requested precision.
inline thread_local bool gpuRetryDisabled = false;
struct SmallFloat {
    mpfr_t x;
    SmallFloat() {
        mpfr_init2(x, 192);
        mpfr_set_zero(x, 1);
    }
    ~SmallFloat() {
        mpfr_clear(x);
    }
};
struct SmallDisk {
    SmallFloat re, im;
    upper::Positive radius;
};
struct SmallOps {
    SmallFloat t, u, re, im, n;
    mpfr_exp_t largest = 0;
    void reset() {
        largest = std::numeric_limits<mpfr_exp_t>::min();
    }
    void track(mpfr_srcptr v) {
        if (!mpfr_number_p(v))
            throw std::runtime_error("GPU derivative range");
        if (!mpfr_zero_p(v))
            largest = std::max(largest, mpfr_get_exp(v));
    }
    upper::Positive rounding() {
        return largest == std::numeric_limits<mpfr_exp_t>::min()
                   ? upper::Positive{}
                   : upper::Positive::power(largest - 192 + 4);
    }
    upper::Positive norm(const SmallDisk &a) {
        mpfr_hypot(n.x, a.re.x, a.im.x, MPFR_RNDU);
        return upper::Positive::abs(n.x);
    }
    void mul(SmallDisk &out, const SmallDisk &a, const SmallDisk &b) {
        auto error =
            upper::add(upper::add(upper::mul(norm(a), b.radius), upper::mul(norm(b), a.radius)),
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
    void addOne(SmallDisk &out) {
        reset();
        mpfr_add_ui(out.re.x, out.re.x, 1, MPFR_RNDN);
        track(out.re.x);
        out.radius = upper::add(out.radius, rounding());
    }
};
inline bool gpuParts(const C &c,
                     uint64_t period,
                     const v4::Control &control,
                     const std::vector<upper::Positive> &radii,
                     std::vector<Part> &parts) {
    if (gpuRetryDisabled || control.stopped())
        return false;
    const char *mode = std::getenv("STMS_GPU_VERIFY");
    bool force = mode && std::string(mode) == "force";
    if ((mode && std::string(mode) == "off") || v4::precision < 256 ||
        (!force && v4::precision < 262144) || v4::precision > 15999000)
        return false;
    const auto exponentLimit = stms_limits::exponentLimit(v4::precision);
    if (mpfr_get_emin() >= -4 * exponentLimit || mpfr_get_emax() <= 4 * exponentLimit)
        return false;
    auto valid = [&](mpfr_srcptr x) {
        return mpfr_number_p(x) && (mpfr_zero_p(x) || std::abs(mpfr_get_exp(x)) <= exponentLimit);
    };
    if (!valid(c.re.x) || !valid(c.im.x))
        return false;
    for (const auto &seed: quality::checkpoints)
        if (!valid(seed.value.re.x) || !valid(seed.value.im.x))
            return false;
    static std::unique_ptr<OrbitGPU> gpu;
    static bool unavailable = false;
    if (unavailable)
        return false;
    const auto begin = v4::Clock::now();
    try {
        if (!gpu)
            gpu = std::make_unique<OrbitGPU>(
                (std::filesystem::path(__FILE__).parent_path().parent_path() / "orbit.spv")
                    .string()
                    .c_str());
        const unsigned F = v4::precision + 128,
                       n = std::max(512u, std::bit_ceil(2 * ((F + 47) / 16)));
        const unsigned tileSize =
            std::max(1u, std::min(8u, unsigned((768ull * 1024 * 1024 - 16ull * n) / (180ull * n))));
        const auto &seeds = quality::checkpoints;
        std::vector<Z> x(seeds.size()), y(seeds.size());
        Z cr, ci;
        auto convert = [&](mpz_ptr out, mpfr_srcptr in) {
            auto e = mpfr_get_z_2exp(out, in);
            auto shift = int64_t(e) + F;
            if (shift >= 0)
                mpz_mul_2exp(out, out, shift);
            else
                mpz_tdiv_q_2exp(out, out, -shift);
        };
        convert(cr.x, c.re.x);
        convert(ci.x, c.im.x);
        const auto unit = upper::Positive::power(-int64_t(F));
        std::vector<uint64_t> current(seeds.size()), end(seeds.size());
        auto A = std::make_unique<SmallDisk[]>(seeds.size()),
             D = std::make_unique<SmallDisk[]>(seeds.size());
        for (size_t j = 0; j < seeds.size(); ++j) {
            convert(x[j].x, seeds[j].value.re.x);
            convert(y[j].x, seeds[j].value.im.x);
            current[j] = seeds[j].index;
            end[j] = j + 1 < seeds.size() ? seeds[j + 1].index : period;
            parts[j].radius = upper::add(radii[j], upper::twice(unit));
            parts[j].critical = upper::Positive::power(0);
            mpfr_set_ui(A[j].re.x, 1, MPFR_RNDN);
        }
        SmallOps ops;
        SmallDisk z, coefficient;
        Z mantissa;
        double gpuSeconds = 0;
        uint64_t iterations = 0;
        for (size_t first = 0; first < seeds.size(); first += tileSize) {
            for (;;) {
                if (control.stopped())
                    throw std::runtime_error("GPU verifier cancelled");
                std::vector<size_t> indices;
                uint64_t remaining = UINT64_MAX;
                for (size_t j = first; j < std::min(seeds.size(), first + tileSize); ++j)
                    if (current[j] < end[j]) {
                        indices.push_back(j);
                        remaining = std::min(remaining, end[j] - current[j]);
                    }
                const unsigned steps = std::min<uint64_t>(32, remaining);
                if (indices.empty())
                    break;
                std::vector<OrbitGPU::Segment> batch;
                for (auto j: indices)
                    batch.push_back({x[j].x, y[j].x, cr.x, ci.x, x[j].x, y[j].x});
                NTTGPU::Timing timing;
                if (gpu->flags & 8) {
                    const bool last = remaining == steps;
                    if (!gpu->pending)
                        gpu->submitOrbit(F, steps, batch, !gpu->resident, last, [&] {
                            return control.stopped();
                        });
                    timing = gpu->finishOrbit([&] {
                        return control.stopped();
                    });
                    // trace was copied out; next device batch overlaps the CPU disk updates below.
                    if (!last) {
                        const unsigned next = std::min<uint64_t>(32, remaining - steps);
                        gpu->submitOrbit(F, next, batch, false, remaining == steps + next, [&] {
                            return control.stopped();
                        });
                    }
                } else
                    timing = gpu->orbit(F, steps, batch, [&] {
                        return control.stopped();
                    });
                gpuSeconds += timing.gpu;
                for (unsigned step = 0; step < steps; ++step)
                    for (size_t k = 0; k < indices.size(); ++k) {
                        const auto j = indices[k];
                        auto &part = parts[j];
                        upper::Positive traceError;
                        for (unsigned comp = 0; comp < 2; ++comp) {
                            auto *p =
                                gpu->trace.data() + ((step * indices.size() + k) * 2 + comp) * 10;
                            uint16_t digits[8];
                            for (unsigned q = 0; q < 8; ++q)
                                digits[q] = p[q];
                            mpz_import(mantissa.x, 8, -1, 2, 0, 0, digits);
                            if (p[9])
                                mpz_neg(mantissa.x, mantissa.x);
                            auto target = comp ? z.im.x : z.re.x;
                            mpfr_set_z(target, mantissa.x, MPFR_RNDN);
                            mpfr_mul_2si(target, target, int64_t(p[8]) * 16 - F, MPFR_RNDN);
                            if (p[8])
                                traceError = upper::add(
                                    traceError, upper::Positive::power(int64_t(p[8]) * 16 - F));
                        }
                        const auto norm = upper::add(ops.norm(z), traceError);
                        mpfr_mul_2ui(coefficient.re.x, z.re.x, 1, MPFR_RNDN);
                        mpfr_mul_2ui(coefficient.im.x, z.im.x, 1, MPFR_RNDN);
                        coefficient.radius = upper::twice(upper::add(part.radius, traceError));
                        ops.mul(D[j], coefficient, D[j]);
                        ops.addOne(D[j]);
                        if (current[j]) {
                            ops.mul(A[j], coefficient, A[j]);
                            part.critical = upper::mul(part.critical,
                                                       upper::twice(upper::add(norm, part.radius)));
                        }
                        // Parameter conversion and component truncation contribute <4*2^-F.
                        part.radius = upper::add(
                            upper::mul(part.radius, upper::add(upper::twice(norm), part.radius)),
                            upper::twice(upper::twice(unit)));
                        if (part.radius.m &&
                            (part.radius.e > -1 || (part.radius.e == -1 && part.radius.m > .5)))
                            throw std::runtime_error("GPU enclosure too wide");
                        ++current[j];
                        ++part.iterations;
                        ++iterations;
                    }
            }
        }
        for (size_t j = 0; j < seeds.size(); ++j) {
            auto &part = parts[j];
            for (unsigned comp = 0; comp < 2; ++comp) {
                auto target = comp ? part.endpoint.im.x : part.endpoint.re.x;
                int rounded = mpfr_set_z(target, comp ? y[j].x : x[j].x, MPFR_RNDN);
                mpfr_mul_2si(target, target, -int64_t(F), MPFR_RNDN);
                if (rounded)
                    part.radius = upper::add(
                        part.radius, upper::Positive::power(mpfr_get_exp(target) - v4::precision));
            }
            for (auto pair: {std::pair{&part.a, &A[j]}, std::pair{&part.b, &D[j]}}) {
                mpfr_set(pair.first->value.re.x, pair.second->re.x, MPFR_RNDN);
                mpfr_set(pair.first->value.im.x, pair.second->im.x, MPFR_RNDN);
                pair.first->radius = pair.second->radius;
            }
            part.status = "COMPLETE";
        }
        std::clog << "GPU_VERIFY iterations=" << iterations << " bits=" << F << " tile=" << tileSize
                  << " gpu_s=" << gpuSeconds
                  << " total_s=" << std::chrono::duration<double>(v4::Clock::now() - begin).count()
                  << std::endl;
        std::ofstream evidence(std::filesystem::path(__FILE__).parent_path().parent_path() /
                                   "gpu_execution.log",
                               std::ios::app);
        evidence << "GPU_VERIFY period=" << period << " iterations=" << iterations << " bits=" << F
                 << " tile=" << tileSize << " gpu_s=" << gpuSeconds
                 << " total_s=" << std::chrono::duration<double>(v4::Clock::now() - begin).count()
                 << '\n';
        return true;
    } catch (const std::exception &e) {
        // Drain submitted work even when CPU-side cancellation interrupts the pipeline.
        if (gpu && gpu->pending) {
            try {
                gpu->finishOrbit([] {
                    return true;
                });
            } catch (...) {
            }
        }
        if (gpu)
            gpu->resident = false;
        std::clog << "GPU_VERIFY_FALLBACK reason=" << e.what() << std::endl;
        for (auto &part: parts)
            part = Part{};
        if (!control.stopped()) {
            gpu.reset();
            unavailable = true;
        }
        return false;
    }
}
