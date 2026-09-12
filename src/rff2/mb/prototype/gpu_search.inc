// Created by GPT-6 on 2026-09-11
// Accelerates the existing low-cost-derivative final search sweep only.
// Higher-precision derivative sweeps retain the CPU kernel.
// Modified by GPT-6 on 2026-09-12
inline bool gpuFastBlocks(
    const C &c, std::vector<Block> &blocks, int goalbits, bool special, unsigned threads) {
    const char *mode = std::getenv("STMS_GPU_SEARCH");
    const bool force = mode && std::string(mode) == "force";
    const unsigned bits = mpf_get_prec(c.r), F = bits + 64;
    if ((mode && std::string(mode) == "off") || (!force && bits < 262144) || F > 16000000 ||
        blocks.empty())
        return false;
    static std::unique_ptr<OrbitGPU> gpu;
    static bool unavailable = false;
    if (unavailable)
        return false;
    try {
        if (!gpu)
            gpu = std::make_unique<OrbitGPU>(
                (std::filesystem::path(__FILE__).parent_path().parent_path() / "orbit.spv")
                    .string()
                    .c_str());
        const unsigned n = std::max(512u, std::bit_ceil(2 * ((F + 47) / 16)));
        const unsigned tileSize =
            std::max(1u, std::min(8u, unsigned((768ull * 1024 * 1024 - 16ull * n) / (180ull * n))));
        Z cr, ci;
        ::F tmp(F + 64);
        mpf_mul_2exp(tmp, c.r, F);
        mpz_set_f(cr.x, tmp);
        mpf_mul_2exp(tmp, c.i, F);
        mpz_set_f(ci.x, tmp);
        auto run = [&](const std::vector<size_t> &indices, unsigned cap, bool commit) {
            const auto begin = std::chrono::steady_clock::now();
            const size_t count = indices.size();
            std::vector<Z> x(count), y(count);
            std::vector<EComplex> A(count, EComplex{1, 0, 0}), B(count);
            std::vector<double> logA(count);
            std::vector<uint64_t> current(count), end(count);
            for (size_t k = 0; k < count; ++k) {
                auto &q = blocks[indices[k]];
                mpf_mul_2exp(tmp, q.z.r, F);
                mpz_set_f(x[k].x, tmp);
                mpf_mul_2exp(tmp, q.z.i, F);
                mpz_set_f(y[k].x, tmp);
                current[k] = q.start;
                end[k] = cap ? std::min<uint64_t>(q.end, q.start + cap) : q.end;
                logA[k] = q.offset;
            }
            Z mr, mi;
            for (size_t first = 0; first < count; first += tileSize)
                for (;;) {
                    checkCancellation();
                    std::vector<size_t> active;
                    unsigned steps = 32;
                    for (size_t k = first; k < std::min(count, first + tileSize); ++k)
                        if (current[k] < end[k]) {
                            active.push_back(k);
                            steps = std::min<uint64_t>(steps, end[k] - current[k]);
                        }
                    if (active.empty())
                        break;
                    std::vector<OrbitGPU::Segment> batch;
                    for (auto k: active)
                        batch.push_back({x[k].x, y[k].x, cr.x, ci.x, x[k].x, y[k].x});
                    gpu->orbit(F, steps, batch, [] {
                        return cancelled();
                    });
                    for (unsigned step = 0; step < steps; ++step)
                        for (size_t k = 0; k < active.size(); ++k) {
                            auto j = active[k];
                            double values[2];
                            long exponents[2];
                            for (unsigned comp = 0; comp < 2; ++comp) {
                                auto *p = gpu->trace.data() +
                                          ((step * active.size() + k) * 2 + comp) * 10;
                                uint16_t digits[8];
                                for (unsigned d = 0; d < 8; ++d)
                                    digits[d] = p[d];
                                auto m = comp ? mi.x : mr.x;
                                mpz_import(m, 8, -1, 2, 0, 0, digits);
                                if (p[9])
                                    mpz_neg(m, m);
                                values[comp] = mpz_get_d_2exp(&exponents[comp], m);
                                exponents[comp] = values[comp]
                                                      ? exponents[comp] + int64_t(p[8]) * 16 - F
                                                      : -1000000000;
                            }
                            const auto e = std::max(exponents[0], exponents[1]);
                            EComplex z =
                                e < -999999999
                                    ? EComplex{}
                                    : normalize(
                                          {std::ldexp(
                                               values[0],
                                               int(std::clamp(exponents[0] - e, -2000L, 2000L))),
                                           std::ldexp(
                                               values[1],
                                               int(std::clamp(exponents[1] - e, -2000L, 2000L))),
                                           e});
                            if (current[j]) {
                                logA[j] += 1 + elg(z);
                                A[j] = emul2(z, A[j]);
                            }
                            B[j] = tangent(z, B[j]);
                            ++current[j];
                        }
                }
            if (commit)
                for (size_t k = 0; k < count; ++k) {
                    auto &q = blocks[indices[k]];
                    mpf_set_z(q.out.r, x[k].x);
                    mpf_div_2exp(q.out.r, q.out.r, F);
                    mpf_set_z(q.out.i, y[k].x);
                    mpf_div_2exp(q.out.i, q.out.i, F);
                    to_f(q.fa, A[k]);
                    to_f(q.fb, B[k]);
                    q.logA = logA[k];
                }
            return std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        };
        if (!force) {
            std::vector<size_t> sample;
            for (size_t k = 0; k < std::min<size_t>(15, blocks.size()); ++k)
                sample.push_back(k * blocks.size() / std::min<size_t>(15, blocks.size()));
            const auto begin = std::chrono::steady_clock::now();
#pragma omp parallel for num_threads(threads) schedule(static)
            for (size_t k = 0; k < sample.size(); ++k) {
                auto &q = blocks[sample[k]];
                C out(bits), a(bits), b(bits);
                fastblock(c,
                          q.z,
                          q.start,
                          std::min<unsigned long>(32, q.end - q.start),
                          goalbits,
                          q.offset,
                          special,
                          out,
                          a,
                          b);
            }
            const double cpu =
                             std::chrono::duration<double>(std::chrono::steady_clock::now() - begin)
                                 .count(),
                         device = run(sample, 32, false);
            const bool use = device < cpu * .8;
            std::clog << "GPU_SEARCH_SELECT bits=" << F << " cpu_s=" << cpu << " gpu_s=" << device
                      << " choice=" << (use ? "GPU" : "CPU") << std::endl;
            if (!use)
                return false;
        }
        std::vector<size_t> all(blocks.size());
        for (size_t j = 0; j < all.size(); ++j)
            all[j] = j;
        const auto seconds = run(all, 0, true);
        std::clog << "GPU_SEARCH blocks=" << all.size() << " bits=" << F << " seconds=" << seconds
                  << std::endl;
        return true;
    } catch (const std::exception &e) {
        std::clog << "GPU_SEARCH_FALLBACK reason=" << e.what() << std::endl;
        if (!cancelled()) {
            gpu.reset();
            unavailable = true;
        }
        return false;
    }
}
