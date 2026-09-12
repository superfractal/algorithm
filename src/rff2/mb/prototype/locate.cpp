// Modified by GPT-6 on 2026-09-10, 2026-09-11, 2026-09-12
/*
 * STMS-RJ: sensitivity-tapered multiple shooting with cached return jets.
 * Experimental nucleus locator, produced 2026-09-10.
 * This implementation does NOT reproduce RFF's output zoom convention.
 * No expected result, filename, known period, or saved orbit is used to solve.
 * See ALGORITHM.md and REPORT_JA.md before integrating into a renderer.
 */
// Integration note: the preserved block above describes the original archive.
// RFF supplies its automatic FPG period; the bridge verifies the final scale.
// Current behavior: ../../../../About_the_Algorithm.md. Attachment provenance
// and unconfirmed terms: ../../../../SOURCES_AND_REFERENCES.md and LICENSE
// in that same directory. The archive's report paths are historical paths.
#include <cstring>
#include <limits>
#include <sstream>
#include "../STMSLimits.hpp"
#include "numerics.inc"
#include "partition.inc"
#include "return_jets.inc"
#ifdef STMS_GPU_VERIFY_ENABLED
#include "../orbit_gpu.hpp"
#include <filesystem>
#include <bit>
#include "gpu_search.inc"
#endif

struct Input {
    float zoom;
    uint64_t field;
    std::string real, imag;
};
uint64_t read_le(std::istream &in, int bytes) {
    uint64_t value = 0;
    for (int j = 0; j < bytes; j++) {
        int ch = in.get();
        if (ch == EOF)
            throw std::runtime_error("truncated RFL");
        value |= uint64_t(static_cast<unsigned char>(ch)) << (8 * j);
    }
    return value;
}
std::string read_string(std::istream &in) {
    const auto n = read_le(in, 8);
    if (n == 0 || n > uint64_t(stms_limits::maxVerificationBits))
        throw std::runtime_error("invalid RFL coordinate length");
    std::string s(n, '\0');
    if (!in.read(s.data(), n))
        throw std::runtime_error("truncated RFL coordinate");
    for (char x: s)
        if (!((x >= '0' && x <= '9') || x == '-' || x == '+' || x == '.'))
            throw std::runtime_error("only plain decimal RFL coordinates are supported");
    return s;
}
Input read_input(const char *path) {
    std::ifstream in(path, std::ios::binary);
    if (!in)
        throw std::runtime_error("cannot open input RFL");
    Input o;
    uint32_t word = static_cast<uint32_t>(read_le(in, 4));
    static_assert(sizeof(float) == 4, "requires binary32 float");
    std::memcpy(&o.zoom, &word, 4);
    if (!std::isfinite(o.zoom) || o.zoom < 0 || o.zoom > stms_limits::maxInputLogZoom)
        throw std::runtime_error("input logZoom outside STMS precision range");
    o.field = read_le(in, 8);
    o.real = read_string(in);
    o.imag = read_string(in);
    if (in.peek() != EOF)
        throw std::runtime_error("unexpected trailing RFL bytes");
    return o;
}
struct StageDiagnostic {
    int iteration;
    double goal, target, logStep, maxState;
    int stateBits, derivativeBits, mode;
    double kernelSeconds, compositionSeconds;
};
struct Result {
    std::vector<StageDiagnostic> stages;
    double pilotSeconds = 0, kernelSeconds = 0, compositionSeconds = 0, jetBuildSeconds = 0,
           jetReplaySeconds = 0, dispatchFallbackSeconds = 0;
    size_t pilotPoints = 0;
    std::vector<Block> checkpoints;
    unsigned long period = 0;
    int iterations = 0;
    double logAB = 0, phase = 0, logCorrection = 0;
    std::string method;
};

// One serial jet build. Later Newton steps reuse those maps with an error test
// on both the starting-state expansion and linearized parameter dependence.
Result solve_jets(C &c, double zoom, bool trace, unsigned long knownPeriod = 0) {
    const unsigned bits = mpf_get_prec(c.r);
    C base(bits), z(bits), A(bits), D(bits), step(bits);
    C seedBase(bits), seedDelta(bits), seedAdjustment(bits);
    std::vector<Block> buildSeeds;
    set(base, c);
    Ops op(bits);
    unsigned long p = 0;
    double mapBuildGoal = zoom + 55;
    double goal = zoom + 55, target = 2 * zoom + 18, logstep = -zoom, curv = zoom + 2;
    Result result;
    result.method = "cached-return-jets";
    for (int it = 0; it < 15; ++it) {
        checkCancellation();
        if (it)
            goal = std::min(target + 8, std::max(zoom + 20, -4 * logstep - 3 * curv + 16));
        std::vector<Block> captured;
        const auto orbitBegin = std::chrono::steady_clock::now();
        C displacement(bits);
        sub(displacement, c, base);
        const bool rebase = it && (goal > -2 * norm10(displacement) - 16 ||
                                   (goal >= target && goal > mapBuildGoal + 32));
        if (rebase) {
            set(base, c);
            persistent_maps.clear();
            mapBuildGoal = goal;
        }
        if (!it || rebase)
            set(seedBase, c);
        auto o = (it && !rebase)
                     ? cached_orbit(base, c, goal, p, z, A, D, &captured)
                     : jet_orbit(c, goal, zoom, it ? p : knownPeriod, z, A, D, &buildSeeds);
        const auto orbitSeconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - orbitBegin).count();
        if (it)
            result.jetReplaySeconds += orbitSeconds;
        else
            result.jetBuildSeconds += orbitSeconds;
        checkCancellation();
        p = o.p;
        if (it == 0)
            target = std::max(zoom + 25, o.logAB + 18);
        op.div(step, z, D);
        sub(c, c, step);
        double prev = logstep;
        logstep = norm10(step);
        if (trace)
            std::cerr << "jet it=" << it << " p=" << p << " goal=" << goal
                      << " logCorrection=" << logstep << " raw=" << o.raw << " jumps=" << o.jumps
                      << " seconds=" << o.seconds << "\n";
        result.period = p;
        result.iterations = it + 1;
        result.logAB = o.logAB;
        result.phase = o.phase;
        result.logCorrection = logstep;
        result.stages.push_back({it, goal, target, logstep, 0, 0, 0, 2, orbitSeconds, 0});
        if (logstep < -target && goal >= target && !captured.empty()) {
            if (captured.size() < 2 && buildSeeds.size() >= 2) {
                sub(seedDelta, c, seedBase);
                for (auto &seed: buildSeeds) {
                    op.mul(seedAdjustment, seed.fb, seedDelta);
                    add(seed.z, seed.z, seedAdjustment);
                }
                captured = std::move(buildSeeds);
            }
            for (size_t j = 0; j < captured.size(); ++j)
                captured[j].end = j + 1 < captured.size() ? captured[j + 1].start : p;
            result.checkpoints = std::move(captured);
            persistent_maps.clear();
            return result;
        }
        if (it)
            curv = std::max(logstep - 2 * prev + 2, zoom - 30);
    }
    throw std::runtime_error("cached-return Newton did not converge");
}

Result solve(C &c,
             double zoom,
             int threads,
             int requested_blocks,
             bool trace,
             unsigned long knownPeriod = 0) {
    omp_set_num_threads(threads);
    const unsigned bits = mpf_get_prec(c.r);
    C s(bits), z(bits), A(bits), D(bits), step(bits), t(bits), temp(bits), res(bits);
    Ops ops(bits);
    mpf_add_ui(t.r, c.r, 2);
    bool special = lg(t.r) < -100 && lg(c.i) < -100;
    std::vector<Point> points;
    points.emplace_back(bits, 0, 0);
    double goal = zoom + 20;
    // Cheap, data-dependent dispatch. No file-number-specific selection.
    unsigned long cap =
        zoom < 500 ? (knownPeriod && knownPeriod <= 8000000 ? knownPeriod : 100000) : 100000000;
    const auto pilotBegin = std::chrono::steady_clock::now();
    auto pilot =
        fastpilot(c, int(goal * std::log2(10)) + 96, special, zoom, z, D, points, cap, knownPeriod);
    auto pilotSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - pilotBegin).count();
    double dispatchFallbackSeconds = 0;
    checkCancellation();
    if (pilot.period == 0) {
        if (zoom >= 500)
            throw std::runtime_error("pilot exceeded 100 million iterations");
        const auto pilotPoints = points.size();
        points.clear();
        C initial = c;
        const auto jetBegin = std::chrono::steady_clock::now();
        try {
            struct BudgetGuard {
                double previous = adaptiveJetBudgetSeconds;
                BudgetGuard(bool bounded) {
                    adaptiveJetBudgetSeconds = bounded ? 1.5 : 0;
                }
                ~BudgetGuard() {
                    adaptiveJetBudgetSeconds = previous;
                }
            } guard(knownPeriod != 0);
            auto result = solve_jets(c, zoom, trace, knownPeriod);
            result.pilotSeconds = pilotSeconds;
            result.pilotPoints = pilotPoints;
            return result;
        } catch (const JetBudgetExceeded &) {
            checkCancellation();
            dispatchFallbackSeconds =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - jetBegin).count();
            set(c, initial);
            persistent_maps.clear();
            points.clear();
            points.emplace_back(bits, 0, 0);
            const auto retryBegin = std::chrono::steady_clock::now();
            pilot = fastpilot(c,
                              int(goal * std::log2(10)) + 96,
                              special,
                              zoom,
                              z,
                              D,
                              points,
                              knownPeriod,
                              knownPeriod);
            pilotSeconds +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - retryBegin)
                    .count();
            checkCancellation();
            if (!pilot.period)
                throw std::runtime_error("direct pilot after jet budget failed");
            if (trace)
                std::clog << "JET_BUDGET_FALLBACK seconds=" << dispatchFallbackSeconds
                          << " period=" << knownPeriod << std::endl;
        }
    }
    unsigned long p = pilot.period;
    double logAB = (pilot.logA + la(D)) / std::log2(10.), target = std::max(zoom + 25, logAB + 18);
    ops.div(step, z, D);
    mpf_neg(step.r, step.r);
    mpf_neg(step.i, step.i);
    double logstep = la(step) / std::log2(10.);
    add(c, c, step);
    const int nblocks =
        requested_blocks ? requested_blocks : (p >= 100000 ? threads * 12 : threads * 8);
    std::vector<Block> blocks;
    size_t m = std::min<size_t>(nblocks, points.size());
    const auto selected =
        sensitivity_partition::boundaries(points, m, p, bits, (target + 12) * std::log2(10.) + 112);
    for (size_t j = 0; j < m; j++) {
        auto &pt = points[selected[j]];
        blocks.emplace_back(bits);
        auto &q = blocks.back();
        q.start = pt.n;
        q.offset = pt.logA;
        ops.mul(temp, pt.D, step);
        add(q.z, pt.z, temp);
    }
    for (size_t j = 0; j < m; j++)
        blocks[j].end = j + 1 < m ? blocks[j + 1].start : p;
    if (trace)
        std::cerr << "pilot p=" << p << " blocks=" << m << " logCorrection=" << logstep << "\n";
    // Compose the connection equations as delta_s[j] = tv[j] + uv[j]*delta_c.
    // Interior states and the parameter must be corrected together.
    std::vector<C> tv, uv;
    for (size_t j = 0; j <= m; j++) {
        tv.emplace_back(bits);
        uv.emplace_back(bits);
    }
    double prevstep = logstep, curvature = zoom + 2;
    // This is the final proposal sweep. Bounded acceptance runs in STMSBridge.
    bool finalProposalSweep = false;
    Result result;
    result.method = "tapered-multiple-shooting";
    result.period = p;
    result.dispatchFallbackSeconds = dispatchFallbackSeconds;
    result.pilotSeconds = pilotSeconds;
    result.pilotPoints = points.size();
    for (int it = 1; it < 20; ++it) {
        checkCancellation();
        logstep = std::max(logstep, -target - 100);
        goal = finalProposalSweep
                   ? target + 12
                   : std::min(target + 12, std::max(zoom + 20, -4 * logstep - 3 * curvature + 16));
        int gbits = int(goal * std::log2(10.)) + 112;
        int abits =
            int((finalProposalSweep ? 48 : std::max(32., -2 * logstep - 2 * curvature + 20)) *
                std::log2(10.));
        const auto kernelBegin = std::chrono::steady_clock::now();
#ifdef STMS_GPU_VERIFY_ENABLED
        const bool usedGpu =
            finalProposalSweep && gpuFastBlocks(c, blocks, gbits, special, threads);
        if (!usedGpu) {
#endif
#pragma omp parallel for schedule(dynamic, 1)
            for (size_t j = 0; j < m; ++j) {
                auto &q = blocks[j];
                int bbits = std::max(192, abits);
                auto out = finalProposalSweep ? fastblock(c,
                                                          q.z,
                                                          q.start,
                                                          q.end - q.start,
                                                          gbits,
                                                          q.offset,
                                                          special,
                                                          q.out,
                                                          q.fa,
                                                          q.fb)
                                              : kernel(c,
                                                       q.z,
                                                       q.start,
                                                       q.end - q.start,
                                                       gbits,
                                                       abits,
                                                       bbits,
                                                       q.offset,
                                                       special,
                                                       false,
                                                       zoom,
                                                       q.out,
                                                       q.fa,
                                                       q.fb);
                q.logA = out.logA;
            }
#ifdef STMS_GPU_VERIFY_ENABLED
        }
#endif
        const auto compositionBegin = std::chrono::steady_clock::now();
        result.kernelSeconds +=
            std::chrono::duration<double>(compositionBegin - kernelBegin).count();
        checkCancellation();
        zero(tv[0]);
        zero(uv[0]);
        double prodAlog = 0, phase = 0;
        for (size_t j = 0; j < m; j++) {
            auto &q = blocks[j];
            if (j + 1 < m)
                sub(res, q.out, blocks[j + 1].z);
            else
                set(res, q.out);
            ops.mul(temp, q.fa, tv[j]);
            add(tv[j + 1], temp, res);
            ops.mul(temp, q.fa, uv[j]);
            add(uv[j + 1], temp, q.fb);
            prodAlog += la(q.fa);
            phase += angle(q.fa);
        }
        ops.div(step, tv[m], uv[m]);
        mpf_neg(step.r, step.r);
        mpf_neg(step.i, step.i);
        prevstep = logstep;
        logstep = la(step) / std::log2(10.);
        add(c, c, step);
        double maxr = -1e300;
        for (size_t j = 1; j < m; j++) {
            ops.mul(temp, uv[j], step);
            add(temp, temp, tv[j]);
            add(blocks[j].z, blocks[j].z, temp);
            maxr = std::max(maxr, la(temp) / std::log2(10.) - blocks[j].offset / std::log2(10.));
        }
        logAB = (prodAlog + la(uv[m])) / std::log2(10.);
        if (trace)
            std::cerr << "shoot it=" << it << " goal=" << goal << " logCorrection=" << logstep
                      << " stateCorrection=" << maxr << " logAB=" << logAB << "\n";
        result.iterations = it + 1;
        result.logAB = logAB;
        result.phase = std::remainder(phase + angle(uv[m]), 2 * std::acos(-1.));
        result.logCorrection = logstep;
        const auto compositionSeconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - compositionBegin)
                .count();
        result.compositionSeconds += compositionSeconds;
        result.stages.push_back(
            {it,
             goal,
             target,
             logstep,
             maxr,
             gbits,
             abits,
             finalProposalSweep ? 1 : 0,
             std::chrono::duration<double>(compositionBegin - kernelBegin).count(),
             compositionSeconds});
        if (logstep < -target && maxr < -target && goal >= target) {
            result.checkpoints = std::move(blocks);
            return result;
        }
        curvature = std::max(logstep - 2 * prevstep + 2, zoom - 30);
        finalProposalSweep = goal >= target && 2 * logstep + curvature < -target - 6;
    }
    throw std::runtime_error("multiple-shooting Newton did not converge");
}

#ifndef STMS_LIBRARY
int main(int argc, char **argv) {
    try {
        auto begin = std::chrono::steady_clock::now();
        if (argc < 2) {
            std::cerr << "Usage: locate input.rfl [--threads 4] [--blocks N] [--trace]\n";
            return 2;
        }
        int threads = 4, blocks = 0;
        bool trace = false;
        for (int j = 2; j < argc; j++) {
            std::string arg = argv[j];
            if (arg == "--trace")
                trace = true;
            else if ((arg == "--threads" || arg == "--blocks") && j + 1 < argc) {
                int x = std::stoi(argv[++j]);
                if (x < 1 || x > 256)
                    throw std::runtime_error("invalid thread/block count");
                if (arg == "--threads")
                    threads = std::min(x, 64);
                else
                    blocks = x;
            } else
                throw std::runtime_error("unknown or incomplete argument: " + arg);
        }
        Input input = read_input(argv[1]);
        const unsigned bits = std::ceil((2 * double(input.zoom) + 160) * std::log2(10.));
        C c(bits);
        if (mpf_set_str(c.r, input.real.c_str(), 10) || mpf_set_str(c.i, input.imag.c_str(), 10))
            throw std::runtime_error("invalid coordinate");
        auto result = solve(c, input.zoom, threads, blocks, trace);
        char *real = nullptr, *imag = nullptr;
        gmp_asprintf(&real, "%.*Ff", int(bits / std::log2(10.)), (mpf_srcptr) c.r);
        gmp_asprintf(&imag, "%.*Ff", int(bits / std::log2(10.)), (mpf_srcptr) c.i);
        const double seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        std::cout << std::setprecision(17) << "{\n  \"schema\": \"nucleus-prototype-v1\",\n"
                  << "  \"method\": \"" << result.method << "\",\n"
                  << "  \"period\": " << result.period << ",\n"
                  << "  \"center_real\": \"" << real << "\",\n"
                  << "  \"center_imag\": \"" << imag << "\",\n"
                  << "  \"input_log_zoom\": " << input.zoom << ",\n"
                  << "  \"input_u64_field\": " << input.field << ",\n"
                  << "  \"log10_abs_A_times_D\": " << result.logAB << ",\n"
                  << "  \"arg_A_times_D\": " << result.phase << ",\n"
                  << "  \"output_log_zoom\": null,\n"
                  << "  \"zoom_status\": \"RFF-compatible zoom not implemented\",\n"
                  << "  \"iterations_including_pilot\": " << result.iterations << ",\n"
                  << "  \"log10_last_parameter_correction\": " << result.logCorrection << ",\n"
                  << "  \"threads_requested\": " << threads << ",\n"
                  << "  \"internal_seconds_including_input\": " << seconds << "\n}\n";
        free(real);
        free(imag);
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}

#endif
