// Created by GPT-6 on 2026-09-11
// Modified by GPT-6 on 2026-09-12
#pragma once
#include <cstdlib>
#include "ntt_gpu.hpp"
struct OrbitGPU : NTTGPU {
    static unsigned configuredFlags();
    static std::string shaderFile(const char *path, unsigned flags);
    explicit OrbitGPU(const char *path, unsigned options = configuredFlags()) :
        NTTGPU(shaderFile(path, options).c_str()), flags(options) {
        if (flags > 15)
            throw std::runtime_error("orbit options");
        profile = std::getenv("RFF_GPU_PROFILE") != nullptr;
    }
    unsigned flags = 15;
    unsigned products() const {
        return (flags & 2) ? 2 : 3;
    }
    struct Segment {
        mpz_srcptr x, y, cr, ci;
        mpz_ptr resultX, resultY;
    };
    unsigned fraction = 0, iterations = 0, segments = 0;
    // pending: a submitted fence must be drained, including during cancellation.
    // resident: device orbit state can continue in the next same-shape chunk.
    bool uploadState = true, downloadState = true, pending = false, resident = false;
    bool profile = false;
    std::vector<unsigned> phases;
    std::vector<double> phaseSeconds;
    std::chrono::steady_clock::time_point callStart, submitStart;
    std::vector<Segment> pendingJobs;
    // Per step and component: eight 16-bit magnitude digits, discarded-limb count, sign.
    std::vector<uint32_t> trace;
    void
    orbitPlan(unsigned bits, unsigned steps, unsigned count, bool first = true, bool last = true);
    void submitOrbit(
        unsigned bits,
        unsigned steps,
        const std::vector<Segment> &jobs,
        bool first,
        bool last,
        const std::function<bool()> &cancel = [] {
            return false;
        });
    Timing finishOrbit(const std::function<bool()> &cancel = [] {
        return false;
    });
    Timing orbit(
        unsigned bits,
        unsigned steps,
        const std::vector<Segment> &jobs,
        const std::function<bool()> &cancel = [] {
            return false;
        });
};
