// Created by GPT-6 on 2026-09-10
// Modified by GPT-6 on 2026-09-10
#pragma once
#include <cstdint>
#include <optional>
#include <string>
namespace merutilm::rff2 { class ParallelRenderState; }
namespace stms_bridge {
struct Proposal {
    std::string real, imag, method, scaledResidual, searchProfile;
    double logZoomLower=0,logZoomUpper=0,logZoom=0, solveSeconds=0, verifySeconds=0;
    uint64_t period=0;
};
std::optional<Proposal> locate(merutilm::rff2::ParallelRenderState&, const std::string&,
    const std::string&, double zoom, uint64_t expectedPeriod, unsigned threads);
}
