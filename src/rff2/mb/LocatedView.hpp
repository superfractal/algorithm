// Created by GPT-6 on 2026-09-10
// Modified by GPT-6 on 2026-09-12
#pragma once
#include "STMSBridge.h"
#include "MB2RenderData.hpp"
namespace located_view {
    using namespace merutilm::rff2;
    struct View {
        FractalSettings settings;
        dex dcMax;
    };
    inline View fromProposal(const FractalSettings &input,
                             dex inputDcMax,
                             const stms_bridge::Proposal &proposal) {
        View out{input, {}};
        out.settings.general.logZoom = static_cast<float>(proposal.logZoom) + 1.5f;
        const int exponent = Perturbator::logZoomToExp10(out.settings.general.logZoom);
        out.settings.reference.center =
            fixed_point_complex_i1(proposal.real, proposal.imag, exponent);
        out.dcMax =
            inputDcMax * rff_math::exp10(input.general.logZoom - out.settings.general.logZoom);
        return out;
    }
    struct Located {
        View view;
        stms_bridge::Proposal proof;
    };
    inline std::unique_ptr<Located> locate(ParallelRenderState &state,
                                           const MB2RenderDataBase &input) {
        if (state.interruptRequested() || !input.getReference() || !input.getPerturbator())
            return nullptr;
        try {
            const auto &settings = input.fractalSettings;
            auto center = settings.reference.center;
            auto proof = stms_bridge::locate(state,
                                             center.real.to_string(),
                                             center.imag.to_string(),
                                             settings.general.logZoom,
                                             input.getReference()->longestPeriod(),
                                             settings.general.threads);
            if (!proof || state.interruptRequested())
                return nullptr;
            auto view = fromProposal(settings, input.getPerturbator()->dcMax, *proof);
            if (state.interruptRequested())
                return nullptr;
            return std::make_unique<Located>(std::move(view), std::move(*proof));
        } catch (...) {
            if (state.interruptRequested())
                return nullptr;
            throw;
        }
    }
} // namespace located_view
