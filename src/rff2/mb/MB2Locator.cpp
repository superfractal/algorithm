//
// Created by Merutilm on 2025-05-16.
//

// Modified by GPT-6 on 2026-09-10, 2026-09-12
#include "MB2Locator.h"
#include "STMSBridge.h"
#include "LocatedView.hpp"

#include "MB2Reference.h"
#include "MB2RenderData.hpp"
#include "Perturbator.h"


namespace merutilm::rff2 {

    std::unique_ptr<fixed_point_complex_i1>
    MB2Locator::findCenterOffset(const MB2RenderDataBase &data) {
        const int exp10 = Perturbator::logZoomToExp10(data.fractalSettings.general.logZoom);
        const MB2ReferenceBase *reference = data.getReference();
        if (!reference)
            return nullptr;

        fixed_point_complex bn = reference->fpgBn.create_variant(exp10, -exp10 * 2);
        fixed_point_complex z = reference->fpgReference.create_variant(exp10, -exp10 * 2);
        fixed_point_complex::neg(bn);
        fixed_point_complex::div(z, z, bn);
        return std::make_unique<fixed_point_complex_i1>(z.real, z.imag, exp10);
    }

    std::unique_ptr<MB2Locator> MB2Locator::locateMinibrot(
        ParallelRenderState &state,
        const MB2RenderDataBase &data,
        std::unique_ptr<ApproxTableCacheBase> &cache,
        const std::function<void(uint64_t, int)> &actionWhileFindingMinibrotCenter,
        const std::function<void(uint64_t, float)> &actionWhileSeriesApprox,
        const std::function<void(uint64_t, float)> &actionWhileCreatingTable,
        const std::function<void(float)> &actionWhileFindingMinibrotZoom) {
        if (state.interruptRequested() || !data.getReference() || !data.getPerturbator())
            return nullptr;
        const auto period = data.getReference()->longestPeriod();
        actionWhileFindingMinibrotCenter(0, 1);
        auto inputCenter = data.fractalSettings.reference.center;
        std::optional<stms_bridge::Proposal> proposal;
        try {
            proposal = stms_bridge::locate(state,
                                           inputCenter.real.to_string(),
                                           inputCenter.imag.to_string(),
                                           data.fractalSettings.general.logZoom,
                                           period,
                                           data.fractalSettings.general.threads);
        } catch (const std::exception &e) {
            vkh::logger::log_err("Locate Minibrot: {}", e.what());
            return nullptr;
        }
        if (!proposal || state.interruptRequested())
            return nullptr;
        auto view = located_view::fromProposal(
            data.fractalSettings, data.getPerturbator()->dcMax, *proposal);
        actionWhileFindingMinibrotZoom(view.settings.general.logZoom);
        if (state.interruptRequested())
            return nullptr;
        return std::make_unique<MB2Locator>(std::move(view.settings), view.dcMax);
    }

    /**
     * This method moves the data, so the paramed data is no longer available.
     * Use the return value instead of this.
     * @return result table
     */
    std::unique_ptr<MB2RenderDataBase> MB2Locator::findAccurateCenterPerturbator(
        ParallelRenderState &state,
        const MB2RenderDataBase &data,
        std::unique_ptr<ApproxTableCacheBase> &cache,
        const std::function<void(uint64_t, int)> &actionWhileFindingMinibrotCenter,
        const std::function<void(uint64_t, float)> &actionWhileSeriesApprox,
        const std::function<void(uint64_t, float)> &actionWhileCreatingTable) {
        // multiply zoom by 2 and find center offset.
        // set the center to center + centerOffset.

        uint64_t longestPeriod = data.getReference()->longestPeriod();
        uint64_t refLen = data.getReference()->length();

        const float logZoom = data.fractalSettings.general.logZoom;
        const FractalSettings &calc = data.fractalSettings;
        FractalSettings doubledZoomCalc = calc;
        const float doubledLogZoom = logZoom * 2;
        const int doubledExp10 = Perturbator::logZoomToExp10(doubledLogZoom);

        doubledZoomCalc.general.logZoom = doubledLogZoom;
        doubledZoomCalc.perturb.absoluteIterationMode = false;
        doubledZoomCalc.perturb.decimalizeIterationMethod = FrtDecimalizeIterationMethod::NONE;


        dex doubledZoomDcMax = data.getPerturbator()->dcMax / rff_math::exp10(logZoom);


        int centerFixCount = 0;

        std::unique_ptr<MB2RenderDataBase> doubledZoomData = nullptr;

        while (doubledZoomData == nullptr || !doubledZoomData->getPerturbator() ||
               !checkMaxIterationOnly(*doubledZoomData)) {
            if (state.interruptRequested()) {
                return nullptr;
            }

            auto center = doubledZoomCalc.reference.center.create_variant(doubledExp10);
            auto centerOffset =
                findCenterOffset(doubledZoomData == nullptr ? data : *doubledZoomData)
                    ->create_variant(doubledExp10);

            fixed_point_complex::add(center, center, centerOffset);

            if (centerOffset.is_strict_zero()) {
                vkh::logger::log_err(
                    "The center could not be found, or you are already in the center");
                return nullptr;
            }
            doubledZoomCalc.reference.center = center;
            ++centerFixCount;

            if (doubledLogZoom < Constants::Fractal::COMPUTESHADER_ZOOM_THRESHOLD) {
                doubledZoomData = std::make_unique<FloatMB2RenderData>(
                    state,
                    doubledZoomCalc,
                    cache,
                    doubledZoomDcMax,
                    Perturbator::logZoomToExp10(doubledLogZoom),
                    refLen,
                    longestPeriod,
                    [&actionWhileFindingMinibrotCenter, &centerFixCount](const uint64_t p) {
                        actionWhileFindingMinibrotCenter(p, centerFixCount);
                    },
                    actionWhileSeriesApprox,
                    actionWhileCreatingTable);

            } else if (doubledLogZoom < Constants::Fractal::MULTITHREAD_ZOOM_THRESHOLD) {
                doubledZoomData = std::make_unique<DoubleMB2RenderData>(
                    state,
                    doubledZoomCalc,
                    cache,
                    doubledZoomDcMax,
                    Perturbator::logZoomToExp10(doubledLogZoom),
                    refLen,
                    longestPeriod,
                    [&actionWhileFindingMinibrotCenter, &centerFixCount](const uint64_t p) {
                        actionWhileFindingMinibrotCenter(p, centerFixCount);
                    },
                    actionWhileSeriesApprox,
                    actionWhileCreatingTable);

            } else {
                doubledZoomData = std::make_unique<DexMB2RenderData>(
                    state,
                    doubledZoomCalc,
                    cache,
                    doubledZoomDcMax,
                    Perturbator::logZoomToExp10(doubledLogZoom),
                    refLen,
                    longestPeriod,
                    [&actionWhileFindingMinibrotCenter, &centerFixCount](const uint64_t p) {
                        actionWhileFindingMinibrotCenter(p, centerFixCount);
                    },
                    actionWhileSeriesApprox,
                    actionWhileCreatingTable);
            }
        }
        return doubledZoomData;
    }

    bool MB2Locator::checkMaxIterationOnly(const MB2RenderDataBase &renderData) {

        const auto it = static_cast<uint64_t>(renderData.getPerturbator()->iterate(
            {renderData.getPerturbator()->dcMax, renderData.getPerturbator()->dcMax / dex(2)}));

        return it == renderData.fractalSettings.perturb.maxIteration;
    }
} // namespace merutilm::rff2
