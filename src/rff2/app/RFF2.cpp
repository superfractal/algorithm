// Modified by GPT-6 on 2026-09-07, 2026-09-10
//
// Created by Merutilm on 2025-08-08.
//

#include "RFF2.hpp"

#include <ranges>
#include <chrono>
#include "../util/sha256.hpp"

#include "../io/RFFLocationBinary.h"
#include "../mb/MB2Locator.h"
#include "../parallel/ParallelArrayDispatcher.h"
#include "../preset/calc/approx/ClcApproxPresets.hpp"
#include "../preset/calc/compress/ClcCompressPresets.hpp"
#include "../preset/calc/sync/ClcSyncPresets.hpp"
#include "../preset/render/compute/RndComputePresets.hpp"
#include "../preset/render/display/RndDisplayPresets.hpp"
#include "../preset/shader/bloom/ShdBloomPresets.hpp"
#include "../preset/shader/color/ShdColorPresets.hpp"
#include "../preset/shader/fog/ShdFogPresets.hpp"
#include "../preset/shader/palette/ShdPalettePresets.hpp"
#include "../preset/shader/slope/ShdSlopePresets.hpp"
#include "../preset/shader/stripe/ShdStripePresets.hpp"
#include "../vulkan/GPCDownsampleForBlur.hpp"
#include "../vulkan/SharedImageContextIndices.hpp"
#include "../vulkan/desc/SharedDescriptorTemplate.hpp"
#include "FnExplore.hpp"
#include "FnFile.hpp"
#include "FnFractal.hpp"
#include "FnPreset.hpp"
#include "FnRender.hpp"
#include "FnShader.hpp"
#include "FnVideo.hpp"
#include "IOUtilities.h"
#include "Utilities.h"
#include "imgui.h"
#include "nfd.hpp"
#include "opencv2/opencv.hpp"
#include "vulkan_helper/engine/executor/ScopedNewCommandBufferExecutor.hpp"
#include "vulkan_helper/engine/window/PlatformWindow.hpp"
#include "vulkan_helper/util/BarrierUtils.hpp"
#include "vulkan_helper/util/BufferImageContextUtils.hpp"


namespace merutilm::rff2 {


    void RFF2::onStart() {
        initialize();
        applyShaderSettings(settings);
        refreshResizeParams(rootWindowContext->getSwapchain().getSwapchainExtent());
        requests.requestRecompute();
    }

    void RFF2::initialize() {
        cursorManager = std::make_unique<CursorManager>(rootWindowContext->getWindow()->getWindow());
        computeShaderManager = std::make_unique<ComputeShaderRenderManager>(*rootWindowContext);
        NFD::Init();
        initImGui();
    }

    void RFF2::onResize(const VkExtent2D newExtent) {
        if (batchTest.active() && batchTest.phase != BatchTest::Phase::SETUP) batchTest.cancel = true;
        state.cancel();
        recreateContexts(newExtent);
        if (newExtent.width > 0 || newExtent.height > 0) {
            engine->getCore().getLogicalDevice().waitDeviceIdle();
            state.cancel();
            refreshResizeParams(newExtent);
            requests.requestRecompute();
            backgroundThreads.notifyAll();
        }
    }


    void RFF2::onQuit() {
        state.cancel();
        if (batchTest.active() && batchTest.log.is_open()) {
            try { batchTest.log << "cancelled: application closed.\n"; batchTest.log.close(); }
            catch (...) { /* An incomplete log cannot be mistaken for a passing run. */ }
        }
        renderer = nullptr;
        NFD::Quit();
    }


    void RFF2::resolveRequests() {

        if (requests.defaultSettingsRequested) {
            applyDefaultSettings();
            requests.defaultSettingsRequested.exchange(false);
            backgroundThreads.notifyAll();
        }

        if (requests.shaderRequested) {
            applyShaderSettings(settings);
            requests.shaderRequested.exchange(false);
            backgroundThreads.notifyAll();
        }

        {
            std::scoped_lock lock2(requests.resizeMutex);
            if (requests.resizeRequested.exchange(false)) {
                onResize(requests.resizeRequestedExtent);
                backgroundThreads.notifyAll();
            }
        }
        {
            std::scoped_lock lock2(requests.createImageMutex);
            if (requests.createImageRequested.exchange(false)) {
                applyCreateImage();
                backgroundThreads.notifyAll();
            }
        }

        auto expected = ComputeState::REQUESTED;
        if (requests.recomputeRequestedState.compare_exchange_weak(expected, ComputeState::RUNNING)) {
            recomputeThreaded();
            // it is threaded, do not notify
        }
    }
    std::unique_ptr<MB2RenderDataBase>
    RFF2::createAppropriateRenderData(const bool computeShader, const float logZoomTest, const float startTime,
                                      const FractalSettings &frt, const dex dcMax, const int exp10,
                                      const uint64_t refInitialCapacity, const uint64_t forcedStrictFPGPeriod) {
        if (computeShader) {
            if (logZoomTest > Constants::Fractal::COMPUTESHADER_ZOOM_THRESHOLD) {
                return std::make_unique<FexMB2RenderData>(state, frt, approxTableCache, dcMax, exp10,
                                                          refInitialCapacity, forcedStrictFPGPeriod,
                                                          FnExplore::getActionWhileRefCalc(*this, startTime),
                                                          FnExplore::getActionWhileSeriesApprox(*this, startTime),
                                                          FnExplore::getActionWhileCreatingTable(*this, startTime));
            } else {
                return std::make_unique<FloatMB2RenderData>(state, frt, approxTableCache, dcMax, exp10,
                                                            refInitialCapacity, forcedStrictFPGPeriod,
                                                            FnExplore::getActionWhileRefCalc(*this, startTime),
                                                            FnExplore::getActionWhileSeriesApprox(*this, startTime),
                                                            FnExplore::getActionWhileCreatingTable(*this, startTime));
            }
        } else {
            if (logZoomTest > Constants::Fractal::MULTITHREAD_ZOOM_THRESHOLD) {
                return std::make_unique<DexMB2RenderData>(state, frt, approxTableCache, dcMax, exp10,
                                                          refInitialCapacity, forcedStrictFPGPeriod,
                                                          FnExplore::getActionWhileRefCalc(*this, startTime),
                                                          FnExplore::getActionWhileSeriesApprox(*this, startTime),
                                                          FnExplore::getActionWhileCreatingTable(*this, startTime));
            } else {
                return std::make_unique<DoubleMB2RenderData>(state, frt, approxTableCache, dcMax, exp10,
                                                             refInitialCapacity, forcedStrictFPGPeriod,
                                                             FnExplore::getActionWhileRefCalc(*this, startTime),
                                                             FnExplore::getActionWhileSeriesApprox(*this, startTime),
                                                             FnExplore::getActionWhileCreatingTable(*this, startTime));
            }
        }
    }

    void RFF2::updateMouseInteraction() {
        double mdx;
        double mdy;
        glfwGetCursorPos(rootWindowContext->getWindow()->getWindow(), &mdx, &mdy);
        const auto mx = static_cast<int>(mdx);
        const auto my = static_cast<int>(mdy);
        const uint16_t x = getMouseXOnIterationBuffer(mx);
        const uint16_t y = getMouseYOnIterationBuffer(my);
        if (renderer->visibleIterationBufferContext == nullptr || x >= getIterationBufferWidth() ||
            y >= getIterationBufferHeight()) {
            return;
        }
        auto it = static_cast<uint64_t>((*renderer->visibleIterationBufferContext)(x, y));
        setStatusMessage(Constants::Status::ITERATION_STATUS,
                         std::format(std::locale("en_US.UTF-8"), "Iterations : {:L}", it));
    }

    void RFF2::update() {
        updateMouseInteraction();
        updateBatchTest();
        resolveRequests();
        invokeUpdaters();
        renderer->render();
    }


    Settings RFF2::genDefaultSettings() {
#ifndef NDEBUG
        return Settings{
                .fractal =
                        FractalSettings{.general = {.bailout = 2.00001f, .logZoom = 2, .threads = 15},
                                        .reference =
                                                {
                                                        .center = fixed_point_complex_i1(
                                                                "-0.85", "0", Perturbator::logZoomToExp10(2)),
                                                        .useParallelRefCalculation = false,
                                                        .sync = ClcSyncPresets::Fast().genRefSync(),
                                                        .compression = ClcCompressPresets::None().genRefComp(),
                                                        .periodMultiplier = 1,
                                                        .reuse = false,
                                                },
                                        .sa = {.use = false,
                                               .appliedTermsCount = 8,
                                               .validatedTermsCount = 1,
                                               .epsilonPower = -5},
                                        .mpa = ClcApproxPresets::UltraFast().genMPA(),
                                        .perturb = {.maxIteration = 300,
                                                    .decimalizeIterationMethod = FrtDecimalizeIterationMethod::LOG_LOG,
                                                    .autoMaxIteration = true,
                                                    .interiorDetectRadiusPower = 12,
                                                    .autoIterationMultiplier = 100,
                                                    .absoluteIterationMode = false}},
                .render = {.display = RndDisplayPresets::Low().genDisplay(),
                           .computeShader = RndComputePresets::None().genComputeShader()},
                .shader = {.palette = ShdPalettePresets::Classic1().genPalette(),
                           .stripe = ShdStripePresets::Disabled().genStripe(),
                           .slope = ShdSlopePresets::Disabled().genSlope(),
                           .color = ShdColorPresets::Disabled().genColor(),
                           .fog = ShdFogPresets::Disabled().genFog(),
                           .bloom = ShdBloomPresets::Disabled().genBloom(),
                           .noiseReduction = {false, 2, 0.1f},
                           .fractal3D = {false, 85, 0, 1, 0, 10.f}},
                .video = {.data = {.defaultZoomIncrement = 2, .isStatic = false},
                          .animation = {.overZoom = 2, .showText = true, .mps = 1},
                          .exportation = {.fps = 60, .bitrate = 9000}},
                .explore = {.autoMoveCursorToCenter = false}};
#else
        return Settings{
                .fractal =
                        FractalSettings{.general = {.bailout = 2.00001f,
                                                    .logZoom = 2,
                                                    .threads = std::thread::hardware_concurrency() - 1},
                                        .reference =
                                                {
                                                        .center = fixed_point_complex_i1(
                                                                "-0.85", "0", Perturbator::logZoomToExp10(2)),
                                                        .useParallelRefCalculation = false,
                                                        .sync = ClcSyncPresets::Fast().genRefSync(),
                                                        .compression = ClcCompressPresets::None().genRefComp(),
                                                        .periodMultiplier = 1,
                                                        .reuse = false,
                                                },
                                        .sa = {.use = false,
                                               .appliedTermsCount = 8,
                                               .validatedTermsCount = 1,
                                               .epsilonPower = -5},
                                        .mpa = ClcApproxPresets::UltraFast().genMPA(),
                                        .perturb = {.maxIteration = 300,
                                                    .decimalizeIterationMethod = FrtDecimalizeIterationMethod::LOG_LOG,
                                                    .autoMaxIteration = true,
                                                    .interiorDetectRadiusPower = 12,
                                                    .autoIterationMultiplier = 100,
                                                    .absoluteIterationMode = false}},
                .render = {.display = RndDisplayPresets::High().genDisplay(),
                           .computeShader = RndComputePresets::None().genComputeShader()},
                .shader = {.palette = ShdPalettePresets::Classic1().genPalette(),
                           .stripe = ShdStripePresets::Disabled().genStripe(),
                           .slope = ShdSlopePresets::Disabled().genSlope(),
                           .color = ShdColorPresets::Disabled().genColor(),
                           .fog = ShdFogPresets::Disabled().genFog(),
                           .bloom = ShdBloomPresets::Disabled().genBloom(),
                           .noiseReduction = {false, 2, 0.1f},
                           .fractal3D = {false, 85, 0, 1, 0, 10.f}},
                .video = {.data = {.defaultZoomIncrement = 2, .isStatic = false},
                          .animation = {.overZoom = 2, .showText = true, .mps = 1},
                          .exportation = {.fps = 60, .bitrate = 9000}},
                .explore = {.autoMoveCursorToCenter = false}};
#endif
    }

    complex<dex> RFF2::offsetConversion(const Settings &s, const int px, const int py) const {
        const double bufOffX = static_cast<double>(px) - static_cast<double>(getIterationBufferWidth()) / 2.0;
        const double bufOffY = static_cast<double>(py) - static_cast<double>(getIterationBufferHeight()) / 2.0;
        return complex{dex(bufOffX), dex(bufOffY)} / getDivisor(s) / dex(s.render.display.clarityMultiplier);
    }

    std::array<int, 2> RFF2::iterationBufferConversion(const Settings &s, const complex<dex> &offset) const {
        const auto [re, im] =
                static_cast<complex<double>>(offset * dex(s.render.display.clarityMultiplier) * getDivisor(s));

        const auto px = static_cast<int>((re < 0 ? std::round(re) : std::ceil(re)) + getIterationBufferWidth() / 2.0);
        const auto py = static_cast<int>((im < 0 ? std::round(im) : std::ceil(im)) + getIterationBufferHeight() / 2.0);
        return {px, py};
    }

    void RFF2::moveCursor(const int px, const int py) const {
        const auto mx = static_cast<int>(static_cast<float>(px) / settings.render.display.clarityMultiplier);
        const auto my = static_cast<int>(static_cast<float>(py) / settings.render.display.clarityMultiplier);
        glfwSetCursorPos(rootWindowContext->getWindow()->getWindow(), mx,
                         rootWindowContext->getSwapchain().getSwapchainExtent().height - my);
    }

    dex RFF2::getDivisor(const Settings &settings) { return rff_math::exp10(settings.fractal.general.logZoom); }


    uint16_t RFF2::calcIterationBufferWidth(const Settings &s) const {
        const float multiplier = s.render.display.clarityMultiplier;
        return static_cast<uint16_t>(static_cast<float>(rootWindowContext->getSwapchain().getSwapchainExtent().width) *
                                     multiplier);
    }

    uint16_t RFF2::calcIterationBufferHeight(const Settings &s) const {
        const float multiplier = s.render.display.clarityMultiplier;
        return static_cast<uint16_t>(static_cast<float>(rootWindowContext->getSwapchain().getSwapchainExtent().height) *
                                     multiplier);
    }

    uint16_t RFF2::getIterationBufferWidth() const { return renderer->visibleIterationBufferContext->getWidth(); }

    uint16_t RFF2::getIterationBufferHeight() const { return renderer->visibleIterationBufferContext->getHeight(); }


    void RFF2::addListeners() {
        auto &eventSystem = rootWindowContext->getWindow()->eventSystem;

        eventSystem.applicationLifecycle.onUpdate.add([this] { update(); });

        eventSystem.resize.onResize.add([this](const int w, const int h) {
            const auto extent = VkExtent2D{static_cast<uint32_t>(w), static_cast<uint32_t>(h)};
            requests.requestResize(extent);
        });

        eventSystem.applicationLifecycle.onStart.add([this] { onStart(); });

        eventSystem.applicationLifecycle.onQuit.add([this] {
            rootWindowContext->core.getLogicalDevice().waitDeviceIdle();
            onQuit();
        });


        eventSystem.mouse.onMouseEnter.add(
                [this] { glfwSetCursor(cursorManager->window, cursorManager->crosshairCursor); });
        eventSystem.mouse.onMouseExit.add([this] { glfwSetCursor(cursorManager->window, nullptr); });

        eventSystem.mouseDrag.onMouseDrag.add(
                [this](const int mb, const int mx, const int my, const int mdx, const int mdy) {
                    if (batchTest.active()) return;
                    const int16_t x = getMouseXOnIterationBuffer(mx);
                    const int16_t y = getMouseYOnIterationBuffer(my);
                    const auto dx = static_cast<int16_t>(getMouseXOnIterationBuffer(mx - mdx) - x);
                    const auto dy = static_cast<int16_t>(getMouseYOnIterationBuffer(my - mdy) - y);
                    const auto dxr = -static_cast<float>(dx) / static_cast<float>(getIterationBufferWidth());
                    const auto dyr = static_cast<float>(dy) / static_cast<float>(getIterationBufferHeight());
                    const auto dz = pow(10.0f, -zoomAnimationInfo.targetLogZoomOffsetAim);

                    zoomAnimationInfo.aimChanged = true;
                    zoomAnimationInfo.targetMouseDragOffset += glm::vec2{dxr * dz, dyr * dz};

                    if (mb == GLFW_MOUSE_BUTTON_LEFT) {
                        const float m = settings.render.display.clarityMultiplier;
                        const float logZoom = settings.fractal.general.logZoom;
                        const int exp10 = Perturbator::logZoomToExp10(logZoom);

                        fixed_point_complex_i1 &center = settings.fractal.reference.center;
                        center.set_exp10(exp10);
                        const fixed_point_complex_i1 add(dex(static_cast<float>(dx) / m) / getDivisor(settings),
                                                         dex(static_cast<float>(dy) / m) / getDivisor(settings), exp10);
                        fixed_point_complex_i1::add(center, center, add);

                        requests.requestRecompute();
                    }
                });
        eventSystem.mouseWheel.onMouseScroll.add([this](const int value) {
            if (batchTest.active()) return;
            settings.fractal.general.logZoom = std::max(Constants::Fractal::ZOOM_MIN, settings.fractal.general.logZoom);
            double mdx;
            double mdy;
            glfwGetCursorPos(rootWindowContext->getWindow()->getWindow(), &mdx, &mdy);
            const auto mx = static_cast<int>(mdx);
            const auto my = static_cast<int>(mdy);
            const int16_t mix = getMouseXOnIterationBuffer(mx);
            const int16_t miy = getMouseYOnIterationBuffer(my);
            zoom(mix, miy, value > 0 ? Constants::Fractal::ZOOM_INTERVAL : -Constants::Fractal::ZOOM_INTERVAL);
        });
    }


    void RFF2::zoom(const int16_t px, const int16_t py, const float logIncrement) {

        settings.fractal.general.logZoom = std::max(Constants::Fractal::ZOOM_MIN, settings.fractal.general.logZoom);
        const int16_t mix = px;
        const int16_t miy = py;
        const auto mxr = static_cast<float>(mix) / static_cast<float>(getIterationBufferWidth()) - 0.5f;
        const auto myr = static_cast<float>(miy) / static_cast<float>(getIterationBufferHeight()) - 0.5f;
        const auto dz = pow(10.0f, -zoomAnimationInfo.targetLogZoomOffsetAim);

        const auto [re, im] = offsetConversion(settings, mix, miy);
        float &logZoom = settings.fractal.general.logZoom;
        fixed_point_complex_i1 &center = settings.fractal.reference.center;
        const int exp10 = Perturbator::logZoomToExp10(logZoom);
        center.set_exp10(exp10);

        const float mz = pow(10.0f, -logIncrement);
        logZoom += logIncrement;
        const fixed_point_complex_i1 add(re * dex(1 - mz), im * dex(1 - mz), exp10);
        fixed_point_complex_i1::add(center, center, add);

        zoomAnimationInfo.aimChanged = true;
        zoomAnimationInfo.stop();
        zoomAnimationInfo.targetLogZoomOffsetAim += logIncrement;
        zoomAnimationInfo.targetMouseZoomOffsetAim += glm::vec2{mxr * dz * (mz - 1), myr * dz * (1 - mz)};
        requests.requestRecompute();
    }


    void RFF2::applyDefaultSettings() {
        rootWindowContext->core.getLogicalDevice().waitDeviceIdle();
        settings = genDefaultSettings();
    }


    void RFF2::applyCreateImage() {
        const uint32_t frameIndex = renderer->getFrameIndex();
        rootWindowContext->getSyncObject().getFence(frameIndex).wait();

        if (requests.createImageRequestedFilename.empty()) {
            const auto path = IOUtilities::ioFileDialog(Constants::File::DESC_IMAGE, IOUtilities::SAVE_FILE,
                                                        Constants::File::EXT_IMAGE);

            if (path == nullptr)
                return;

            requests.createImageRequestedFilename = path->string();
        }
        const auto &imgCtx = rootWindowContext->getSharedImageContext().getImageContextMF(
                SharedImageContextIndices::MF_MAIN_RENDER_IMAGE_SECONDARY)[frameIndex];

        vkh::BufferContext bufCtx = vkh::BufferContext::createContext(
                rootWindowContext->core,
                {
                        .size = imgCtx.capacity,
                        .usage = VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                        .properties = VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                });
        vkh::BufferContext::mapMemory(rootWindowContext->core, bufCtx);
        // NEW COMMAND BUFFER
        {
            const auto executor =
                    vkh::ScopedNewCommandBufferExecutor(rootWindowContext->core, rootWindowContext->getCommandPool());
            vkh::BarrierUtils::cmdImageMemoryBarrier(
                    executor.getCommandBuffer().getCommandBufferHandle(), imgCtx.image, VK_ACCESS_SHADER_WRITE_BIT,
                    VK_ACCESS_TRANSFER_READ_BIT, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL,
                    VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, 0, 1, VK_IMAGE_ASPECT_COLOR_BIT,
                    VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT);
            vkh::BufferImageContextUtils::cmdCopyImageToBuffer(executor.getCommandBuffer(), imgCtx, bufCtx);
        }
        auto img = cv::Mat(static_cast<int>(imgCtx.extent.height), static_cast<int>(imgCtx.extent.width), CV_16UC4,
                           bufCtx.mappedMemory);
        cv::cvtColor(img, img, cv::COLOR_RGBA2BGRA);
        bool saved = false;
        try {
            saved = cv::imwrite(requests.createImageRequestedFilename, img);
        } catch (...) {
            vkh::BufferContext::unmapMemory(rootWindowContext->core, bufCtx);
            vkh::BufferContext::destroyContext(rootWindowContext->core, bufCtx);
            throw;
        }
        vkh::BufferContext::unmapMemory(rootWindowContext->core, bufCtx);
        vkh::BufferContext::destroyContext(rootWindowContext->core, bufCtx);
        if (!saved) throw std::runtime_error("Cannot save PNG");
    }

    void RFF2::invokeUpdaters() {
        static float time = rootWindowContext->getWindow()->getTime();
        const float t = rootWindowContext->getWindow()->getTime();
        const float dt = t - time;
        time = t;

        if (canShowPreview && !zoomAnimationInfo.aimChanged) {
            renderer->updateStagingBuffer |= renderer->visibleIterationBufferContext->fill();
            renderer->descriptorStorage->iteration->applyMaxIteration();
            zoomAnimationInfo.reset();
        }

        zoomAnimationInfo.update(dt);

        renderer->descriptorStorage->smoothZoom->set(zoomAnimationInfo.targetMouseDragOffset +
                                                             zoomAnimationInfo.targetMouseZoomOffset,
                                                     zoomAnimationInfo.targetLogZoomOffset);
    }

    void RFF2::applyShaderSettings(const Settings &s) const {
        using namespace SharedDescriptorTemplate;
        vkh::Core &core = rootWindowContext->core;
        core.getLogicalDevice().waitDeviceIdle();

        renderer->descriptorStorage->palette->set(s.shader.palette, !(batchTest.active() && !batchTest.locate));
        renderer->descriptorStorage->stripe->set(s.shader.stripe);
        renderer->descriptorStorage->color->set(s.shader.color);
        renderer->descriptorStorage->fog->set(s.shader.fog);
        renderer->descriptorStorage->bloom->set(s.shader.bloom);
        renderer->descriptorStorage->noiseReduction->set(s.shader.noiseReduction);
        renderer->descriptorStorage->camera3d->set(s.shader.fractal3D);
        renderer->descriptorStorage->fractal3d->set(s.shader.fractal3D);
    }

    void RFF2::refreshResizeParams(const VkExtent2D swapchainExtent) const {
        const uint16_t iw = calcIterationBufferWidth(settings);
        const uint16_t ih = calcIterationBufferHeight(settings);
        const auto &[dWidth, dHeight] =
                RendererUtils::getBlurredImageExtent(swapchainExtent, settings.render.display.clarityMultiplier);
        const auto &[sWidth, sHeight] = rootWindowContext->getSwapchain().getSwapchainExtent();

        // shared
        renderer->descriptorStorage->iteration->resetIterationBuffer(iw, ih);
        renderer->descriptorStorage->batchResult->resizeBatchResultBuffer(iw, ih);
        renderer->computeIterateFloat->resizeWriteBuffer(iw, ih);
        renderer->computeIterateFex->resizeWriteBuffer(iw, ih);
        renderer->descriptorStorage->renderMetaIterationVariant->resetIterationBuffer(iw, ih);

        // unique
        renderer->rccDownsample->downsample->setRescaledResolution(GPCDownsampleForBlur::DESC_INDEX_RESAMPLE_IMAGE_FOG,
                                                                   {dWidth, dHeight});
        renderer->rccDownsample->downsample->setRescaledResolution(
                GPCDownsampleForBlur::DESC_INDEX_RESAMPLE_IMAGE_BLOOM, {dWidth, dHeight});
        renderer->rccPresentPrepare->smoothZoom->setRescaledResolution({sWidth, sHeight});


        renderer->rg1->fractal3d->resetPiplineVI(iw, ih);
        renderer->computeIterateFloat->setExtent({iw, ih});
        renderer->computeIterateFex->setExtent({iw, ih});
        renderer->computeIgnoreIsolated->setExtent({iw, ih});

        renderer->visibleIterationBufferContext = std::make_unique<GraphicsMatrixBuffer<double>>(
                rootWindowContext->core, iw, ih, VK_BUFFER_USAGE_TRANSFER_SRC_BIT | VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT |
                        VK_MEMORY_PROPERTY_HOST_CACHED_BIT);
        renderer->updateStagingBuffer = true;
    }

    void RFF2::registerRenderers() {
        renderer = registerRenderer<RFF2Renderer>(*engine, *rootWindowContext, settings, zoomAnimationInfo,
                                                  [this] { renderImGui(); });
        createImGuiContext(renderer->imguiRenderContext);
    }

    void RFF2::initImGui() {

        const ImGuiIO &io = ImGui::GetIO();
        const std::filesystem::path path =
                vkh::ExecutableUtils::getExecutableDirectory() / ".." / "res" / "IBMPlexSansKR-Medium.ttf";
        io.Fonts->AddFontFromFileTTF(path.string().data(), 20.0f, nullptr, io.Fonts->GetGlyphRangesKorean());

        ImGuiStyle &style = ImGui::GetStyle();

        style.Colors[ImGuiCol_WindowBg] = ImVec4(0.08f, 0.09f, 0.11f, .8f);
        style.Colors[ImGuiCol_TitleBg] = ImVec4(0.06f, 0.07f, 0.08f, .8f);
        style.Colors[ImGuiCol_TitleBgActive] = ImVec4(0.12f, 0.14f, 0.18f, .9f);
        style.Colors[ImGuiCol_Border] = ImVec4(0.25f, 0.25f, 0.28f, .8f);
        style.Colors[ImGuiCol_Button] = ImVec4(0.24f, 0.24f, 0.24f, .9f);
        style.Colors[ImGuiCol_ButtonHovered] = ImVec4(0.34f, 0.34f, 0.34f, 1.0f);
        style.Colors[ImGuiCol_ButtonActive] = ImVec4(0.18f, 0.32f, 0.56f, 1.0f);
        style.Colors[ImGuiCol_SliderGrab] = ImVec4(0.24f, 0.50f, 0.95f, 1.0f);
        style.Colors[ImGuiCol_SliderGrabActive] = ImVec4(0.34f, 0.60f, 1.00f, 1.0f);
        style.Colors[ImGuiCol_CheckMark] = ImVec4(0.24f, 0.50f, 0.95f, 1.0f);

        style.Colors[ImGuiCol_Tab] = ImVec4(0.18f, 0.18f, 0.20f, 1.0f);
        style.Colors[ImGuiCol_TabHovered] = ImVec4(0.24f, 0.50f, 0.95f, 1.0f);
        style.Colors[ImGuiCol_TabSelected] = ImVec4(0.24f, 0.50f, 0.95f, 0.85f);
        style.Colors[ImGuiCol_TabDimmed] = ImVec4(0.16f, 0.16f, 0.17f, 1.0f);
        style.Colors[ImGuiCol_TabDimmedSelected] = ImVec4(0.20f, 0.20f, 0.22f, 1.0f);

        style.WindowRounding = 8.0f;
        style.ChildRounding = 4.0f;
        style.FrameRounding = 6.0f;
        style.GrabRounding = 6.0f;
        style.PopupRounding = 8.0f;
        style.TabRounding = 6.0f;
        style.ScrollbarRounding = 8.0f;

        style.FrameBorderSize = 0.0f;
        style.WindowBorderSize = 0.0f;
        style.ChildBorderSize = 0.0f;
    }


    void RFF2::renderImGui() {

        renderControlImGui();
        renderStatusImGui();
    }

    void RFF2::renderControlImGui() {
        ImGui::SetNextWindowSizeConstraints(ImVec2(480, 180), ImVec2(FLT_MAX, FLT_MAX));
        ImGui::Begin("Control");
        const bool testing = batchTest.active();
        ImGui::BeginDisabled(testing);
        if (ImGui::BeginTabBar("Control")) {
            if (ImGui::BeginTabItem("File")) {
                FnFile::saveMap(*this);
                FnFile::saveImage(*this);
                FnFile::saveLocation(*this);
                FnFile::loadMap(*this);
                FnFile::loadLocation(*this);
                if (ImGui::Button("Batch Render SHA-256", ImVec2(-FLT_MIN, 0))) startBatchTest(false);
                if (ImGui::Button("Batch Locate Minibrot", ImVec2(-FLT_MIN, 0))) startBatchTest(true);
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Fractal")) {
                FnFractal::reference(*this);
                FnFractal::iterations(*this);
                FnFractal::sa(*this);
                FnFractal::mpa(*this);
                FnFractal::automaticIterations(*this);
                FnFractal::absoluteIterationMode(*this);
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Render")) {
                FnRender::setResolutionProperties(*this);
                FnRender::setRenderProperties(*this);
                FnRender::setComputeShader(*this);
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Presets")) {
                FnPreset::calculation(*this);
                FnPreset::render(*this);
                FnPreset::resolution(*this);
                FnPreset::shader(*this);

                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Shader")) {
                FnShader::palette(*this);
                FnShader::stripe(*this);
                FnShader::slope(*this);
                FnShader::color(*this);
                FnShader::fog(*this);
                FnShader::bloom(*this);
                FnShader::noiseReduction(*this);
                FnShader::fractal3D(*this);
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Video")) {
                FnVideo::dataSettings(*this);
                FnVideo::animationSettings(*this);
                FnVideo::exportSettings(*this);
                FnVideo::generateVidKeyframes(*this);
                FnVideo::exportZoomVideo(*this);
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Explore")) {
                FnExplore::recompute(*this);
                FnExplore::reset(*this);
                FnExplore::cancelRender(*this);
                FnExplore::moveCursorToCenter(*this);
                FnExplore::reuseReference(*this);
                FnExplore::moveToCenter(*this);
                FnExplore::goToOriginalReference(*this);
                FnExplore::locateCenteredReference(*this);
                FnExplore::locateMinibrot(*this);
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }

        ImGui::EndDisabled();
        if (testing && ImGui::Button("Cancel Test", ImVec2(-FLT_MIN, 0))) {
            batchTest.cancel = true;
            state.interrupt();
        }
        if (!batchTest.message.empty()) ImGui::TextWrapped("%s", batchTest.message.c_str());
        ImGui::End();
    }

    void RFF2::renderStatusImGui() const {

        const float height = ImGui::GetTextLineHeight() + ImGui::GetStyle().WindowPadding.y * 2;
        ImGui::SetNextWindowPos(ImVec2(0, ImGui::GetIO().DisplaySize.y - height));

        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, height));

        ImGui::Begin("StatusBar", nullptr,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                             ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
        ImGui::PushStyleVar(ImGuiStyleVar_CellPadding, ImVec2(0, 0));
        if (ImGui::BeginTable("StatusBarTable", static_cast<int>(statusMessages.size()),
                              ImGuiTableFlags_BordersInner)) {
            for (const auto &statusMessage: statusMessages) {
                ImGui::TableNextColumn();
                ImGui::TextUnformatted(statusMessage.c_str());
            }
            ImGui::EndTable();
        }

        ImGui::PopStyleVar();
        ImGui::End();
    }

    void RFF2::refreshSharedImgContexts(const VkExtent2D extent) {
        using namespace SharedImageContextIndices;
        auto &sharedImg = rootWindowContext->getSharedImageContext();
        sharedImg.cleanupContexts();
        auto iiiGetter = [](const VkExtent2D ex, const VkFormat format, const VkImageUsageFlags usage) {
            return vkh::ImageInitInfo{
                    .imageType = VK_IMAGE_TYPE_2D,
                    .imageViewType = VK_IMAGE_VIEW_TYPE_2D,
                    .imageFormat = format,
                    .extent = {ex.width, ex.height, 1},
                    .useMipmap = VK_FALSE,
                    .arrayLayers = 1,
                    .samples = VK_SAMPLE_COUNT_1_BIT,
                    .imageTiling = VK_IMAGE_TILING_OPTIMAL,
                    .usage = usage,
                    .initialLayout = VK_IMAGE_LAYOUT_UNDEFINED,
                    .properties = VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
            };
        };

        const auto internalImageExtent =
                RendererUtils::getInternalImageExtent(extent, settings.render.display.clarityMultiplier);
        const auto blurredImageExtent =
                RendererUtils::getBlurredImageExtent(extent, settings.render.display.clarityMultiplier);

        sharedImg.appendMultiframeImageContext(MF_MAIN_RENDER_IMAGE_PRIMARY,
                                               iiiGetter(internalImageExtent, VK_FORMAT_R16G16B16A16_UNORM,
                                                         VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT |
                                                                 VK_IMAGE_USAGE_INPUT_ATTACHMENT_BIT |
                                                                 VK_IMAGE_USAGE_SAMPLED_BIT));
        sharedImg.appendMultiframeImageContext(
                MF_MAIN_RENDER_IMAGE_SECONDARY,
                iiiGetter(internalImageExtent, VK_FORMAT_R16G16B16A16_UNORM,
                          VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_INPUT_ATTACHMENT_BIT |
                                  VK_IMAGE_USAGE_TRANSFER_SRC_BIT | VK_IMAGE_USAGE_SAMPLED_BIT));
        sharedImg.appendMultiframeImageContext(
                MF_MAIN_RENDER_IMAGE_DEPTH,
                iiiGetter(internalImageExtent, VK_FORMAT_D32_SFLOAT, VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT));
        sharedImg.appendMultiframeImageContext(MF_MAIN_RENDER_DOWNSAMPLED_IMAGE_PRIMARY,
                                               iiiGetter(blurredImageExtent, VK_FORMAT_R8G8B8A8_UNORM,
                                                         VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT |
                                                                 VK_IMAGE_USAGE_SAMPLED_BIT |
                                                                 VK_IMAGE_USAGE_STORAGE_BIT));
        sharedImg.appendMultiframeImageContext(MF_MAIN_RENDER_DOWNSAMPLED_IMAGE_SECONDARY,
                                               iiiGetter(blurredImageExtent, VK_FORMAT_R8G8B8A8_UNORM,
                                                         VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT |
                                                                 VK_IMAGE_USAGE_SAMPLED_BIT |
                                                                 VK_IMAGE_USAGE_STORAGE_BIT));
    }

    void RFF2::overwriteMatrixFromMap(const RFFDynamicMapBinary &map) const {
        rootWindowContext->core.getLogicalDevice().waitDeviceIdle();
        const uint32_t iw = getIterationBufferWidth();
        const uint32_t ih = getIterationBufferHeight();
        if (iw != map.width || ih != map.height) {
            vkh::logger::log_err("Map size mismatch, {}x{} required but provided {}x{}", iw, ih, map.width, map.height);
            return;
        }

        renderer->descriptorStorage->iteration->setMaxIteration(static_cast<double>(map.maxIteration));
        renderer->descriptorStorage->iteration->applyMaxIteration();
        renderer->visibleIterationBufferContext->fill(map.iterations);
        renderer->updateStagingBuffer = true;
    }

    std::filesystem::path RFF2::getBackupLocationPath() {
        return vkh::ExecutableUtils::getExecutableDirectory() /
               std::format("{}.{}", Constants::File::BACKUP_FILE_NAME, Constants::File::EXT_LOCATION);
    }

    void RFF2::saveBackup() const {
        const auto path = getBackupLocationPath();
        saveCurrentLocation(path);
    }

    void RFF2::saveCurrentLocation(const std::filesystem::path &path) const {
        auto frt = settings.fractal; // clone the settings
        auto &center = frt.reference.center;
        RFFLocationBinary(frt.general.logZoom, center.real.to_string(), center.imag.to_string(),
                          frt.perturb.maxIteration)
                .exportFile(path);
    }

    void RFF2::loadLocation(const std::filesystem::path &path) {
        const RFFLocationBinary location = RFFLocationBinary::read(path);

        settings.fractal.reference.center = fixed_point_complex_i1(location.getReal(), location.getImag(),
                                                                   Perturbator::logZoomToExp10(location.getLogZoom()));
        settings.fractal.general.logZoom = location.getLogZoom();
        settings.fractal.perturb.maxIteration = location.getMaxIteration();
        requests.requestRecompute();
    }

    int16_t RFF2::getMouseXOnIterationBuffer(const int mx) const {
        const float multiplier = settings.render.display.clarityMultiplier;
        return static_cast<int16_t>(static_cast<float>(mx) * multiplier);
    }

    int16_t RFF2::getMouseYOnIterationBuffer(const int my) const {
        const float multiplier = settings.render.display.clarityMultiplier;
        return static_cast<int16_t>(static_cast<float>(getIterationBufferHeight()) -
                                    static_cast<float>(my) * multiplier);
    }

    void RFF2::checkBackupLoad() {
        const auto path = getBackupLocationPath();
        if (std::filesystem::exists(path)) {
            if (!std::filesystem::exists(path))
                return; // user deleted file manually
            loadLocation(path);
        }
    }

    void RFF2::recomputeThreaded() {

        state.createThread([this] {
            const double locateSeconds = pendingLocateSeconds.exchange(-1.0);
            Settings s = this->settings; // clone the settings
            const auto start = rootWindowContext->getWindow()->getTime();
            bool success = false;

            try {
                matchSettingsBeforeCreatingRenderData(s);
                success = prepareRenderData(start, s);
                if (success) matchSettingsAfterCreatingRenderData(s);

                if (success) {
                    beforeIterationFill(s);
                    success = fillIteration(start, s);
                    if (success && locateSeconds >= 0.0) {
                        setStatusMessage(Constants::Status::RENDER_STATUS,
                            std::format("Done | Locate: {:.3f}s", locateSeconds));
                    }
                }
            } catch (allocation_cancelled &) {
                vkh::logger::log("Memory allocation cancelled by user");
            } catch (const std::exception &e) {
                setStatusMessage(Constants::Status::RENDER_STATUS, e.what());
            }

            afterComputeFinally(success);
        });
    }

    void RFF2::moveCursorToCenter() const {
        const std::unique_ptr<fixed_point_complex_i1> off = MB2Locator::findCenterOffset(*renderData);
        auto offDex = static_cast<complex<dex>>(*off);
        if (renderData->getPerturbator()) {
            offDex -= renderData->getPerturbator()->off;
        }
        const auto [width, height] = rootWindowContext->getSwapchain().getSwapchainExtent();

        // multiplying 1.01 to attract reference center to client center
        const std::array<int, 2> ib = iterationBufferConversion(settings, offDex * dex(1.01));

        if (ib[0] >= 0 && ib[1] >= 0 && ib[0] < width && ib[1] < height) {
            moveCursor(ib[0], ib[1]);
        }
    }

    void RFF2::beforeIterationFill(Settings &s) const {
        if (settings.explore.autoMoveCursorToCenter) {
            moveCursorToCenter();
        }

        renderer->descriptorStorage->iteration->setMaxIteration(static_cast<double>(s.fractal.perturb.maxIteration));

        saveBackup();
    }

    void RFF2::matchSettingsBeforeCreatingRenderData(Settings &s) {
        if (s.render.computeShader.use) {
            s.fractal.reference.sync.referenceSynchronizationInterval = 1;
            s.fractal.mpa.useCompress = false;
            s.fractal.reference.compression.compressCriteria = 0;
        }
    }

    void RFF2::matchSettingsAfterCreatingRenderData(Settings &s) const {
        s.fractal = renderData->fractalSettings;
    }

    bool RFF2::prepareRenderData(const float startTime, Settings &s) {

        canShowPreview = false;

        if (state.interruptRequested())
            return false;

        auto &frt = s.fractal;
        const float logZoom = frt.general.logZoom;

        setStatusMessage(Constants::Status::ZOOM_STATUS,
                         std::format("Zoom : {:.06f}E{:d}", pow(10, fmod(logZoom, 1)), static_cast<int>(logZoom)));

        const complex<dex> offset = offsetConversion(s, 0, 0);
        const dex dcMax = offset.norm_approx();

        static uint64_t capacity = 0;
        if (renderData && renderData->getReference()) {
            capacity = renderData->getReference()->length();
        }

        std::function actionPerRefCalcIteration = FnExplore::getActionWhileRefCalc(*this, startTime);
        std::function actionPerSeriesApproxIteration = FnExplore::getActionWhileSeriesApprox(*this, startTime);
        std::function actionPerCreatingTableIteration = FnExplore::getActionWhileCreatingTable(*this, startTime);


        if (state.interruptRequested())
            return false;


        const int exp10 = Perturbator::logZoomToExp10(logZoom);
        if (frt.reference.reuse) {
            if (!renderData || !renderData->getReference() || !renderData->getPerturbator()) {
                vkh::logger::log_err("Do not reuse Reference during reference calculation!!!");
                this->settings.fractal.reference.reuse = false;
                requests.requestRecompute();
                return false;
            }

            fixed_point_complex_i1 center = frt.reference.center.create_variant(exp10);
            const fixed_point_complex_i1 referenceCenter = renderData->getReference()->center.create_variant(exp10);
            fixed_point_complex::sub(center, center, referenceCenter);
            const dex distance = static_cast<complex<dex>>(center).norm_approx();


            renderData->translate(frt.general.logZoom, dcMax + distance, frt.perturb, frt.reference.center,
                                  actionPerSeriesApproxIteration);
        } else {
            renderData = nullptr;


            renderData = createAppropriateRenderData(s.render.computeShader.use, logZoom, startTime, frt, dcMax, exp10,
                                                     capacity, 0);
        }

        // sync settings

        const MB2ReferenceBase *reference = renderData->getReference();
        if (!reference || state.interruptRequested())
            return false;

        size_t refLength = reference->length();
        size_t mpaLen = approxTableCache ? approxTableCache->tableSizeUsed : 0;

        setStatusMessage(Constants::Status::PERIOD_STATUS,
                         std::format("Period : {:L} ({:L}, {:L})", reference->longestPeriod(), refLength, mpaLen));
        if (state.interruptRequested())
            return false;

        return true;
    }

    void RFF2::fillIterationMultithreaded(const float startTime, const Settings &s) {
        std::atomic renderPixelsCount = 0;
        const uint16_t w = getIterationBufferWidth();
        const uint16_t h = getIterationBufferHeight();

        static std::vector<double> actualIterationMatrix(0);
        if (actualIterationMatrix.size() != w * h)
            actualIterationMatrix.resize(w * h);


        uint32_t len = static_cast<uint32_t>(w) * h;

        auto rendered = std::vector<uint8_t>(len);

        auto func = [&s, this, &renderPixelsCount, &rendered](const uint16_t x, const uint16_t y, const uint16_t xRes,
                                                              const uint16_t yRes, float, float, const uint32_t i,
                                                              double) {
            assert(i < rendered.size());
            rendered[i] = true;
            const auto dc = offsetConversion(s, x, y);
            const double iteration = renderData->getPerturbator()->iterate(dc);

            renderer->visibleIterationBufferContext->set(x, y, iteration);

            auto my = static_cast<int16_t>(y + 1);
            while (my < yRes && !rendered[my * xRes + x]) {
                renderer->visibleIterationBufferContext->set(x, my, iteration);
                ++my;
            }

            ++renderPixelsCount;
            return iteration;
        };
        const auto previewer =
                ParallelArrayDispatcher<double>(state, actualIterationMatrix, w, h, s.fractal.general.threads,
                                                s.render.display.pixelRenderPriority, std::move(func));


        auto statusThread = std::jthread([&renderPixelsCount, len, this, startTime](const std::stop_token &stop) {
            static float time = rootWindowContext->getWindow()->getTime();
            while (!stop.stop_requested()) {
                const float elapsed = rootWindowContext->getWindow()->getTime() - time;
                if (elapsed > Constants::Status::UI_REFRESH_INTERVAL) {
                    time = rootWindowContext->getWindow()->getTime();
                    float ratio = static_cast<float>(renderPixelsCount.load()) / static_cast<float>(len) * 100;
                    setStatusMessage(Constants::Status::TIME_STATUS,
                                     std::format("Time : {}", Utilities::formatTime(time - startTime)));
                    setStatusMessage(Constants::Status::RENDER_STATUS, std::format("Calculation : {:.3f}%", ratio));
                }
            }
        });


        canShowPreview = true;
        previewer.dispatch();

        statusThread.request_stop();
        statusThread.join();

        if (state.interruptRequested())
            return;

        const auto syncer = ParallelArrayDispatcher<double>(
                state, actualIterationMatrix, w, h, s.fractal.general.threads, RndPixelRenderPriority::SEQUENTIAL,
                [this](const uint16_t x, const uint16_t y, uint16_t, uint16_t, float, float, uint32_t, const double a) {
                    renderer->visibleIterationBufferContext->set(x, y, a);
                    return 0;
                });

        syncer.dispatch();
    }


    bool RFF2::fillIteration(const float startTime, const Settings &s) {

        if (state.interruptRequested())
            return false;

        const auto floatDat = dynamic_cast<FloatMB2RenderData *>(renderData.get());
        const auto fexDat = dynamic_cast<FexMB2RenderData *>(renderData.get());
        if ((floatDat || fexDat) && s.render.computeShader.use) {
            renderer->rg0->iterationPalette->setPerturbationMainIterator(PerturbationMainIterator::GPU);
            if (floatDat) {
                fillIterationComputeShader<float, fex>(startTime, s);
            } else {
                fillIterationComputeShader<fex, float>(startTime, s);
            }
        } else {
            renderer->rg0->iterationPalette->setPerturbationMainIterator(PerturbationMainIterator::CPU);
            renderer->computeIterateFex->clearMeta(*computeShaderManager->commandPool);
            renderer->computeIterateFloat->clearMeta(*computeShaderManager->commandPool);
            fillIterationMultithreaded(startTime, s);
        }

        if (state.interruptRequested()) {
            return false;
        }

        setStatusMessage(Constants::Status::RENDER_STATUS, "Done");
        return true;
    }

    void RFF2::afterComputeFinally(const bool success) {
        if (!success) {
            // vkh::logger::log("Recompute cancelled.");
        }
        auto expected = ComputeState::RUNNING;
        requests.recomputeRequestedState.compare_exchange_strong(expected, success ? ComputeState::IDLE : ComputeState::CANCELLED);
        backgroundThreads.notifyAll();
    }

    void RFF2::startBatchTest(const bool locate) {
        if (batchTest.active()) return;
        try {
            const std::filesystem::path root = R"(C:\software\RFF_Astra\Debug)";
            const auto folder = root / (locate ? "locate_minibrot_test/location" : "sha256_test");
            std::vector<std::filesystem::path> files;
            for (const auto &entry : std::filesystem::directory_iterator(folder)) {
                if (entry.is_regular_file() && entry.path().extension() == ".rfl") files.push_back(entry.path());
            }
            std::ranges::sort(files);
            if (files.empty()) throw std::runtime_error("No RFL files: " + folder.string());
            state.cancel();
            requests.recomputeRequestedState = ComputeState::IDLE;
            batchTest.savedSettings = settings;
            batchTest.savedExtent = rootWindowContext->getSwapchain().getSwapchainExtent();
            batchTest.files = std::move(files);
            batchTest.locate = locate;
            batchTest.index = 0;
            batchTest.cancel = false;
            batchTest.workerDone = false;
            const auto parent = locate ? folder.parent_path() : folder;
            for (unsigned i = 1;; ++i) {
                batchTest.output = parent / std::format("batch_output_{:04}", i);
                if (std::filesystem::create_directory(batchTest.output)) break;
            }
            batchTest.log.open(batchTest.output / (locate ? "locate_log.txt" : "sha256_log.txt"));
            if (!batchTest.log) throw std::runtime_error("Cannot create test log");
            batchTest.log.exceptions(std::ios::badbit | std::ios::failbit);
            batchTest.phase = BatchTest::Phase::SETUP;
            batchTest.frames = 0;
            batchTest.message = "Preparing test: " + batchTest.output.string();
            if (!locate) {
                settings = genDefaultSettings();
                settings.render.display = RndDisplayPresets::High().genDisplay();
                settings.shader.palette = ShdPalettePresets::Classic1().genPalette();
                settings.shader.palette.animationSpeed = 0;
                settings.fractal.reference.compression = ClcCompressPresets::None().genRefComp();
                settings.fractal.mpa.useCompress = false;
                settings.explore.autoMoveCursorToCenter = false;
                rootWindowContext->getWindow()->setResolution(1280, 720);
            }
            settings.fractal.reference.reuse = false;
            zoomAnimationInfo.reset();
            requests.requestShader();
            requests.requestResize(rootWindowContext->getSwapchain().getSwapchainExtent());
            if (locate) {
                batchTest.log << "Locate Minibrot; steady_clock seconds; preparation/render excluded\n"
                              << "Coordinates are saved as RFL. Distance is log10(|result - suppose|).\n"
                              << "Proximity is a comparison aid, not a proof of minibrot identity.\n";
            } else {
                batchTest.log << "index\tsource\treuse\tlogZoom\tmaxIteration\trfm_sha256\tpng_sha256\n";
            }
            batchTest.log.flush();
            std::ofstream info(batchTest.output / "run_info.txt");
            info.exceptions(std::ios::badbit | std::ios::failbit);
            const auto &gpu = rootWindowContext->core.getPhysicalDeviceLoader().getPhysicalDeviceProperties();
            info << "GPU=" << gpu.deviceName << "\ndriver=" << gpu.driverVersion
                 << "\nthreads=" << settings.fractal.general.threads
                 << "\ncomputeShader=" << settings.render.computeShader.use
                 << "\ninput=" << folder.string() << '\n';
            if (!locate) info << "resolution=1280x720\npalette=Classic 1\nlinearInterpolation=off\nreferenceCompression=off\nreferenceReuse=alternating 0,1\n";
            info.close();
        } catch (const std::exception &e) {
            if (batchTest.active()) finishBatchTest(std::string("ERROR: ") + e.what());
            else batchTest.message = std::string("ERROR: ") + e.what();
        }
    }

    void RFF2::finishBatchTest(const std::string &message) {
        state.cancel();
        requests.recomputeRequestedState = ComputeState::IDLE;
        try {
            if (batchTest.log.is_open()) {
                batchTest.log << message << '\n';
                batchTest.log.close();
            }
        } catch (...) { batchTest.log.exceptions(std::ios::goodbit); batchTest.log.close(); }
        batchTest.log.clear();
        batchTest.phase = BatchTest::Phase::IDLE;
        batchTest.message = message + " " + batchTest.output.string();
        settings = *batchTest.savedSettings;
        // The previous reference was replaced by the test; the restored view needs a fresh one.
        settings.fractal.reference.reuse = false;
        rootWindowContext->getWindow()->setResolution(batchTest.savedExtent.width, batchTest.savedExtent.height);
        requests.requestShader();
        requests.requestResize(rootWindowContext->getSwapchain().getSwapchainExtent());
        requests.requestRecompute();
    }

    void RFF2::updateBatchTest() {
        using Phase = BatchTest::Phase;
        if (!batchTest.active()) return;
        try {
            if (batchTest.cancel) { finishBatchTest("cancelled."); return; }
            if (batchTest.phase == Phase::SETUP) {
                // Let GLFW resize events and the resulting render settle before loading entry one.
                if (requests.resizeRequested || requests.recomputeRequestedState == ComputeState::RUNNING ||
                    requests.recomputeRequestedState == ComputeState::REQUESTED) return;
                if (!batchTest.locate && (getIterationBufferWidth() != 1280 || getIterationBufferHeight() != 720))
                    throw std::runtime_error("Fixed 1280x720 test resolution is unavailable");
                batchTest.phase = Phase::NEXT;
            }
            if (batchTest.phase == Phase::NEXT) {
                if (batchTest.index == batchTest.files.size()) { finishBatchTest("done."); return; }
                state.cancel();
                const auto &path = batchTest.files[batchTest.index];
                const auto location = RFFLocationBinary::read(path);
                if (!location.hasData()) throw std::runtime_error("Cannot read: " + path.string());
                settings.fractal.reference.reuse = !batchTest.locate && batchTest.index % 2 == 1;
                if (!settings.fractal.reference.reuse) {
                    renderData.reset();
                    approxTableCache.reset();
                }
                batchTest.message = std::format("{}/{}: {}", batchTest.index + 1, batchTest.files.size(), path.filename().string());
                loadLocation(path);
                batchTest.phase = Phase::RENDER;
                return;
            }
            if (batchTest.phase == Phase::RENDER) {
                const auto compute = requests.recomputeRequestedState.load();
                if (compute == ComputeState::CANCELLED) throw std::runtime_error("Render cancelled or failed");
                if (compute != ComputeState::IDLE) return;
                state.cancel(); // Join before reading renderData or starting the locator.
                if (batchTest.locate) {
                    batchTest.phase = Phase::LOCATE;
                    batchTest.workerDone = false;
                    batchTest.workerSuccess = false;
                    batchTest.workerError.clear();
                    state.createThread([this] {
                        try {
                            const auto start = rootWindowContext->getWindow()->getTime();
                            const auto timer = std::chrono::steady_clock::now();
                            auto locator = MB2Locator::locateMinibrot(state, *renderData, approxTableCache,
                                FnExplore::getActionWhileFindingMBCenter(*this, renderData->getReference()->longestPeriod(), start),
                                FnExplore::getActionWhileSeriesApprox(*this, start),
                                FnExplore::getActionWhileCreatingTable(*this, start),
                                FnExplore::getActionWhileFindingZoom(*this, start));
                            batchTest.locateSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - timer).count();
                            if (locator && !state.interruptRequested()) {
                                settings.fractal.reference.center = locator->settings.reference.center;
                                settings.fractal.general.logZoom = locator->settings.general.logZoom - MB2Locator::MINIBROT_LOG_ZOOM_OFFSET;
                                batchTest.workerSuccess = true;
                            }
                        } catch (const std::exception &e) { batchTest.workerError = e.what(); }
                          catch (...) { batchTest.workerError = "Locator failed"; }
                        batchTest.workerDone = true;
                    });
                } else {
                    batchTest.phase = Phase::CAPTURE;
                    batchTest.frames = rootWindowContext->core.getPhysicalDeviceLoader().getMaxFramesInFlight() + 2;
                }
                return;
            }
            if (batchTest.phase == Phase::LOCATE) {
                if (!batchTest.workerDone) return;
                state.cancel();
                if (!batchTest.workerSuccess) throw std::runtime_error("Locate failed/cancelled: " + batchTest.workerError);
                const auto &file = batchTest.files[batchTest.index];
                saveCurrentLocation(batchTest.output / file.filename());
                // Ensure export succeeded even though the legacy RFL writer reports errors only via logger.
                (void)SHA256::file(batchTest.output / file.filename());
                auto &center = settings.fractal.reference.center;
                batchTest.log << std::format("\n{}\tseconds={:.9f}\tlogZoom={:.9g}\n", file.filename().string(), batchTest.locateSeconds, settings.fractal.general.logZoom)
                              << "result.real=" << center.real.to_string() << '\n'
                              << "result.imag=" << center.imag.to_string() << '\n';
                const auto expectedPath = file.parent_path().parent_path() / "suppose" / file.filename();
                const auto expected = RFFLocationBinary::read(expectedPath);
                if (expected.hasData()) {
                    const auto exp10 = Perturbator::logZoomToExp10(std::max(settings.fractal.general.logZoom, expected.getLogZoom()));
                    auto delta = center.create_variant(exp10);
                    const fixed_point_complex_i1 expectedCenter(expected.getReal(), expected.getImag(), exp10);
                    fixed_point_complex::sub(delta, delta, expectedCenter);
                    const auto distance = static_cast<complex<dex>>(delta).norm_approx();
                    batchTest.log << "suppose.real=" << expected.getReal() << '\n'
                                  << "suppose.imag=" << expected.getImag() << '\n'
                                  << std::format("suppose.logZoom={:.9g}\ndelta.logZoom={:.9g}\nlog10.distance={}\n",
                                      expected.getLogZoom(), settings.fractal.general.logZoom - expected.getLogZoom(),
                                      distance == dex(0) ? "-inf (exact match)" : std::format("{}", rff_math::log10(distance)));
                } else batchTest.log << "suppose=MISSING_OR_UNREADABLE\n";
                setStatusMessage(Constants::Status::TIME_STATUS, std::format("Time : {:.6f}s", batchTest.locateSeconds));
                batchTest.log.flush();
                ++batchTest.index;
                batchTest.phase = Phase::NEXT;
                return;
            }
            if (batchTest.phase == Phase::CAPTURE) {
                if (requests.resizeRequested || requests.recomputeRequestedState != ComputeState::IDLE)
                    throw std::runtime_error("Render or resolution changed during capture");
                if (batchTest.frames-- > 0) return;
                const auto stem = std::format("{:04}", batchTest.index + 1);
                const auto mapPath = batchTest.output / (stem + ".rfm");
                const auto pngPath = batchTest.output / (stem + ".png");
                generateMap().exportFile(mapPath);
                requests.createImageRequestedFilename = pngPath.string();
                applyCreateImage();
                const auto &frt = renderData->fractalSettings;
                batchTest.log << std::format("{}\t{}\t{}\t{:.9g}\t{}\t{}\t{}\n", batchTest.index + 1,
                    batchTest.files[batchTest.index].filename().string(), batchTest.index % 2,
                    frt.general.logZoom, frt.perturb.maxIteration, SHA256::file(mapPath), SHA256::file(pngPath));
                batchTest.log.flush();
                ++batchTest.index;
                batchTest.phase = Phase::NEXT;
            }
        } catch (const std::exception &e) { finishBatchTest(std::string("ERROR: ") + e.what()); }
    }

} // namespace merutilm::rff2
