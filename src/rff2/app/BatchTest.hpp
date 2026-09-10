// created by GPT-6 on 2026-09-07
#pragma once
#include <atomic>
#include <filesystem>
#include <fstream>
#include <optional>
#include <vector>
#include "../settings/Settings.h"

namespace merutilm::rff2 {
    struct BatchTest {
        enum class Phase { IDLE, SETUP, NEXT, RENDER, LOCATE, CAPTURE };
        Phase phase = Phase::IDLE;
        bool locate = false;
        bool cancel = false;
        size_t index = 0;
        unsigned frames = 0;
        std::vector<std::filesystem::path> files;
        std::filesystem::path output;
        std::ofstream log;
        std::optional<Settings> savedSettings;
        VkExtent2D savedExtent{};
        std::atomic<bool> workerDone = false;
        bool workerSuccess = false;
        double locateSeconds = 0;
        std::string workerError;
        std::string message;
        [[nodiscard]] bool active() const { return phase != Phase::IDLE; }
    };
}
