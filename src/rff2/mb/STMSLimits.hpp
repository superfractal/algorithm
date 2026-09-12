// Created by GPT-6 on 2026-09-10
// Modified by GPT-6 on 2026-09-12
#pragma once
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
namespace stms_limits {
    // The kernels use signed int bit counts. Reserve headroom for intermediate
    // counts and exponent arithmetic instead of imposing a decimal zoom ceiling.
    inline constexpr unsigned maxVerificationBits = std::numeric_limits<int>::max() / 32;
    inline const double maxInputLogZoom =
        std::floor((maxVerificationBits / std::log2(10.) - 160) / 2);
    inline unsigned bitsForZoom(double zoom) {
        if (!std::isfinite(zoom) || zoom < 0 || zoom > maxInputLogZoom)
            throw std::runtime_error("STMS requested precision exceeds kernel integer capacity");
        return static_cast<unsigned>(std::ceil((2 * zoom + 160) * std::log2(10.)));
    }
    inline int64_t exponentLimit(int64_t bits) {
        return std::max<int64_t>(1000000, 2 * bits + 1024);
    }
} // namespace stms_limits
