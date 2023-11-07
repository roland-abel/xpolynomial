/// @file chebyshev_polynomial.h
///
/// @author Roland Abel
/// @date 08.10.2023

#ifndef UTILS_H_
#define UTILS_H_

#include <ranges>
#include <limits>
#include <cmath>

namespace xmath {

    template<typename T, typename FP = T>
    inline bool nearly_zero(T a, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return std::abs(a) < epsilon;
    }

    template<typename T, typename FP = T>
    inline bool nearly_equal(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return nearly_zero(a - b, epsilon);
    }

    template<typename T, typename FP = T>
    inline bool greater_than(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return a - epsilon > b;
    }

    template<typename T, typename FP = T>
    inline bool greater_than_or_equal(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return a + epsilon > b;
    }

    template<typename T, typename FP = T>
    inline bool less_than(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return a + epsilon < b;
    }

    template<typename T, typename FP = T>
    inline bool less_than_or_equal(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return a - epsilon < b;
    }

    inline bool is_even(long a) {
        return a % 2 == 0;
    }

    inline bool is_odd(long a) {
        return a % 2 == 1;
    }

    /// @brief
    /// @tparam R
    /// @param r
    /// @return
    template<std::ranges::range R>
    auto to_vector(R &&r) {
        auto r_common = r | std::views::common;
        return std::vector(r_common.begin(), r_common.end());
    }
}

#endif // UTILS_H_