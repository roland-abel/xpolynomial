/// @file utils.h
/// @brief Provides some numeric helper functions.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#ifndef UTILS_H_
#define UTILS_H_

#include <concepts>
#include <ranges>
#include <limits>
#include <complex>
#include <cmath>
#include <vector>

namespace xmath {

    /// @brief Rounds a floating-point value to the nearest integer, with halves away from zero.
    /// @tparam T A floating-point type.
    /// @param x The value to round.
    /// @return The rounded value.
    template<std::floating_point T>
    constexpr T round_half_away_from_zero(T x) {
        if (x != x) {
            return x;
        }

        if (x < 0) {
            return -round_half_away_from_zero(-x);
        }

        const auto half_rounded_up = x + static_cast<T>(0.5);

        if (half_rounded_up >= static_cast<T>(std::numeric_limits<long long>::max())) {
            return x;
        }

        return static_cast<T>(static_cast<long long>(half_rounded_up));
    }

    /// @brief Checks if a value is nearly zero within a specified epsilon.
    /// @tparam T The data type of the value.
    /// @tparam FP The floating-point precision type.
    /// @param a The value to be checked.
    /// @param epsilon The epsilon value for comparison (default is machine epsilon).
    /// @return True if the value is nearly zero; otherwise, false.
    template<typename T, typename FP = T>
    constexpr bool nearly_zero(T a, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        if constexpr (std::floating_point<T>) {
            return (a < 0 ? -a : a) < epsilon;
        } else {
            return std::abs(a) < epsilon;
        }
    }

    /// @brief Checks if two values are nearly equal within a specified epsilon.
    /// @tparam T The data type of the values.
    /// @tparam FP The floating-point precision type.
    /// @param a The first value.
    /// @param b The second value.
    /// @param epsilon The epsilon value for comparison (default is machine epsilon).
    /// @return True if the values are nearly equal; otherwise, false.
    template<typename T, typename FP = T>
    constexpr bool nearly_equal(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return nearly_zero(a - b, epsilon);
    }

    /// @brief Checks if a value is greater than another within a specified epsilon.
    /// @tparam T The data type of the values.
    /// @tparam FP The floating-point precision type.
    /// @param a The first value.
    /// @param b The second value.
    /// @param epsilon The epsilon value for comparison (default is machine epsilon).
    /// @return True if the first value is greater than the second; otherwise, false.
    template<typename T, typename FP = T>
    constexpr bool greater_than(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return a - epsilon > b;
    }

    /// @brief Checks if a value is greater than or equal to another within a specified epsilon.
    /// @tparam T The data type of the values.
    /// @tparam FP The floating-point precision type.
    /// @param a The first value.
    /// @param b The second value.
    /// @param epsilon The epsilon value for comparison (default is machine epsilon).
    /// @return True if the first value is greater than or equal to the second; otherwise, false.
    template<typename T, typename FP = T>
    constexpr bool greater_than_or_equal(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return a + epsilon > b;
    }

    /// @brief Checks if a value is less than another within a specified epsilon.
    /// @tparam T The data type of the values.
    /// @tparam FP The floating-point precision type.
    /// @param a The first value.
    /// @param b The second value.
    /// @param epsilon The epsilon value for comparison (default is machine epsilon).
    /// @return True if the first value is less than the second; otherwise, false.
    template<typename T, typename FP = T>
    constexpr bool less_than(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return a + epsilon < b;
    }

    /// @brief Checks if a value is less than or equal to another within a specified epsilon.
    /// @tparam T The data type of the values.
    /// @tparam FP The floating-point precision type.
    /// @param a The first value.
    /// @param b The second value.
    /// @param epsilon The epsilon value for comparison (default is machine epsilon).
    /// @return True if the first value is less than or equal to the second; otherwise, false.
    template<typename T, typename FP = T>
    constexpr bool less_than_or_equal(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return a - epsilon < b;
    }

    /// @brief Checks if a given integer is even.
    /// @tparam T An integer type.
    /// @param a The integer value to check.
    /// @return True if the integer is even; otherwise, false.
    template<std::integral T>
    constexpr bool is_even(T a) {
        return a % 2 == 0;
    }

    /// @brief Checks if a given integer is odd.
    /// @tparam T An integer type.
    /// @param a The integer value to check.
    /// @return True if the integer is odd; otherwise, false.
    template<std::integral T>
    constexpr bool is_odd(T a) {
        return a % 2 != 0;
    }

    /// @brief Converts a range to a vector.
    /// @tparam R The type of the input range.
    /// @param r The input range.
    /// @return A vector containing the elements of the input range.
    template<std::ranges::range R>
    constexpr auto to_vector(R &&r) {
        auto r_common = r | std::views::common;
        return std::vector(r_common.begin(), r_common.end());
    }
}

#endif // UTILS_H_