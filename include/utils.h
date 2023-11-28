/// @file utils.h
/// @brief
///
/// @author Roland Abel
/// @date October 8, 2023
///
/// Copyright (c) 2023 Roland Abel
///
/// This software is released under the MIT License.
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
/// THE SOFTWARE.

#ifndef UTILS_H_
#define UTILS_H_

#include <ranges>
#include <limits>
#include <cmath>

namespace xmath {

    /// @brief Checks if a value is nearly zero within a specified epsilon.
    /// @tparam T The data type of the value.
    /// @tparam FP The floating-point precision type.
    /// @param a The value to be checked.
    /// @param epsilon The epsilon value for comparison (default is machine epsilon).
    /// @return True if the value is nearly zero; otherwise, false.
    template<typename T, typename FP = T>
    inline bool nearly_zero(T a, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return std::abs(a) < epsilon;
    }

    /// @brief Checks if two values are nearly equal within a specified epsilon.
    /// @tparam T The data type of the values.
    /// @tparam FP The floating-point precision type.
    /// @param a The first value.
    /// @param b The second value.
    /// @param epsilon The epsilon value for comparison (default is machine epsilon).
    /// @return True if the values are nearly equal; otherwise, false.
    template<typename T, typename FP = T>
    inline bool nearly_equal(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
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
    inline bool greater_than(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
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
    inline bool greater_than_or_equal(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
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
    inline bool less_than(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
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
    inline bool less_than_or_equal(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return a - epsilon < b;
    }

    /// @brief Checks if a given integer is even.
    /// @param a The integer value to check.
    /// @return True if the integer is even; otherwise, false.
    inline bool is_even(long a) {
        return a % 2 == 0;
    }

    /// @brief Checks if a given integer is odd.
    /// @param a The integer value to check.
    /// @return True if the integer is odd; otherwise, false.
    inline bool is_odd(long a) {
        return a % 2 != 0;
    }

    /// @brief Converts a range to a vector.
    /// @tparam R The type of the input range.
    /// @param r The input range.
    /// @return A vector containing the elements of the input range.
    template<std::ranges::range R>
    auto to_vector(R &&r) {
        auto r_common = r | std::views::common;
        return std::vector(r_common.begin(), r_common.end());
    }
}

#endif // UTILS_H_