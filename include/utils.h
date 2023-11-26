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