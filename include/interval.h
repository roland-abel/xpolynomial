/// @file interval.h
/// @brief Interval over real numbers.
///
/// @author Roland Abel
/// @date November 7, 2023
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

#pragma once

#include <cmath>
#include <functional>
#include "utils.h"

namespace xmath {

    template<typename T>
    struct floating_point_specification {
        using floating_point_type = T;
        static constexpr floating_point_type one = 1.0;
        static constexpr floating_point_type zero = 0.0;
    };

    template<>
    struct floating_point_specification<double> {
        using floating_point_type = double;
        static constexpr floating_point_type one = 1.0;
        static constexpr floating_point_type zero = 0.0;
        static constexpr floating_point_type epsilon = 1e-9;
    };

    template<>
    struct floating_point_specification<float> {
        using floating_point_type = float;
        static constexpr floating_point_type one = 1.0f;
        static constexpr floating_point_type zero = 0.0f;
        static constexpr floating_point_type epsilon = 1e-5f;
    };

    /// @brief Represents a real interval with specified boundary conditions.
    /// @tparam T The data type of the interval endpoints.
    /// @tparam FT The floating-point specification type.
    template<typename T, typename FT = floating_point_specification<T>>
    class interval : std::pair<T, T> {
    public:
        using spec = floating_point_specification<T>;
        using floating_point_type = typename spec::floating_point_type;
        static constexpr floating_point_type epsilon = spec::epsilon;

        /// Indicates whether the interval boundaries are open or closed.
        enum interval_bounds {
            /// Indicated that the endpoint is
            /// to be excluded the interval.
            opened,

            /// Indicated that the endpoint is
            /// to be included the interval.
            closed
        };

    public:
        /// @brief Creates the closed interval [0., 1.].
        interval()
                : interval(spec::zero, spec::one) {
        }

        /// Constructor that creates a closed/opened interval with the boundaries `a` and `b`.
        /// @param a The lower endpoint of the interval.
        /// @param b The upper endpoint of the interval.
        /// @param lower The lower boundary closed or opened (default opened).
        /// @param upper The lower boundary closed or opened (default closed).
        interval(T a, T b, interval_bounds lower = opened, interval_bounds upper = closed)
                : std::pair<T, T>(a, b) {
            lower_ = lower;
            upper_ = upper;
        }

        /// @brief Gets the lower point of the interval.
        /// @return The lower point.
        T lower() const { return std::pair<T, T>::first; }

        /// @brief Gets the end point of the interval.
        /// @return The end point.
        T upper() const { return std::pair<T, T>::second; }

        /// @brief Gets the length of the interval.
        /// @return The length.
        T length() const { return upper() - lower(); }

        /// @brief Gets a value indicated whether the interval is empty.
        /// @return True if the interval ist empty, otherwise false.
        [[nodiscard]] bool is_empty() const {
            return greater_than(lower(), upper()) || (is_degenerate() && !is_closed());
        }

        /// @brief Gets a values indicated whether the interval is closed.
        /// @return True if the interval is closed; otherwise false.
        [[nodiscard]] bool is_closed() const { return lower_ == closed && upper_ == closed; }

        /// @brief Gets a values indicated whether the interval is opened.
        /// @return True if the interval is opened; otherwise false.
        [[nodiscard]] bool is_opened() const { return lower_ == opened && upper_ == opened; }

        /// @brief Gets a values indicated whether the lower boundary of the interval is open.
        /// @return True if the lower boundary is opened; otherwise false.
        [[nodiscard]] bool is_lower_open() const { return lower_ == opened; }

        /// @brief Gets a values indicated whether the upper boundary of the interval is open.
        /// @return True if the upper boundary is opened; otherwise false.
        [[nodiscard]] bool is_upper_open() const { return upper_ == opened; }

        /// @brief Gets a values indicated whether the lower boundary of the interval is closed.
        /// @return True if the lower boundary is closed; otherwise false.
        [[nodiscard]] bool is_lower_closed() const { return lower_ == closed; }

        /// @brief Gets a values indicated whether the upper boundary of the interval is closed.
        /// @return True if the upper boundary is closed; otherwise false.
        [[nodiscard]] bool is_upper_closed() const { return upper_ == closed; }

        /// @brief Gets a values indicated whether the interval is half open.
        /// @return True if the interval is half open; otherwise false.
        [[nodiscard]] bool is_half_open() const { return lower_ != upper_; }

        /// @brief Gets a values indicated whether the interval is degenerate, e.g. the intervals boundary are equals.
        /// @return True if the interval is degenerate; otherwise false.
        [[nodiscard]] bool is_degenerate() const { return nearly_equal(lower(), upper(), epsilon); }

        /// Gets a tuple of two intervals created by the current interval by bisection.
        /// @return A tuple of two intervals.
        std::pair<interval<T>, interval<T>>
        bisect(interval_bounds lower_bounds = opened, interval_bounds upper_bounds = closed) const {
            const auto c = (lower() + upper()) / static_cast<T>(2.);
            return std::make_pair(
                    interval(lower(), c, lower_bounds, upper_bounds),
                    interval(c, upper(), lower_bounds, upper_bounds));
        }

        /// Gets a function that map linear the interval to the given interval.
        /// @param I The interval.
        /// @return The linear mapping function.
        std::function<T(const T &x)> linear_transform(const interval<T> &I) const {
            const auto a = lower(), b = upper(), alpha = I.lower(), beta = I.upper();
            const auto m = (beta - alpha) / (b - a), c = (alpha * b - beta * a) / (b - a);

            auto map = [m, c](const T &t) {
                return m * t + c;
            };
            return map;
        }

    private:
        interval_bounds lower_;
        interval_bounds upper_;
    };
}
