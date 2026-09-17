/// @file interval.h
/// @brief Interval over real numbers.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#ifndef INTERVAL_H_
#define INTERVAL_H_

#include <cmath>
#include <functional>
#include "utils.h"

namespace xmath {

    /// @brief Specifies the floating-point type and its constants.
    /// @tparam T The floating-point type.
    template<typename T>
    struct floating_point_specification {
        using floating_point_type = T;
        static constexpr floating_point_type one = 1.0;
        static constexpr floating_point_type zero = 0.0;
    };

    /// @brief Floating-point specification for `double`.
    template<>
    struct floating_point_specification<double> {
        using floating_point_type = double;
        static constexpr floating_point_type one = 1.0;
        static constexpr floating_point_type zero = 0.0;
        static constexpr floating_point_type epsilon = 1e-9;
    };

    /// @brief Floating-point specification for `float`.
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
    class interval {
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
        constexpr interval()
                : interval(spec::zero, spec::one) {
        }

        /// Constructor that creates a closed/opened interval with the boundaries `lower` and `upper`.
        /// @param lower The lower endpoint of the interval.
        /// @param upper The upper endpoint of the interval.
        /// @param lower_bounds The lower boundary closed or opened (default opened).
        /// @param upper_bounds The upper boundary closed or opened (default closed).
        constexpr interval(T lower, T upper,
                           interval_bounds lower_bounds = opened, interval_bounds upper_bounds = closed)
                : lower_(lower), upper_(upper), lower_bounds_(lower_bounds), upper_bounds_(upper_bounds) {
        }

        /// @brief Gets the lower point of the interval.
        /// @return The lower point.
        constexpr T lower() const noexcept { return lower_; }

        /// @brief Gets the end point of the interval.
        /// @return The end point.
        constexpr T upper() const noexcept { return upper_; }

        /// @brief Gets the length of the interval.
        /// @return The length.
        constexpr T length() const noexcept { return upper() - lower(); }

        /// @brief Gets a value indicated whether the interval is empty.
        /// @return True if the interval is empty, otherwise false.
        [[nodiscard]] constexpr bool is_empty() const noexcept {
            return greater_than(lower(), upper()) || (is_degenerate() && !is_closed());
        }

        /// @brief Gets a value indicated whether the interval is closed.
        /// @return True if the interval is closed; otherwise false.
        [[nodiscard]] constexpr bool is_closed() const noexcept { return lower_bounds_ == closed && upper_bounds_ == closed; }

        /// @brief Gets a value indicated whether the interval is opened.
        /// @return True if the interval is opened; otherwise false.
        [[nodiscard]] constexpr bool is_opened() const noexcept { return lower_bounds_ == opened && upper_bounds_ == opened; }

        /// @brief Gets a value indicated whether the lower boundary of the interval is open.
        /// @return True if the lower boundary is opened; otherwise false.
        [[nodiscard]] constexpr bool is_lower_open() const noexcept { return lower_bounds_ == opened; }

        /// @brief Gets a value indicated whether the upper boundary of the interval is open.
        /// @return True if the upper boundary is opened; otherwise false.
        [[nodiscard]] constexpr bool is_upper_open() const noexcept { return upper_bounds_ == opened; }

        /// @brief Gets a value indicated whether the lower boundary of the interval is closed.
        /// @return True if the lower boundary is closed; otherwise false.
        [[nodiscard]] constexpr bool is_lower_closed() const noexcept { return lower_bounds_ == closed; }

        /// @brief Gets a value indicated whether the upper boundary of the interval is closed.
        /// @return True if the upper boundary is closed; otherwise false.
        [[nodiscard]] constexpr bool is_upper_closed() const noexcept { return upper_bounds_ == closed; }

        /// @brief Gets a value indicated whether the interval is half open.
        /// @return True if the interval is half open; otherwise false.
        [[nodiscard]] constexpr bool is_half_open() const noexcept { return lower_bounds_ != upper_bounds_; }

        /// @brief Gets a value indicated whether the interval is degenerate, e.g. the intervals boundary are equals.
        /// @return True if the interval is degenerate; otherwise false.
        [[nodiscard]] constexpr bool is_degenerate() const noexcept { return nearly_equal(lower(), upper(), epsilon); }

        /// Gets a tuple of two intervals created by the current interval by bisection.
        /// @return A tuple of two intervals.
        constexpr std::pair<interval<T>, interval<T>>
        bisect(interval_bounds lower_bounds = opened, interval_bounds upper_bounds = closed) const noexcept {
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
        /// The lower endpoint of the interval.
        T lower_;
        /// The upper endpoint of the interval.
        T upper_;
        /// The lower boundary type of the interval.
        interval_bounds lower_bounds_;
        /// The upper boundary type of the interval.
        interval_bounds upper_bounds_;
    };
}

#endif // INTERVAL_H_
