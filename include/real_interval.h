/// @file real_interval.h
///
/// @author Roland Abel
/// @date 07.11.2023

#ifndef REAL_INTERVAL_H_
#define REAL_INTERVAL_H_

#include <limits>
#include <cmath>
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

    template<typename T, typename FT = floating_point_specification<T>>
    class real_interval : std::pair<T, T> {
    public:
        using spec = floating_point_specification<T>;
        using floating_point_type = spec::floating_point_type;
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
        /// Creates the closed interval [0., 1.].
        real_interval()
                : real_interval(spec::zero, spec::one) {
        }

        /// Constructor that creates a closed/opened interval with the boundaries `a` and `b`.
        /// @param a The left endpoint of the interval.
        /// @param a The right endpoint of the interval.
        /// @param left The left boundary closed or opened.
        /// @param right The left boundary closed or opened.
        real_interval(T a, T b, interval_bounds left = closed, interval_bounds right = closed)
                : std::pair<T, T>(a, b) {
            left_ = left;
            right_ = right;
        }

        /// Gets the start point of the interval.
        /// @return The start point.
        inline T start() const { return std::pair<T, T>::first; }

        /// Gets the end point of the interval.
        /// @return The end point.
        inline T end() const { return std::pair<T, T>::second; }

        /// Gets the length of the interval.
        /// @return The length.
        inline T length() const { return end() - start(); }

        /// Gets a value indicated whether the interval is empty.
        /// @return True if the interval ist empty, otherwise false.
        inline bool is_empty() const {
            return greater_than(start(), end()) || (is_degenerate() && !is_closed());
        }

        /// Gets a values indicated whether the interval is closed.
        /// @return True if the interval is closed; otherwise false.
        inline bool is_closed() const { return left_ == closed && right_ == closed; }

        /// Gets a values indicated whether the interval is opened.
        /// @return True if the interval is opened; otherwise false.
        inline bool is_opened() const { return left_ == opened && right_ == opened; }

        inline bool is_left_open() const { return left_ == opened; }

        inline bool is_right_open() const { return right_ == opened; }

        inline bool is_left_closed() const { return left_ == closed; }

        inline bool is_right_closed() const { return right_ == closed; }

        inline bool is_half_open() const { return left_ != right_; }

        inline bool is_degenerate() const { return nearly_equal(start(), end(), epsilon); }

        /// Gets a tuple of two intervals created by the current interval by bisection.
        /// @return A tuple of two intervals.
        std::pair<real_interval<T>, real_interval<T>> bisect() const {
            return std::make_pair(
                    real_interval(start(), (start() + end()) / static_cast<T>(2.)),
                    real_interval((start() + end()) / static_cast<T>(2.), end()));
        }

    private:
        interval_bounds left_;
        interval_bounds right_;
    };
}

#endif // REAL_INTERVAL_H_
