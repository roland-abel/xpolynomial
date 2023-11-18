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
        /// @brief Creates the closed interval [0., 1.].
        real_interval()
                : real_interval(spec::zero, spec::one) {
        }

        /// Constructor that creates a closed/opened interval with the boundaries `a` and `b`.
        /// @param a The lower endpoint of the interval.
        /// @param a The upper endpoint of the interval.
        /// @param lower The lower boundary closed or opened.
        /// @param upper The lower boundary closed or opened.
        real_interval(T a, T b, interval_bounds lower = closed, interval_bounds upper = closed)
                : std::pair<T, T>(a, b) {
            lower_ = lower;
            upper_ = upper;
        }

        /// @brief Gets the lower point of the interval.
        /// @return The lower point.
        inline T lower() const { return std::pair<T, T>::first; }

        /// @brief Gets the end point of the interval.
        /// @return The end point.
        inline T upper() const { return std::pair<T, T>::second; }

        /// @brief Gets the length of the interval.
        /// @return The length.
        inline T length() const { return upper() - lower(); }

        /// @brief Gets a value indicated whether the interval is empty.
        /// @return True if the interval ist empty, otherwise false.
        inline bool is_empty() const {
            return greater_than(lower(), upper()) || (is_degenerate() && !is_closed());
        }

        /// @brief Gets a values indicated whether the interval is closed.
        /// @return True if the interval is closed; otherwise false.
        inline bool is_closed() const { return lower_ == closed && upper_ == closed; }

        /// @brief Gets a values indicated whether the interval is opened.
        /// @return True if the interval is opened; otherwise false.
        inline bool is_opened() const { return lower_ == opened && upper_ == opened; }

        /// @brief Gets a values indicated whether the lower boundary of the interval is open.
        /// @return True if the lower boundary is opened; otherwise false.
        inline bool is_lower_open() const { return lower_ == opened; }

        /// @brief Gets a values indicated whether the upper boundary of the interval is open.
        /// @return True if the upper boundary is opened; otherwise false.
        inline bool is_upper_open() const { return upper_ == opened; }

        /// @brief Gets a values indicated whether the lower boundary of the interval is closed.
        /// @return True if the lower boundary is closed; otherwise false.
        inline bool is_lower_closed() const { return lower_ == closed; }

        /// @brief Gets a values indicated whether the upper boundary of the interval is closed.
        /// @return True if the upper boundary is closed; otherwise false.
        inline bool is_upper_closed() const { return upper_ == closed; }

        /// @brief Gets a values indicated whether the interval is half open.
        /// @return True if the interval is half open; otherwise false.
        inline bool is_half_open() const { return lower_ != upper_; }

        /// @brief Gets a values indicated whether the interval is degenerate, e.g. the intervals boundary are equals.
        /// @return True if the interval is degenerate; otherwise false.
        inline bool is_degenerate() const { return nearly_equal(lower(), upper(), epsilon); }

        /// Gets a tuple of two intervals created by the current interval by bisection.
        /// @return A tuple of two intervals.
        std::pair<real_interval<T>, real_interval<T>> bisect() const {
            return std::make_pair(
                    real_interval(lower(), (lower() + upper()) / static_cast<T>(2.)),
                    real_interval((lower() + upper()) / static_cast<T>(2.), upper()));
        }

    private:
        interval_bounds lower_;
        interval_bounds upper_;
    };
}

#endif // REAL_INTERVAL_H_
