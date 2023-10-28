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

    template<typename T>
    class interval : std::pair<T, T> {
    public:
        interval(T a, T b) : std::pair<T, T>(a, b) {
        }

        inline T start() const { return std::pair<T, T>::first; }

        inline T end() const { return std::pair<T, T>::second; }

        inline T length() const { return end() - start(); }

        inline bool is_empty() const { return start() > end(); }

        std::pair<interval<T>, interval<T>> bisect() const {
            return std::make_pair(
                    interval(start(), (start() + end()) / 2.),
                    interval((start() + end()) / 2., end()));
        }
    };

    template<typename T, typename FP = T>
    inline bool nearly_zero(T a, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return std::abs(a) < epsilon;
    }

    template<typename T, typename FP = T>
    inline bool nearly_equal(T a, T b, FP epsilon = std::numeric_limits<FP>::epsilon()) {
        return nearly_zero(a - b, epsilon);
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
    template <std::ranges::range R>
    auto to_vector(R&& r) {
        auto r_common = r | std::views::common;
        return std::vector(r_common.begin(), r_common.end());
    }
}

#endif // UTILS_H_