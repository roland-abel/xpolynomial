//
// Created by abel on 30.08.2023.
//

#ifndef UTILS_H_
#define UTILS_H_

#include <limits>
#include <cmath>

namespace xmath {

    template<typename T>
    class interval : std::pair<T, T> {
    public:
        interval(T a, T b) : std::pair<T, T>(a, b) {
        }

        inline T start() const { return std::pair<T, T>::first; }

        inline T end() const { return std::pair<T, T>::second;; }

        inline T length() const { return end() - start(); }

        inline bool is_empty() const { return start() > end(); }

        std::pair<interval<T>, interval<T>> bisect() const {
            return std::make_pair(
                    interval(start(), (start() + end()) / 2.),
                    interval((start() + end()) / 2., end()));
        }
    };

    template<typename T>
    inline bool nearly_zero(T a, T tolerance = std::numeric_limits<T>::epsilon()) {
        return std::abs(a) < tolerance;
    }

    template<typename T>
    inline bool nearly_equal(T a, T b, T tolerance = std::numeric_limits<T>::epsilon()) {
        return nearly_zero(a - b, tolerance);
    }

    inline bool is_even(long a) {
        return a % 2 == 0;
    }

    inline bool is_odd(long a) {
        return a % 2 == 1;
    }

    /// @brief Gets the number of sign changes of the given sequence.
    /// @param sequence The sequence for which the number of sign changes are to determined.
    /// @param tolerance
    /// @return The number of sign changes.
    template<typename T>
    size_t sign_changes(const std::vector<T> &sequence, T tolerance = 1e-5) {
        auto changes = 0;
        auto size = sequence.size();

        if (size <= 1) {
            return 0; // A sequence with only one value hasn't a sign change.
        }

        int prev_sign = (sequence[0] >= 0) ? 1 : -1;

        for (int i = 1; i < size; ++i) {
            if (nearly_zero<T>(sequence[i], tolerance)) {
                continue;
            }

            int current_sign = (sequence[i] >= 0) ? 1 : -1;
            if (current_sign != prev_sign) {
                changes++;
                prev_sign = current_sign;
            }
        }
        return changes;
    }
}

#endif // UTILS_H_