/// @file real_interval.h
///
/// @author Roland Abel
/// @date 07.11.2023

#ifndef REAL_INTERVAL_H_
#define REAL_INTERVAL_H_

#include <limits>
#include <cmath>

namespace xmath {

    template<typename T>
    class real_interval : std::pair<T, T> {
    public:
        real_interval(T a, T b) : std::pair<T, T>(a, b) {
        }

        inline T start() const { return std::pair<T, T>::first; }

        inline T end() const { return std::pair<T, T>::second; }

        inline T length() const { return end() - start(); }

        inline bool is_empty() const { return start() > end(); }

        std::pair<real_interval<T>, real_interval<T>> bisect() const {
            return std::make_pair(
                    real_interval(start(), (start() + end()) / 2.),
                    real_interval((start() + end()) / 2., end()));
        }
    };
}

#endif // REAL_INTERVAL_H_
