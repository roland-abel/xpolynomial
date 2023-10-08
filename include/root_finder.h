/// @file root_finder.h
///
/// @author Roland Abel
/// @date 08.10.2023

#ifndef ROOT_FINDER_H_
#define ROOT_FINDER_H_

#include <utility>
#include <iostream>
#include <vector>

namespace xmath {

    /// @brief A class for finding roots using various numerical methods.
    /// @tparam T
    template<typename T>
    class root_finder {
    public:
        using value_type = T;

        /// @brief Perform the bisection method to find a zero point of the given function within the specified interval.
        /// @param func The function for which the zero point needs to be found.
        /// @param I The interval.
        /// @param tolerance The desired accuracy for the root approximation (default is 1e-15).
        /// @return The approximate zero point of the function within the given interval.
        static value_type bisection(
                const std::function<value_type(value_type)> &func,
                const interval<value_type> &I,
                value_type tolerance = 1e-5);

        /// @brief Perform the regula falsi method to find a zero point of the given function within the specified interval.
        /// @param func The function for which the zero point needs to be found.
        /// @param I The interval.
        /// @param tolerance The desired accuracy for the root approximation (default is 1e-15).
        /// @return The approximate zero point of the function within the given interval.
        static value_type regula_falsi(
                const std::function<value_type(value_type)> &func,
                const interval<value_type> &I,
                value_type tolerance = 1e-5);

        /// @brief
        /// @param func
        /// @param derive
        /// @param initial
        /// @param max_iterations
        /// @param tolerance
        /// @return
        static value_type newton_raphson(
                const std::function<value_type(value_type)> &func,
                const std::function<value_type(value_type)> &derive,
                value_type initial,
                int max_iterations = 100,
                value_type tolerance = 1e-5);
    };
}

#include "root_finder.tpp"

#endif // ROOT_FINDER_H_