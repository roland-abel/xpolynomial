/// @file root_finder.h
///
/// @author Roland Abel
/// @date 08.10.2023

#ifndef ROOT_FINDER_H_
#define ROOT_FINDER_H_

#include <utility>
#include <iostream>
#include <vector>
#include <functional>

namespace xmath {

    /// @brief A class for finding roots using various numerical methods.
    /// @tparam T
    template<typename T>
    class root_finder {
    public:
        using value_type = T;

        /// @brief Perform the bisection method to find a zero point of the given function within the specified real_interval.
        /// @param func The function for which the zero point needs to be found.
        /// @param I The real interval.
        /// @param tolerance The desired accuracy for the root approximation (default is 1e-15).
        /// @return The approximate zero point of the function within the given real_interval,
        /// or NaN if the Intermediate Value Theorem conditions are not met (The function has the same sign at both
        /// real_interval boundaries).
        static value_type bisection(
                const std::function<value_type(value_type)> &func,
                const real_interval<value_type> &I,
                value_type tolerance = 1e-5);

        /// @brief Perform the regula falsi method to find a zero point of the given function within the specified real_interval.
        /// @param func The function for which the zero point needs to be found.
        /// @param I The real_interval.
        /// @param epsilon The desired accuracy for the root approximation (default is 1e-15).
        /// @return The approximate zero point of the function within the given real_interval.
        ///
        /// @note If no solution is found, an undefined value may be returned (NaN).
        /// It is recommended to check the validity of the result.
        static value_type regula_falsi(
                const std::function<value_type(value_type)> &func,
                const real_interval<value_type> &I,
                value_type epsilon = 1e-5);

        /// @brief Computes an approximation of a root for an equation using the Newton-Raphson method.
        /// @param func The function for which to find the root.
        /// @param derive The derivative of the func function.
        /// @param initial The initial value for the iteration.
        /// @param max_iterations The maximum number of iterations (default: 100).
        /// @param tolerance The tolerance for the approximation to the root (default: 1e-5).
        /// @return The approximation of the solution for the equation.
        ///
        /// @note If no solution is found, an undefined value may be returned (NaN).
        /// It is recommended to check the validity of the result.
        static value_type newton_raphson(
                const std::function<value_type(const value_type &)> &func,
                const std::function<value_type(const value_type &)> &derive,
                value_type initial,
                int max_iterations = 100,
                value_type tolerance = 1e-5);
    };
}

#include "root_finder.tpp"

#endif // ROOT_FINDER_H_