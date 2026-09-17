/// @file root_finder.h
/// @brief Root finder class using various numerical methods.
///
/// @author Roland Abel
/// @date October 8, 2023
///
/// Copyright (c) 2026 Roland Abel

#ifndef ROOT_FINDER_H_
#define ROOT_FINDER_H_

#include <optional>
#include "interval.h"

namespace xmath {

    /// @brief A class for finding roots using various numerical methods.
    /// @tparam T The data type of the coefficients in the polynomials.
    template<typename T>
    class root_finder {
    public:
        using value_type = T;

        /// @brief Perform the bisection method to find a zero point of the given function within the specified interval.
        /// This method failed if the Intermediate Value Theorem conditions are not met (The function has the same sign at both interval boundaries).
        /// @tparam F The type of the function for which the zero point needs to be found.
        /// @param func The function for which the zero point needs to be found.
        /// @param I The real interval.
        /// @param epsilon The desired accuracy for the root approximation (default is 1e-15).
        /// @return The approximate zero point of the function within the given interval.
        template<typename F>
        static std::optional<T> bisection(
                F &&func,
                const interval<value_type> &I,
                value_type epsilon = 1e-15);

        /// @brief Perform the regula falsi method to find a zero point of the given function within the specified interval.
        /// @tparam F The type of the function for which the zero point needs to be found.
        /// @param func The function for which the zero point needs to be found.
        /// @param I The interval.
        /// @param epsilon The desired accuracy for the root approximation (default is 1e-15).
        /// @return The approximate zero point of the function within the given interval.
        ///
        /// @note If no solution is found, the returned optional<> has not a value.
        template<typename F>
        static std::optional<T> regula_falsi(
                F &&func,
                const interval<value_type> &I,
                value_type epsilon = 1e-15);

        /// @brief Computes an approximation of a root for an equation using the Newton-Raphson method.
        /// @tparam F The type of the function for which to find the root.
        /// @tparam G The type of the derivative function.
        /// @param func The function for which to find the root.
        /// @param derivative The derivative of the func function.
        /// @param initial The initial value for the iteration.
        /// @param max_iterations The maximum number of iterations (default: 100).
        /// @param epsilon The epsilon for the approximation to the root (default: 1e-15).
        /// @return The approximation of the solution for the equation.
        ///
        /// @note If no solution is found, the returned optional<> has not a value.
        template<typename F, typename G>
        static std::optional<T> newton_raphson(
                F &&func,
                G &&derivative,
                value_type initial,
                int max_iterations = 100,
                value_type epsilon = 1e-15);
    };
}

#include "root_finder.tpp"

#endif // ROOT_FINDER_H_