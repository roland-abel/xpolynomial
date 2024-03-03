/// @file root_finder.h
/// @brief Root finder class using various numerical methods.
///
/// @author Roland Abel
/// @date October 8, 2023
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

#ifndef ROOT_FINDER_H_
#define ROOT_FINDER_H_

#include <functional>
#include <optional>

namespace xmath {

    /// @brief A class for finding roots using various numerical methods.
    /// @tparam T The data type of the coefficients in the polynomials.
    template<typename T>
    class root_finder {
    public:
        using value_type = T;

        /// @brief Perform the bisection method to find a zero point of the given function within the specified interval.
        /// This method failed if the Intermediate Value Theorem conditions are not met (The function has the same sign at both interval boundaries).
        /// @param func The function for which the zero point needs to be found.
        /// @param I The real interval.
        /// @param epsilon The desired accuracy for the root approximation (default is 1e-15).
        /// @return The approximate zero point of the function within the given interval.
        static std::optional<T> bisection(
                const std::function<value_type(value_type)> &func,
                const interval<value_type> &I,
                value_type epsilon = 1e-15);

        /// @brief Perform the regula falsi method to find a zero point of the given function within the specified interval.
        /// @param func The function for which the zero point needs to be found.
        /// @param I The interval.
        /// @param epsilon The desired accuracy for the root approximation (default is 1e-15).
        /// @return The approximate zero point of the function within the given interval.
        ///
        /// @note If no solution is found, the returned optional<> has not a value.
        static std::optional<T> regula_falsi(
                const std::function<value_type(value_type)> &func,
                const interval<value_type> &I,
                value_type epsilon = 1e-15);

        /// @brief Computes an approximation of a root for an equation using the Newton-Raphson method.
        /// @param func The function for which to find the root.
        /// @param derive The derivative of the func function.
        /// @param initial The initial value for the iteration.
        /// @param max_iterations The maximum number of iterations (default: 100).
        /// @param epsilon The epsilon for the approximation to the root (default: 1e-15).
        /// @return The approximation of the solution for the equation.
        ///
        /// @note If no solution is found, the returned optional<> has not a value.
        static std::optional<T> newton_raphson(
                const std::function<value_type(const value_type &)> &func,
                const std::function<value_type(const value_type &)> &derive,
                value_type initial,
                int max_iterations = 100,
                value_type epsilon = 1e-15);
    };
}

#include "root_finder.tpp"

#endif // ROOT_FINDER_H_