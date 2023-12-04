/// @file chebyshev_polynomial.h
/// @brief
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

#ifndef CHEBYSHEV_POLYNOMIAL_H_
#define CHEBYSHEV_POLYNOMIAL_H_

#include <vector>
#include <functional>
#include <cstdint>
#include "real_interval.h"
#include "polynomial.h"

namespace xmath {

    /// @brief Provides methods useful for dealing with Chebyshev polynomials of the first and second kind.
    /// @tparam T The data type of the coefficients.
    template<typename T>
    class chebyshev_polynomial {
    public:
        using polynomial_sequence = std::vector<polynomial<T>>;
        using value_type = polynomial<T>::value_type;
        using values_type = polynomial<T>::values_type;

        static_assert(std::is_floating_point<T>::value, "The type parameter must be a floating point type.");

        /// @brief Creates first kind Chebyshev polynomial T_n for the given N.
        /// @param N The N of the polynomial.
        /// @return The Chebyshev polynomial of the first kind of the N.
        static polynomial<T> create_1st_kind(size_t N);

        /// Creates first kind Chebyshev polynomial T_n for the given N.
        /// @param N The N of the polynomial.
        /// @param chebyshev_cache The Chebyshev polynomial cache is used to avoid repeated calculations.
        /// @return The Chebyshev polynomial of the first kind of the N.
        static polynomial<T> create_1st_kind(size_t N, chebyshev_polynomial<T>::polynomial_sequence &chebyshev_cache);

        /// @brief Computes the Chebyshev nodes for a given N within a specified interval.
        /// @param N  The N of the Chebyshev nodes.
        /// @return A vector containing the Chebyshev nodes for the specified N within the interval.
        static values_type chebyshev_nodes(size_t N, const real_interval<T> &I = real_interval<value_type>(-1., 1.));

        /// Calculates the Chebyshev series of order n for the given point and the coefficients alphas.
        /// @param alphas A vector containing the coefficients for the Chebyshev series polynomial.
        /// @return The Chebyshev series of order n for the given coefficients alphas as polynomial.
        static polynomial<T> chebyshev_series(const values_type &alphas);

        /// Evaluates the Chebyshev series for the given coefficients and input value by using Clenshaw algorithm.
        /// @param alphas A vector containing the coefficients for the Chebyshev series polynomial.
        /// @param x The input x value.
        /// @return The value of the Chebyshev series of order n for the given coefficients alphas at x.
        static value_type clenshaw(const values_type &alphas, const value_type &x);

        /// @brief Compute the numerical integral for a function defined in [-1., 1.] with the weight function
        /// 1/(1-x^2) by using the Chebyshev-Gauss quadrature rule.
        /// @param func The function f defined in the interval [-1., 1.].
        /// @param N The number of iteration.
        /// @return The approximation of the integral ∑_i=1^N w_i f(x_i) where the weights are the constant
        /// w_i = π/N and x_i's are the Chebyshev nodes.
        static value_type chebyshev_quadrature(std::function<value_type(value_type)> func, uint32_t N = 5);

    private:
        static polynomial_sequence chebyshev_1st_kind_polynomials_;
    };
}

#include "chebyshev_polynomial.tpp"

#endif // CHEBYSHEV_POLYNOMIAL_H_