/// @file polynomial_interpolation.h
/// @brief
///
/// @author Roland Abel
/// @date November 25, 2023
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

#ifndef POLYNOMIAL_INTERPOLATION_H_
#define POLYNOMIAL_INTERPOLATION_H_

#include "polynomial.h"

namespace xmath {

    /// @brief Provides methods for polynomial interpolation using Lagrange basis.
    /// @tparam T The data type of the coefficients in the polynomials.
    template<typename T>
    class polynomial_interpolation {
    public:
        using polynomial_sequence = std::vector<polynomial<T>>;
        using value_type = polynomial<T>::value_type;
        using values_type = polynomial<T>::values_type;

        static_assert(std::is_floating_point<T>::value, "The type parameter must be a floating point type.");

        /// Gets the Lagrange basis for polynomials of degree k = xs.size().
        /// @param xs The distinct x points x_0, ..., x_k.
        /// @return The Lagrange basis for polynomials. For each polynomial l_j is l_j(x_m) = 0 if j != m and l_j(x_j) = 1.
        static polynomial_sequence lagrange_basis(const std::vector<T> &xs);

        /// Performs Lagrange interpolation to compute the polynomial that passes through a set of given points.
        /// @param x_values The vector containing the x coordinates of the interpolation points.
        /// @param y_values The vector containing the corresponding y coordinates of the interpolation points.
        /// @return A polynomial representing the Lagrange interpolation polynomial.
        static polynomial<T> lagrange_interpolation(const values_type &x_values, const values_type &y_values);
    };
}

#include "polynomial_interpolation.tpp"

#endif // POLYNOMIAL_INTERPOLATION_H_