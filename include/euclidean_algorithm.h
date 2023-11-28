/// @file euclidean_algorithm.h
/// @brief
///
/// @author Roland Abel
/// @date October 20, 2023
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

#ifndef EUCLIDEAN_ALGORITHM_H_
#define EUCLIDEAN_ALGORITHM_H_

#include <vector>
#include "polynomial.h"

namespace xmath {

    /// @brief
    template<typename T>
    class euclidean_algorithm {
    public:
        /// @brief Computes the greatest common divisor (gcd) of two polynomials. A greatest common divisor
        /// of p and q is a polynomial d that divides p and q, and such that every common divisor of p and q
        /// also divides d.
        /// @param p The left polynomial.
        /// @param q The right polynomial.
        /// @return Returns the gcd of p and q.
        static polynomial<T> euclidean(const polynomial<T> &p, const polynomial<T> &q);

        /// @brief Computes the greatest common divisor (gcd) of two polynomials and two polynomials s and t such that
        /// gcd == s * p + t * q by using the extended Euclidean algorithm.
        /// Example usage:
        /// @code{.cpp}
        ///
        /// auto X = polynomial<double>::monomial(1, 1.0);
        /// auto p = X.pow(4) - 2 * X.pow(3) - 6 * X.pow(2) + 12 * X + 15;
        /// auto q = X.pow(3) + X.pow(2) - 4 * X - 4;
        /// auto [s, t, g] = EuclideanAlgorithm::extended_euclidean(p, q);
        /// s * p + t * q == g
        ///
        /// @endcode
        /// @param p The first polynomial for the extended Euclidean algorithm.
        /// @param q The second polynomial for the extended Euclidean algorithm.
        /// @return Returns polynomials s, t, g such that g = gcd(p, q) = s*p + t*q.
        static std::tuple<polynomial<T>, polynomial<T>, polynomial<T>> extended_euclidean(
                const polynomial<T> &p,
                const polynomial<T> &q);
    };
}

#include "euclidean_algorithm.tpp"

#endif // EUCLIDEAN_ALGORITHM_H_