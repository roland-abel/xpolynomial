/// @file square_free_decomposition.h
/// @brief Implements Yun's algorithm for for square-free decomposition.
/// @see https://en.wikipedia.org/wiki/Square-free_polynomial
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

#ifndef SQUARE_FREE_DECOMPOSITION_H_
#define SQUARE_FREE_DECOMPOSITION_H_

#include <vector>
#include "polynomial.h"

namespace xmath {

    /// @brief A square-free decomposition or square-free factorization of a polynomial is a factorization
    /// into powers of square-free polynomials p = q_1 q_2^2 q_3^3 ... qn^n where those of the q_i are
    /// non-constant and pairwise coprime square-free polynomials.
    /// See: https://en.wikipedia.org/wiki/Square-free_polynomial
    /// @tparam T The data type of the coefficients in the polynomials.
    template<typename T>
    class square_free_decomposition {
    public:
        using value_type = typename polynomial<T>::value_type;
        using polynomial_sequence = std::vector<polynomial<T>>;

        /// @brief Checks if the given polynomial is square-free. A polynomial p is square-free if and only
        /// if the greatest common division of the polynomial p and its derivative p' is constant.
        /// @return True if the polynomial p is square-free; otherwise false.
        static bool is_square_free(const polynomial<T> &p);

        /// @brief Computes the content of a polynomial.
        /// This method requires a polynomial whose coefficients are given exactly as integers.
        /// @param p The polynomial.
        /// @return The content of the polynomial, which is the greatest common divisor (gcd) of its coefficients,
        /// or optional has not a value if p is not integral.
        static std::optional<T> content(const polynomial<T> &p);

        /// Calculates the primitive part of a polynomial.
        /// @param p The polynomial.
        /// @return The primitive part of the polynomial or optional has not a value if p is not integral.
        static std::optional<polynomial<T>> primitive_part(const polynomial<T> &p);

        /// @brief Gets the square-free decomposition of the polynomial p.
        /// This method requires a polynomial whose coefficients are given exactly as integers.
        /// @param p The polynomial.
        /// @return The square-free decomposition of p.
        static std::optional<polynomial_sequence> yun_algorithm(const polynomial<T> &p);

        /// @brief Gets the polynomial from the square-free decomposition.
        /// This method requires a polynomial whose coefficients are given exactly as integers.
        /// @param square_free_seq The square-free decomposition.
        /// @return The polynomial which square-free decomposition is given.
        static polynomial<T> from_square_free_decomposition(const polynomial_sequence &square_free_seq);
    };
}

#include "square_free_decomposition.tpp"

#endif // SQUARE_FREE_DECOMPOSITION_H_
