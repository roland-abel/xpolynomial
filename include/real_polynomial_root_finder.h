/// @file real_polynomial_root_finder.h
/// @brief Root finder class of polynomials with real coefficients.
///
/// @author Roland Abel
/// @date August 19, 2023
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

#ifndef REAL_POLYNOMIAL_ROOT_FINDER_H_
#define REAL_POLYNOMIAL_ROOT_FINDER_H_

#include <utility>
#include <iostream>
#include <vector>
#include "interval.h"
#include "polynomial.h"

namespace xmath {

    /// @brief A class for finding roots of polynomials with real coefficients using various numerical methods.
    /// @tparam T The data type of the coefficients in the polynomials.
    template<typename T>
    class real_polynomial_root_finder {
    public:
        using value_type = typename polynomial<T>::value_type;
        using polynomial_sequence = std::vector<polynomial<T>>;
        using roots_type = std::vector<value_type>;
        using multiplicities_type = std::vector<unsigned short>;

        /// @brief Computes the roots of a quadratic polynomial using the quadratic formula.
        /// @param p The quadratic polynomial for which the roots are to be computed.
        /// @return std::pair contains the two roots of the quadratic polynomial,
        /// or no value if the polynomial has no real roots or is not quadratic.
        static std::optional<std::tuple<T, T>> quadratic_roots(const polynomial<T> &p);

        /// @brief Checks if a given polynomial has a cubic normal form. A cubic normal form is given if the polynomial
        /// can be expressed in the following from: p = X^3 + aX + b.
        /// @param p The polynomial to be checked for a cubic normal form.
        /// @return True if the given polynomial has a cubic normal form; false otherwise.
        static bool has_cubic_normal_form(const polynomial<T> &p);

        /// @brief Finds the real roots of the cubic polynomial given in the normal form p = X^3 + aX + b.
        /// @param p The cubic polynomial in normal form.
        /// @return A tuple containing the real roots of the polynomial,
        /// or the tuple (NaN, NaN, NaN) if no valid real roots are found.
        static std::vector<T> cubic_normal_form_roots(const polynomial<T> &p);

        /// @brief Finds the real roots of a cubic polynomial.
        /// @param p The cubic polynomial to find the roots of.
        /// @return A tuple containing the real roots of the polynomial.
        static std::vector<T> cubic_roots(const polynomial<T> &p);

        /// @brief Find a root of a polynomial using Newton's method, starting from a given initial guess.
        /// @param p The polynomial for which the root needs to be found.
        /// @param initial The initial guess for the root.
        /// @param tolerance The desired accuracy for the root approximation (default is 1e-15).
        /// @param max_iterations The maximum number of iterations for the algorithm (default is 100).
        /// @return The approximate root of the polynomial.
        static std::optional<T> newton_raphson(
                const polynomial<T> &p,
                value_type initial,
                int max_iterations = 100,
                value_type tolerance = polynomial<T>::epsilon);

        /// @brief Gets the number of sign changes of the coefficients of the given polynomial.
        /// @param p The polynomial for which the number of sign changes are to determined.
        /// @return The number of sign changes.
        static size_t sign_changes(const polynomial<T> &p);

        /// @brief Calculates the Cauchy's bounds for the real roots of a polynomial.
        /// @see https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots
        /// @param p The polynomial for which to determine the Cauchy's bounds.
        /// @return The calculated Cauchy's bounds for the real roots of p.
        static std::optional<T> cauchy_bounds(const polynomial<T> &p);

        /// @brief Calculates the Lagrange's bound for the real roots of a polynomial.
        /// @see https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots
        /// @param p The polynomial for which to determine the Lagrange's bounds.
        /// @return The calculated Lagrange's bounds for the real roots of p.
        static std::optional<T> lagrange_bounds(const polynomial<T> &p);

        /// @brief Gets the Sturm's polynomial sequence.
        /// @see https://en.wikipedia.org/wiki/Sturm%27s_theorem
        /// @param p The polynomial for which the Sturm's sequence are to determine.
        /// @return The canonical sequence.
        static polynomial_sequence sturm_sequence(const polynomial<T> &p);

        /// @brief Gets the sign variations at x for the given sequence of polynomials.
        /// The sign variations at a point of the Sturm sequence of p is the number of sign changes.
        /// @see https://en.wikipedia.org/wiki/Sturm%27s_theorem
        /// @param seq The sequence of polynomials.
        /// @param x The value for which the sign variations are to determined.
        /// @return Sequence of sign variations.
        static std::vector<short> sign_variations(const polynomial_sequence &seq, const value_type &x);

        /// @brief Counts the number of distinct real roots of a polynomial within a specified half-open interval,
        /// if the given polynomial is square-free polynomial.
        /// @param p The square-free polynomial for which the number of roots are to determined.
        /// @param I The real half-open interval (a, b] in which to find the roots of p.
        /// @return The number of distinct real roots within the specified interval,
        /// or no value if the polynomial p is not square-free.
        static std::optional<size_t> number_distinct_roots(const polynomial<T> &p, const interval<T> &I);

        /// @brief Counts the number of distinct real roots of a square-free polynomial.
        /// @param p The square-free polynomial for which the number of roots are to determined.
        /// @return The number of distinct real roots, or no value if the polynomial p is not square-free.
        static std::optional<size_t> number_distinct_roots(const polynomial<T> &p);

        /// @brief Isolates the real roots of a polynomial within specified intervals.
        /// @param p The polynomial for which to isolate the real roots.
        /// @return A vector of intervals containing the isolated real roots
        /// or an empty vector if no roots are found or if polynomial is not square-free.
        static std::vector<interval<T>> root_isolation(const polynomial<T> &p);

        /// @brief Finds the roots of a polynomial. This method requires a polynomial whose coefficients
        /// are given exactly as integers.
        /// @param p The polynomial for which to find the roots.
        /// @param epsilon The desired precision (default 1e-15).
        /// @return A tuple containing the roots and their multiplicities, or empty lists if the coefficients of
        /// the given polynomial are not all integers.
        static std::tuple<roots_type, multiplicities_type> find_roots(const polynomial<T> &p, T epsilon = 1e-15);
    };
}

#include "real_polynomial_root_finder.tpp"

#endif // REAL_POLYNOMIAL_ROOT_FINDER_H_