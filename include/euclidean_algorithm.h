/// @file euclidean_algorithm.h
/// @brief Implements the Euclidean algorithm for polynomials.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#ifndef EUCLIDEAN_ALGORITHM_H_
#define EUCLIDEAN_ALGORITHM_H_

#include <tuple>
#include "polynomial.h"

namespace xmath {

    /// @brief Provides methods for performing the Euclidean algorithm on polynomials.
    /// @tparam T The data type of the coefficients in the polynomials.
    template<typename T>
    class euclidean_algorithm {
    public:
        /// @brief Computes the greatest common divisor (gcd) of two polynomials.
        /// A greatest common divisor of p and q is a polynomial d that divides p and q,
        /// and such that every common divisor of p and q also divides d.
        /// @param p The left polynomial.
        /// @param q The right polynomial.
        /// @return Returns the gcd of p and q.
        static polynomial<T> euclidean(const polynomial<T> &p, const polynomial<T> &q);

        /// @brief Computes the greatest common divisor (gcd) of two polynomials
        /// and two polynomials s and t such that gcd == s * p + t * q using the
        /// extended Euclidean algorithm.
        ///
        /// Example usage:
        /// @code{.cpp}
        /// auto X = polynomial<double>::monomial(1, 1.0);
        /// auto p = X.pow(4) - 2 * X.pow(3) - 6 * X.pow(2) + 12 * X + 15;
        /// auto q = X.pow(3) + X.pow(2) - 4 * X - 4;
        /// auto [s, t, g] = euclidean_algorithm<double>::extended_euclidean(p, q);
        /// assert(s * p + t * q == g);
        /// @endcode
        ///
        /// @param p The first polynomial for the extended Euclidean algorithm.
        /// @param q The second polynomial for the extended Euclidean algorithm.
        /// @return Returns a tuple containing polynomials s, t, and g such that g = gcd(p, q) = s * p + t * q.
        static std::tuple<polynomial<T>, polynomial<T>, polynomial<T>> extended_euclidean(
                const polynomial<T> &p,
                const polynomial<T> &q);
    };
}

#include "euclidean_algorithm.tpp"

#endif // EUCLIDEAN_ALGORITHM_H_