/// @file euclidean_algorithm.h
///
/// @author Roland Abel
/// @date 08.10.2023

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

        /// @brief Computes the extended Euclidean algorithm for two polynomials.
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