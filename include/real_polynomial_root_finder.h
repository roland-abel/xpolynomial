//
// Created by abel on 26.08.2023.
//

#ifndef REAL_POLYNOMIAL_ROOT_FINDER_H_
#define REAL_POLYNOMIAL_ROOT_FINDER_H_

#include <utility>
#include <iostream>
#include <vector>
#include "polynomial.h"

namespace xmath {

    /// @brief A class for finding roots of polynomials with real coefficients using various numerical methods.
    template<typename T>
    class real_polynomial_root_finder {
    public:
        using value_type = polynomial<T>::value_type;
        using polynomial_sequence = std::vector<polynomial<T>>;

        /// @brief Computes the roots of a quadratic polynomial using the quadratic formula.
        /// @param p The quadratic polynomial for which the roots are to be computed.
        /// @return std::pair contains the two roots of the quadratic polynomial, or
        ///         std::pair(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN())
        ///         if the polynomial has no real roots or is not quadratic.
        static std::tuple<value_type, value_type> quadratic_roots(const polynomial<T> &p);

        /// @brief
        /// @param p
        /// @return
        static bool has_cubic_normal_form(const polynomial<T> &p);

        /// @brief Finds the real roots of the cubic polynomial given in the normal form p = X^3 + aX + b.
        /// @param p The cubic polynomial in normal form.
        /// @return The real roots of the polynomial.
        static std::tuple<value_type, value_type, value_type> cubic_normal_form_roots(const polynomial<T> &p);

        /// @brief
        /// @param p
        /// @return
        static std::tuple<value_type, value_type, value_type> cubic_roots(const polynomial<T> &p);

        /// @brief Find a root of a polynomial using Newton's method, starting from a given initial guess.
        /// @param p The polynomial for which the root needs to be found.
        /// @param initial The initial guess for the root.
        /// @param tolerance The desired accuracy for the root approximation (default is 1e-15).
        /// @param max_iterations The maximum number of iterations for the algorithm (default is 100).
        /// @return The approximate root of the polynomial.
        static value_type newton_raphson(
                const polynomial<T> &p,
                value_type initial,
                int max_iterations = 100,
                value_type tolerance = polynomial<T>::tolerance);

        /// @brief Gets the number of sign changes of the coefficients of the given polynomial.
        /// @param p The polynomial for which the number of sign changes are to determined.
        /// @return The number of sign changes.
        static size_t sign_changes(const polynomial<T> &p);

        /// @brief
        /// @param p
        /// @return
        static value_type cauchy_bound(const polynomial<T> &p);

        /// @brief Gets the Sturm's polynomial sequence.
        /// @param p The polynomial for which the Sturm's sequence are to determine.
        /// @return The canonical sequence.
        static polynomial_sequence sturm_sequence(const polynomial<T> &p);

        /// @brief Gets the sign variations at x for the given sequence of polynomials.
        /// @param seq The sequence of polynomials.
        /// @param x The value for which the sign variations are to determined.
        /// @return Sequence of sign variations.
        static std::vector<short> sign_variations(const polynomial_sequence &seq, const value_type &x);

        /// @brief Gets the number of distinct real roots of p within the interval I, if the given square-free polynomial p.
        /// @return The number of distinct real roots of p,
        /// @param p The square-free polynomial for which the number of roots are to determined.
        /// @param I The real interval in which to find the roots of p.
        /// @return
        static int number_distinct_roots(const polynomial<T> &p, const interval<T> &I);

        /// @brief Gets the number of distinct real roots of p, if the given square-free polynomial p.
        /// @param p The square-free polynomial for which the number of roots are to determined.
        /// @return The number of distinct real roots of p,
        static int number_distinct_roots(const polynomial<T> &p);

        /// @brief Gets disjoint intervals, such that each one contains extract one real root, and together
        /// they contains all's the roots of the  given polynomial.
        /// @param p The polynomial.
        /// @return The list of the disjoint intervals.
        static std::vector<interval<T>> root_isolation(const polynomial<T> &p);

        /// @brief
        /// @param p
        /// @param tolerance
        /// @return
        static std::vector<value_type> find_roots(const polynomial<T> &p, value_type tolerance = polynomial<T>::tolerance);
    };
}

#include "real_polynomial_root_finder.tpp"

#endif // REAL_POLYNOMIAL_ROOT_FINDER_H_