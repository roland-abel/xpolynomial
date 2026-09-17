/// @file complex_polynomial_root_finder.h
/// @brief Root finder class for polynomials with complex coefficients.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#ifndef COMPLEX_POLYNOMIAL_ROOT_FINDER_H_
#define COMPLEX_POLYNOMIAL_ROOT_FINDER_H_

#include <vector>
#include <complex>
#include "complex_polynomial.h"

namespace xmath {

    /// @brief Provides methods for finding roots of complex polynomials using various algorithms.
    /// @tparam T The data type of the coefficients in the complex polynomial.
    template<typename T>
    class complex_polynomial_root_finder {
    public:
        /// @brief Gets a vector of the n-th roots of unity.
        /// @param n A positive integer specifying the order of the roots.
        /// @return A vector containing the n-th roots of unity.
        static std::vector<std::complex<T>> nth_roots_of_unity(int n);

        /// @brief Finds roots of a complex polynomial using the Durand-Kerner method.
        /// @param p The complex polynomial for which roots are to be found.
        /// @param initial_points The initial points to start the iteration.
        /// @param max_iterations The maximum number of iterations (default: 100).
        /// @return A vector containing the roots of the complex polynomial.
        static std::vector<std::complex<T>> durand_kerner_method(
                const complex_polynomial<T> &p,
                const std::vector<std::complex<T>> &initial_points,
                size_t max_iterations = 100);

        /// @brief Finds roots of a complex polynomial using the Aberth-Ehrlich method.
        /// @param p The complex polynomial for which roots are to be found.
        /// @param initial_points The initial points to start the iteration.
        /// @param max_iterations The maximum number of iterations (default: 100).
        /// @return A vector containing the roots of the complex polynomial.
        static std::vector<std::complex<T>> aberth_ehrlich_method(
                const complex_polynomial<T> &p,
                const std::vector<std::complex<T>> &initial_points,
                size_t max_iterations = 100);
    };
}

#include "complex_polynomial_root_finder.tpp"

#endif // COMPLEX_POLYNOMIAL_ROOT_FINDER_H_
