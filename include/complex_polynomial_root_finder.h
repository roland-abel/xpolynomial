/// @file complex_polynomial_root_finder.h
/// @see Jenkins–Traub algorithm - https://en.wikipedia.org/wiki/Jenkins%E2%80%93Traub_algorithm
/// @see https://arxiv.org/pdf/1806.06280.pdf
///
/// @author Roland Abel
/// @date 20.10.2023

#include <vector>
#include <complex>
#include "complex_polynomial.h"

#ifndef COMPLEX_POLYNOMIAL_ROOT_FINDER_H_
#define COMPLEX_POLYNOMIAL_ROOT_FINDER_H_

namespace xmath {

    /// @brief
    /// @tparam T
    template<typename T>
    class complex_polynomial_root_finder {
    public:
        using complex_numbers = std::vector<std::complex<T>>;

        /// @brief Gets a vector of th n-th roots of unity.
        /// @param n Positive integer.
        /// @return The n-th roots of unity.
        static std::vector<std::complex<T>> nth_roots_of_unity(int n);

        /// @brief
        /// @param p
        /// @param initial_points
        /// @return
        static std::vector<std::complex<T>> durand_kerner_method(
                const complex_polynomial<T> &p,
                const std::vector<std::complex<T>> &initial_points,
                size_t max_iterations = 100);
    };
}

#include "complex_polynomial_root_finder.tpp"

#endif // COMPLEX_POLYNOMIAL_ROOT_FINDER_H_
