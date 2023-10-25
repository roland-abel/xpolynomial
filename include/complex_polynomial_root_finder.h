/// @file complex_polynomial_root_finder.h
/// @see Jenkins–Traub algorithm - https://en.wikipedia.org/wiki/Jenkins%E2%80%93Traub_algorithm
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
    };
}

#include "complex_polynomial_root_finder.tpp"

#endif // COMPLEX_POLYNOMIAL_ROOT_FINDER_H_
