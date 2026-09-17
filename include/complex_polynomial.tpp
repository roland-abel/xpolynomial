/// @file complex_polynomial.tpp
/// @brief Defines template classes for polynomials with complex coefficients.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#ifndef COMPLEX_POLYNOMIAL_TPP_H_
#define COMPLEX_POLYNOMIAL_TPP_H_

#include "complex_polynomial.h"

namespace xmath {
    template<typename T>
    using real_values_type = typename polynomial<T>::values_type;

    template<typename T>
    auto separate(const complex_polynomial<T> &polynomial) -> decltype(auto) {
        auto real_coeffs = real_values_type<T>();
        auto imag_coeffs = real_values_type<T>();

        for (auto z: polynomial.coefficients()) {
            real_coeffs.push_back(z.real());
            imag_coeffs.push_back(z.imag());
        }
        return std::pair<real_polynomial<T>, real_polynomial<T> >(real_coeffs, imag_coeffs);
    }
}

#endif // COMPLEX_POLYNOMIAL_TPP_H_
