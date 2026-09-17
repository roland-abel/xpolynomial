/// @file polynomial_interpolation.tpp
/// @brief Lagrange polynomial interpolation.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#ifndef POLYNOMIAL_INTERPOLATION_TPP_H_
#define POLYNOMIAL_INTERPOLATION_TPP_H_

#include <cstdint>
#include "polynomial_interpolation.h"

namespace xmath {
    template<typename T>
    using polynomial_sequence = typename polynomial_interpolation<T>::polynomial_sequence;

    template<typename T>
    typename polynomial_interpolation<T>::polynomial_sequence
    polynomial_interpolation<T>::lagrange_basis(const std::vector<T> &xs) {
        const auto &one = polynomial<T>::one();
        const auto &X = polynomial<T>::monomial(1);

        const size_t N = xs.size();
        auto lagrange = [&](const uint16_t j) {
            auto b = one;
            for (size_t i = 0; i < N; ++i) {
                if (i == j) {
                    continue;
                }
                b *= (X - xs[i]) / (xs[j] - xs[i]);
            }
            return b;
        };

        polynomial_sequence basis;
        for (size_t i = 0; i < N; ++i) {
            basis.push_back(lagrange(i));
        }
        return basis;
    }

    template<typename T>
    std::optional<polynomial<T>> polynomial_interpolation<T>::lagrange_interpolation(
            const values_type &x_values,
            const values_type &y_values) {
        if (x_values.size() != y_values.size() || x_values.empty()) {
            return {};
        }

        const auto basis_polynomials = lagrange_basis(x_values);

        polynomial<T> p;
        for (size_t i = 0; i < x_values.size(); ++i) {
            p += basis_polynomials[i] * y_values[i];
        }
        return p;
    }
}

#endif // POLYNOMIAL_INTERPOLATION_TPP_H_
