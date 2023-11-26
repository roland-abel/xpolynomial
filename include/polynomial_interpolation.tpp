/// @file polynomial_interpolation.tpp
///
/// @author Roland Abel
/// @date 25.11.2023

#include "polynomial_interpolation.h"

namespace xmath {

    namespace {
        template<typename T>
        using polynomial_sequence = polynomial_interpolation<T>::polynomial_sequence;
    }

    template<typename T>
    polynomial_interpolation<T>::polynomial_sequence
    polynomial_interpolation<T>::lagrange_basis(const std::vector<T> &xs) {
        const auto &one = polynomial<T>::one();
        const auto &X = polynomial<T>::monomial(1);

        const uint16_t N = xs.size();
        auto lagrange = [&](const uint16_t j) {
            auto b = one;
            for (int i = 0; i < N; ++i) {
                if (i == j) {
                    continue;
                }
                b *= (X - xs[i]) / (xs[j] - xs[i]);
            }
            return b;
        };

        polynomial_sequence basis;
        for (int i = 0; i < N; ++i) {
            basis.push_back(lagrange(i));
        }
        return basis;
    }

    template<typename T>
    polynomial<T> polynomial_interpolation<T>::lagrange_interpolation(
            const values_type &x_values,
            const values_type &y_values) {

        if (x_values.size() != y_values.size() || x_values.empty()) {
            throw std::invalid_argument("Input vectors must have the same non-zero size.");
        }

        const auto basis_polynomials = lagrange_basis(x_values);

        polynomial<T> p;
        for (size_t i = 0; i < x_values.size(); ++i) {
            p += basis_polynomials[i] * y_values[i];
        }
        return p;
    }
}