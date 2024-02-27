/// @file polynomial_interpolation.tpp
/// @brief Lagrange polynomial interpolation.
///
/// @author Roland Abel
/// @date November 25, 2023
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

#pragma once

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
    std::optional<polynomial<T> > polynomial_interpolation<T>::lagrange_interpolation(
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
