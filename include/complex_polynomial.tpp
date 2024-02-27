/// @file complex_polynomial_root_finder.tpp
/// @brief Defines template classes for polynomials with complex coefficients.
///
/// @author Roland Abel
/// @date October 20, 2023
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
