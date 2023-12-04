/// @file square_free_decomposition.tpp
/// @brief
///
/// @author Roland Abel
/// @date October 8, 2023
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

#include <vector>
#include <numeric>
#include <ranges>
#include "euclidean_algorithm.h"
#include "square_free_decomposition.h"

namespace xmath {

    template<typename T>
    bool square_free_decomposition<T>::is_square_free(const polynomial<T> &p) {
        return euclidean_algorithm<T>::euclidean(p, p.derive()).is_constant();
    }

    template<typename T>
    std::optional<T> square_free_decomposition<T>::content(const polynomial<T> &p) {
        if (!p.is_integer()) {
            return {};
        }

        auto integer_coeffs = p.coefficients() | std::views::transform([](double coefficient) {
            return std::abs(static_cast<long>(coefficient));
        });

        long gcd = 0;
        for (auto c: integer_coeffs) {
            gcd = std::gcd(gcd, c);
        }
        return gcd;
    }

    template<typename T>
    std::optional<polynomial<T>> square_free_decomposition<T>::primitive_part(const polynomial<T> &p) {
        return content(p).transform([&p](auto c) { return p / c; });
    }

    template<typename T>
    polynomial<T> square_free_decomposition<T>::from_square_free_decomposition(
            const square_free_decomposition::polynomial_sequence &square_free_seq) {

        auto q = polynomial<T>::one();
        for (int k = 0; k < square_free_seq.size(); ++k) {
            q *= square_free_seq[k].pow(k + 1);
        }
        return q;
    }

    template<typename T>
    std::optional<typename square_free_decomposition<T>::polynomial_sequence>
    square_free_decomposition<T>::yun_algorithm(const polynomial<T> &p) {
        auto square_free_seq = xmath::square_free_decomposition<T>::polynomial_sequence();
        if (is_square_free(p)) {
            square_free_seq.push_back(p);
            return square_free_seq;
        }

        return primitive_part(p).transform([&](auto pp) {
            auto p_prim = pp.derive();
            auto q = euclidean_algorithm<T>::euclidean(pp, p_prim).to_integer();
            auto b = pp / q;
            auto c = p_prim / q;
            auto d = c - b.derive();

            q = euclidean_algorithm<T>::euclidean(b, d);
            square_free_seq.push_back(q);

            while (!d.is_zero()) {
                b = b / q;
                c = d / q;
                d = c - b.derive();

                q = euclidean_algorithm<T>::euclidean(b, d);
                square_free_seq.push_back(q);
            }
            return square_free_seq;
        });
    }
}

