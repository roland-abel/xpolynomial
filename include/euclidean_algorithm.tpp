/// @file euclidean_algorithm.tpp
/// @brief Implements the Euclidean algorithm for polynomials.
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

#include "euclidean_algorithm.h"

namespace xmath {

    template<typename T>
    polynomial<T> euclidean_algorithm<T>::euclidean(const polynomial<T> &p, const polynomial<T> &q) {
        auto a = p;
        auto b = q;

        while (!b.is_zero()) {
            auto r = a % b; // remainder
            a = b;
            b = r;
        }
        return a.normalize();
    }

    template<typename T>
    std::tuple<polynomial<T>, polynomial<T>, polynomial<T>>
    euclidean_algorithm<T>::extended_euclidean(const polynomial<T> &p, const polynomial<T> &q) {
        auto zero = polynomial<T>::zero();
        auto one = polynomial<T>::one();

        auto a = p;
        auto b = q;

        auto a1 = one;
        auto a2 = zero;

        auto b1 = zero;
        auto b2 = one;

        while (!b.is_zero()) {
            auto [quotient, remainder] = a.divide(b);

            a = b;
            b = remainder;

            auto r1 = a1 - quotient * b1;
            auto r2 = a2 - quotient * b2;

            a1 = b1;
            a2 = b2;

            b1 = r1;
            b2 = r2;
        }
        return std::make_tuple(a1, a2, a);
    }
}