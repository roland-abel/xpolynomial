/// @file nth_roots_of_unity.cpp
/// @brief Sample to demonstrate finding nth roots of unity and checking polynomial roots.
///
/// This program calculates the nth roots of unity and checks if a corresponding polynomial has those roots.
/// It then prints whether the polynomial has roots and the calculated n-th roots.
///
/// @author Roland Abel
/// @date October 28, 2023
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

#include <iostream>
#include "complex_polynomial_root_finder.h"

using namespace std;
using namespace xmath;

namespace {
    using ComplexPolynomial = complex_polynomial<double>;
    using RealPolynomial = polynomial<double, polynomial_specification<double>>;
    using RootFinder = complex_polynomial_root_finder<double>;

    auto X = RealPolynomial::monomial(1, 1.0);
    auto Y = RealPolynomial::monomial(1, 1.0);
    auto Z = ComplexPolynomial::monomial(1, 1.0);
}

int main() {
    auto n = 13;
    auto p = Z.pow(n) - 1.; // p(Z) = Z^n - 1
    auto roots = RootFinder::nth_roots_of_unity(n);

    cout << "Has roots: " << (p.has_roots(roots) ? "true" : "false") << endl;
    cout << n << "th-roots:" << endl;
    for (auto z: roots) {
        cout << z << endl;
    }
    return 0;
}