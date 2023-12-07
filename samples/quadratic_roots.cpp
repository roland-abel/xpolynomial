/// @file quadratic_roots.cpp
/// @brief Example to demonstrating quadratic polynomial creation and root finding.
///
/// This program creates a quadratic polynomial and finds its roots using the `real_polynomial_root_finder<>` class.
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

#include "real_polynomial_root_finder.h"

using namespace std;
using namespace xmath;

namespace {
    using RootFinder = real_polynomial_root_finder<double>;
    auto X = polynomial<double>::monomial(1, 1.0);
}

auto main() -> int {
    // Create quadratic polynomial
    auto p = .25 * X.pow(2) - 1.5 * X - 1;

    // Find root for the polynomial p
    auto [r1, r2] = RootFinder::quadratic_roots(p).value();

    cout << "Polynomial: " << p << endl
         << "Roots:" << endl
         << "r1 = " << r1 << endl
         << "r2 = " << r2 << endl;

    return 0;
}
