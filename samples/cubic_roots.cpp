/// @file quadratic_roots.cpp
/// @brief Example to demonstrate polynomial creation and root finding.
///
/// This program creates a cubic polynomial and finds its roots using the RootFinder class.
/// It then prints the polynomial and the calculated roots.
///
/// @author Roland Abel
/// @date October 9, 2023
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

int main() {
    // Create cubic polynomial
    auto p = 3.5 * (X - 7.125).pow(2) * (X - 4.5);

    // Find root for the polynomial p
    auto [r1, r2, r3] = RootFinder::cubic_roots(p);

    cout << "Polynomial: " << p << endl
         << "Roots:" << endl
         << "r1 = " << r1 << endl
         << "r2 = " << r2 << endl
         << "r3 = " << r3 << endl;

    return 0;
}
