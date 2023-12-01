/// @file chebyshev_1st.cpp
/// @brief Example to demonstrate the computation of Chebyshev Polynomials of the 1st kind.
///
/// This program calculates and prints Chebyshev Polynomials of the 1st kind up to a specified maximum order.
/// It uses the `chebyshev_polynomial<>` class to generate the polynomials.
///
/// @author Roland Abel
/// @date October 14, 2023
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

#include "chebyshev_polynomial.h"
#include "real_polynomial_root_finder.h"

using namespace std;
using namespace xmath;

namespace {
    using ChebyshevPolynomial = chebyshev_polynomial<double>;
    using RootFinder = real_polynomial_root_finder<double>;
}

int main() {
    const auto max_order = 10;

    cout << "Chebyshev Polynomials of 1st kind: " << endl;
    for (int n = 0; n <= max_order; ++n) {
        auto T_n = ChebyshevPolynomial::create_1st_kind(n);
        cout << "T_" << n << ": " << T_n << endl;
    }
    cout << endl;

    cout << "Roots of T_10:" << endl;
    auto [roots, multiplicities] = RootFinder::find_roots(ChebyshevPolynomial::create_1st_kind(10));
    for (int k = 0; k < roots.size(); ++k) {
        cout << "r[" << k << "] = " << roots[k] << endl;
    }

    return 0;
}