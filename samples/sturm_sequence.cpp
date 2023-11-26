/// @file sturm_sequence.cpp
/// @brief Sample to demonstrate polynomial creation, Sturm's' sequence, and counting distinct roots.
///
/// This program creates a polynomial with 5 distinct real roots and performs operations
/// using the `real_polynomial_root_finder<>` class, including generating the Sturm sequence and
/// counting the number of distinct roots.
///
/// @author Roland Abel
/// @date October 11, 2023
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
    // Create a polynomial with 5 distinct read roots.
    auto p = (X + 2.3) * (X + 1.25) * (X - 0.75) * (X - 1.45) * (X - 2.85);

    auto canonical_seq = RootFinder::sturm_sequence(p);
    auto number_roots = RootFinder::number_distinct_roots(p);

    cout << "Polynomial: " << p << endl << endl
         << "Number of roots: " << number_roots << endl
         << "Canonical polynomial sequence:" << endl;

    for (const auto& q: canonical_seq) {
        cout << q << endl;
    }
    return 0;
}