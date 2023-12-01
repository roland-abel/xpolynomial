/// @file square_free.cpp
/// @brief Sample for the `square-free decomposition` class.
///
/// @author Roland Abel
/// @date November 28, 2023
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
#include "square_free_decomposition.h"

using namespace std;
using namespace xmath;

namespace {
    auto X = polynomial<double>::monomial(1, 1.0);
    using SquareFree = square_free_decomposition<double>;
}

int main() {
    // Create a polynomial
    auto p = X * (X - 4) * (X + 3.).pow(2) * (X.pow(2) + X).pow(3) * (X.pow(2) + X).pow(4);
    cout << "p = " << p << endl;

    if (!p.is_integer()) {
        cout << "The coefficients of the polynomial must be integers." << endl;
        return -1;
    }

    auto square_free_seq = SquareFree::yun_algorithm(p);
    for (int k = 0; k < square_free_seq.size(); ++k) {
        cout << "q" << k << " = " << square_free_seq[k] << endl;
    }

    cout << "p = " << SquareFree::from_square_free_decomposition(square_free_seq) << endl << endl;
    return 0;
}