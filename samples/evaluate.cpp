/// @file evaluate.cpp
/// @brief Example to demonstrating polynomial creation and evaluation.
///
/// This program creates a 4th-degree polynomial and evaluates it for a range of values.
/// It then prints the result of the polynomial evaluation for each value.
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

#include <numeric>
#include <vector>
#include <iostream>
#include "polynomial.h"

using namespace std;
using namespace xmath;

namespace {
    auto X = polynomial<double>::monomial(1, 1.0);
}

int main() {
    // Create a 4th degree polynomial
    auto p = 3 * X.pow(4) - 2.5 * X.pow(3) + X.pow(2) - X + 1;

    // Evaluate the polynomial for a range of values
    vector<double> values(10);
    iota(values.begin(), values.end(), 1);

    for (auto x: values) {
        cout << "p(" << x << ") = " << p(x) << endl;
    }
    return 0;
}
