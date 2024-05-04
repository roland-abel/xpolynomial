/// @file interpolation.cpp
/// @brief Example to demonstrating polynomial interpolation.
///
/// This program creates a 4th-degree polynomial and evaluates it for a range of values.
/// It then prints the result of the polynomial evaluation for each value.
///
/// @author Roland Abel
/// @date December 4, 2023
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
#include "polynomial_interpolation.h"

using namespace std;
using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using Interpolation = polynomial_interpolation<double>;
}

auto main() -> int {
    const auto x_values = {-3., -2., -1., 0., 1., 2., 3.};
    const auto y_values = {-2.4, -1.5, 1.1, 2.5, -3.6, -1.25, -2.1};

    const auto p = Interpolation::lagrange_interpolation(x_values, y_values).value();
    cout << "p(x) = " << p << endl;

    for (auto x: x_values) {
        cout << "p(" << x << ") = " << p(x) << endl;
    }
    return 0;
}
