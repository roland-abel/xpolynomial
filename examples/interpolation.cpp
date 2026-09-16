/// @file interpolation.cpp
/// @brief Example to demonstrating polynomial interpolation.
///
/// This program creates a 4th-degree polynomial and evaluates it for a range of values.
/// It then prints the result of the polynomial evaluation for each value.
///
/// @author Roland Abel
/// @date December 4, 2023
///
/// Copyright (c) 2026 Roland Abel

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
