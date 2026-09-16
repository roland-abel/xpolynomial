/// @file quadratic_roots.cpp
/// @brief Example to demonstrating quadratic polynomial creation and root finding.
///
/// This program creates a quadratic polynomial and finds its roots using the `real_polynomial_root_finder<>` class.
///
/// @author Roland Abel
/// @date October 28, 2023
///
/// Copyright (c) 2026 Roland Abel

#include "real_polynomial_root_finder.h"

using namespace std;
using namespace xmath;

namespace {
    using RootFinder = real_polynomial_root_finder<double>;
    const auto X = polynomial<double>::monomial(1, 1.0);
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
