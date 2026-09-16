/// @file quadratic_roots.cpp
/// @brief Example to demonstrating polynomial creation and root finding.
///
/// This program creates a cubic polynomial and finds its roots using the RootFinder class.
/// It then prints the polynomial and the calculated roots.
///
/// @author Roland Abel
/// @date October 9, 2023
///
/// Copyright (c) 2026 Roland Abel

#include <xpolynomial.h>

using namespace std;
using namespace xmath;

namespace {
    using RootFinder = real_polynomial_root_finder<double>;
    const auto X = polynomial<double>::monomial(1, 1.0);
}

auto main() -> int {
    // Create cubic polynomial
    auto p = 3.5 * (X - 7.125).pow(2) * (X - 4.5);
    cout << "p(x) = " << p << endl;

    // Find root for the polynomial p
    cout << "Roots of p" << endl;
    for (auto r: RootFinder::cubic_roots(p)) {
        cout << r << endl;
    }
    return 0;
}
