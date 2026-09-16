/// @file sturm_sequence.cpp
/// @brief Example to demonstrating polynomial creation and finding roots with multiplicities.
///
/// This program creates a polynomial and finds its roots along with their multiplicities using
/// the `real_polynomial_root_finder<>` class.
///
/// @author Roland Abel
/// @date October 11, 2023
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
    // Create a polynomial
    auto p = (X + 3.).pow(3) * (X - 1.) * (X.pow(2) + X - 2).pow(3);
    auto [roots, multiplicities] = RootFinder::find_roots(p);

    cout << "Polynomial: " << p << endl << endl;
    for (int k = 0; k < roots.size(); ++k) {
        cout << "Root: r[" << k << "] = " << roots[k] << ", Multiplicity: "
             << multiplicities[k] << endl;
    }
    return 0;
}