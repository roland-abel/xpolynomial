/// @file quadratic_roots.cpp
/// @brief Finds the roots of a cubic polynomial.
///
/// @author Roland Abel
/// @date 09.10.2023

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
