/// @file quadratic_roots.cpp
/// @brief Finds the roots of a quadratic polynomial.
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
    // Create quadratic polynomial
    auto p = .25 * X.pow(2) - 1.5 * X - 1;

    // Find root for the polynomial p
    auto [r1, r2] = RootFinder::quadratic_roots(p);

    cout << "Polynomial: " << p << endl
         << "Roots:" << endl
         << "r1 = " << r1 << endl
         << "r2 = " << r2 << endl;

    return 0;
}
