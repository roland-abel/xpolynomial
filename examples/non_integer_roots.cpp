/// @file non_integer_roots.cpp
/// @brief Example to demonstrating root finding with non-integer coefficients.
///
/// This program shows that `real_polynomial_root_finder<>::find_roots` returns
/// an empty result for polynomials whose coefficients are not all integers.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#include <iostream>
#include <xpolynomial.h>

using namespace std;
using namespace xmath;

namespace {
    using RootFinder = real_polynomial_root_finder<double>;
    const auto X = polynomial<double>::monomial(1, 1.0);
}

auto main() -> int {
    auto p = 1.5 * (X + 2.) * (X - 1.);

    cout << "p(x) = " << p << endl
         << "is_integer = " << boolalpha << p.is_integer() << endl << endl;

    auto [roots, multiplicities] = RootFinder::find_roots(p);
    if (roots.empty()) {
        cout << "find_roots requires integer coefficients and returned no roots." << endl;
    } else {
        for (std::size_t k = 0; k < roots.size(); ++k) {
            cout << "root r[" << k << "] = " << roots[k]
                 << ", multiplicity = " << multiplicities[k] << endl;
        }
    }
    return 0;
}