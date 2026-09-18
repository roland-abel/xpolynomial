/// @file nth_roots_of_unity.cpp
/// @brief Example to demonstrating   finding nth roots of unity and checking polynomial roots.
///
/// This program calculates the nth roots of unity and checks if a corresponding polynomial has those roots.
/// It then prints whether the polynomial has roots and the calculated n-th roots.
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
    using ComplexPolynomial = complex_polynomial<double>;
    using RootFinder = complex_polynomial_root_finder<double>;

    const auto Z = ComplexPolynomial::monomial(1, 1.0);
}

auto main() -> int {
    auto n = 13;
    auto p = Z.pow(n) - 1.; // p(Z) = Z^n - 1
    auto roots = RootFinder::nth_roots_of_unity(n);

    cout << "Has roots: " << (p.has_roots(roots) ? "true" : "false") << endl;
    cout << n << "th-roots:" << endl;
    for (auto z: roots) {
        cout << z << endl;
    }
    return 0;
}