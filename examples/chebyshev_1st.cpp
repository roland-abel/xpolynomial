/// @file chebyshev_1st.cpp
/// @brief Example to demonstrating the computation of Chebyshev Polynomials of the 1st kind.
///
/// This program calculates and prints Chebyshev Polynomials of the 1st kind up to a specified maximum order.
/// It uses the `chebyshev_polynomial<>` class to generate the polynomials.
///
/// @author Roland Abel
/// @date October 14, 2023
///
/// Copyright (c) 2026 Roland Abel

#include "chebyshev_polynomial.h"
#include "real_polynomial_root_finder.h"

using namespace std;
using namespace xmath;

namespace {
    using ChebyshevPolynomial = chebyshev_polynomial<double>;
    using RootFinder = real_polynomial_root_finder<double>;
}

auto main() -> int {
    const auto max_order = 10;

    cout << "Chebyshev Polynomials of 1st kind: " << endl;
    for (int n = 0; n <= max_order; ++n) {
        auto T_n = ChebyshevPolynomial::create_1st_kind(n);
        cout << "T_" << n << ": " << T_n << endl;
    }
    cout << endl;

    cout << "Roots of T_10:" << endl;
    auto [roots, multiplicities] = RootFinder::find_roots(ChebyshevPolynomial::create_1st_kind(10));
    for (std::size_t k = 0; k < roots.size(); ++k) {
        cout << "r[" << k << "] = " << roots[k] << endl;
    }

    return 0;
}