/// @file legendre.cpp
/// @brief Example to demonstrating the computation of Legendre polynomials.
///
/// This program calculates and prints the Legendre polynomials P_n for a range
/// of orders using the `legendre_polynomial<>` class.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#include <iostream>
#include <xpolynomial.h>

using namespace std;
using namespace xmath;

auto main() -> int {
    const auto max_order = 8;

    cout << "Legendre polynomials P_n:" << endl;
    for (int n = 0; n <= max_order; ++n) {
        cout << "P_" << n << ": " << legendre_polynomial<double>::create(n) << endl;
    }
    return 0;
}