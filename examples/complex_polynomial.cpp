/// @file complex_polynomial.cpp
/// @brief Example to demonstrating polynomials with complex coefficients.
///
/// This program creates a complex polynomial, evaluates it, and separates it
/// into its real and imaginary parts using the `complex_polynomial<>` class.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#include <complex>
#include <iostream>
#include <xpolynomial.h>

using namespace std;
using namespace xmath;

namespace {
    using ComplexPolynomial = complex_polynomial<double>;
    constexpr auto i = std::complex(0., 1.);
}

auto main() -> int {
    auto p = ComplexPolynomial({2. - i, 3. + 2. * i, -1. + i});

    cout << "p(x) = " << p << endl
         << "p(1) = " << p(1.) << endl << endl;

    auto [re, im] = separate(p);
    cout << "Re(p) = " << re << endl
         << "Im(p) = " << im << endl;

    return 0;
}