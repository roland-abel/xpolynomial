/// @file euclidean.cpp
/// @brief Computes the greatest common divisor (gcd) of two polynomials and two polynomials s and t such that
/// the gcd is given by s * p + t * q.
///
/// @author Roland Abel
/// @date November 28, 2023
///
/// Copyright (c) 2026 Roland Abel

#include <iostream>
#include "euclidean_algorithm.h"

using namespace std;
using namespace xmath;

namespace {
    using Euclidean = euclidean_algorithm<double>;
    const auto X = polynomial<double>::monomial(1, 1.0);
}

auto main() -> int {
    auto p = X.pow(4) - 2 * X.pow(3) - 6 * X.pow(2) + 12 * X + 15;
    auto q = X.pow(3) + X.pow(2) - 4 * X - 4;

    // s, t, g such that g = gcd(p, q) = s*p + t*q
    auto [s, t, g] = Euclidean::extended_euclidean(p, q);

    cout << "p = " << p.to_string() << endl
         << "q = " << q.to_string() << endl << endl
         << "g = gcd(p, q) = " << g << endl
         << "s = " << s << endl
         << "t = " << t << endl
         << "g = s*p + t*q = " << s * p + t * q << endl;

    return 0;
}