/// @file sturm_sequence.cpp
/// @brief Example to demonstrating polynomial creation, Sturm's' sequence, and counting distinct roots.
///
/// This program creates a polynomial with 5 distinct real roots and performs operations
/// using the `real_polynomial_root_finder<>` class, including generating the Sturm sequence and
/// counting the number of distinct roots.
///
/// @author Roland Abel
/// @date September 17, 2026
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
    // Create a polynomial with 5 distinct read roots.
    auto p = (X + 2.3) * (X + 1.25) * (X - 0.75) * (X - 1.45) * (X - 2.85);

    auto canonical_seq = RootFinder::sturm_sequence(p);
    auto number_roots = RootFinder::number_distinct_roots(p).value();

    cout << "p(x) = " << p << endl << endl
         << "Number of roots: " << number_roots << endl
         << "Canonical polynomial sequence:" << endl;

    for (const auto& q: canonical_seq) {
        cout << q << endl;
    }
    return 0;
}