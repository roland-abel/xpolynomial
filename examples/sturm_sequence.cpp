/// @file sturm_sequence.cpp
/// @brief Determines the number of distinct real roots by using the Sturm's methode.
///
/// @author Roland Abel
/// @date 11.10.2023

#include "real_polynomial_root_finder.h"

using namespace std;
using namespace xmath;

namespace {
    using RootFinder = real_polynomial_root_finder<double>;
    auto X = polynomial<double>::monomial(1, 1.0);
}

int main() {
    // Create a polynomial with 5 distinct read roots.
    auto p = (X + 2.3) * (X + 1.25) * (X - 0.75) * (X - 1.45) * (X - 2.85);

    auto canonical_seq = RootFinder::sturm_sequence(p);
    auto number_roots = RootFinder::number_distinct_roots(p);

    cout << "Polynomial: " << p
         << "Number of roots: " << number_roots << endl
         << "Canonical polynomial sequence:" << endl;

    for (auto q: canonical_seq) {
        cout << q;
    }
    return 0;
}