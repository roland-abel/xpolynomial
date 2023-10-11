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
    // Create a polynomial
    auto p = (X + 1 / 3.).pow(3) * (X - 1.5) * (X.pow(2) + X - 2).pow(2);
    auto [roots, multiplicities] = RootFinder::find_roots(p);

    cout << "Polynomial: " << p << endl << endl;
    for (int k = 0; k < roots.size(); ++k) {
        cout << "Root: r[" << k << "] = " << roots[k]
             << ", Multiplicity: " << multiplicities[k] << endl;
    }
    return 0;
}