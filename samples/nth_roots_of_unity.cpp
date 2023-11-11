/// @file nth_roots_of_unity.cpp.cpp
/// @brief
///
/// @author Roland Abel
/// @date 28.10.2023.

#include <iostream>
#include "complex_polynomial_root_finder.h"

using namespace std;
using namespace xmath;

namespace {
    using ComplexPolynomial = complex_polynomial<double>;
    using RealPolynomial = polynomial<double, polynomial_specification<double>>;
    using RootFinder = complex_polynomial_root_finder<double>;

    auto X = RealPolynomial::monomial(1, 1.0);
    auto Y = RealPolynomial::monomial(1, 1.0);
    auto Z = ComplexPolynomial::monomial(1, 1.0);
}

int main() {
    auto n = 13;
    auto p = Z.pow(n) - 1.;
    auto roots = RootFinder::nth_roots_of_unity(n);

    cout << "Has roots: " << p.has_roots(roots) << endl;
    cout << n << "th-roots:" << endl;
    for (auto z: roots) {
        cout << z << endl;
    }
    return 0;
}