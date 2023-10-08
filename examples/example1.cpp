/// @file example1.cpp
/// @brief
///
/// @author Roland Abel
/// @date 08.10.2023

#include <numeric>
#include <vector>
#include <iostream>
#include "polynomial.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    Polynomial X = Polynomial::monomial(1, 1.0);
}

int main() {

    // Create a 4th degree polynomial p = 3X^4 - 2.5X^3 + X^2 - X +1
    auto p = 3 * X.pow(4) - 2.5 * X.pow(3) + X.pow(2) - X + 1;

    // Evaluate the polynomial for a range of values
    std::vector<double> values(10);
    std::iota(values.begin(), values.end(), 1);

    for (auto x: values) {
        std::cout << "p(" << x << ") = " << p(x) << std::endl;
    }
    return 0;
}
