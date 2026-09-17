/// @file evaluate.cpp
/// @brief Example to demonstrating polynomial creation and evaluation.
///
/// This program creates a 4th-degree polynomial and evaluates it for a range of values.
/// It then prints the result of the polynomial evaluation for each value.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#include <numeric>
#include <vector>
#include <iostream>
#include <xpolynomial.h>

using namespace std;
using namespace xmath;

namespace {
    constexpr double p_at(double x) {
        const auto X = polynomial<double>::monomial(1, 1.0);
        const auto p = 3 * X.pow(4) - 2.5 * X.pow(3) + X.pow(2) - X + 1;
        return p(x);
    }

    static_assert(p_at(0.0) == 1.0);
    static_assert(p_at(1.0) == 1.5);
}

auto main() -> int {
    // Evaluate the polynomial for a range of values.
    vector<double> values(10);
    std::iota(values.begin(), values.end(), 1);

    for (auto x: values) {
        cout << "p(" << x << ") = " << p_at(x) << endl;
    }
    return 0;
}
