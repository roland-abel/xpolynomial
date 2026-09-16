/// @file square_free.cpp
/// @brief Example to demonstrating a square-free decomposition of a polynomial.
///
/// @author Roland Abel
/// @date November 28, 2023
///
/// Copyright (c) 2026 Roland Abel

#include <iostream>
#include "square_free_decomposition.h"

using namespace std;
using namespace xmath;

namespace {
    using SquareFree = square_free_decomposition<double>;
    const auto X = polynomial<double>::monomial(1, 1.0);
}

auto main() -> int {
    auto p = X * (X - 4) * (X + 3.).pow(2) * (X.pow(2) + X).pow(3) * (X.pow(2) + X).pow(4);
    cout << "p = " << p << endl;

    auto square_free_seq = SquareFree::yun_algorithm(p).value();
    for (auto k = 0; k < square_free_seq.size(); ++k) {
        cout << "q" << k << " = " << square_free_seq[k] << endl;
    }
    return 0;
}
