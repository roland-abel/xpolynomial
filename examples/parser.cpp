/// @file parser.tpp
/// @brief
///
/// @author Roland Abel
/// @date March 02, 2023
///
/// Copyright (c) 2026 Roland Abel

#include <iostream>
#include "polynomial_parser.h"

using namespace std;
using namespace xmath;
using namespace xmath::parser;

auto main() -> int {
    parse_polynomial("-(X^3 - 5*X^2 + 4*X)^2 + 6*X^2").transform([](const polynomial_t &p) {
        cout << p << endl;
    });
    return 0;
}

