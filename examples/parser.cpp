/// @file parser.cpp
/// @brief
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#include <iostream>
#include <xpolynomial.h>

using namespace std;
using namespace xmath;
using namespace xmath::parser;

auto main() -> int {
    parse_polynomial("-(X^3 - 5*X^2 + 4*X)^2 + 6*X^2").transform([](const polynomial_t &p) {
        cout << p << endl;
    });
    return 0;
}

