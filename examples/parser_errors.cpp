/// @file parser_errors.cpp
/// @brief Example to demonstrating the error handling of the polynomial parser.
///
/// This program parses several expressions, including invalid ones, and prints
/// the descriptive `error_t` reported by `parser::parse_polynomial`.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#include <iostream>
#include <xpolynomial.h>

using namespace std;
using namespace xmath;

namespace {
    auto error_message(xmath::parser::error_t error) -> const char * {
        switch (error) {
            case xmath::parser::error_t::UNEXPECTED_END:
                return "unexpected end";
            case xmath::parser::error_t::EMPTY_EXPRESSION:
                return "empty expression";
            case xmath::parser::error_t::INVALID_VARIABLE:
                return "invalid variable";
            case xmath::parser::error_t::INVALID_TOKEN:
                return "invalid token";
            case xmath::parser::error_t::INVALID_OPERATOR:
                return "invalid operator";
            case xmath::parser::error_t::INVALID_POWER_EXPONENT:
                return "invalid power exponent";
            case xmath::parser::error_t::INVALID_NUMBER:
                return "invalid number";
            case xmath::parser::error_t::DIVISION_BY_ZERO:
                return "division by zero";
            case xmath::parser::error_t::OPERAND_EXPECTED:
                return "operand expected";
        }
        return "unknown error";
    }
}

auto main() -> int {
    for (const auto *expression: {"2*X^3 - 3*X + 1", "", "X^^2", "1/(X - X)", "a + 1", "3/0"}) {
        auto result = xmath::parser::parse_polynomial(expression);

        if (result.has_value()) {
            cout << "'" << expression << "' -> p(x) = " << result.value() << endl;
        } else {
            cout << "'" << expression << "' -> error: " << error_message(result.error()) << endl;
        }
    }
    return 0;
}