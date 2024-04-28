/// @file polynomial_parser_tests.cpp
/// @brief
///
/// @author Roland Abel
/// @date March 02, 2023
///
/// Copyright (c) 2024 Roland Abel
///
/// This software is released under the MIT License.
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
/// THE SOFTWARE.

#include <gtest/gtest.h>
#include <vector>

#include "polynomial_parser.h"

using namespace xmath;
using namespace xmath::parser;

namespace {
    constexpr double epsilon = 1e-9;
    auto P = [](double_t value) -> polynomial_t { return polynomial_t{value}; };
}

TEST(PolynomialParserTests, IsOperatorTest) {
    ASSERT_TRUE(is_operator('+'));
    ASSERT_TRUE(is_operator('-'));
    ASSERT_TRUE(is_operator('*'));
    ASSERT_TRUE(is_operator('/'));
    ASSERT_TRUE(is_operator('^'));

    ASSERT_FALSE(is_operator('.'));
    ASSERT_FALSE(is_operator('('));
    ASSERT_FALSE(is_operator(')'));
}

TEST(PolynomialParserTests, GetNextCharacterTest) {
    const auto result = get_next_character("12+34+X", 3);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value(), '3');
}

TEST(PolynomialParserTests, ScanDigitTest) {
    const auto state = scan_number("7");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<number_t>(token).value, 7.0);
}

TEST(PolynomialParserTests, ScanIntegerNumberTest) {
    const auto state = scan_number("57");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 2);
    EXPECT_EQ(std::get<number_t>(token).value, 57.0);
}

TEST(PolynomialParserTests, ScanFloatingPointNumberTest) {
    const auto state = scan_number("X + 12.45 + 3", 4);
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 9);
    EXPECT_EQ(std::get<number_t>(token).value, 12.45);
}

TEST(PolynomialParserTests, ScanIdentifierTest) {
    const auto state = scan_identifier("X");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    // EXPECT_EQ(std::get<variable_t>(token), identifier_marker);
}

TEST(PolynomialParserTests, ScanPlusOperatorTest) {
    const auto state = scan_operator("+");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<operator_t>(token), operator_t::PLUS);
}

TEST(PolynomialParserTests, ScanMinusOperatorTest) {
    const auto state = scan_operator("-");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<operator_t>(token), operator_t::MINUS);
}

TEST(PolynomialParserTests, ScanMultiplyOperatorTest) {
    const auto state = scan_operator("*");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<operator_t>(token), operator_t::MULTIPLY);
}

TEST(PolynomialParserTests, ScanDivideOperatorTest) {
    const auto state = scan_operator("/");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<operator_t>(token), operator_t::DIVIDE);
}

TEST(PolynomialParserTests, ScanOpenParenthesisTest) {
    const auto state = scan_parenthesis("(");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<parenthesis_t>(token), parenthesis_t::OPENED);
}

TEST(PolynomialParserTests, ScanCloseParenthesisTest) {
    const auto state = scan_parenthesis(")");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<parenthesis_t>(token), parenthesis_t::CLOSED);
}

TEST(PolynomialParserTests, ScanTokenEmptyExpressionTest) {
    const auto state = scan_token("");
    ASSERT_FALSE(state.has_value());

    EXPECT_EQ(state.error(), error_t::EMPTY_EXPRESSION);
}

TEST(PolynomialParserTests, ScanTokenInvalidPositionTest) {
    constexpr auto invalid_pos = 25;
    const auto state = scan_token("3 + X", invalid_pos);

    EXPECT_FALSE(state.has_value());
    EXPECT_EQ(state.error(), error_t::UNEXPECTED_END);
}

TEST(PolynomialParserTests, ScanTokenFloatingNumberTest) {
    const auto state = scan_token("3.14");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 4);
    EXPECT_EQ(std::get<number_t>(token).value, 3.14);
}

TEST(PolynomialParserTests, ScanTokenIdentifierTest) {
    const auto state = scan_token("X");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_TRUE(std::get<variable_t>(token).value == 'X');
}

TEST(PolynomialParserTests, ScanTokenCloseParenthesisTest) {
    const auto state = scan_token(")");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<parenthesis_t>(token), parenthesis_t::CLOSED);
}

TEST(PolynomialParserTests, ScanTokenOpenParenthesisTest) {
    const auto state = scan_token("(");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<parenthesis_t>(token), parenthesis_t::OPENED);
}

TEST(PolynomialParserTests, ScanTokenPlusOperatorTest) {
    const auto state = scan_token("+");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<operator_t>(token), operator_t::PLUS);
}

TEST(PolynomialParserTests, ScanTokenMinusOperatorTest) {
    const auto state = scan_token("-");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<operator_t>(token), operator_t::MINUS);
}

TEST(PolynomialParserTests, ScanTokenMultiplyOperatorTest) {
    const auto state = scan_token("*");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<operator_t>(token), operator_t::MULTIPLY);
}

TEST(PolynomialParserTests, ScanTokenDivideOperatorTest) {
    const auto state = scan_token("/");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<operator_t>(token), operator_t::DIVIDE);
}

TEST(PolynomialParserTests, ScanTokenPowerOperatorTest) {
    const auto state = scan_token("^");
    ASSERT_TRUE(state.has_value());

    auto [token, pos] = state.value();

    EXPECT_EQ(pos, 1);
    EXPECT_EQ(std::get<operator_t>(token), operator_t::POWER);
}

TEST(PolynomialParserTests, ScanTokenInvalidTokenTest) {
    const auto state = scan_token("$");
    ASSERT_FALSE(state.has_value());

    EXPECT_EQ(state.error(), error_t::INVALID_TOKEN);
}

TEST(PolynomialParserTests, TokenizeEmptyExpressionTest) {
    const auto result = tokenize("");
    ASSERT_FALSE(result.has_value());

    EXPECT_EQ(result.error(), error_t::EMPTY_EXPRESSION);
}

TEST(PolynomialParserTests, TokenizeInvaildExpressionTest) {
    const auto result = tokenize("$ + 6");
    ASSERT_FALSE(result.has_value());

    EXPECT_EQ(result.error(), error_t::INVALID_TOKEN);
}

TEST(PolynomialParserTests, TokenizeExpressionTest) {
    const auto result = tokenize("3 * (5 + 2.8) - X^5");
    ASSERT_TRUE(result.has_value());

    const auto &tokens = result.value();
    const tokens_t expected_tokens = {
            {number_t{3.0}},
            {operator_t::MULTIPLY},
            {parenthesis_t::OPENED},
            {number_t{5.0}},
            {operator_t::PLUS},
            {number_t{2.8}},
            {parenthesis_t::CLOSED},
            {operator_t::MINUS},
            {variable_t{'X'}},
            {operator_t::POWER},
            {number_t{5.0}},
            {END}
    };

    ASSERT_EQ(tokens.size(), 12);

    for (int i = 0; i < tokens.size(); ++i) {
        EXPECT_EQ(tokens[i], expected_tokens[i]);
    }
}

TEST(PolynomialParserTests, ToPostfixTokensTests) {
    const auto result = tokenize("3 * (5 + 2.8) - X^4").and_then(to_postfix_tokens);
    ASSERT_TRUE(result.has_value());

    const auto &tokens = result.value();
    const tokens_t expected_tokens = {
            {number_t{3.0}},
            {number_t{5.0}},
            {number_t{2.8}},
            {operator_t::PLUS},
            {operator_t::MULTIPLY},
            {variable_t{'X'}},
            {number_t{4.0}},
            {operator_t::POWER},
            {operator_t::MINUS}
    };

    ASSERT_EQ(tokens.size(), expected_tokens.size());

    for (int i = 0; i < tokens.size(); ++i) {
        EXPECT_EQ(tokens[i], expected_tokens[i]);
    }
}

TEST(PolynomialParserTests, ComputePlusTest) {
    const auto p = X.pow(2) + 3;
    const auto q = 3 * X.pow(3) - 4;
    const auto result = compute(operator_t::PLUS, p, q);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value(), p + q);
}

TEST(PolynomialParserTests, ComputeMinusTest) {
    const auto p = X.pow(2) + 3;
    const auto q = 3 * X.pow(3) - 4;
    const auto result = compute(operator_t::MINUS, p, q);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value(), p - q);
}

TEST(PolynomialParserTests, ComputeMultiplyTest) {
    const auto p = X.pow(2) + 3;
    const auto q = 3 * X.pow(3) - 4;
    const auto result = compute(operator_t::MULTIPLY, p, q);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value(), p * q);
}

TEST(PolynomialParserTests, ComputeDivTest) {
    const auto p = (X.pow(2) + 3).pow(2) * (X.pow(2) + 3);
    const auto q = X.pow(2) + 3;
    const auto result = compute(operator_t::DIVIDE, p, q);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value(), p / q);
}

TEST(PolynomialParserTests, ComputeDivByZeroTest) {
    const auto p = (X.pow(2) + 3).pow(2) * (X.pow(2) + 3);
    const auto q = P(0);
    const auto result = compute(operator_t::DIVIDE, p, q);

    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error(), error_t::DIVISION_BY_ZERO);
}

TEST(PolynomialParserTests, ComputePowTest) {
    const auto p = X.pow(2) + 3;
    const auto result = compute(operator_t::POWER, p, P(3));

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value(), p.pow(3));
}

TEST(PolynomialParserTests, ComputeInvaildPowerExpoentTest) {
    const auto p = X.pow(2) + 3;
    const auto result = compute(operator_t::POWER, p, P(3.5));

    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error(), error_t::INVALID_POWER_EXPONENT);
}

TEST(PolynomialParserTests, ComputeRpnEmptyTest) {
    items_t empty = items_t{};
    const auto result = evaluate_rpn(empty);

    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error(), error_t::EMPTY_EXPRESSION);
}

TEST(PolynomialParserTests, ComputeRpnPlusOperatorTest) {
    // "3  X +" == "3 + X"
    items_t input = items_t{P(3.0), X, operator_t::PLUS};
    const auto result = evaluate_rpn(input);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value(), X + 3.0);
}

TEST(PolynomialParserTests, ComputeRpnMinusOperatorTest) {
    // "3  X -" == "3 - X"
    items_t input = items_t{P(3.0), X, operator_t::MINUS};
    const auto result = evaluate_rpn(input);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value(), P(3.0) - X);
}

TEST(PolynomialParserTests, ComputeRpnMultiOperatorTest) {
    // "5  X *" == "5*X"
    items_t input = items_t{P(5.0), X, operator_t::MULTIPLY};
    const auto result = evaluate_rpn(input);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value(), P(5.0) * X);
}

TEST(PolynomialParserTests, ComputeRpnDivOperatorTest) {
    // "X  2. /" == "X/2"
    items_t input = items_t{X, P(2.0), operator_t::DIVIDE};
    const auto result = evaluate_rpn(input);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value(), X / 2.);
}

TEST(PolynomialParserTests, ComputeRpnPowerOperatorTest) {
    // "X  2. ^" == "X^2"
    items_t input = items_t{X, P(2.0), operator_t::POWER};
    const auto result = evaluate_rpn(input);

    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result.value(), X.pow(2));
}

TEST(PolynomialParserTests, ComputeRpnTest) {
    const auto plus = operator_t::PLUS;
    const auto minus = operator_t::MINUS;
    const auto multi = operator_t::MULTIPLY;
    const auto power = operator_t::POWER;

    const std::vector<std::pair<items_t, polynomial_t >> rpn_input = {
            {items_t{P(3.), X, plus},                                                     3. + X},
            {items_t{P(2.), P(3.), X, plus, multi},                                       2. * (3. + X)},
            {items_t{X, P(2.), power},                                                    X.pow(2)},
            {items_t{X, P(2.), power, P(4.), plus, P(3.), X, P(4.), power, minus, multi}, (X.pow(2) + 4.) * (P(3.) - X.pow(4))}
    };

    for (const auto &input: rpn_input) {
        const auto result = evaluate_rpn(input.first);
        ASSERT_TRUE(result.has_value());
        EXPECT_EQ(result, input.second);
    }
}

TEST(PolynomialParserTests, ParseTest) {
    using expression_polynomial_t = std::pair<std::string, polynomial_t>;
    const std::vector<expression_polynomial_t> expected_values = {
            {"0",                               P(0)},
            {"1",                               P(1)},
            {"3.9",                             P(3.9)},
            {"5 - 6",                           -P(1)},
            {"5 * 6",                           P(30)},
            {"X",                               X},
            {"X - X",                           P(0)},
            {"0 * X",                           P(0)},
            {"1 * X",                           X},
            {"2*X",                             2 * X},
            {"2*2*X",                           4 * X},
            {"2*2*X*X * X *X",                  4 * X.pow(4)},
            {"X + 3",                           X + 3},
            {"X - 2.5",                         X - 2.5},
            {"3*X",                             3 * X},
            {"3*X + 5",                         3 * X + 5},
            {"3*X + 3",                         3 * X + 3},
            {"(4 + X^4) * (X^3 - 2)",           (4 + X.pow(4)) * (X.pow(3) - 2)},
            {"(4 + X^4)^2 * ((X^3 - 2)^2) - 1", (4 + X.pow(4)).pow(2) * ((X.pow(3) - 2).pow(2)) - 1},
            {"X^9 / X^2",                       X.pow(7)},
            // {"-1", -P(1)},
    };

    for (const auto &[expr, p]: expected_values) {
        const auto &result = parse(expr);

        ASSERT_TRUE(result.has_value());
        EXPECT_EQ(result.value(), p);
    }
}
