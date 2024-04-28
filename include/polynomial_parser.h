/// @file polynomial_parser.h
/// @brief
///
/// @author Roland Abel
/// @date March 02, 2024
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

#pragma once

#include <expected>
#include <variant>
#include <stack>
#include <regex>
#include <ranges>
#include <map>
#include <functional>

#include "polynomial.h"

namespace xmath::parser {
    template<class... Ts>
    struct overloaded : Ts ... {
        using Ts::operator()...;
    };

    enum class parenthesis_t {
        OPENED,
        CLOSED
    };

    enum class error_t {
        UNEXPECTED_TOKEN,
        UNEXPECTED_END,
        EMPTY_EXPRESSION,
        INVALID_IDENTIFIER,
        INVALID_TOKEN,
        INVALID_OPERATOR,
        INVALID_POWER_EXPONENT,
        INVALID_NUMBER,
        DIVISION_BY_ZERO,
        POLYNOMIAL_EXPECTED
    };

    enum class operator_t {
        PLUS,
        MINUS,
        MULTIPLY,
        DIVIDE,
        POWER
    };

    struct end_marker_t {
        bool operator==(end_marker_t const &) const = default;

        bool operator!=(end_marker_t const &) const = default;
    };

    static end_marker_t END{};

    struct variable_t {
        char8_t value;

        bool operator==(variable_t const &) const = default;

        bool operator!=(variable_t const &) const = default;
    };

    struct number_t {
        double_t value;

        bool operator==(number_t const &) const = default;

        bool operator!=(number_t const &) const = default;
    };

    using token_t = std::variant<operator_t, parenthesis_t, variable_t, number_t, end_marker_t>;

    struct scan_state_t {
        const token_t token;
        const uint16_t position{};
    };

    using tokens_t = std::vector<token_t>;
    using polynomial_t = xmath::polynomial<double_t>;
    using item_t = std::variant<polynomial_t, operator_t>;
    using items_t = std::vector<item_t>;
    using stack_t = std::stack<polynomial_t>;

    using character_result_t = std::expected<char8_t, error_t>;
    using scan_result_t = std::expected<scan_state_t, error_t>;
    using tokenize_result_t = std::expected<tokens_t, error_t>;
    using items_result_t = std::expected<items_t, error_t>;
    using polynomial_result_t = std::expected<polynomial_t, error_t>;

    const auto zero = polynomial_t::zero();
    const auto one = polynomial_t::one();
    const auto X = polynomial_t::monomial(1, 1.0);

    static const auto operator_map = std::map<char8_t, operator_t>{
            {'+', operator_t::PLUS},
            {'-', operator_t::MINUS},
            {'*', operator_t::MULTIPLY},
            {'/', operator_t::DIVIDE},
            {'^', operator_t::POWER}
    };

    static const auto parenthesis_map = std::map<char8_t, parenthesis_t>{
            {'(', parenthesis_t::OPENED},
            {')', parenthesis_t::CLOSED}
    };

    inline const auto is_operator = [](const char8_t &character) noexcept {
        return operator_map.contains(character);
    };

    inline const auto is_parenthesis = [](const char8_t &character) noexcept {
        return parenthesis_map.contains(character);
    };

    inline const auto to_operator = [](const char8_t &character) {
        return operator_map.at(character);
    };

    inline const auto to_parenthesis = [](const char8_t &character) {
        return parenthesis_map.at(character);
    };

    const auto get_next_character = [](const std::string &expression, const uint16_t pos = 0) noexcept -> character_result_t {
        return pos < expression.size()
               ? character_result_t{expression[pos]}
               : std::unexpected<error_t>(error_t::UNEXPECTED_END);
    };

    scan_result_t scan_number(const std::string &expression, const uint16_t pos = 0) noexcept {
        using string_result_t = std::expected<std::string, error_t>;
        static std::regex number_regex("^[0-9]*(\\.[0-9]*)?");
        auto get_substring = [&]() {
            return pos < expression.size()
                   ? string_result_t{expression.substr(pos)}
                   : std::unexpected<error_t>(error_t::UNEXPECTED_END);
        };
        auto read_number = [&](const std::string &str) {
            std::smatch base_match;
            return std::regex_search(str, base_match, number_regex)
                   ? string_result_t{base_match[0]}
                   : std::unexpected<error_t>(error_t::INVALID_NUMBER);
        };
        auto make_state = [&](const std::string &number) -> scan_result_t {
            return scan_state_t(number_t{std::stod(number)}, number.length() + pos);
        };
        return get_substring().and_then(read_number)
                .and_then(make_state);
    }

///
/// @param expression
/// @param pos
/// @return
    scan_result_t scan_identifier(const std::string &expression, const uint16_t &pos = 0) noexcept {
        auto is_identifier = [](const auto &character) {
            return std::isalnum(character);
        };
        auto make_state = [&](const auto &character) -> scan_result_t {
            return is_identifier(character)
                   ? scan_result_t{scan_state_t(variable_t{character}, pos + 1)}
                   : std::unexpected<error_t>(error_t::UNEXPECTED_END);
        };
        return get_next_character(expression, pos)
                .and_then(make_state);
    }

    inline scan_result_t scan_operator(const std::string &expression, const uint16_t &pos = 0) noexcept {
        auto make_state = [&](const auto &character) -> scan_result_t {
            return is_operator(character)
                   ? scan_result_t{scan_state_t(to_operator(character), pos + 1)}
                   : std::unexpected{error_t::INVALID_TOKEN};
        };
        return get_next_character(expression, pos)
                .and_then(make_state);
    }

///
/// @param expression
/// @param pos
/// @return
    scan_result_t scan_parenthesis(const std::string &expression, const uint16_t &pos = 0) noexcept {
        auto make_state = [&](const auto &character) -> scan_result_t {
            return is_parenthesis(character)
                   ? scan_result_t{scan_state_t(to_parenthesis(character), pos + 1)}
                   : std::unexpected{error_t::INVALID_TOKEN};
        };

        return get_next_character(expression, pos).and_then(make_state);
    }

///
/// @param expression
/// @param pos
/// @return
    scan_result_t scan_token(const std::string &expression, const uint16_t &pos = 0) noexcept {
        if (expression.empty()) {
            return std::unexpected{error_t::EMPTY_EXPRESSION};
        }
        if (pos > expression.size()) {
            return std::unexpected{error_t::UNEXPECTED_END};
        }
        if (pos == expression.size()) {
            return scan_state_t(END, pos);
        }

        auto is_digit = [](const char8_t &character) noexcept { return isdigit(character); };
        auto is_alpha = [](const char8_t &character) { return isalpha(character); };
        auto is_space = [](const char8_t &character) { return std::isspace(character); };

        auto scan_character = [&](const char8_t &character) -> scan_result_t {
            if (is_space(character)) {
                return scan_token(expression, pos + 1);
            }
            if (is_digit(character)) {
                return scan_number(expression, pos);
            }
            if (is_alpha(character)) {
                return scan_identifier(expression, pos);
            }
            if (is_parenthesis(character)) {
                return scan_parenthesis(expression, pos);
            }
            if (is_operator(character)) {
                return scan_operator(expression, pos);
            }
            return std::unexpected{error_t::INVALID_TOKEN};
        };

        return get_next_character(expression, pos).and_then(scan_character);
    }

///
/// @param expression The expression.
/// @return A vector of the all parsed token.
    tokenize_result_t tokenize(const std::string &expression) noexcept {
        const auto is_end_token = [](const token_t &token) noexcept {
            return std::holds_alternative<end_marker_t>(token);
        };

        tokens_t tokens{};
        const std::function<tokenize_result_t(scan_result_t)> tokenize_ = [&](const auto &result) {
            return result.and_then([&](const auto &s) {
                tokens.push_back(s.token);
                if (is_end_token(s.token)) {
                    return tokenize_result_t{tokens};
                }
                return tokenize_(scan_token(expression, s.position));
            }).or_else([](const auto &error) {
                return tokenize_result_t{std::unexpected{error}};
            });
        };
        return tokenize_(scan_token(expression));
    }

/// Converts the given list of infix token to a list of post fix token by using the shunting-yard algorithm.
/// @param infix_tokens The given collection of infix tokens.
/// @return The list of post fix tokens.
    tokenize_result_t to_postfix_tokens(const tokens_t &infix_tokens) {
        tokens_t postfix_tokens = {};
        std::stack<token_t> operator_stack = {};

        auto precedence = [](const operator_t op) -> uint8_t {
            static const auto precedence_map = std::map<operator_t, uint8_t>{
                    {operator_t::PLUS,     1},
                    {operator_t::MINUS,    1},
                    {operator_t::MULTIPLY, 2},
                    {operator_t::DIVIDE,   2},
                    {operator_t::POWER,    3}
            };
            return precedence_map.contains(op) ? precedence_map.at(op) : 0;
        };

        auto top_precedence_greater_or_equal = [&](const operator_t op) -> bool {
            const auto top_token = operator_stack.top();
            if (!std::holds_alternative<operator_t>(top_token)) {
                return false;
            }

            const auto top_operator = std::get<operator_t>(top_token);
            return precedence(top_operator) >= precedence(op);
        };

        auto top_is_parenthesis = [&]() -> bool {
            return std::holds_alternative<parenthesis_t>(operator_stack.top());
        };

        auto apply_operator = [&](const operator_t &op) {
            while (!operator_stack.empty() && top_precedence_greater_or_equal(op)) {
                postfix_tokens.emplace_back(operator_stack.top());
                operator_stack.pop();
            }
            operator_stack.push(op);
        };

        auto apply_parenthesis = [&](const parenthesis_t &parenthesis) {
            switch (parenthesis) {
                case parenthesis_t::OPENED:
                    operator_stack.push(parenthesis);
                    break;
                case parenthesis_t::CLOSED:
                    while (!operator_stack.empty() && !top_is_parenthesis()) {
                        postfix_tokens.emplace_back(operator_stack.top());
                        operator_stack.pop();
                    }
                    operator_stack.pop();
                    break;
            }
        };

        auto apply_variable = [&](const variable_t &variable) { postfix_tokens.emplace_back(variable); };
        auto apply_number = [&](const number_t &number) { postfix_tokens.emplace_back(number); };
        auto apply_end = [&](const end_marker_t &end_token) {};

        auto apply_for_token = overloaded{
                apply_operator,
                apply_parenthesis,
                apply_variable,
                apply_number,
                apply_end
        };

        for (const auto &token: infix_tokens) {
            std::visit(apply_for_token, token);
        }

        while (!operator_stack.empty()) {
            postfix_tokens.push_back(operator_stack.top());
            operator_stack.pop();
        }
        return postfix_tokens;
    }

///
/// @param op
/// @param left
/// @param right
/// @return
    polynomial_result_t compute(const operator_t op, const polynomial_t &left, const polynomial_t &right) {
        auto plus = [](const auto &p, const auto &q) {
            return p + q;
        };
        auto minus = [](const auto &p, const auto &q) {
            return p - q;
        };
        auto multiply = [](const auto &p, const auto &q) {
            return p * q;
        };
        auto divide = [](const auto &p, const auto &q) {
            return !q.is_zero()
                   ? polynomial_result_t{p / q}
                   : std::unexpected{error_t::DIVISION_BY_ZERO};
        };
        auto pow = [](const auto &p, const auto &q) {
            auto is_non_negative_integer = [](const auto &p) {
                return p.is_constant() && p.is_integer() && p.at(0) >= 0;
            };
            return is_non_negative_integer(q)
                   ? polynomial_result_t{p.pow(static_cast<uint32_t>(q.at(0)))}
                   : std::unexpected{error_t::INVALID_POWER_EXPONENT};
        };

        using func_t = std::function<polynomial_result_t(const polynomial_t &, const polynomial_t &)>;
        static const auto compute_map = std::map<operator_t, func_t>{
                {operator_t::PLUS,     plus},
                {operator_t::MINUS,    minus},
                {operator_t::MULTIPLY, multiply},
                {operator_t::DIVIDE,   divide},
                {operator_t::POWER,    pow}
        };
        return compute_map.at(op)(left, right);
    }

///
/// @param items
/// @return
    polynomial_result_t evaluate_rpn(const items_t &items) {
        using item_result_t = std::expected<item_t, error_t>;
        if (items.empty()) {
            return std::unexpected{error_t::EMPTY_EXPRESSION};
        }

        auto stack = stack_t{};
        auto pop = [&stack]() -> item_result_t {
            if (stack.empty()) {
                return std::unexpected{error_t::EMPTY_EXPRESSION};
            }

            const auto item = stack.top();
            stack.pop();

            return item_t{item};
        };

        auto pop_polynomial = [&]() -> polynomial_result_t {
            return pop().and_then([&](const item_t &item) {
                return std::holds_alternative<polynomial_t>(item)
                       ? polynomial_result_t{std::get<polynomial_t>(item)}
                       : std::unexpected{error_t::POLYNOMIAL_EXPECTED};
            });
        };

        auto push = [&stack](const polynomial_t &polynomial) {
            stack.push(polynomial);
        };

        auto type_of_item = overloaded{
                [&](const polynomial_t &polynomial) { push(polynomial); },
                [&](const operator_t &op) {
                    auto _ = pop_polynomial().transform([&](const auto &right) {
                        pop_polynomial().transform([&](const auto &left) {
                            return compute(op, left, right).transform(push);
                        });
                    });
                }
        };

        for (const item_t &item: items) {
            std::visit(type_of_item, item);
        }
        return pop_polynomial();
    }

///
/// @param tokens
/// @return
    items_result_t to_items(const tokens_t &tokens) {
        auto to_item = [&](const token_t &token) -> item_t {
            return std::visit(overloaded{
                    [](const operator_t &op) { return item_t{op}; },
                    [](const variable_t &) { return item_t{X}; },
                    [](const number_t &number) { return item_t{polynomial_t{number.value}}; },
                    [](auto x) { return item_t{}; }
            }, token);
        };

        return tokens | transform(to_item) | std::ranges::to<std::vector>();
    }

///
/// @param expression
/// @return
    polynomial_result_t parse(const std::string &expression) {
        return tokenize(expression)
                .and_then(to_postfix_tokens)
                .and_then(to_items)
                .and_then(evaluate_rpn);
    }
}
