/// @file polynomial_parser.h
/// @brief Provides a parser for polynomial expressions given as strings.
///
/// @author Roland Abel
/// @date March 02, 2024
///
/// Copyright (c) 2026 Roland Abel

#ifndef POLYNOMIAL_PARSER_H_
#define POLYNOMIAL_PARSER_H_

#include <expected>
#include <variant>
#include <stack>
#include <regex>
#include <ranges>
#include <map>
#include <functional>

#include "polynomial.h"

namespace xmath::parser {
    using std::views::transform;

    /// @brief Helper for creating an overloaded visitor from multiple lambdas.
    /// @tparam Ts The types of the callables.
    template<class... Ts>
    struct overloaded : Ts ... {
        using Ts::operator()...;
    };

    /// @brief The supported binary and unary operators.
    enum class operator_t {
        PLUS,
        MINUS,
        MULTIPLY,
        DIVIDE,
        POWER,
        SIGN_MINUS,
        SIGN_PLUS
    };

    /// @brief The opening and closing parenthesis.
    enum class parenthesis_t {
        OPENED,
        CLOSED
    };

    /// @brief The possible parser errors.
    enum class error_t {
        UNEXPECTED_END,
        EMPTY_EXPRESSION,
        INVALID_VARIABLE,
        INVALID_TOKEN,
        INVALID_OPERATOR,
        INVALID_POWER_EXPONENT,
        INVALID_NUMBER,
        DIVISION_BY_ZERO,
        OPERAND_EXPECTED
    };

    /// @brief Marks the end of the token list.
    struct end_marker_t {
        bool operator==(end_marker_t const &) const = default;

        bool operator!=(end_marker_t const &) const = default;
    };

    /// The end marker token.
    inline constexpr end_marker_t END{};

    /// @brief A variable token, e.g. 'X'.
    struct variable_t {
        char8_t value;

        bool operator==(variable_t const &) const = default;

        bool operator!=(variable_t const &) const = default;
    };

    /// @brief A number token.
    struct number_t {
        double_t value;

        bool operator==(number_t const &) const = default;

        bool operator!=(number_t const &) const = default;
    };

    /// @brief The type of a single token.
    using token_t = std::variant<operator_t, parenthesis_t, variable_t, number_t, end_marker_t>;

    /// @brief Pairs a scanned token with the position of the next character.
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

    /// The monomial X used by the parser.
    inline const auto X = polynomial_t::monomial(1, 1.0);

    /// Maps a character to the corresponding operator.
    inline const auto operator_map = std::map<char8_t, operator_t>{
            {'+', operator_t::PLUS},
            {'-', operator_t::MINUS},
            {'*', operator_t::MULTIPLY},
            {'/', operator_t::DIVIDE},
            {'^', operator_t::POWER}
    };

    /// Maps a character to the corresponding parenthesis.
    inline const auto parenthesis_map = std::map<char8_t, parenthesis_t>{
            {'(', parenthesis_t::OPENED},
            {')', parenthesis_t::CLOSED}
    };

    /// Maps a plus/minus operator to the corresponding sign operator.
    inline const auto sign_operator_map = std::map<operator_t, operator_t>{
            {operator_t::PLUS,  operator_t::SIGN_PLUS},
            {operator_t::MINUS, operator_t::SIGN_MINUS}
    };

    /// Checks if the given character is an operator.
    inline const auto is_operator = [](const char8_t &ch) noexcept { return operator_map.contains(ch); };
    /// Checks if the given character is a parenthesis.
    inline const auto is_parenthesis = [](const char8_t &ch) noexcept { return parenthesis_map.contains(ch); };
    /// Converts the given character to the corresponding operator.
    inline const auto to_operator = [](const char8_t &ch) { return operator_map.at(ch); };
    /// Converts the given character to the corresponding parenthesis.
    inline const auto to_parenthesis = [](const char8_t &ch) { return parenthesis_map.at(ch); };

    /// Gets the next character of the expression at the given position.
    /// @param expression The expression to read from.
    /// @param pos The position of the character. Defaults to 0.
    /// @return The character or an error if the end of the expression is reached.
    const auto get_next_character = [](const std::string &expression, const uint16_t pos = 0) noexcept {
        return pos < expression.size()
               ? character_result_t{expression[pos]}
               : std::unexpected<error_t>(error_t::UNEXPECTED_END);
    };

    /// Scans a number from the given expression starting from the specified position.
    /// @param expression The expression to scan.
    /// @param pos The starting position from which to scan the number. Defaults to 0.
    /// @return The result of scanning the number.
    inline scan_result_t scan_number(const std::string &expression, const uint16_t pos = 0) noexcept {
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
            return !number.empty()
                   ? scan_result_t{scan_state_t(number_t{std::stod(number)}, number.length() + pos)}
                   : std::unexpected<error_t>(error_t::INVALID_NUMBER);
        };
        return get_substring().and_then(read_number)
                .and_then(make_state);
    }

    /// Scans a variable from the given expression starting from the specified position.
    /// @param expression The expression to scan.
    /// @param pos The starting position from which to scan the variable. Defaults to 0.
    /// @param variable The variable used in the expression. Defaults to 'X'.
    /// @return The result of scanning the variable.
    inline scan_result_t scan_variable(const std::string &expression, const uint16_t &pos = 0, const char8_t &variable = 'X') noexcept {
        auto is_variable = [&](const auto &character) {
            return std::isalnum(character) && variable == character;
        };
        auto make_state = [&](const auto &character) -> scan_result_t {
            return is_variable(character)
                   ? scan_result_t{scan_state_t(variable_t{variable}, pos + 1)}
                   : std::unexpected<error_t>(error_t::INVALID_VARIABLE);
        };
        return get_next_character(expression, pos)
                .and_then(make_state);
    }

    /// Scans an operator from the given expression starting from the specified position.
    /// @param expression The expression to scan.
    /// @param pos The starting position from which to scan the operator. Defaults to 0.
    /// @return The result of scanning the operator.
    inline scan_result_t scan_operator(const std::string &expression, const uint16_t &pos = 0) noexcept {
        auto make_state = [&](const auto &character) -> scan_result_t {
            return is_operator(character)
                   ? scan_result_t{scan_state_t(to_operator(character), pos + 1)}
                   : std::unexpected{error_t::INVALID_TOKEN};
        };
        return get_next_character(expression, pos)
                .and_then(make_state);
    }

    /// Scans a parenthesis from the given expression starting from the specified position.
    /// @param expression The expression to scan.
    /// @param pos The starting position from which to scan the parenthesis. Defaults to 0.
    /// @return The result of scanning the parenthesis.
    inline scan_result_t scan_parenthesis(const std::string &expression, const uint16_t &pos = 0) noexcept {
        auto make_state = [&](const auto &character) -> scan_result_t {
            return is_parenthesis(character)
                   ? scan_result_t{scan_state_t(to_parenthesis(character), pos + 1)}
                   : std::unexpected{error_t::INVALID_TOKEN};
        };

        return get_next_character(expression, pos).and_then(make_state);
    }

    /// Scans the token from the given expression starting from the specified position.
    /// @param expression The expression to scan.
    /// @param pos The starting position from which to scan the token. Defaults to 0.
    /// @param variable The variable used in the expression. Defaults to 'X'.
    /// @return The result of scanning the token.
    inline scan_result_t scan_token(const std::string &expression, const uint16_t &pos = 0, const char8_t &variable = 'X') noexcept {
        if (expression.empty()) {
            return std::unexpected{error_t::EMPTY_EXPRESSION};
        }
        if (pos > expression.size()) {
            return std::unexpected{error_t::UNEXPECTED_END};
        }
        if (pos == expression.size()) {
            return scan_state_t(END, pos);
        }

        auto is_digit = [](const char8_t &character) { return isdigit(character); };
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
                return scan_variable(expression, pos, variable);
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

    /// Tokenizes the given expression.
    /// @param expression The expression to tokenize.
    /// @param variable The variable used in the expression. Defaults to 'X'.
    /// @return The result of tokenizing the expression.
    inline tokenize_result_t tokenize(const std::string &expression, const char8_t &variable = 'X') noexcept {
        const auto is_end_token = [](const token_t &token) noexcept {
            return std::holds_alternative<end_marker_t>(token);
        };

        tokens_t tokens{};
        const std::function<tokenize_result_t(scan_result_t)> tokenize_ = [&](const scan_result_t &result) {
            return result.and_then([&](const auto &s) {
                tokens.emplace_back(s.token);
                if (is_end_token(s.token)) {
                    return tokenize_result_t{tokens};
                }
                return tokenize_(scan_token(expression, s.position, variable));
            }).or_else([](const auto &error) {
                return tokenize_result_t{std::unexpected{error}};
            });
        };

        return tokenize_(scan_token(expression));
    }

    /// Converts the tokens so that the PLUS and the MINUS operators are converted to sign operators.
    /// @param tokens The tokens to transform.
    /// @return The transformed tokens.
    inline tokenize_result_t convert_tokens_with_signs(const tokens_t &tokens) noexcept {
        const auto is_operator_ = [](const token_t &token, const operator_t op) noexcept {
            return std::holds_alternative<operator_t>(token) && (std::get<operator_t>(token) == op);
        };
        const auto is_open_parenthesis = [](const token_t &token) noexcept {
            return std::holds_alternative<parenthesis_t>(token) && std::get<parenthesis_t>(token) == parenthesis_t::OPENED;
        };
        const auto is_plus_or_minus_operator = [&](const token_t &token) noexcept {
            return is_operator_(token, operator_t::PLUS) || is_operator_(token, operator_t::MINUS);
        };

        bool last_token_was_operator_or_open_parenthesis = true;

        auto to_sign_token = [&](const token_t &token) {
            if (is_plus_or_minus_operator(token) && last_token_was_operator_or_open_parenthesis) {
                const auto op = std::get<operator_t>(token);
                return token_t{operator_t{sign_operator_map.at(op)}};
            }
            return token;
        };

        auto to_sign_operator_ = [&](const token_t &token) {
            const auto transformed_token = to_sign_token(token);

            last_token_was_operator_or_open_parenthesis = is_plus_or_minus_operator(token) || is_open_parenthesis(token);
            return transformed_token;
        };

        return tokens | transform(to_sign_operator_) | std::ranges::to<std::vector>();
    }

    /// Converts the given infix token list to a post fix token list by using the shunting-yard algorithm.
    /// @param infix The infix collection to convert.
    /// @return The list of post fix tokens.
    inline tokenize_result_t convert_to_postfix(const tokens_t &infix) {
        tokens_t postfix = {};
        std::stack<token_t> operator_stack = {};

        // Gets the precedence of the given operator. The higher the number, the more priority it has.
        auto precedence = [](const operator_t op) -> uint8_t {
            static const auto precedence_map = std::map<operator_t, uint8_t>{
                    {operator_t::PLUS,       1},
                    {operator_t::MINUS,      1},
                    {operator_t::MULTIPLY,   2},
                    {operator_t::DIVIDE,     2},
                    {operator_t::SIGN_MINUS, 3},
                    {operator_t::SIGN_PLUS,  3},
                    {operator_t::POWER,      4}
            };
            return precedence_map.contains(op) ? precedence_map.at(op) : 0;
        };

        auto is_operator_ = [&](const token_t &token) { return std::holds_alternative<operator_t>(token); };
        auto is_parenthesis_ = [&](const token_t &token) { return std::holds_alternative<parenthesis_t>(token); };

        /// Gets `true` if the top item is a parenthesis; otherwise `false`. Only the "(" can be on the top of the stack.
        auto top_is_parenthesis = [&]() -> bool { return is_parenthesis_(operator_stack.top()); };

        // Gets `true` if the top item is an operator and its precedence is greater than that of the given operators.
        auto top_precedence_greater_or_equal = [&](const operator_t op) -> bool {
            const auto top_token = operator_stack.top();
            if (!is_operator_(top_token)) {
                return false;
            }
            const auto top_operator = std::get<operator_t>(top_token);
            return precedence(top_operator) >= precedence(op);
        };

        auto move_operator_stack_to_postfix = [&](const auto &condition) {
            // Moves the elements from the top of the stack to the post fiy queue until the condition is fulfilled.
            while (condition()) {
                postfix.emplace_back(operator_stack.top());
                operator_stack.pop();
            }
        };

        auto process_operator = [&](const operator_t &op) {
            move_operator_stack_to_postfix([&]() {
                return !operator_stack.empty() && top_precedence_greater_or_equal(op);
            });
            operator_stack.emplace(op);
        };

        auto process_parenthesis = [&](const parenthesis_t &parenthesis) {
            switch (parenthesis) {
                case parenthesis_t::OPENED:
                    operator_stack.emplace(parenthesis);
                    break;
                case parenthesis_t::CLOSED:
                    move_operator_stack_to_postfix([&]() {
                        return !operator_stack.empty() && !top_is_parenthesis();
                    });
                    operator_stack.pop();
                    break;
            }
        };

        auto process_variable = [&](const variable_t &variable) { postfix.emplace_back(variable); };
        auto process_number = [&](const number_t &number) { postfix.emplace_back(number); };
        auto process_end = [&](const end_marker_t &) {};

        auto process_token = [&](const token_t &token) {
            std::visit(overloaded{
                    process_operator,
                    process_parenthesis,
                    process_variable,
                    process_number,
                    process_end
            }, token);
        };

        std::for_each(infix.begin(), infix.end(), process_token);
        move_operator_stack_to_postfix([&]() {
            return !operator_stack.empty();
        });

        return postfix;
    }

    /// Applies a binary operator to two polynomials.
    /// @param op The binary operator to apply.
    /// @param left_operand The left operand polynomial.
    /// @param right_operand The right operand polynomial.
    /// @return The result of applying the binary operator to the polynomials.
    inline polynomial_result_t apply_binary_operator(
            const operator_t op,
            const polynomial_t &left_operand,
            const polynomial_t &right_operand) {

        auto plus = [](const auto &p, const auto &q) { return p + q; };
        auto minus = [](const auto &p, const auto &q) { return p - q; };
        auto multiply = [](const auto &p, const auto &q) { return p * q; };
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

        return compute_map.at(op)(left_operand, right_operand);
    }

    /// Applies a unary operator to a polynomial.
    /// @param op The unary operator to apply.
    /// @param operand The operand polynomial.
    /// @return The result of applying the unary operator to the polynomial.
    inline polynomial_result_t apply_unary_operator(const operator_t op, const polynomial_t &operand) {
        auto sign_plus = [](const auto &p) { return p; };
        auto sign_minus = [](const auto &p) { return (-1) * p; };

        using func_t = std::function<polynomial_result_t(const polynomial_t &)>;
        static const auto compute_map = std::map<operator_t, func_t>{
                {operator_t::SIGN_PLUS,  sign_plus},
                {operator_t::SIGN_MINUS, sign_minus}
        };

        const auto is_sign_operator = [&](const operator_t &op) noexcept {
            return op == operator_t::SIGN_MINUS || op == operator_t::SIGN_PLUS;
        };

        return is_sign_operator(op)
               ? compute_map.at(op)(operand)
               : std::unexpected{error_t::INVALID_OPERATOR};
    }

    /// Evaluates a postfix expression represented by items.
    /// @param postfix The postfix expression to evaluate.
    /// @return The result of evaluating the postfix expression.
    inline polynomial_result_t evaluate(const items_t &postfix) {
        using item_result_t = std::expected<item_t, error_t>;

        if (postfix.empty()) {
            return std::unexpected{error_t::EMPTY_EXPRESSION};
        }

        auto stack = stack_t{};

        auto pop = [&]() -> item_result_t {
            if (stack.empty()) {
                return std::unexpected{error_t::UNEXPECTED_END};
            }
            const auto item = stack.top();
            stack.pop();

            return item_t{item};
        };
        auto pop_polynomial = [&]() {
            auto is_polynomial = [](const item_t &item) {
                return std::holds_alternative<polynomial_t>(item);
            };
            return pop().and_then([&](const item_t &item) {
                return is_polynomial(item)
                       ? polynomial_result_t{std::get<polynomial_t>(item)}
                       : std::unexpected{error_t::OPERAND_EXPECTED};
            }).or_else([](const auto &) {
                return polynomial_result_t{std::unexpected{error_t::OPERAND_EXPECTED}};
            });
        };
        auto push_polynomial = [&](const polynomial_t &polynomial) -> polynomial_result_t {
            stack.push(polynomial);
            return polynomial;
        };
        auto is_unary_operator = [](const operator_t &op) {
            return op == operator_t::SIGN_PLUS || op == operator_t::SIGN_MINUS;
        };
        auto apply_operator = [&](const operator_t &op) {
            return pop_polynomial().and_then([&](const auto &right) {
                if (is_unary_operator(op)) {
                    return apply_unary_operator(op, right)
                            .and_then(push_polynomial);
                } else {
                    return pop_polynomial().and_then([&](const auto &left) {
                        return apply_binary_operator(op, left, right)
                                .and_then(push_polynomial);
                    });
                }
            });
        };

        auto apply_item_on_stack = [&](const item_t &item) -> polynomial_result_t {
            return std::visit(overloaded{
                    [&](const polynomial_t &p) { return push_polynomial(p); },
                    [&](const operator_t &op) { return apply_operator(op); }
            }, item);
        };

        for (const auto &result: postfix | transform(apply_item_on_stack)) {
            if (!result.has_value()) {
                return std::unexpected<error_t>{result.error()};
            }
        }

        return pop_polynomial();
    }

    /// Transforms each token form the given list to an item which contains either a polynomial or an operator.
    /// @param tokens The token vector.
    /// @return The list of items (either a polynomial or an operator).
    inline items_result_t convert_to_items(const tokens_t &tokens) {
        auto to_item = [&](const token_t &token) -> item_t {
            return std::visit(overloaded{
                    [](const operator_t &op) { return item_t{op}; },
                    [](const variable_t &) { return item_t{X}; },
                    [](const number_t &number) { return item_t{polynomial_t{number.value}}; },
                    [](auto) { return item_t{}; }
            }, token);
        };

        return tokens | transform(to_item) | std::ranges::to<std::vector>();
    }

    /// Parses a polynomial expression given as a string.
    /// @param expression The polynomial expression to parse.
    /// @param variable The variable used in the polynomial expression.
    /// @return The result of the polynomial expression.
    inline polynomial_result_t parse_polynomial(const std::string &expression, const char8_t &variable = 'X') {
        return tokenize(expression, variable)
                .and_then(convert_tokens_with_signs)
                .and_then(convert_to_postfix)
                .and_then(convert_to_items)
                .and_then(evaluate);
    }
}

#endif // POLYNOMIAL_PARSER_H_
