/// @file polynomial_tests.cpp
/// @brief
///
/// @author Roland Abel
/// @date October 8, 2023
///
/// Copyright (c) 2023 Roland Abel
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
#include <cmath>
#include <numbers>
#include "polynomial.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using values_type = Polynomial::values_type;

    constexpr double epsilon = Polynomial::epsilon;
    const auto zero = Polynomial::zero();
    const auto one = Polynomial::one();
    const auto X = Polynomial::monomial(1, 1.0);

    void check_polynomial_division(
        const Polynomial &p,
        const Polynomial &q,
        const Polynomial &s,
        const Polynomial &r) {
        auto [s_prim, r_prim] = p.divide(q);
        EXPECT_EQ(s_prim, s);
        EXPECT_EQ(r_prim, r);
        EXPECT_EQ(p, s_prim * q + r);
    }
}

// Tests for the default constructor.
TEST(PolynomialTests, DefaultConstructor) {
    const auto p = Polynomial();
    EXPECT_EQ(p.degree(), 0);
}

// Tests for the copy constructor.
TEST(PolynomialTests, CopyConstructor) {
    const auto p = Polynomial({1, 2, 3});
    const auto q = Polynomial(p);

    EXPECT_TRUE(p == q);
}

// Tests for the constructor with coefficients.
TEST(PolynomialTests, ConstructorWithCoefficients) {
    EXPECT_EQ(Polynomial({0}).degree(), 0);
    EXPECT_EQ(Polynomial({0, 0, 0}).degree(), 0);
    EXPECT_EQ(Polynomial({0, 0, 1}).degree(), 2); // x^2
    EXPECT_EQ(Polynomial({0, 1, -1}).degree(), 2); // x - x^2
    EXPECT_EQ(Polynomial({1, 2, 1}).degree(), 2); // 1 + 2x + x^2
}

// Tests for the zero polynomial.
TEST(PolynomialTests, ZeroPloynomialTest) {
    EXPECT_TRUE(Polynomial::zero().is_zero());
    EXPECT_TRUE(Polynomial().is_zero());
    EXPECT_TRUE(Polynomial({0, 0, 0, 0}).is_zero());
    EXPECT_FALSE(Polynomial({1}).is_zero());

    EXPECT_EQ(Polynomial::zero(), zero);
    EXPECT_EQ(Polynomial::zero(), Polynomial());
    EXPECT_EQ(Polynomial::zero(), Polynomial({0, 0, 0, 0}));
    EXPECT_EQ(Polynomial::zero().degree(), 0);
}

// Tests for the constant 1 polynomial.
TEST(PolynomialTests, OnePloynomialTest) {
    EXPECT_TRUE(Polynomial::one().is_one());
    EXPECT_TRUE(Polynomial::one().is_constant());
    EXPECT_TRUE(Polynomial::one().is_linear());
    EXPECT_EQ(Polynomial::one(), one);
    EXPECT_EQ(Polynomial::one(), Polynomial({1}));
    EXPECT_EQ(Polynomial::one().degree(), 0);
}

// Tests for constant polynomials (degree = 0).
TEST(PolynomialTests, ConstantPolynomialTest) {
    EXPECT_TRUE(Polynomial({0}).is_constant());
    EXPECT_TRUE(Polynomial({1}).is_constant());
    EXPECT_TRUE(Polynomial({2, 0, 0}).is_constant());

    EXPECT_EQ(Polynomial({2}).degree(), 0);
}

// Tests for linear polynomials (degree = 1).
TEST(PolynomialTests, LinearPolynomialTest) {
    EXPECT_TRUE(Polynomial({0, 1}).is_linear());
    EXPECT_TRUE(Polynomial({-1, 1}).is_linear());

    EXPECT_FALSE(Polynomial({1, 1, 2}).is_linear());
}

// Tests for quadratic polynomials (degree = 2).
TEST(PolynomialTests, QuadraticPolynomialTest) {
    EXPECT_TRUE(Polynomial({0, 0, 1}).is_quadratic()); // x^2
    EXPECT_TRUE(Polynomial({0, 1, 1}).is_quadratic()); // x + x^2
    EXPECT_TRUE(Polynomial({1, 2, -1}).is_quadratic()); // 1 + 2x - x^2

    EXPECT_FALSE(Polynomial({0, 0, 1, 2, -1}).is_quadratic());
    EXPECT_FALSE(Polynomial({3, 0, 0, 1, 2, -1}).is_quadratic());
    EXPECT_FALSE(Polynomial({2, -1}).is_quadratic());
}

TEST(PolynomialTests, CubicPolynomialTest) {
    EXPECT_TRUE(Polynomial({0, 0, 1, 1}).is_cubic()); // x^3 + x^2
    EXPECT_FALSE(Polynomial({0, 0, 1, 2, -1}).is_cubic());
}

// Tests for the leading coefficient.
TEST(PolynomialTests, LeadingCoefficient) {
    EXPECT_NEAR(Polynomial({1, 0, 0}).leading_coefficient(), 1.0, epsilon);
    EXPECT_NEAR(Polynomial({3.0, 1.0, 2.0}).leading_coefficient(), 2.0, epsilon);
    EXPECT_NEAR(Polynomial({0, 0, 3.5, 1.0, 4.0}).leading_coefficient(), 4.0, epsilon);
}

// Tests for monomial.
TEST(PolynomialTests, MonomialTest) {
    EXPECT_EQ(Polynomial::monomial(1, 1.0), X);
    EXPECT_EQ(Polynomial::monomial(1, 2.0), 2 * X);
    EXPECT_EQ(Polynomial::monomial(4, 3.5), 3.5 * X.pow(4));
    EXPECT_EQ(Polynomial::monomial(10, -2.5), -2.5 * X.pow(10));
    EXPECT_EQ(Polynomial::monomial(10, 3.5).degree(), 10);
}

TEST(PolynomialTests, ToStringTest) {
    EXPECT_EQ(zero.to_string(), "0");
    EXPECT_EQ(one.to_string(), "1");
    EXPECT_EQ((-one).to_string(), "-1");
    EXPECT_EQ(X.to_string(), "x");
    EXPECT_EQ((-X).to_string(), "-x");
    EXPECT_EQ((-X.pow(12) - 1).to_string(), "-x^12 - 1");
    EXPECT_EQ((-X.pow(12) + X - 1).to_string(), "-x^12 + x - 1");
    EXPECT_EQ((X.pow(2)).to_string(), "x^2");
    EXPECT_EQ((-X.pow(3)).to_string(), "-x^3");
    EXPECT_EQ((-X.pow(3) + X.pow(2)).to_string(), "-x^3 + x^2");
    EXPECT_EQ((-X.pow(3) - 2.4 * X.pow(2)).to_string(), "-x^3 - 2.4x^2");
    EXPECT_EQ((-3.2 * X.pow(6) - 1.4 * X.pow(2) - 1).to_string(), "-3.2x^6 - 1.4x^2 - 1");
    EXPECT_EQ((-(1. / 3.) * X.pow(6) - 1.456 * X.pow(2) - (1. / 4.)).to_string(), "-0.333333x^6 - 1.456x^2 - 0.25");
}

// Tests for the equal operator.
TEST(PolynomialTests, EqualOperatorTest) {
    const auto p = Polynomial({1, 2, 3});

    EXPECT_TRUE(p == Polynomial({1, 2, 3}));
    EXPECT_TRUE(p == Polynomial({1, 2, 3, 0, 0}));
    EXPECT_FALSE(p == Polynomial({0, 2, 3, 0, 0}));
}

// Tests for the not equal operator.
TEST(PolynomialTests, NotEqualOperatorTest) {
    const auto p = Polynomial({1, 2, 3, 4});

    EXPECT_FALSE(p != Polynomial({1, 2, 3, 4}));
    EXPECT_FALSE(p != Polynomial({1, 2, 3, 4, 0, 0}));
    EXPECT_TRUE(p != Polynomial({0, 2, 3, 4}));
}

TEST(PolynomialTests, UnitaryPlusOperatorTest) {
    const auto p = Polynomial({1, -2, 3, -4});
    EXPECT_EQ(+p, p);
}

TEST(PolynomialTests, UnitaryMinusOperatorTest) {
    const auto p = Polynomial({-1, 2, -3, 4});
    EXPECT_EQ(-p, Polynomial({1, -2, 3, -4}));
    EXPECT_EQ(-(-p), p);
}

// Tests for the index operator.
TEST(PolynomialTests, IndexOperatorTest) {
    const auto p = Polynomial({0.5, 2.3, 4, -1, 3.6});

    EXPECT_NEAR(p[4], 3.6, epsilon);
    EXPECT_NEAR(p[3], -1, epsilon);
    EXPECT_NEAR(p[2], 4, epsilon);
    EXPECT_NEAR(p[1], 2.3, epsilon);
    EXPECT_NEAR(p[0], 0.5, epsilon);
}

// Tests for the power function.
TEST(PolynomialTests, PowerFunctionTest) {
    const auto p = Polynomial({0, 1});

    EXPECT_EQ(p.pow(0), one);
    EXPECT_EQ(p.pow(1), p);
    EXPECT_EQ(p.pow(2), X * X);
    EXPECT_EQ(p.pow(3), X * X * X);
    EXPECT_EQ(p.pow(4), X * X * X * X);

    EXPECT_EQ((1 + X).pow(2), Polynomial({1, 2, 1}));
    EXPECT_EQ(1 + 2 * X + X.pow(2), Polynomial({1, 2, 1}));
    EXPECT_EQ((2 + X.pow(2)).pow(3), Polynomial({8, 0, 12, 0, 6, 0, 1}));
    EXPECT_EQ((2 + X.pow(2)).pow(3), X.pow(6) + 6 * X.pow(4) + 12 * X.pow(2) + 8);
}

// Tests for the index operator returning zero when out of range.
TEST(PolynomialTests, IndexOperatorGetZeroWhenOutOfRangeTest) {
    const auto p = Polynomial({0.5, 2.3, 4, -1, 3.6});
    EXPECT_NEAR(p[10], 0.0, epsilon);
}

// Tests polynomial addition with a scalar.
TEST(PolynomialTests, AdditionWithScalarTest) {
    const auto p = Polynomial({1, 2, 3, 4}); // 1 + 2x + 3x^2 + 4x^3
    EXPECT_EQ((p + 3.5), Polynomial({4.5, 2, 3, 4}));

    const auto q = Polynomial({1, 2, 3, 0});
    EXPECT_EQ(q + 3.5, Polynomial({4.5, 2, 3}));
    EXPECT_EQ(q, Polynomial({1, 2, 3, 0}));
}

// Tests polynomial subtraction with a scalar.
TEST(PolynomialTests, SubstrationWithScalarTest) {
    const auto p = Polynomial({1, 2, 3, 4}); // 1 + 2x + 3x^2 + 4x^3
    EXPECT_EQ(p - 3.5, Polynomial({-2.5, 2, 3, 4}));

    const auto q = Polynomial({1, 2, 3, 0});
    EXPECT_EQ((q - 2.5), Polynomial({-1.5, 2, 3}));
    EXPECT_EQ(q, Polynomial({1, 2, 3, 0}));
}

// Tests for multiplication with a scalar.
TEST(PolynomialTests, MultiplicationWithScalarTest) {
    const auto p = Polynomial({1, 2, 3, 4});
    const auto q = p * 2.5;

    EXPECT_EQ(q, 2.5 * p);
    EXPECT_EQ(q, Polynomial({2.5, 5, 7.5, 10.0}));
    EXPECT_EQ(p, Polynomial({1, 2, 3, 4}));
}

// Tests for division with a scalar.
TEST(PolynomialTests, DivisionWithScalarTest) {
    const auto p = Polynomial({1, 2, 3, 4});
    const auto q = p / 2.5;

    EXPECT_EQ(q, p / 2.5);
    EXPECT_EQ(q, Polynomial({1. / 2.5, 2. / 2.5, 3. / 2.5, 4. / 2.5}));
    EXPECT_EQ(p, Polynomial({1, 2, 3, 4}));
}

// Tests polynomial addition.
TEST(PolynomialTests, PolynomialAdditionTest) {
    const auto p = Polynomial({1, 2, 4, 6});
    const auto q = Polynomial({2, 2, 1});
    EXPECT_EQ(p + q, Polynomial({3, 4, 5, 6}));
}

TEST(PolynomialTests, CompoundAssignmentOperatorScalarAdditionTest) {
    auto p = 3 * X.pow(4) - 2 * X.pow(3) + X.pow(2) + 2;
    EXPECT_EQ(p += 2.5, 3 * X.pow(4) - 2 * X.pow(3) + X.pow(2) + 4.5);
}

TEST(PolynomialTests, CompoundAssignmentOperatorPolynomialAdditionTest) {
    auto p = Polynomial({1, 2, 3, 4});
    const auto q = Polynomial({2, 2, 0});

    p += q;
    EXPECT_EQ(p, Polynomial({3, 4, 3, 4}));
}

// Tests polynomial subtraction.
TEST(PolynomialTests, PolynomialSubstractionTest) {
    const auto p = Polynomial({1, 2, 3, 4});
    const auto q = Polynomial({2, 2, 0});

    EXPECT_EQ(p - q, Polynomial({-1, 0, 3, 4}));
    EXPECT_TRUE((p - p).is_zero());
}

TEST(PolynomialTests, CompoundAssignmentOperatorScalarSubstractionTest) {
    auto p = 3 * X.pow(4) - 2 * X.pow(3) + X.pow(2) + 2;
    EXPECT_EQ(p -= 2.5, 3 * X.pow(4) - 2 * X.pow(3) + X.pow(2) - .5);
}

TEST(PolynomialTests, CompoundAssignmentOperatorPolynomialSubstractionTest) {
    auto p = Polynomial({1, 2, 3, 4});
    const auto q = Polynomial({2, 2, 0});

    p -= q;
    EXPECT_EQ(p, Polynomial({-1, 0, 3, 4}));
}

// Tests polynomial multiplication.
TEST(PolynomialTests, PolynomialMultiplicationTest) {
    const auto p = Polynomial({1, 2, 3, 4}); // 1 + 2x + 3x^2 + 4x^3
    const auto q = Polynomial({2, 2, 0}); // 2 + 2x

    EXPECT_EQ(p * q, Polynomial({2, 6, 10, 14, 8, 0})); // 2 + 6x + 10x^2 + 14x^3 + 8x^4
}

// Tests multiplication with a zero polynomial.
TEST(PolynomialTests, MultiplicationWithZeroPolynomialTest) {
    const auto p = Polynomial({0, 0, 1, 2, 3, 4});
    EXPECT_TRUE((p * zero).is_zero());
}

TEST(PolynomialTests, CompoundAssignmentScalarMultiplicationOperatorTest) {
    auto p = 3 * X.pow(4) - 2 * X.pow(3) + X.pow(2) + 2;
    EXPECT_EQ(p *= 2, 6 * X.pow(4) - 4 * X.pow(3) + 2 * X.pow(2) + 4);
}

TEST(PolynomialTests, CompoundAssignmentPolynomialMultiplicationOperatorTest) {
    auto p = Polynomial({1, 2, 3, 4}); // 1 + 2x + 3x^2 + 4x^3
    const auto q = Polynomial({2, 2, 0}); // 2 + 2x

    p *= q;
    EXPECT_EQ(p, Polynomial({2, 6, 10, 14, 8, 0})); // 2 + 6x + 10x^2 + 14x^3 + 8x^4
}

TEST(PolynomialTests, IsNormalizeTest) {
    const auto p = 1 + 2 * X + 3 * X.pow(2) + X.pow(3);
    EXPECT_TRUE(p.is_normalized());
    EXPECT_TRUE(!(4.3 * p).is_normalized());
}

TEST(PolynomialTests, NormalizeTest) {
    const auto p = 1 + 2 * X + 3 * X.pow(2) + 4 * X.pow(3);
    const auto q = p.normalize();

    EXPECT_EQ(q.leading_coefficient(), 1.0);
    EXPECT_TRUE(q.is_normalized());
    EXPECT_EQ(q, (1. / 4.) + (1. / 2.) * X + (3. / 4.) * X.pow(2) + X.pow(3));
}

TEST(PolynomialTests, IsIntegerTest) {
    ASSERT_TRUE(Polynomial({}).is_integer());
    ASSERT_TRUE(Polynomial({1.0, 2.0, 3.0}).is_integer());
    ASSERT_TRUE(Polynomial({-1.0, -2.0, 3.0}).is_integer());
    ASSERT_TRUE(Polynomial({3. / 3., -2.0, 3.0}).is_integer());
    ASSERT_FALSE(Polynomial({1.0, 2.5, 3.0}).is_integer());
    ASSERT_FALSE(Polynomial({1.0, 1. / 3., 3.0}).is_integer());
}

// Tests the derivative of a first-order polynomial.
TEST(PolynomialTests, DerivativeOfFirstOrderTest) {
    const auto p = one - 6 * X.pow(2) + 2 * X.pow(3) + 3 * X.pow(4) + 4 * X.pow(5);
    EXPECT_EQ(p.derive(), -12 * X + 6 * X.pow(2) + 12 * X.pow(3) + 20 * X.pow(4));

    EXPECT_EQ((X.pow(4) + 1).derive(), 4 * X.pow(3));
    EXPECT_EQ((X.pow(3) + 1).derive(), 3 * X.pow(2));
    EXPECT_EQ((X.pow(2) + 1).derive(), 2 * X);
    EXPECT_EQ((X + 1).derive(), one);

    EXPECT_TRUE(one.derive().is_zero());
    EXPECT_TRUE(zero.derive().is_zero());
}

// Tests the indefinite integral of the polynomial.
TEST(PolynomialTests, IntegrateTest) {
    EXPECT_EQ(one.integrate(), X);
    EXPECT_EQ((3 * one).integrate(), 3 * X);
    EXPECT_EQ(X.integrate(), .5 * X.pow(2));
    EXPECT_EQ((X.pow(3) + 4 * X.pow(2)).integrate(), 0.25 * X.pow(4) + (4. / 3) * X.pow(3));
    EXPECT_EQ((X.pow(2)).integrate(), (1. / 3) * X.pow(3));
}

// Tests polynomial evaluation at different points.
TEST(PolynomialTests, EvaluateTest) {
    const auto p1 = Polynomial({1, 1, 2}); // 1 + x + 2x^2

    EXPECT_NEAR(p1.evaluate(-2), 7.0, epsilon);
    EXPECT_NEAR(p1.evaluate(-1), 2.0, epsilon);
    EXPECT_NEAR(p1.evaluate(0), 1.0, epsilon);
    EXPECT_NEAR(p1.evaluate(1), 4.0, epsilon);
    EXPECT_NEAR(p1.evaluate(2), 11.0, epsilon);

    const auto p2 = Polynomial({0, -1, 1}); // -x + x^2
    EXPECT_NEAR(p2.evaluate(-1), 2., epsilon);
    EXPECT_NEAR(p2.evaluate(0), 0., epsilon);
    EXPECT_NEAR(p2.evaluate(1), 0., epsilon);

    const auto p3 = Polynomial({1, -2, 1});
    EXPECT_NEAR(p3.evaluate(1), 0., epsilon);
    EXPECT_NEAR(p3.evaluate(-1), 4., epsilon);
}

// Tests polynomial evaluation using the call operator.
TEST(PolynomialTests, EvaluateWithOperatorTest) {
    auto f = [](const double x) { return -3 + 2 * x - std::pow(x, 2) + pow(x, 3) + 2 * pow(x, 4); };
    const auto p = -3 + 2 * X - X.pow(2) + X.pow(3) + 2 * X.pow(4);

    for (const auto x: {-2.0, -1.0, 0.0, 1.0, 2.0}) {
        EXPECT_NEAR(p(x), f(x), epsilon);
    }
}

TEST(PolynomialTests, IsRootTest) {
    const auto p = X.pow(2) - 2.;
    ASSERT_NEAR(std::sqrt(2) * std::sqrt(2) - 2., 0., epsilon);
    ASSERT_NEAR(p(std::sqrt(2)), 0., epsilon);

    EXPECT_TRUE(p.is_root(std::sqrt(2)));
    EXPECT_TRUE(p.is_root(-std::sqrt(2)));
    EXPECT_FALSE(p.is_root(1.0));
}

TEST(PolynomialTests, FromRootsTest) {
    auto p = Polynomial::from_roots({});
    EXPECT_EQ(p, one);

    values_type roots = {-3};
    p = Polynomial::from_roots(roots);
    EXPECT_TRUE(p.has_roots(roots));
    EXPECT_TRUE(p.is_normalized());
    EXPECT_EQ(p.degree(), 1);

    roots = {-1, 0, 1};
    p = Polynomial::from_roots(roots);
    EXPECT_TRUE(p.has_roots(roots));
    EXPECT_TRUE(p.is_normalized());
    EXPECT_EQ(p.degree(), 3);

    roots = {1.2, 1.5, -2.4, 6.3};
    p = Polynomial::from_roots(roots);
    EXPECT_TRUE(p.has_roots(roots));
    EXPECT_TRUE(p.is_normalized());
    EXPECT_EQ(p.degree(), 4);

    roots = {-1.3, 0.2, 1.0, 4.1, 3.1, 8.12};
    p = Polynomial::from_roots(roots);
    EXPECT_TRUE(p.has_roots(roots));
    EXPECT_TRUE((1.2 * p).has_roots(roots));
    EXPECT_TRUE(p.is_normalized());
    EXPECT_EQ(p.degree(), 6);

    constexpr auto pi = std::numbers::pi;
    roots = {-pi, pi};

    p = Polynomial::from_roots(roots);
    EXPECT_TRUE(p.has_roots(roots));
    EXPECT_NEAR(p(pi), 0, epsilon);
    EXPECT_NEAR(p(-pi), 0, epsilon);
}

// Tests polynomial divide function.
TEST(PolynomialTests, DivideTest) {
    check_polynomial_division(X, one, X, zero);
    check_polynomial_division(X.pow(2), 4. * one, .25 * X.pow(2), zero);

    check_polynomial_division(X.pow(2) - 2 * X + 1., X - 1., X - 1, zero);
    check_polynomial_division(X.pow(2) - 2 * X + 1., X + 2., X - 4., 9. * one);

    check_polynomial_division(
        3 * X.pow(3) + X.pow(2) + X + 5, 5 * X.pow(2) + -3 * X + 1,
        3. / 5 * X + (14. / 25), (52. / 25) * X + (111. / 25));

    check_polynomial_division(
        3 * X.pow(3) + X.pow(2) + X + 5, 5 * X.pow(2) + -3 * X + 1,
        3. / 5 * X + (14. / 25), (52. / 25) * X + (111. / 25));
}

TEST(PolynomialTests, DivideWithLargeLeadingCoefficientTest) {
    constexpr auto lc = 1e10;
    constexpr auto eps = 1e-13;

    EXPECT_FALSE(nearly_equal(lc + 1., lc, eps));
    EXPECT_FALSE(nearly_equal(lc - 1., lc, eps));
    EXPECT_FALSE(nearly_zero(1. / lc, eps));

    check_polynomial_division(
        lc * X.pow(4) - 1,
        X - 1.,
        lc * (X.pow(3) + X.pow(2) + X.pow(1) + 1.), -one + lc);
}

// Tests polynomial divide function.
TEST(PolynomialTests, DivisionOperatorTest) {
    const auto quotient = (3 * X.pow(3) + X.pow(2) + X + 5) / (5 * X.pow(2) + -3 * X + 1);
    EXPECT_EQ(quotient, 3. / 5 * X + (14. / 25));
}

TEST(PolynomialTests, ModuloOperatorTest) {
    const auto remainder = (3 * X.pow(3) + X.pow(2) + X + 5) % (5 * X.pow(2) + -3 * X + 1);
    EXPECT_EQ(remainder, 52. / 25 * X + (111. / 25));
}

TEST(PolynomialTests, CompoundAssignmentScalarDivisionOperatorTest) {
    auto p = 3 * X.pow(4) - 2 * X.pow(3) + X.pow(2) + 2;
    EXPECT_EQ(p /= 2, (3. / 2) * X.pow(4) - X.pow(3) + .5 * X.pow(2) + 1);
}

TEST(PolynomialTests, CompoundAssignmentPolynomialDivisionOperatorTest) {
    auto p = 3 * X.pow(3) + X.pow(2) + X + 5;
    const auto q = 5 * X.pow(2) + -3 * X + 1;

    p /= q;
    EXPECT_EQ(p, (3. / 5) * X + (14. / 25));
}

TEST(PolynomialTests, CompoundAssignmentPolynomialModuloOperatorTest) {
    auto p = 3 * X.pow(3) + X.pow(2) + X + 5;
    const auto q = 5 * X.pow(2) + -3 * X + 1;

    p %= q;
    EXPECT_EQ(p, (52. / 25) * X + (111. / 25));
}

// Test for the compose function.
TEST(PolynomialTests, ComposeTest) {
    const Polynomial r = X;
    const Polynomial p = 1 + 2 * X;
    const Polynomial q = 3 + 4 * X;
    const Polynomial s = -1 + 2 * X + 4 * X.pow(2);

    EXPECT_EQ(p.compose(r), p);
    EXPECT_EQ(p.compose(q), 1 + 2 * q); // p(q(x))
    EXPECT_EQ(q.compose(p), 3 + 4 * p); // q(p(x))
    EXPECT_EQ(s.compose(p), -1 + 2 * (1 + 2 * X) + 4 * (1 + 2 * X).pow(2)); // s(p(x))
    EXPECT_EQ(s.compose(q), -1 + 2 * (3 + 4 * X) + 4 * (3 + 4 * X).pow(2)); // s(q(x))
}
