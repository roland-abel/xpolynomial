/// @file polynomial_test.cpp
/// @brief This file contains unit tests for the Polynomial class and related functions.
///
/// @author Roland Abel
/// @date 19.08.2023

#include <gtest/gtest.h>
#include <cmath>
#include "polynomial.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using values_type = Polynomial::values_type;

    constexpr double tolerance = Polynomial::tolerance;
    Polynomial zero = Polynomial::zero();
    Polynomial one = Polynomial::one();
    Polynomial X = Polynomial::monomial(1, 1.0);
}

// Tests for the default constructor.
TEST(PolynomialTests, DefaultConstructur) {
    auto p = Polynomial();
    EXPECT_EQ(p.degree(), 0);
}

// Tests for the copy constructor.
TEST(PolynomialTests, CopyConstructur) {
    auto p = Polynomial({1, 2, 3});
    auto q = Polynomial(p);

    EXPECT_TRUE(p == q);
}

// Tests for the constructor with coefficients.
TEST(PolynomialTests, ConstructurWithCoefficients) {
    EXPECT_EQ(Polynomial({0}).degree(), 0);
    EXPECT_EQ(Polynomial({0, 0, 0}).degree(), 0);
    EXPECT_EQ(Polynomial({0, 0, 1}).degree(), 2);   // x^2
    EXPECT_EQ(Polynomial({0, 1, -1}).degree(), 2);  // x - x^2
    EXPECT_EQ(Polynomial({1, 2, 1}).degree(), 2);   // 1 + 2x + x^2
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
    EXPECT_TRUE(Polynomial({0, 0, 1}).is_quadratic());// x^2
    EXPECT_TRUE(Polynomial({0, 1, 1}).is_quadratic());// x + x^2
    EXPECT_TRUE(Polynomial({1, 2, -1}).is_quadratic());// 1 + 2x - x^2

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
    EXPECT_NEAR(Polynomial({1, 0, 0}).leading_coefficient(), 1.0, tolerance);
    EXPECT_NEAR(Polynomial({3.0, 1.0, 2.0}).leading_coefficient(), 2.0, tolerance);
    EXPECT_NEAR(Polynomial({0, 0, 3.5, 1.0, 4.0}).leading_coefficient(), 4.0, tolerance);
}

// Tests for monomial.
TEST(PolynomialTests, MonomialTest) {
    EXPECT_EQ(Polynomial::monomial(1, 1.0), X);
    EXPECT_EQ(Polynomial::monomial(1, 2.0), 2 * X);
    EXPECT_EQ(Polynomial::monomial(4, 3.5), 3.5 * X.pow(4));
    EXPECT_EQ(Polynomial::monomial(10, -2.5), -2.5 * X.pow(10));
    EXPECT_EQ(Polynomial::monomial(10, 3.5).degree(), 10);
}

// Tests for the equal operator.
TEST(PolynomialTests, EqualOperatorTest) {
    auto p = Polynomial({1, 2, 3});

    EXPECT_TRUE(p == Polynomial({1, 2, 3}));
    EXPECT_TRUE(p == Polynomial({1, 2, 3, 0, 0}));
    EXPECT_FALSE(p == Polynomial({0, 2, 3, 0, 0}));
}

// Tests for the not equal operator.
TEST(PolynomialTests, NotEqualOperatorTest) {
    auto p = Polynomial({1, 2, 3, 4});

    EXPECT_FALSE(p != Polynomial({1, 2, 3, 4}));
    EXPECT_FALSE(p != Polynomial({1, 2, 3, 4, 0, 0}));
    EXPECT_TRUE(p != Polynomial({0, 2, 3, 4}));
}

TEST(PolynomialTests, UnitaryPlusOperatorTest) {
    auto p = Polynomial({1, -2, 3, -4});
    EXPECT_EQ(+p, p);
}

TEST(PolynomialTests, UnitaryMinusOperatorTest) {
    auto p = Polynomial({-1, 2, -3, 4});
    EXPECT_EQ(-p, Polynomial({1, -2, 3, -4}));
    EXPECT_EQ(-(-p), p);
}

// Tests for the index operator.
TEST(PolynomialTests, IndexOperatorTest) {
    auto p = Polynomial({0.5, 2.3, 4, -1, 3.6});

    EXPECT_NEAR(p[4], 3.6, tolerance);
    EXPECT_NEAR(p[3], -1, tolerance);
    EXPECT_NEAR(p[2], 4, tolerance);
    EXPECT_NEAR(p[1], 2.3, tolerance);
    EXPECT_NEAR(p[0], 0.5, tolerance);
}

// Tests for the power function.
TEST(PolynomialTests, PowerFunctionTest) {
    Polynomial p = Polynomial({0, 1});

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
    EXPECT_NEAR(p[10], 0.0, tolerance);
}

// Tests polynomial addition with a scalar.
TEST(PolynomialTests, AdditionWithScalarTest) {
    auto p = Polynomial({1, 2, 3, 4});  // 1 + 2x + 3x^2 + 4x^3
    EXPECT_EQ((p + 3.5), Polynomial({4.5, 2, 3, 4}));

    p = Polynomial({1, 2, 3, 0});
    EXPECT_EQ((p + 3.5), Polynomial({4.5, 2, 3}));
    EXPECT_EQ(p, Polynomial({1, 2, 3, 0}));
}

// Tests polynomial subtraction with a scalar.
TEST(PolynomialTests, SubstrationWithScalarTest) {
    auto p = Polynomial({1, 2, 3, 4});  // 1 + 2x + 3x^2 + 4x^3
    EXPECT_EQ((p - 3.5), Polynomial({-2.5, 2, 3, 4}));

    p = Polynomial({1, 2, 3, 0});
    EXPECT_EQ((p - 2.5), Polynomial({-1.5, 2, 3}));
    EXPECT_EQ(p, Polynomial({1, 2, 3, 0}));
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
    const auto p = Polynomial({1, 2, 3, 4});
    const auto q = Polynomial({2, 2, 0});
    EXPECT_EQ(p + q, Polynomial({3, 4, 3, 4}));
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
    auto p = Polynomial({1, 2, 3, 4});  // 1 + 2x + 3x^2 + 4x^3
    auto q = Polynomial({2, 2, 0});     // 2 + 2x

    EXPECT_EQ(p * q, Polynomial({2, 6, 10, 14, 8, 0})); // 2 + 6x + 10x^2 + 14x^3 + 8x^4
}

// Tests multiplication with a zero polynomial.
TEST(PolynomialTests, MultiplicationWithZeroPolynomialTest) {
    auto p = Polynomial({0, 0, 1, 2, 3, 4});
    EXPECT_TRUE((p * zero).is_zero());
}

TEST(PolynomialTests, CompoundAssignmentScalarMultiplicationOperatorTest) {
    auto p = 3 * X.pow(4) - 2 * X.pow(3) + X.pow(2) + 2;
    EXPECT_EQ(p *= 2, 6 * X.pow(4) - 4 * X.pow(3) + 2 * X.pow(2) + 4);
}

TEST(PolynomialTests, CompoundAssignmentPolynomialMultiplicationOperatorTest) {
    auto p = Polynomial({1, 2, 3, 4});          // 1 + 2x + 3x^2 + 4x^3
    const auto q = Polynomial({2, 2, 0});       // 2 + 2x

    p *= q;
    EXPECT_EQ(p, Polynomial({2, 6, 10, 14, 8, 0})); // 2 + 6x + 10x^2 + 14x^3 + 8x^4
}

TEST(PolynomialTests, NormalizeTest) {
    auto p = 1 + 2 * X + 3 * X.pow(2) + 4 * X.pow(3);
    const auto q = p.normalize();

    EXPECT_EQ(q.leading_coefficient(), 1.0);
    EXPECT_EQ(q, (1. / 4.) + (1. / 2.) * X + (3. / 4.) * X.pow(2) + X.pow(3));
}

// Tests the derivative of a first-order polynomial.
TEST(PolynomialTests, DerivativeOfFirstOrderTest) {
    auto p = one - 6 * X.pow(2) + 2 * X.pow(3) + 3 * X.pow(4) + 4 * X.pow(5);
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
    EXPECT_EQ(one.integrate(), zero);
    EXPECT_EQ(X.integrate(), .5 * X.pow(2));
    EXPECT_EQ((X.pow(3) + 4 * X.pow(2)).integrate(), 0.25 * X.pow(4) + (4. / 3) * X.pow(3));
    EXPECT_EQ((X.pow(2)).integrate(), (1. / 3) * X.pow(3));
}

// Tests polynomial evaluation at different points.
TEST(PolynomialTests, EvaluateTest) {
    auto p = Polynomial({1, 1, 2});    // 1 + x + 2x^2

    EXPECT_NEAR(p.evaluate(-2), 7.0, tolerance);
    EXPECT_NEAR(p.evaluate(-1), 2.0, tolerance);
    EXPECT_NEAR(p.evaluate(0), 1.0, tolerance);
    EXPECT_NEAR(p.evaluate(1), 4.0, tolerance);
    EXPECT_NEAR(p.evaluate(2), 11.0, tolerance);

    p = Polynomial({0, -1, 1}); // -x + x^2
    EXPECT_NEAR(p.evaluate(-1), 2., tolerance);
    EXPECT_NEAR(p.evaluate(0), 0., tolerance);
    EXPECT_NEAR(p.evaluate(1), 0., tolerance);
}

// Tests polynomial evaluation using the call operator.
TEST(PolynomialTests, EvaluateWithOperatorTest) {
    auto f = [](double x) { return -3 + 2 * x - std::pow(x, 2) + pow(x, 3) + 2 * pow(x, 4); };
    auto p = -3 + 2 * X - X.pow(2) + X.pow(3) + 2 * X.pow(4);

    auto xs = {-2.0, -1.0, 0.0, 1.0, 2.0};
    for (auto x: xs) {
        EXPECT_NEAR(p(x), f(x), tolerance);
    }
}

TEST(PolynomialTests, HasRootTest) {
    auto p = X - 1.47;
    EXPECT_TRUE(p.has_root(1.47));

    p = X.pow(2) - 2;
    EXPECT_TRUE(p.has_root(std::sqrt(2)));
    EXPECT_TRUE(p.has_root(-std::sqrt(2)));
    EXPECT_FALSE(p.has_root(1.0));
}

// Tests polynomial creation from roots.
TEST(PolynomialTests, FromRootsTest) {
    auto p = Polynomial::from_roots({});
    EXPECT_EQ(p, one);

    values_type roots = {-3};
    p = Polynomial::from_roots(roots);
    EXPECT_TRUE(p.has_roots(roots));
    EXPECT_EQ(p.degree(), 1);

    roots = {-1, 0, 1};
    p = Polynomial::from_roots(roots);
    EXPECT_TRUE(p.has_roots(roots));
    EXPECT_EQ(p.degree(), 3);

    roots = {1.2, 1.5, -2.4, 6.3};
    p = Polynomial::from_roots(roots);
    EXPECT_TRUE(p.has_roots(roots));
    EXPECT_EQ(p.degree(), 4);

    roots = {-1.3, 0.2, 1.0, 4.1, 3.1, 8.12};
    p = Polynomial::from_roots(roots);
    EXPECT_TRUE(p.has_roots(roots));
    EXPECT_EQ(p.degree(), 6);
}

// Tests polynomial divide function.
TEST(PolynomialTests, DivideTest) {
    auto [quotient, remainder] = X.divide(one);    // X/1 = X remainder 0

    EXPECT_EQ(quotient, X);
    EXPECT_EQ(remainder, zero);

    auto [q1, r1] = X.pow(2).divide(4. * one); // X^2/4 = 0.25*X^2
    EXPECT_EQ(q1, .25 * X.pow(2));
    EXPECT_EQ(r1, zero);

    auto [q2, r2] = (3 * X.pow(3) + X.pow(2) + X + 5).divide(5 * X.pow(2) + -3 * X + 1);
    EXPECT_EQ(q2, (3. / 5) * X + (14. / 25));
    EXPECT_EQ(r2, (52. / 25) * X + (111. / 25));
}

// Tests polynomial divide function.
TEST(PolynomialTests, DivisionOperatorTest) {
    auto quotient = (3 * X.pow(3) + X.pow(2) + X + 5) / (5 * X.pow(2) + -3 * X + 1);
    EXPECT_EQ(quotient, (3. / 5) * X + (14. / 25));
}

TEST(PolynomialTests, ModuloOperatorTest) {
    auto remainder = (3 * X.pow(3) + X.pow(2) + X + 5) % (5 * X.pow(2) + -3 * X + 1);
    EXPECT_EQ(remainder, (52. / 25) * X + (111. / 25));
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
    Polynomial r = X;
    Polynomial p = 1 + 2 * X;
    Polynomial q = 3 + 4 * X;
    Polynomial s = -1 + 2 * X + 4 * X.pow(2);

    EXPECT_EQ(p.compose(r), p);
    EXPECT_EQ(p.compose(q), 1 + 2 * q); // p(q(x))
    EXPECT_EQ(q.compose(p), 3 + 4 * p); // q(p(x))
    EXPECT_EQ(s.compose(p), -1 + 2 * (1 + 2 * X) + 4 * (1 + 2 * X).pow(2)); // s(p(x))
    EXPECT_EQ(s.compose(q), -1 + 2 * (3 + 4 * X) + 4 * (3 + 4 * X).pow(2)); // s(q(x))
}
