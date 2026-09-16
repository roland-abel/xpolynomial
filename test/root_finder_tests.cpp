/// @file root_finder_tests.cpp
/// @brief
///
/// @author Roland Abel
/// @date October 8, 2023
///
/// Copyright (c) 2026 Roland Abel

#include <gtest/gtest.h>
#include <cmath>
#include <numbers>
#include "polynomial.h"
#include "real_polynomial_root_finder.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using RootFinder = root_finder<double>;

    constexpr auto epsilon = 1e-5;
    const auto X = Polynomial::monomial(1u);
}

TEST(RootFinderTests, BisectionWithIncorrectsEndpointsTest) {
    const auto I = interval(3., 4.);
    const auto root = RootFinder::bisection(X, I);

    EXPECT_FALSE(root.has_value());
}

TEST(RootFinderTests, BisectionTest) {
    const auto p = 4 * X.pow(2) + .5 * X - 4;

    EXPECT_NEAR(RootFinder::bisection(p, interval(.0, 2.)).value(), (std::sqrt(257) - 1.) / 16., epsilon);
    EXPECT_NEAR(RootFinder::bisection(p, interval(-2., .0)).value(), (-1. - std::sqrt(257)) / 16., epsilon);
}

TEST(RootFinderTests, Bisection2Test) {
    const auto p = Polynomial::from_roots({-2, 0, -1, 1});
    const auto root = RootFinder::bisection(p, interval(-3., -1.5)).value();

    ASSERT_NEAR(root, -2., epsilon);
    ASSERT_NEAR(p(root), 0., epsilon);

    EXPECT_TRUE(p.is_root(root));
}

TEST(RootFinderTests, Bisection3Test) {
    const auto p = X.pow(3) - .75 * X;
    const auto root = RootFinder::bisection(p, interval(-.875, -.4375)).value();

    EXPECT_NEAR(root, -std::sqrt(.75), epsilon);
    EXPECT_NEAR(p(root), 0., epsilon);
    EXPECT_TRUE(p.is_root(root));
}

TEST(RootFinderTests, NewtonRaphsonForQudraticPolynomialTest) {
    const auto p1 = 4 * X.pow(2) + .5 * X - 4;
    const auto q1 = p1.derive();

    EXPECT_NEAR(RootFinder::newton_raphson(p1, q1, 0.5).value(), 0.939451, epsilon);
    EXPECT_NEAR(RootFinder::newton_raphson(p1, q1, -0.5).value(), -1.06445, epsilon);

    const auto p2 = X.pow(3) - 3;
    const auto q2 = p2.derive();
    EXPECT_NEAR(RootFinder::newton_raphson(X.pow(3) - 3, q2, 1.).value(), 1.44224, epsilon);
}

TEST(RootFinderTests, NewtonRaphsonForCubicPolynomialTest) {
    auto p = 2 * X.pow(3) - 3 * X - 1;
    auto q = p.derive();

    EXPECT_NEAR(RootFinder::newton_raphson(p, q, -0.9).value(), -1, epsilon);
    EXPECT_NEAR(RootFinder::newton_raphson(p, q, -0.4).value(), -0.36602, epsilon);
    EXPECT_NEAR(RootFinder::newton_raphson(p, q, 1.3).value(), 1.36602, epsilon);
}

TEST(RootFinderTests, NewtonRaphsonCosTest) {
    auto func = [](const auto &x) { return std::cos(x); };
    auto f_prim = [](auto x) { return -std::sin(x); };

    constexpr auto pi_half = std::numbers::pi / 2.;
    EXPECT_NEAR(RootFinder::newton_raphson(func, f_prim, 1.1, 100).value(), pi_half, epsilon);
}

TEST(RootFinderTests, NewtonRaphsonFailTest) {
    auto p = X.pow(3) - 3 * X - 1;
    auto q = p.derive();

    EXPECT_FALSE(RootFinder::newton_raphson(p, q, -1.).has_value());
}

TEST(RootFinderTests, RegulaFalsiTest) {
    const auto func = [](const double x) { return 2 * std::cos(x); };
    const auto I = interval(0.25, std::numbers::pi);
    const auto zero_point = RootFinder::regula_falsi(func, I).value();

    EXPECT_NEAR(zero_point, .5 * std::numbers::pi, epsilon);
}

TEST(RootFinderTests, RegulaFalsi2Test) {
    const auto p = X.pow(5) - 10 * X.pow(4) + 40 * X.pow(3) - 80 * X.pow(2) + 80 * X - 30;
    const auto I = interval(0., 2.);
    const auto root = RootFinder::regula_falsi(p, I).value();

    EXPECT_NEAR(root, 0.85130254011, epsilon);
}