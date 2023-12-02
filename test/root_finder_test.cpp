/// @file root_finder_tests.cpp
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
#include "polynomial.h"
#include "real_polynomial_root_finder.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using RootFinder = root_finder<double>;

    constexpr auto epsilon = 1e-5;
    auto X = Polynomial::monomial(1);
}

TEST(RootFinderTests, BisectionWithIncorrectsEndpointsTest) {
    auto I = real_interval(3., 4.);
    auto root = RootFinder::bisection(X, I);

    EXPECT_TRUE(std::isnan(root));
}

TEST(RootFinderTests, BisectionTest) {
    auto p = 4 * X.pow(2) + .5 * X - 4;

    EXPECT_NEAR(RootFinder::bisection(p, real_interval(.0, 2.)), (std::sqrt(257) - 1.) / 16., epsilon);
    EXPECT_NEAR(RootFinder::bisection(p, real_interval(-2., .0)), (-1. - std::sqrt(257)) / 16., epsilon);
}

TEST(RootFinderTests, Bisection2Test) {
    auto p = Polynomial::from_roots({-2, 0, -1, 1});
    auto root = RootFinder::bisection(p, real_interval(-3., -1.5));

    ASSERT_NEAR(root, -2., epsilon);
    ASSERT_NEAR(p(root), 0., epsilon);

    EXPECT_TRUE(p.is_root(root));
}

TEST(RootFinderTests, Bisection3Test) {
    auto p = X.pow(3) - .75 * X;
    auto root = RootFinder::bisection(p, real_interval(-.875, -.4375));

    EXPECT_NEAR(root, -std::sqrt(.75), epsilon);
    EXPECT_NEAR(p(root), 0., epsilon);
    EXPECT_TRUE(p.is_root(root));
}

TEST(RootFinderTests, NewtonRaphsonForQudraticPolynomialTest) {
    auto p = 4 * X.pow(2) + .5 * X - 4;
    auto q = p.derive();

    EXPECT_NEAR(RootFinder::newton_raphson(p, q, 0.5), 0.939451, epsilon);
    EXPECT_NEAR(RootFinder::newton_raphson(p, q, -0.5), -1.06445, epsilon);

    p = X.pow(3) - 3;
    q = p.derive();
    EXPECT_NEAR(RootFinder::newton_raphson(p, q, 1.), 1.44224, epsilon);
}

TEST(RootFinderTests, NewtonRaphsonForCubicPolynomialTest) {
    auto p = 2 * X.pow(3) - 3 * X - 1;
    auto q = p.derive();

    EXPECT_NEAR(RootFinder::newton_raphson(p, q, -0.9), -1, epsilon);
    EXPECT_NEAR(RootFinder::newton_raphson(p, q, -0.4), -0.36602, epsilon);
    EXPECT_NEAR(RootFinder::newton_raphson(p, q, 1.3), 1.36602, epsilon);
}

TEST(RootFinderTests, NewtonRaphsonCosTest) {
    auto func = [](const auto &x) { return std::cos(x); };
    auto f_prim = [](auto x) { return -std::sin(x); };

    const auto pi_half = std::numbers::pi / 2.;
    EXPECT_NEAR(RootFinder::newton_raphson(func, f_prim, 1.1, 100), pi_half, epsilon);
}

TEST(RootFinderTests, NewtonRaphsonFailTest) {
    auto p = X.pow(3) - 3 * X - 1;
    auto q = p.derive();

    EXPECT_TRUE(std::isnan(RootFinder::newton_raphson(p, q, -1.)));
}

TEST(RootFinderTests, RegulaFalsiTest) {
    auto func = [](double x) { return 2 * std::cos(x); };
    auto I = real_interval(0.25, std::numbers::pi);
    auto zero_point = RootFinder::regula_falsi(func, I);

    EXPECT_NEAR(zero_point, .5 * std::numbers::pi, epsilon);
}

TEST(RootFinderTests, RegulaFalsi2Test) {
    auto p = X.pow(5) - 10 * X.pow(4) + 40 * X.pow(3) - 80 * X.pow(2) + 80 * X - 30;
    auto I = real_interval(0., 2.);
    auto root = RootFinder::regula_falsi(p, I);

    EXPECT_NEAR(root, 0.85130254011, epsilon);
}