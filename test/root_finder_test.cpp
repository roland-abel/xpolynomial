/// @file root_finder_tests.h
///
/// @author Roland Abel
/// @date 08.10.2023

#include <gtest/gtest.h>
#include <cmath>
#include "polynomial.h"
#include "real_polynomial_root_finder.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using RootFinder = root_finder<double>;
    using value_type = Polynomial::value_type;

    constexpr auto epsilon = Polynomial::epsilon;
    auto zero = Polynomial::zero();
    auto one = Polynomial::one();
    auto X = Polynomial::monomial(1);
}

TEST(RootFinderTests, BisectionWithIncorrectsEndpointsTest) {
    auto I = real_interval(3., 4.);
    auto root = RootFinder::bisection(X, I, epsilon);

    EXPECT_TRUE(std::isnan(root));
}

TEST(RootFinderTests, BisectionTest) {
    auto p = 4 * X.pow(2) + .5 * X - 4;

    EXPECT_NEAR(RootFinder::bisection(p, real_interval(.0, 2.), epsilon), 0.9394510000000000, epsilon);
    EXPECT_NEAR(RootFinder::bisection(p, real_interval(-2., .0), epsilon), -1.06445, epsilon);
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

TEST(RootFinderTests, NewtonRaphsonFailTest) {
    auto p = X.pow(3) - 3 * X - 1;
    auto q = p.derive();

    EXPECT_TRUE(std::isnan(RootFinder::newton_raphson(p, q, -1.)));
}

TEST(RootFinderTests, RegulaFalsi1Test) {
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