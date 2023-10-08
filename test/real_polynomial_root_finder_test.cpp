/// @file polynomial_test.cpp
/// @brief This file contains unit tests for the Polynomial class and related functions.
///
/// @author Roland Abel
/// @date 19.08.2023

#include <gtest/gtest.h>
#include <cmath>
#include "polynomial.h"
#include "real_polynomial_root_finder.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using RootFinder = real_polynomial_root_finder<double>;
    using value_type = Polynomial::value_type;

    constexpr auto tolerance = Polynomial::tolerance;
    auto zero = Polynomial::zero();
    auto one = Polynomial::one();
    auto X = Polynomial::monomial(1);
}

TEST(RealPolynomialRootFinderTests, NotQuadraticPolynomialTest) {
    auto roots = RootFinder::quadratic_roots(X + 1);
    EXPECT_TRUE(std::isnan(get<0>(roots)));
    EXPECT_TRUE(std::isnan(get<1>(roots)));

    roots = RootFinder::quadratic_roots(zero);
    EXPECT_TRUE(std::isnan(get<0>(roots)));
    EXPECT_TRUE(std::isnan(get<1>(roots)));

    roots = RootFinder::quadratic_roots(X.pow(3));
    EXPECT_TRUE(std::isnan(get<0>(roots)));
    EXPECT_TRUE(std::isnan(get<1>(roots)));
}

TEST(RealPolynomialRootFinderTests, QuadraticPolynomialWhitoutRealRootsTest) {
    auto p = X.pow(2) + 1;
    auto roots = RootFinder::quadratic_roots(p);

    EXPECT_TRUE(std::isnan(get<0>(roots)));
    EXPECT_TRUE(std::isnan(get<1>(roots)));
}

TEST(RealPolynomialRootFinderTests, QuadraticPolynomialWithTwoRootsTest) {
    auto p = 4 * X.pow(2) + .5 * X - 4;
    auto roots = RootFinder::quadratic_roots(p);

    EXPECT_NEAR(get<0>(roots), (-.5 + std::sqrt(64.25)) / 8., tolerance);
    EXPECT_NEAR(get<1>(roots), (-.5 - std::sqrt(64.25)) / 8., tolerance);
}

TEST(RealPolynomialRootFinderTests, QuadraticPolynomialWithOneRootTest) {
    auto p = (X - 0.3).pow(2);
    auto roots = RootFinder::quadratic_roots(p);

    EXPECT_NEAR(get<0>(roots), .3, tolerance);
    EXPECT_NEAR(get<1>(roots), .3, tolerance);
}

TEST(RealPolynomialRootFinderTests, CubicPolynomialWithOneRootsTest) {
    auto p = X.pow(3) - 5.;
    auto roots = RootFinder::cubic_roots(p);

    EXPECT_NEAR(get<0>(roots), std::pow(5., 1. / 3.), tolerance);
    EXPECT_TRUE(std::isnan(get<1>(roots)));
    EXPECT_TRUE(std::isnan(get<2>(roots)));
}

TEST(RealPolynomialRootFinderTests, CubicNormalFormPolynomialWithThreeUnequalRootsTest) {
    auto p = X.pow(3) - 2 * X + 1;

    EXPECT_EQ((X - 1) * (X.pow(2) + X - 1), p);
    EXPECT_TRUE(p.has_roots({
                                    1.,
                                    -.5 + std::sqrt(5. / 4.),
                                    -.5 - std::sqrt(5. / 4.)
                            }));

    auto [r1, r2, r3] = RootFinder::cubic_roots(p);
    EXPECT_NEAR(r1, 1., tolerance);
    EXPECT_NEAR(r2, -.5 - std::sqrt(5. / 4.), tolerance);
    EXPECT_NEAR(r3, -.5 + std::sqrt(5. / 4.), tolerance);
}

TEST(RealPolynomialRootFinderTests, CubicPolynomialWithThreeUnequalRootsTest) {
    auto p = (X - 2.) * (X + 4.) * (X - 5.);
    EXPECT_TRUE(p.has_roots({5., 2., -4}));

    auto [r1, r2, r3] = RootFinder::cubic_roots(p);
    EXPECT_NEAR(r1, 5., tolerance);
    EXPECT_NEAR(r2, -4., tolerance);
    EXPECT_NEAR(r3, 2., tolerance);
}

TEST(RealPolynomialRootFinderTests, CubicPolynomialWithTwoEqualRootsTest) {
    auto p = 3.5 * (X - 7.).pow(2) * (X - 4.);
    EXPECT_TRUE(p.has_roots({7., 4.}));

    auto [r1, r2, r3] = RootFinder::cubic_roots(p);
    EXPECT_NEAR(r1, 4., tolerance);
    EXPECT_NEAR(r2, 7., tolerance);
    EXPECT_NEAR(r3, 7., tolerance);
}

TEST(RealPolynomialRootFinderTests, CubicPolynomialWithThreeEqualRootsTest) {
    auto p = 2.25 * X.pow(3);

    auto [r1, r2, r3] = RootFinder::cubic_roots(p);
    EXPECT_NEAR(r1, 0., tolerance);
    EXPECT_NEAR(r2, 0., tolerance);
    EXPECT_NEAR(r3, 0., tolerance);

    p = -1.5 * (X + 7.5).pow(3);
    EXPECT_TRUE(p.has_roots({-7.5}));

    auto roots = RootFinder::cubic_roots(p);
    EXPECT_NEAR(get<0>(roots), -7.5, tolerance);
    EXPECT_NEAR(get<1>(roots), -7.5, tolerance);
    EXPECT_NEAR(get<2>(roots), -7.5, tolerance);
}

TEST(RealPolynomialRootFinderTests, SignChangesOfCoefficients) {
    EXPECT_EQ(RootFinder::sign_changes(zero), 0);
    EXPECT_EQ(RootFinder::sign_changes(5 * one), 0);
    EXPECT_EQ(RootFinder::sign_changes(-5 * X - 1), 0);
    EXPECT_EQ(RootFinder::sign_changes(5 * X + 1), 0);
    EXPECT_EQ(RootFinder::sign_changes(-5 * X + 1), 1);
    EXPECT_EQ(RootFinder::sign_changes(5 * X - 1), 1);
    EXPECT_EQ(RootFinder::sign_changes(5 * X.pow(5) - 1 * X.pow(4) - X.pow(3) - X.pow(2) + X - 1), 3);
}

TEST(RealPolynomialRootFinderTests, CauchyBoundTest) {
    EXPECT_EQ(RootFinder::cauchy_bound(3 * X.pow(4) - 6 * X.pow(3) - 2 * X.pow(2) - 9), 4);
}

TEST(RealPolynomialRootFinderTests, SturmSequenceTest) {
    auto p = X.pow(4) + X.pow(3) - X - 1;
    auto seq = RootFinder::sturm_sequence(p);

    EXPECT_EQ(seq.size(), 5);
    EXPECT_EQ(seq[0], p);
    EXPECT_EQ(seq[1], p.derive());
    EXPECT_EQ(seq[2], (3. / 16.) * X.pow(2) + (3. / 4.) * X + (15. / 16.));
    EXPECT_EQ(seq[3], -32 * X + -64);
    EXPECT_EQ(seq[4], -(3. / 16.) * one);

    EXPECT_TRUE(seq.back().is_constant());
}

TEST(RealPolynomialRootFinderTests, SignVariationsTest) {
    auto p = X.pow(4) + X.pow(3) - X - 1;
    auto seq = RootFinder::sturm_sequence(p);

    auto variations = RootFinder::sign_variations(seq, -10);
    EXPECT_EQ(variations, std::vector<short>({+1, -1, +1, +1, -1}));

    variations = RootFinder::sign_variations(seq, 10);
    EXPECT_EQ(variations, std::vector<short>({+1, +1, +1, -1, -1}));
}

TEST(RealPolynomialRootFinderTests, NumberDistinctRootsTest) {
    auto p = X.pow(4) + X.pow(3) - X - 1;
    EXPECT_EQ(RootFinder::number_distinct_roots(p), 2);

    p = Polynomial::from_roots({-2.5, -1.15, 0., .5, 1.25, 4.125, 6.5}); // 7 roots
    EXPECT_EQ(RootFinder::number_distinct_roots(p), 7);
}

TEST(RealPolynomialRootFinderTests, NonSquareFreePolynomialTest) {
    auto p = (X - 2).pow(3); // non square-free polynomial
    EXPECT_EQ(RootFinder::number_distinct_roots(p), std::numeric_limits<int>::quiet_NaN());
}

TEST(RealPolynomialRootFinderTests, RootIsolationTest) {
    auto p = Polynomial::from_roots({-2.5, -1.15, 0., .5, 1.25, 4.125, 6.5}); // 7 roots (is square-free polynomial)
    auto intervals = RootFinder::root_isolation(p);

    EXPECT_EQ(intervals.size(), 7);
    for (auto I: intervals) {
        EXPECT_EQ(RootFinder::number_distinct_roots(p, I), 1);
    }
}

TEST(RealPolynomialRootFinderTests, FindRootsTest) {
    auto p = Polynomial::from_roots({-2.5, -1.15, 0., .5, 1.25, 4.125, 6.5});
    auto roots = RootFinder::find_roots(p, 1e-5);

    EXPECT_EQ(roots.size(), 7);
    EXPECT_TRUE(p.has_roots(roots));
}