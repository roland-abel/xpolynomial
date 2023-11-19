/// @file real_polynomial_root_finder_tests.h
///
/// @author Roland Abel
/// @date 08.10.2023

#include <gtest/gtest.h>
#include <cmath>
#include "polynomial.h"
#include "real_polynomial_root_finder.h"
#include "chebyshev_polynomial.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using RootFinder = real_polynomial_root_finder<double>;
    using ChebyshevPolynomial = chebyshev_polynomial<double>;
    using value_type = Polynomial::spec::value_type;

    constexpr auto epsilon = Polynomial::spec::epsilon;
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

    EXPECT_NEAR(get<0>(roots), (-.5 + std::sqrt(64.25)) / 8., epsilon);
    EXPECT_NEAR(get<1>(roots), (-.5 - std::sqrt(64.25)) / 8., epsilon);
}

TEST(RealPolynomialRootFinderTests, QuadraticPolynomialWithOneRootTest) {
    auto p = (X - 0.3).pow(2);
    auto roots = RootFinder::quadratic_roots(p);

    EXPECT_NEAR(get<0>(roots), .3, epsilon);
    EXPECT_NEAR(get<1>(roots), .3, epsilon);
}

TEST(RealPolynomialRootFinderTests, CubicPolynomialWithOneRootsTest) {
    auto p = X.pow(3) - 5.;
    auto roots = RootFinder::cubic_roots(p);

    EXPECT_NEAR(get<0>(roots), std::pow(5., 1. / 3.), epsilon);
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
    EXPECT_NEAR(r1, 1., epsilon);
    EXPECT_NEAR(r2, -.5 - std::sqrt(5. / 4.), epsilon);
    EXPECT_NEAR(r3, -.5 + std::sqrt(5. / 4.), epsilon);
}

TEST(RealPolynomialRootFinderTests, CubicPolynomialWithThreeUnequalRootsTest) {
    auto p = (X - 2.) * (X + 4.) * (X - 5.);
    EXPECT_TRUE(p.has_roots({5., 2., -4}));

    auto [r1, r2, r3] = RootFinder::cubic_roots(p);
    EXPECT_NEAR(r1, 5., epsilon);
    EXPECT_NEAR(r2, -4., epsilon);
    EXPECT_NEAR(r3, 2., epsilon);
}

TEST(RealPolynomialRootFinderTests, CubicPolynomialWithTwoEqualRootsTest) {
    auto p = 3.5 * (X - 7.).pow(2) * (X - 4.);
    EXPECT_TRUE(p.has_roots({7., 4.}));

    auto [r1, r2, r3] = RootFinder::cubic_roots(p);
    EXPECT_NEAR(r1, 4., epsilon);
    EXPECT_NEAR(r2, 7., epsilon);
    EXPECT_NEAR(r3, 7., epsilon);
}

TEST(RealPolynomialRootFinderTests, CubicPolynomialWithThreeEqualRootsTest) {
    auto p = 2.25 * X.pow(3);

    auto [r1, r2, r3] = RootFinder::cubic_roots(p);
    EXPECT_NEAR(r1, 0., epsilon);
    EXPECT_NEAR(r2, 0., epsilon);
    EXPECT_NEAR(r3, 0., epsilon);

    p = -1.5 * (X + 7.5).pow(3);
    EXPECT_TRUE(p.has_roots({-7.5}));

    auto roots = RootFinder::cubic_roots(p);
    EXPECT_NEAR(get<0>(roots), -7.5, epsilon);
    EXPECT_NEAR(get<1>(roots), -7.5, epsilon);
    EXPECT_NEAR(get<2>(roots), -7.5, epsilon);
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

TEST(RealPolynomialRootFinderTests, CauchysBoundsTest) {
    EXPECT_NEAR(RootFinder::cauchy_bounds(3 * X.pow(4) - 6 * X.pow(3) - 2 * X.pow(2) - 9), 4., epsilon);
}

TEST(RealPolynomialRootFinderTests, LagrangesBoundsTest) {
    EXPECT_NEAR(RootFinder::lagrange_bounds(3 * X.pow(4) - 6 * X.pow(3) - 2 * X.pow(2) - 9), 17. / 3., epsilon);
    EXPECT_NEAR(RootFinder::lagrange_bounds(.1 * X.pow(4) - 6 * X.pow(3) - 2 * X.pow(2) - 9), 170., epsilon);
    EXPECT_NEAR(RootFinder::lagrange_bounds(100 * X.pow(4) - 6 * X.pow(3) - 2 * X.pow(2) - 9), 1., epsilon);
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

TEST(RealPolynomialRootFinderTests, NonSquareRootIsolationTest) {
    auto p = (X - 2).pow(3); // non square-free polynomial
    auto intervals = RootFinder::root_isolation(p);

    EXPECT_TRUE(intervals.empty());
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
    auto [roots, _] = RootFinder::find_roots(p);

    EXPECT_EQ(roots.size(), 7);
    EXPECT_TRUE(p.has_roots(roots));
}

TEST(RealPolynomialRootFinderTests, FindRootsForPolynomialWithMultipleRootsTest) {
    auto p = Polynomial::from_roots({-2.5, -1.15, 1.15, .5, .5, .5});
    auto [roots, multiplicities] = RootFinder::find_roots(p);

    EXPECT_EQ(roots.size(), 4);
    EXPECT_TRUE(p.has_roots(roots));
    EXPECT_EQ(multiplicities, std::vector<unsigned short>({1, 1, 1, 3}));
}

TEST(RealPolynomialRootFinderTests, NewtonRaphsonTest) {
    auto p = X.pow(2) - 2;
    EXPECT_NEAR(RootFinder::newton_raphson(p, 1.1), std::sqrt(2), epsilon);
}

TEST(RealPolynomialRootFinderTests, FindChebyshevRootsTest) {
    const auto n = 17;
    const auto precision = 1e-11;

    const auto T_n = ChebyshevPolynomial::create_1st_kind(n);
    const auto [roots, multiplicities] = RootFinder::find_roots(T_n, precision);
    const auto nodes = ChebyshevPolynomial::chebyshev_nodes(n);

    EXPECT_EQ(roots.size(), n);
    EXPECT_TRUE(T_n.has_roots(roots));
}