/// @file real_polynomial_root_finder_tests.cpp
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
#include "test_utilities.h"
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
    EXPECT_FALSE(RootFinder::quadratic_roots(X + 1).has_value());
    EXPECT_FALSE(RootFinder::quadratic_roots(zero).has_value());
    EXPECT_FALSE(RootFinder::quadratic_roots(X.pow(3)).has_value());
}

TEST(RealPolynomialRootFinderTests, QuadraticPolynomialWhitoutRealRootsTest) {
    EXPECT_FALSE(RootFinder::quadratic_roots(X.pow(2) + 1).has_value());
}

TEST(RealPolynomialRootFinderTests, QuadraticPolynomialWithTwoRootsTest) {
    auto p = 4 * X.pow(2) + .5 * X - 4;
    auto roots = RootFinder::quadratic_roots(p).value();

    EXPECT_NEAR(get<0>(roots), (-.5 + std::sqrt(64.25)) / 8., epsilon);
    EXPECT_NEAR(get<1>(roots), (-.5 - std::sqrt(64.25)) / 8., epsilon);
}

TEST(RealPolynomialRootFinderTests, QuadraticPolynomialWithOneRootTest) {
    auto p = (X - 0.3).pow(2);
    auto roots = RootFinder::quadratic_roots(p).value();

    EXPECT_NEAR(get<0>(roots), .3, epsilon);
    EXPECT_NEAR(get<1>(roots), .3, epsilon);
}

TEST(RealPolynomialRootFinderTests, CubicPolynomialWithOneRootsTest) {
    auto p = X.pow(3) - 5.;
    auto roots = RootFinder::cubic_roots(p);

    EXPECT_EQ(roots.size(), 1);
    EXPECT_NEAR(roots[0], std::pow(5., 1. / 3.), epsilon);
}

TEST(RealPolynomialRootFinderTests, CubicNormalFormPolynomialWithThreeUnequalRootsTest) {
    auto p = X.pow(3) - 2 * X + 1;

    EXPECT_EQ((X - 1) * (X.pow(2) + X - 1), p);
    EXPECT_TRUE(p.has_roots({1., -.5 + std::sqrt(5. / 4.), -.5 - std::sqrt(5. / 4.)}));

    auto roots = RootFinder::cubic_roots(p);
    EXPECT_NEAR(roots[0], 1., epsilon);
    EXPECT_NEAR(roots[1], -.5 - std::sqrt(5. / 4.), epsilon);
    EXPECT_NEAR(roots[2], -.5 + std::sqrt(5. / 4.), epsilon);
}

TEST(RealPolynomialRootFinderTests, CubicPolynomialWithThreeUnequalRootsTest) {
    auto p = (X - 2.) * (X + 4.) * (X - 5.);
    EXPECT_TRUE(p.has_roots({5., 2., -4}));

    auto roots = RootFinder::cubic_roots(p);
    EXPECT_NEAR(roots[0], 5., epsilon);
    EXPECT_NEAR(roots[1], -4., epsilon);
    EXPECT_NEAR(roots[2], 2., epsilon);
}

TEST(RealPolynomialRootFinderTests, CubicPolynomialWithTwoEqualRootsTest) {
    auto p = 3.5 * (X - 7.).pow(2) * (X - 4.);
    EXPECT_TRUE(p.has_roots({7., 4.}));

    auto roots = RootFinder::cubic_roots(p);
    EXPECT_NEAR(roots[0], 4., epsilon);
    EXPECT_NEAR(roots[1], 7., epsilon);
    EXPECT_NEAR(roots[2], 7., epsilon);
}

TEST(RealPolynomialRootFinderTests, CubicPolynomialWithThreeEqualRootsTest) {
    auto p = 2.25 * X.pow(3);

    auto roots = RootFinder::cubic_roots(p);
    EXPECT_NEAR(roots[0], 0., epsilon);
    EXPECT_NEAR(roots[0], 0., epsilon);
    EXPECT_NEAR(roots[0], 0., epsilon);

    p = -1.5 * (X + 7.5).pow(3);
    EXPECT_TRUE(p.has_roots({-7.5}));

    roots = RootFinder::cubic_roots(p);
    EXPECT_NEAR(roots[0], -7.5, epsilon);
    EXPECT_NEAR(roots[0], -7.5, epsilon);
    EXPECT_NEAR(roots[0], -7.5, epsilon);
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
    EXPECT_NEAR(RootFinder::cauchy_bounds(3 * X.pow(4) - 6 * X.pow(3) - 2 * X.pow(2) - 9).value(), 4., epsilon);
}

TEST(RealPolynomialRootFinderTests, LagrangesBoundsTest) {
    EXPECT_NEAR(RootFinder::lagrange_bounds(3 * X.pow(4) - 6 * X.pow(3) - 2 * X.pow(2) - 9).value(), 17. / 3., epsilon);
    EXPECT_NEAR(RootFinder::lagrange_bounds(.1 * X.pow(4) - 6 * X.pow(3) - 2 * X.pow(2) - 9).value(), 170., epsilon);
    EXPECT_NEAR(RootFinder::lagrange_bounds(100 * X.pow(4) - 6 * X.pow(3) - 2 * X.pow(2) - 9).value(), 1., epsilon);
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
    EXPECT_EQ(RootFinder::number_distinct_roots(p).value(), 2);

    p = Polynomial::from_roots({-2.5, -1.15, 0., .5, 1.25, 4.125, 6.5}); // 7 roots
    EXPECT_EQ(RootFinder::number_distinct_roots(p), 7);
}

TEST(RealPolynomialRootFinderTests, NonSquareFreePolynomialTest) {
    auto p = (X - 2).pow(3); // non square-free polynomial
    EXPECT_FALSE(RootFinder::number_distinct_roots(p).has_value());
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
    auto p = Polynomial::from_roots({-2, 0, 0, -1, 1});
    ASSERT_TRUE(p.is_integer());

    auto [roots, _] = RootFinder::find_roots(p);

    EXPECT_EQ(roots.size(), 4);
    EXPECT_TRUE(p.has_roots(roots));
    EXPECT_UNIQUE(roots, epsilon);
}

TEST(RealPolynomialRootFinderTests, FindRootsForPolynomialWithMultipleRootsTest) {
    auto p = Polynomial::from_roots({-2, -1, 1, 5, 5, 5});
    ASSERT_TRUE(p.is_integer());

    auto [roots, multiplicities] = RootFinder::find_roots(p);

    EXPECT_EQ(roots.size(), 4);
    EXPECT_TRUE(p.has_roots(roots));
    EXPECT_EQ(multiplicities, std::vector<unsigned short>({1, 1, 1, 3}));
}

TEST(RealPolynomialRootFinderTests, NewtonRaphsonTest) {
    auto p = X.pow(2) - 2;
    EXPECT_NEAR(RootFinder::newton_raphson(p, 1.1).value(), std::sqrt(2), epsilon);
}

TEST(RealPolynomialRootFinderTests, FindChebyshevRootsTest) {
    const auto n = 17;

    const auto T_n = ChebyshevPolynomial::create_1st_kind(n);
    const auto [roots, multiplicities] = RootFinder::find_roots(T_n);
    const auto nodes = ChebyshevPolynomial::chebyshev_nodes(n);

    EXPECT_EQ(roots.size(), n);
    EXPECT_TRUE(T_n.has_roots(roots));
    EXPECT_UNIQUE(roots, epsilon);
}