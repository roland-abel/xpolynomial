/// @file chebyshev_polynomial_tests.h
/// @brief
///
/// @author Roland Abel
/// @date October 14, 2023
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
#include "chebyshev_polynomial.h"

using namespace xmath;

namespace {
    using std::numbers::pi;
    using Polynomial = polynomial<double>;
    using Interval = real_interval<double>;
    using ChebyshevPolynomial = chebyshev_polynomial<double>;

    constexpr auto epsilon = Polynomial::epsilon;
    auto zero = Polynomial::zero();
    auto one = Polynomial::one();
    auto X = Polynomial::monomial(1, 1.0);
}

TEST(ChebyshvPolynomialTest, FirstKindChebyshvPolynomialsTest) {
    auto T_0 = ChebyshevPolynomial::create_1st_kind(0);
    EXPECT_EQ(T_0, one);

    auto T_1 = ChebyshevPolynomial::create_1st_kind(1);
    EXPECT_EQ(T_1, X);

    auto T_2 = ChebyshevPolynomial::create_1st_kind(2);
    EXPECT_EQ(T_2, 2 * X.pow(2) - 1);

    auto T_3 = ChebyshevPolynomial::create_1st_kind(3);
    EXPECT_EQ(T_3, 4 * X.pow(3) - 3 * X);

    auto T_4 = ChebyshevPolynomial::create_1st_kind(4);
    EXPECT_EQ(T_4, 8 * X.pow(4) - 8 * X.pow(2) + 1);

    auto T_5 = ChebyshevPolynomial::create_1st_kind(5);
    EXPECT_EQ(T_5, 16 * X.pow(5) - 20 * X.pow(3) + 5 * X);

    auto T_6 = ChebyshevPolynomial::create_1st_kind(6);
    EXPECT_EQ(T_6, 32 * X.pow(6) - 48 * X.pow(4) + 18 * X.pow(2) - 1);

    auto T_7 = ChebyshevPolynomial::create_1st_kind(7);
    EXPECT_EQ(T_7, 64 * X.pow(7) - 112 * X.pow(5) + 56 * X.pow(3) - 7 * X);

    auto T_8 = ChebyshevPolynomial::create_1st_kind(8);
    EXPECT_EQ(T_8, 128 * X.pow(8) - 256 * X.pow(6) + 160 * X.pow(4) - 32 * X.pow(2) + 1);

    auto T_9 = ChebyshevPolynomial::create_1st_kind(9);
    EXPECT_EQ(T_9, 256 * X.pow(9) - 576 * X.pow(7) + 432 * X.pow(5) - 120 * X.pow(3) + 9 * X);

    auto T_10 = ChebyshevPolynomial::create_1st_kind(10);
    EXPECT_EQ(T_10, 512 * X.pow(10) - 1280 * X.pow(8) + 1120 * X.pow(6) - 400 * X.pow(4) + 50 * X.pow(2) - 1);
}

TEST(ChebyshvPolynomialTest, FirstKindChebyshvPolynomialsWithEmptyCacheTest) {
    ChebyshevPolynomial::polynomial_sequence chebyshev_cache;

    auto T_0 = ChebyshevPolynomial::create_1st_kind(0, chebyshev_cache);
    EXPECT_EQ(T_0, one);

    auto T_1 = ChebyshevPolynomial::create_1st_kind(1, chebyshev_cache);
    EXPECT_EQ(T_1, X);

    auto T_2 = ChebyshevPolynomial::create_1st_kind(2, chebyshev_cache);
    EXPECT_EQ(T_2, 2 * X.pow(2) - 1);

    auto T_9 = ChebyshevPolynomial::create_1st_kind(9, chebyshev_cache);
    EXPECT_EQ(T_9, 256 * X.pow(9) - 576 * X.pow(7) + 432 * X.pow(5) - 120 * X.pow(3) + 9 * X);

    auto T_10 = ChebyshevPolynomial::create_1st_kind(10, chebyshev_cache);
    EXPECT_EQ(T_10, 512 * X.pow(10) - 1280 * X.pow(8) + 1120 * X.pow(6) - 400 * X.pow(4) + 50 * X.pow(2) - 1);

    EXPECT_EQ(11, chebyshev_cache.size());
}

TEST(ChebyshvPolynomialTest, Roots1stKindZeroOrderTest) {
    auto roots = ChebyshevPolynomial::chebyshev_nodes(0, Interval(-1., 1.));
    EXPECT_EQ(roots.size(), 0);
}

TEST(ChebyshvPolynomialTest, ChebyshevNodesTest) {
    size_t n = 15;
    auto nodes = ChebyshevPolynomial::chebyshev_nodes(n, Interval(-1., 1.));
    auto T_n = ChebyshevPolynomial::create_1st_kind(n);

    EXPECT_EQ(nodes.size(), n);
    EXPECT_TRUE(T_n.has_roots(nodes));
}

TEST(ChebyshvPolynomialTest, ChebyshevNodesForIntervalTest) {
    const auto n = 5;
    auto I = Interval(-2.0, 3.0);

    const auto A = .5 * (I.lower() + I.upper());
    const auto B = .5 * (I.upper() - I.lower());

    auto chebyshev_nodes = ChebyshevPolynomial::chebyshev_nodes(n, Interval(-1., 1.));
    const auto nodes = ChebyshevPolynomial::chebyshev_nodes(n, I);

    ASSERT_EQ(n, nodes.size());
    for (size_t i = 0; i < n; ++i) {
        EXPECT_NEAR(nodes[i], A + B * chebyshev_nodes[i], epsilon);
    }
}

TEST(ChebyshvPolynomialTest, CheckChebyshv1stKindFormularTest) {
    // Check: T_n(cos(x)) == cos(n * x)
    auto x = 3.2;
    auto n = 20;

    auto T_n = ChebyshevPolynomial::create_1st_kind(n);
    EXPECT_NEAR(T_n(std::cos(x)), std::cos(n * x), epsilon);
}

TEST(ChebyshvPolynomialTest, ClenshawTest) {
    auto alphas = std::vector<double>{-1., 1.25, 2.5, -3.5, 4.2};
    auto chebyshev = ChebyshevPolynomial::chebyshev_series(alphas);
    auto xs = std::vector<double>{-1., -.5, 0., .5, 1.};

    for (auto const &x: xs) {
        EXPECT_NEAR(chebyshev(x), ChebyshevPolynomial::clenshaw(alphas, x), epsilon);
    }
}

TEST(ChebyshvPolynomialTest, ChebyshevGaussQuadratureForMonomialsTest) {
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature(one), pi, epsilon);
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature(X), 0., epsilon);
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature(X.pow(2)), pi / 2., epsilon);
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature(X.pow(3)), 0., epsilon);
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature(X.pow(4)), 3. * pi / 8., epsilon);
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature(X.pow(5)), 0., epsilon);
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature(X.pow(6)), 5. * pi / 16., epsilon);
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature(X.pow(7)), 0., epsilon);
}

TEST(ChebyshvPolynomialTest, ChebyshevGaussQuadratureTest) {
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature([](auto x) { return std::sin(x - 1.); }), -2.02285, epsilon);
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature([](auto x) { return std::exp(x); }), 3.97746, epsilon);
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature([](auto x) { return std::exp(x) - 1.; }), 0.835871, epsilon);
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature([](auto x) { return std::log(x + 4.); }), 4.30489, epsilon);
}

TEST(ChebyshvPolynomialTest, ChebyshevGaussQuadratureOverIntervalTest) {
    EXPECT_NEAR(ChebyshevPolynomial::chebyshev_quadrature(X.pow(2), Interval(-0.5, 0.5)), pi / 2., epsilon);
}
