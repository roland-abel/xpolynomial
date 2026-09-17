/// @file legendre_polynomial_tests.cpp
/// @brief
///
/// @author Roland Abel
/// @date October 8, 2023
///
/// Copyright (c) 2026 Roland Abel

#include <gtest/gtest.h>
#include <cmath>
#include <xpolynomial.h>

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using LegendrePolynomial = legendre_polynomial<double>;

    constexpr auto epsilon = Polynomial::epsilon;
    const auto zero = Polynomial::zero();
    const auto one = Polynomial::one();
    const auto X = Polynomial::monomial(1, 1.0);
}

TEST(LegendrePolynomialTests, LegendrePolynomials) {
    auto p0 = LegendrePolynomial::create(0);
    EXPECT_EQ(p0, one);

    auto p1 = LegendrePolynomial::create(1);
    EXPECT_EQ(p1, X);

    auto p2 = LegendrePolynomial::create(2);
    EXPECT_EQ(p2, .5 * (3 * X.pow(2) - 1));

    auto p3 = LegendrePolynomial::create(3);
    EXPECT_EQ(p3, .5 * (5 * X.pow(3) - 3 * X));

    auto p4 = LegendrePolynomial::create(4);
    EXPECT_EQ(p4, 1. / 8 * (35. * X.pow(4) - 30 * X.pow(2) + 3));

    auto p5 = LegendrePolynomial::create(5);
    EXPECT_EQ(p5, 1. / 8 * (63 * X.pow(5) - (70 * X.pow(3)) + (15 * X)));

    auto p6 = LegendrePolynomial::create(6);
    EXPECT_EQ(p6, 1. / 16 * (231 * X.pow(6) - 315 * X.pow(4) + 105 * X.pow(2) - 5));

    auto p7 = LegendrePolynomial::create(7);
    EXPECT_EQ(p7, 1. / 16 * (429 * X.pow(7) - 693 * X.pow(5) + 315 * X.pow(3) - 35 * X));

    auto p8 = LegendrePolynomial::create(8);
    EXPECT_EQ(p8, (1. / 128) * (6435 * X.pow(8) - 12012 * X.pow(6) + 6930 * X.pow(4) - 1260 * X.pow(2) + 35));
}

TEST(LegendrePolynomialTests, LegendrePolynomialsSkipUncachedOrders) {
    const auto p = LegendrePolynomial::create(18);

    EXPECT_EQ(p.degree(), 18);
    for (const auto x : {-0.75, -0.25, 0.0, 0.25, 0.75}) {
        EXPECT_NEAR(p(x), std::legendre(18, x), epsilon);
    }

    for (size_t order = 9; order <= 18; ++order) {
        const auto cached = LegendrePolynomial::create(order);

        EXPECT_EQ(cached.degree(), order);
        for (const auto x : {-0.75, -0.25, 0.0, 0.25, 0.75}) {
            EXPECT_NEAR(cached(x), std::legendre(static_cast<unsigned int>(order), x), epsilon);
        }
    }

    EXPECT_EQ(LegendrePolynomial::create(18), p);
}

TEST(LegendrePolynomialTests, LegendrePolynomialsBeyondPrecomputedValues) {
    for (size_t order = 9; order <= 15; ++order) {
        const auto p = LegendrePolynomial::create(order);

        EXPECT_NEAR(p(1.0), 1.0, epsilon);
        EXPECT_NEAR(p(-1.0), (order % 2 == 0) ? 1.0 : -1.0, epsilon);

        if (order % 2 == 1) {
            EXPECT_NEAR(p(0.0), 0.0, epsilon);
        }
    }
}

TEST(LegendrePolynomialTests, LegendrePolynomialsSatisfyRecurrence) {
    const auto x = 0.3;
    for (size_t order = 10; order <= 15; ++order) {
        const auto p_n = LegendrePolynomial::create(order);
        const auto p_nm1 = LegendrePolynomial::create(order - 1);
        const auto p_nm2 = LegendrePolynomial::create(order - 2);

        const auto n = static_cast<double>(order);
        const auto expected = (2. * n - 1.) / n * x * p_nm1(x) - (n - 1.) / n * p_nm2(x);

        EXPECT_NEAR(p_n(x), expected, epsilon);
    }
}
