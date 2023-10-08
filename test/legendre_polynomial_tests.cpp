//
// Created by abel on 19.08.2023.
//

#include <gtest/gtest.h>
#include "legendre_polynomial.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using LegendrePolynomial = legendre_polynomial<double>;

    constexpr auto tolerance = Polynomial::tolerance;
    auto zero = Polynomial::zero();
    auto one = Polynomial::one();
    auto X = Polynomial::monomial(1, 1.0);
}

TEST(LegendrePolynomialTests, TheFirstLegendrePolynomials) {
    auto p0 = LegendrePolynomial(0);
    EXPECT_EQ(p0, one);

    auto p1 = LegendrePolynomial(1);
    EXPECT_EQ(p1, X);

    auto p2 = LegendrePolynomial(2);
    EXPECT_EQ(p2, .5 * (3 * X.pow(2) - 1));

    auto p3 = LegendrePolynomial(3);
    EXPECT_EQ(p3, .5 * (5 * X.pow(3) - 3 * X));

    auto p4 = LegendrePolynomial(4);
    EXPECT_EQ(p4, (1. / 8) * (35. * X.pow(4) - 30 * X.pow(2) + 3));

    auto p5 = LegendrePolynomial(5);
    EXPECT_EQ(p5, (1. / 8) * (63 * X.pow(5) - (70 * X.pow(3)) + (15 * X)));

    auto p6 = LegendrePolynomial(6);
    EXPECT_EQ(p6, (1. / 16) * (231 * X.pow(6) - 315 * X.pow(4) + 105 * X.pow(2) - 5));

    auto p7 = LegendrePolynomial(7);
    EXPECT_EQ(p7, (1. / 16) * (429 * X.pow(7) - 693 * X.pow(5) + 315 * X.pow(3) - 35 * X));

    auto p8 = LegendrePolynomial(8);
    EXPECT_EQ(p8, (1. / 128) * (6435 * X.pow(8) - 12012 * X.pow(6) + 6930 * X.pow(4) - 1260 * X.pow(2) + 35));
}

TEST(LegendrePolynomialTests, EvaluationLegendrePolynomial) {
}