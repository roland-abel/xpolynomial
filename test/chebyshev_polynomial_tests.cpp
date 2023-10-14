/// @file legendre_polynomial_tests.h
///
/// @author Roland Abel
/// @date 14.10.2023

#include <gtest/gtest.h>
#include "chebyshev_polynomial.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using ChebyshvPolynomial = chebyshev_polynomial<double>;

    constexpr auto epsilon = Polynomial::tolerance;
    auto zero = Polynomial::zero();
    auto one = Polynomial::one();
    auto X = Polynomial::monomial(1, 1.0);
}

TEST(ChebyshvPolynomialTest, FirstKindChebyshvPolynomials) {
    auto T_0 = ChebyshvPolynomial::create_1st_kind(0);
    EXPECT_EQ(T_0, one);

    auto T_1 = ChebyshvPolynomial::create_1st_kind(1);
    EXPECT_EQ(T_1, X);

    auto T_2 = ChebyshvPolynomial::create_1st_kind(2);
    EXPECT_EQ(T_2, 2 * X.pow(2) - 1);

    auto T_3 = ChebyshvPolynomial::create_1st_kind(3);
    EXPECT_EQ(T_3, 4 * X.pow(3) - 3 * X);

    auto T_4 = ChebyshvPolynomial::create_1st_kind(4);
    EXPECT_EQ(T_4, 8 * X.pow(4) - 8 * X.pow(2) + 1);

    auto T_5 = ChebyshvPolynomial::create_1st_kind(5);
    EXPECT_EQ(T_5, 16 * X.pow(5) - 20 * X.pow(3) + 5 * X);

    auto T_6 = ChebyshvPolynomial::create_1st_kind(6);
    EXPECT_EQ(T_6, 32 * X.pow(6) - 48 * X.pow(4) + 18 * X.pow(2) - 1);

    auto T_7 = ChebyshvPolynomial::create_1st_kind(7);
    EXPECT_EQ(T_7, 64 * X.pow(7) - 112 * X.pow(5) + 56 * X.pow(3) - 7 * X);

    auto T_8 = ChebyshvPolynomial::create_1st_kind(8);
    EXPECT_EQ(T_8, 128 * X.pow(8) - 256 * X.pow(6) + 160 * X.pow(4) - 32 * X.pow(2) + 1);

    auto T_9 = ChebyshvPolynomial::create_1st_kind(9);
    EXPECT_EQ(T_9, 256 * X.pow(9) - 576 * X.pow(7) + 432 * X.pow(5) - 120 * X.pow(3) + 9 * X);

    auto T_10 = ChebyshvPolynomial::create_1st_kind(10);
    EXPECT_EQ(T_10, 512 * X.pow(10) - 1280 * X.pow(8) + 1120 * X.pow(6) - 400 * X.pow(4) + 50 * X.pow(2) - 1);
}
