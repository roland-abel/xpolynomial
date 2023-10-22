/// @file polynomial_tests.h
///
/// @author Roland Abel
/// @date 08.10.2023

#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include "complex_polynomial.h"

using namespace xmath;

namespace {
    using ComplexPolynomial = complex_polynomial<double>;
    using RealPolynomial = polynomial<double, polynomial_specification<double>>;

    constexpr auto I = ComplexPolynomial::value_type(0., 1.);

    constexpr double epsilon = ComplexPolynomial::epsilon;
    auto zero = ComplexPolynomial::zero();
    auto one = ComplexPolynomial::one();

    auto X = RealPolynomial::monomial(1, 1.0);
    auto Y = RealPolynomial::monomial(1, 1.0);
    auto Z = ComplexPolynomial::monomial(1, 1.0);
}

TEST(ComplexPolynomialTests, ValueTypeCheck) {
    auto z = ComplexPolynomial::value_type(2. + 4.4 * I);

    EXPECT_NEAR(z.real(), std::complex<double>(2., 4.4).real(), epsilon);
    EXPECT_NEAR(z.imag(), std::complex<double>(2., 4.4).imag(), epsilon);
}

TEST(ComplexPolynomialTests, DefaultConstructor) {
    auto p = ComplexPolynomial();
    EXPECT_EQ(p.degree(), 0);
    EXPECT_EQ(p, zero);
}

TEST(ComplexPolynomialTests, ZeroPloynomialTest) {
    EXPECT_TRUE(ComplexPolynomial::zero().is_zero());
    EXPECT_TRUE(ComplexPolynomial().is_zero());
    EXPECT_TRUE(ComplexPolynomial({0, 0, 0, 0}).is_zero());

    EXPECT_EQ(ComplexPolynomial::zero(), zero);
    EXPECT_EQ(ComplexPolynomial::zero(), ComplexPolynomial());
    EXPECT_EQ(ComplexPolynomial::zero(), ComplexPolynomial({0, 0, 0, 0}));
    EXPECT_EQ(ComplexPolynomial::zero().degree(), 0);
}

TEST(ComplexPolynomialTests, OnePloynomialTest) {
    EXPECT_TRUE(ComplexPolynomial::one().is_constant());
    EXPECT_TRUE(ComplexPolynomial::one().is_linear());
    EXPECT_EQ(ComplexPolynomial::one(), one);
    EXPECT_EQ(ComplexPolynomial::one(), ComplexPolynomial({std::complex<double>(1, 0)}));
    EXPECT_EQ(ComplexPolynomial::one().degree(), 0);
}

TEST(ComplexPolynomialTests, ConstructorWithCoefficients) {
    auto p = ComplexPolynomial({1. + I, 2. - 3. * I, 4. - 6. * I});

    EXPECT_EQ(p.degree(), 2);
    EXPECT_NEAR(p[0].real(), 1., epsilon);
    EXPECT_NEAR(p[0].imag(), 1., epsilon);
    EXPECT_NEAR(p[1].real(), 2., epsilon);
    EXPECT_NEAR(p[1].imag(), -3., epsilon);
    EXPECT_NEAR(p[2].real(), 4., epsilon);
    EXPECT_NEAR(p[2].imag(), -6., epsilon);
}

TEST(ComplexPolynomialTests, SeparateComplexPolynomialTest) {
    auto p = (Z - I).pow(3);
    EXPECT_EQ(p.degree(), 3);

    auto [real, imag] = separate(p);
    EXPECT_EQ(real, X.pow(3) - 3 * X);
    EXPECT_EQ(imag, -3 * Y.pow(2) + 1);

    EXPECT_EQ(p, real + I * imag);
}