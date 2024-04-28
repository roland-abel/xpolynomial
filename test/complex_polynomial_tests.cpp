/// @file polynomial_tests.cpp
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
#include <complex>
#include "test_utilities.h"
#include "complex_polynomial.h"

using namespace xmath;

namespace {
    using ComplexPolynomial = complex_polynomial<double>;
    using RealPolynomial = polynomial<double>;

    constexpr auto i = std::complex(0., 1.);
    constexpr auto I = ComplexPolynomial::value_type(i);

    constexpr double epsilon = ComplexPolynomial::epsilon;
    const auto zero = ComplexPolynomial::zero();
    const auto one = ComplexPolynomial::one();

    const auto X = RealPolynomial::monomial(1, 1.0);
    const auto Y = RealPolynomial::monomial(1, 1.0);
    const auto Z = ComplexPolynomial::monomial(1, 1.0);
}

TEST(ComplexPolynomialTests, ValueTypeCheck) {
    auto z = ComplexPolynomial::value_type(2. + 4.4 * I);

    EXPECT_NEAR(z.real(), std::complex(2., 4.4).real(), epsilon);
    EXPECT_NEAR(z.imag(), std::complex(2., 4.4).imag(), epsilon);
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

TEST(ComplexPolynomialTests, EvaluateTest) {
    const auto p = Z.pow(2) - one;
    EXPECT_COMPLEX_NEAR(p(1), 0. + 0. * i, epsilon);
    EXPECT_COMPLEX_NEAR(p(-1.), 0. + 0. * i, epsilon);
}
