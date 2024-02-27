/// @file polynomial_interpolation_tests.cpp
/// @brief
///
/// @author Roland Abel
/// @date November 25, 2023
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
#include "polynomial_interpolation.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using Interpolation = polynomial_interpolation<double>;

    constexpr double epsilon = Polynomial::epsilon;
}

TEST(PolynomialInterpolation, LagrangeBasisFunctionsTest) {
    const std::vector<double> xs = {1.0, 2.5, 3., 5.};
    const auto basis = Interpolation::lagrange_basis(xs);

    EXPECT_EQ(basis.size(), 4);

    EXPECT_EQ(basis[0].degree(), 3);
    EXPECT_EQ(basis[1].degree(), 3);
    EXPECT_EQ(basis[2].degree(), 3);
    EXPECT_EQ(basis[3].degree(), 3);

    EXPECT_TRUE(basis[0].has_roots({2.5, 3., 5.}));
    EXPECT_TRUE(basis[1].has_roots({1.0, 3., 5.}));
    EXPECT_TRUE(basis[2].has_roots({1.0, 2.5, 5.}));
    EXPECT_TRUE(basis[3].has_roots({1.0, 2.5, 3.}));

    EXPECT_NEAR(basis[0](xs[0]), 1.0, epsilon);
    EXPECT_NEAR(basis[1](xs[1]), 1.0, epsilon);
    EXPECT_NEAR(basis[2](xs[2]), 1.0, epsilon);
    EXPECT_NEAR(basis[3](xs[3]), 1.0, epsilon);
}

TEST(PolynomialInterpolation, LagrangeInterpolationTest) {
    const std::vector x_values = {1.0, 2.0, 3.0};
    const std::vector y_values = {2.0, 4.0, 6.0};
    const auto p = Interpolation::lagrange_interpolation(x_values, y_values).value();

    EXPECT_DOUBLE_EQ(p(1.0), 2.0);
    EXPECT_DOUBLE_EQ(p(2.0), 4.0);
    EXPECT_DOUBLE_EQ(p(3.0), 6.0);
}
