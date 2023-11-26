/// @file matrix_tests.h
///
/// @author Roland Abel
/// @date 25.11.2023

#include <gtest/gtest.h>
#include "polynomial_interpolation.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using values_type = Polynomial::values_type;
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
    const std::vector<double> x_values = {1.0, 2.0, 3.0};
    const std::vector<double> y_values = {2.0, 4.0, 6.0};
    const auto p = Interpolation::lagrange_interpolation(x_values, y_values);

    EXPECT_DOUBLE_EQ(p(1.0), 2.0);
    EXPECT_DOUBLE_EQ(p(2.0), 4.0);
    EXPECT_DOUBLE_EQ(p(3.0), 6.0);
}