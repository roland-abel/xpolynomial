/// @file matrix_tests.h
///
/// @author Roland Abel
/// @date 08.10.2023

#include <gtest/gtest.h>
#include "matrix.h"

using namespace xmath;

namespace {
    using Matrix = matrix<double>;
    using value_type = Matrix::value_type;
    using values_type = Matrix::values_type;

    constexpr value_type epsilon = Matrix::epsilon;
}

TEST(MatrixTests, DefaultConstructor) {
    auto M = Matrix(2, 2);
    EXPECT_TRUE(M.is_zero());
}

TEST(MatrixTests, CheckDimensions) {
    auto M = Matrix(2, 4);

    // expect number of rows/columns
    EXPECT_EQ(2, M.rows());
    EXPECT_EQ(4, M.cols());
}

TEST(MatrixTests, ConstructorWithConstValues) {
    const auto value = 2.1f;
    auto M = Matrix(2, 2, value);

    // expect initial value
    EXPECT_EQ(M(0, 0), value);
    EXPECT_EQ(M(1, 0), value);
    EXPECT_EQ(M(0, 1), value);
    EXPECT_EQ(M(1, 1), value);
}

TEST(MatrixTests, WrongInitialValueList) {
    const auto coeffs = {1.1, 2.4, -0.7, 1.0, 2.2};
    auto action = [&coeffs]() { return Matrix(2, 3, coeffs); };

    // throws exception
    EXPECT_ANY_THROW(action());
}

TEST(MatrixTests, InitialWithValuesList) {
    const auto coeffs = {1.1, 2.4, -0.7, 1.0, 2.2, -5.2};
    auto M = Matrix(2, 3, coeffs);

    // expect number of rows/columns
    EXPECT_EQ(2, M.rows());
    EXPECT_EQ(3, M.cols());

    // expect right initial values
    EXPECT_EQ(M(0, 0), 1.1);
    EXPECT_EQ(M(0, 1), 2.4);
    EXPECT_EQ(M(0, 2), -0.7);
    EXPECT_EQ(M(1, 0), 1.0);
    EXPECT_EQ(M(1, 1), 2.2);
    EXPECT_EQ(M(1, 2), -5.2);
}

TEST(MatrixTests, InitialListConstructor) {
    auto M = Matrix({{1.1, 2.4, -0.7f},
                     {1.0, 2.2, -5.2, 6.4f}});

    // expect number of rows/columns
    EXPECT_EQ(2, M.rows());
    EXPECT_EQ(4, M.cols());
}

TEST(MatrixTests, MatrixProxyTest) {
    auto M = Matrix(3, 2, {1., 2., 3., 4., 5., 6.});

    EXPECT_NEAR(M[0][0], 1., epsilon);
    EXPECT_NEAR(M[0][1], 2., epsilon);
    EXPECT_NEAR(M[1][1], 4., epsilon);
    EXPECT_NEAR(M[1][0], 3., epsilon);
    EXPECT_NEAR(M[2][0], 5., epsilon);
    EXPECT_NEAR(M[2][1], 6., epsilon);
}

TEST(MatrixTests, ZeroMatrixTest) {
    EXPECT_TRUE(Matrix(2, 2, {.0, .0, .0, .0}).is_zero());
    EXPECT_FALSE(Matrix(2, 2, {.0, .1, .0, .0}).is_zero());
}

TEST(MatrixTests, CheckEmpty) {
    // expect an empty Matrix
    EXPECT_TRUE(Matrix().is_empty());

    // expect a non empty Matrix
    EXPECT_FALSE(Matrix(2, 2).is_empty());
}

TEST(MatrixTests, CheckSquareMatrix) {
    // expect square Matrix
    EXPECT_TRUE(Matrix(4, 4).is_square());

    // expect non square Matrix
    EXPECT_FALSE(Matrix(4, 5).is_square());
}

TEST(MatrixTests, CheckIndex) {
    auto M = Matrix(2, 3);

    // expect right index
    EXPECT_EQ(0, M.index(0, 0));
    EXPECT_EQ(1, M.index(0, 1));
    EXPECT_EQ(2, M.index(0, 2));
    EXPECT_EQ(3, M.index(1, 0));
    EXPECT_EQ(4, M.index(1, 1));
    EXPECT_EQ(5, M.index(1, 2));
}

TEST(MatrixTests, TestSymmetricalMatrix) {
    auto M = Matrix({{1.1,  2.4,  -0.7},
                     {2.4,  2.2,  -5.2},
                     {-0.7, -5.2, -5.2}});
    EXPECT_TRUE(M.is_symmetrical());

    M = Matrix({{1.1,       1.0 / 3.0},
                {0.3333333, 2.2},});
    EXPECT_TRUE(M.is_symmetrical());

    M = Matrix({{1.1, 2.4},
                {3.4, 2.2},});
    EXPECT_FALSE(M.is_symmetrical());

    // is not a square matrix
    M = Matrix(2, 3, {1.1, 2.4, -0.7, 1.0, 2.2, -5.2});
    EXPECT_FALSE(M.is_symmetrical());
}

TEST(MatrixTests, TestTranspose) {
    auto A = Matrix(2, 4, {1., 2., 3., 4., 5., 6., 7., 8.});
    auto A_transposed = Matrix(4, 2, {1., 5., 2., 6., 3., 7., 4., 8.});

    EXPECT_TRUE(A.transpose() == A_transposed);
    EXPECT_TRUE(A.transpose().transpose() == A);
}

TEST(MatrixTests, MatrixAdditionTest) {
    const auto A = Matrix(2, 4, {1., 2., 3., 4., 5., 6., 7., 2.5});
    const auto B = Matrix(2, 4, {8., 7., 6., 5., -4., 3., 2., 1.});

    EXPECT_EQ(A + B, Matrix(2, 4, {9., 9., 9., 9., 1., 9., 9., 3.5}));
}

TEST(MatrixTests, MatrixSubtractionTest) {
    const auto A = Matrix(2, 4, {1., 2., 3., 4., 5., 6., 7., 2.5});
    const auto B = Matrix(2, 4, {8., 7., 6., 5., -4., 3., 2., 1.});

    EXPECT_EQ(A - B, Matrix(2, 4, {-7., -5., -3., -1., 9., 3., 5., 1.5}));
    EXPECT_EQ(B - A, Matrix(2, 4, {7., 5., 3., 1., -9., -3., -5., -1.5}));
}
