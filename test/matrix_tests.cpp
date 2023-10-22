/// @file matrix_tests.h
///
/// @author Roland Abel
/// @date 08.10.2023

#include <gtest/gtest.h>
#include <numeric>
#include "matrix.h"

using namespace xmath;

namespace {
    using Matrix = matrix<double>;
    using coeff_type = Matrix::coeff_type;
    using coeffs_type = Matrix::coeffs_type;
    constexpr coeff_type epsilon = Matrix::tolerance;
}

Matrix CreateTestMatrix(size_t num_rows, size_t num_cols, const coeff_type start_value = 1.0) {
    auto values = coeffs_type(num_rows * num_cols);
    std::iota(values.begin(), values.end(), start_value);

    return {num_rows, num_cols, values};
}

TEST(MatrixTests, CheckDimensions) {
    auto m = Matrix(2, 4);

    // expect number of rows/columns
    EXPECT_EQ(2, m.rows());
    EXPECT_EQ(4, m.cols());
}

TEST(MatrixTests, ConstructorWithConstValues) {
    const auto value = 2.1f;
    auto m = Matrix(2, 2, value);

    // expect initial value
    EXPECT_EQ(m(0, 0), value);
    EXPECT_EQ(m(1, 0), value);
    EXPECT_EQ(m(0, 1), value);
    EXPECT_EQ(m(1, 1), value);
}

TEST(MatrixTests, WrongInitialValueList) {
    const coeffs_type values = {1.1, 2.4, -0.7, 1.0, 2.2};
    auto action = [&values]() { return Matrix(2, 3, values); };

    // throws exception
    EXPECT_ANY_THROW(action());
}

TEST(MatrixTests, InitialWithValuesList) {
    const coeffs_type values = {1.1, 2.4, -0.7, 1.0, 2.2, -5.2};
    auto m = Matrix(2, 3, values);

    // expect number of rows/columns
    EXPECT_EQ(2, m.rows());
    EXPECT_EQ(3, m.cols());

    // expect right initial values
    EXPECT_EQ(m(0, 0), 1.1);
    EXPECT_EQ(m(0, 1), 2.4);
    EXPECT_EQ(m(0, 2), -0.7);
    EXPECT_EQ(m(1, 0), 1.0);
    EXPECT_EQ(m(1, 1), 2.2);
    EXPECT_EQ(m(1, 2), -5.2);
}

TEST(MatrixTests, InitialListConstructor) {
    auto m = Matrix
            ({
                     {1.1, 2.4, -0.7f},
                     {1.0, 2.2, -5.2, 6.4f}
             });

    // expect number of rows/columns
    EXPECT_EQ(2, m.rows());
    EXPECT_EQ(4, m.cols());
}

TEST(MatrixTests, MatrixProxyTest) {
    auto m = CreateTestMatrix(3, 2);

    //
    EXPECT_NEAR(m[0][0], 1., epsilon);
    EXPECT_NEAR(m[0][1], 2., epsilon);
    EXPECT_NEAR(m[1][1], 4., epsilon);
    EXPECT_NEAR(m[1][0], 3., epsilon);
    EXPECT_NEAR(m[2][0], 5., epsilon);
    EXPECT_NEAR(m[2][1], 6., epsilon);
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
    auto m = Matrix(2, 3);

    // expect right index
    EXPECT_EQ(0, m.index(0, 0));
    EXPECT_EQ(1, m.index(0, 1));
    EXPECT_EQ(2, m.index(0, 2));
    EXPECT_EQ(3, m.index(1, 0));
    EXPECT_EQ(4, m.index(1, 1));
    EXPECT_EQ(5, m.index(1, 2));
}

TEST(MatrixTests, TestSymmetricalMatrix) {
    auto m = Matrix
            ({
                     {1.1,  2.4,  -0.7},
                     {2.4,  2.2,  -5.2},
                     {-0.7, -5.2, -5.2}
             });
    EXPECT_TRUE(m.is_symmetrical());

    m = Matrix
            ({
                     {1.1,               1.0 / 3.0},
                     {0.333333333333333, 2.2},
             });
    EXPECT_TRUE(m.is_symmetrical());

    m = Matrix
            ({
                     {1.1, 2.4},
                     {3.4, 2.2},
             });
    EXPECT_FALSE(m.is_symmetrical());

    // is not a square matrix
    m = Matrix
            (2, 3,
             {1.1, 2.4, -0.7,
              1.0, 2.2, -5.2
             });
    EXPECT_FALSE(m.is_symmetrical());
}

TEST(MatrixTests, TestTranspose) {
    auto m = CreateTestMatrix(4, 4);
    auto transposed = Matrix
            ({
                     {1.0, 5.0, 9.0,  13.0},
                     {2.0, 6.0, 10.0, 14.0},
                     {3.0, 7.0, 11.0, 15.0},
                     {4.0, 8.0, 12.0, 16.0}
             });

    EXPECT_TRUE(m.transpose() == transposed);
    EXPECT_TRUE(m.transpose().transpose() == m);
}