/// @file square_free_decomposition_tests.cpp
/// @brief Unit test for `square_free_decomposition_tests`.
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
#include "polynomial.h"
#include "square_free_decomposition.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using SquareFree = square_free_decomposition<double>;
    using value_type = Polynomial::value_type;

    constexpr auto epsilon = Polynomial::epsilon;
    auto zero = Polynomial::zero();
    auto one = Polynomial::one();
    auto X = Polynomial::monomial(1);
}

TEST(SquareFreeDecompositionTests, IsSquareFreeTest) {
    EXPECT_TRUE(SquareFree::is_square_free(X - 1));
    EXPECT_TRUE(SquareFree::is_square_free((X - 1) * (X - 2)));
    EXPECT_TRUE(SquareFree::is_square_free(X.pow(4) + 1));
    EXPECT_TRUE(SquareFree::is_square_free((X.pow(4) + 1) * (X.pow(2) + 1)));

    EXPECT_FALSE(SquareFree::is_square_free((X - 2).pow(2)));
    EXPECT_FALSE(SquareFree::is_square_free((X - 1) * (X - 2).pow(2)));
}

TEST(SquareFreeDecompositionTests, ContentTest) {
    ASSERT_EQ(SquareFree::content({2.0, 3.0, 5.0, 7.0}), 1);
    ASSERT_EQ(SquareFree::content({1.0, 2.0, 3.0}), 1);
    ASSERT_EQ(SquareFree::content({2.0, 4.0, 6.0}), 2);
    ASSERT_EQ(SquareFree::content({8.0, 12.0, 16.0}), 4);
    ASSERT_EQ(SquareFree::content({8.0, -12.0, 16.0}), 4);
    ASSERT_EQ(SquareFree::content({15.0, -5.0, -30.0}), 5);
}

TEST(SquareFreeDecompositionTests, PrimitivePartTest) {
    ASSERT_EQ(SquareFree::primitive_part({2., 3., 7., 11.}), Polynomial({2., 3., 7., 11.}));
    ASSERT_EQ(SquareFree::primitive_part({2., 4., 6.}), Polynomial({1., 2., 3.}));
    ASSERT_EQ(SquareFree::primitive_part({15., 5., 30.}), Polynomial({3., 1., 6.}));
    ASSERT_EQ(SquareFree::primitive_part({8.0, -12.0, 16.0}), Polynomial({2.0, -3.0, 4.0}));
}

TEST(SquareFreeDecompositionTests, YunAlgorithmNotIntegralPolynomialTest) {
    auto p = Polynomial ::from_roots({{1., 2., 1. / 4.}});
    ASSERT_TRUE(SquareFree::is_square_free(p) && !p.is_integer());
    EXPECT_EQ(SquareFree::yun_algorithm(p).size(), 1);

    auto q = Polynomial ::from_roots({{1., 2., 1., 1./2., 3.}});
    ASSERT_TRUE(!SquareFree::is_square_free(q) && !q.is_integer());
    EXPECT_EQ(SquareFree::yun_algorithm(q).size(), 0);
}

TEST(SquareFreeDecompositionTests, YunAlgorithm1Test) {
    auto p = (X - 3.).pow(3) * (X - 2.).pow(2) * (X - 1.);
    auto seq = SquareFree::yun_algorithm(p);

    EXPECT_EQ(seq.size(), 3);
    EXPECT_EQ(p, SquareFree::from_square_free_decomposition(seq));
}

TEST(SquareFreeDecompositionTests, YunAlgorithm2Test) {
    auto p = X.pow(2) * (X.pow(2) + 2).pow(3);
    auto seq = SquareFree::yun_algorithm(p);

    EXPECT_EQ(seq.size(), 3);
    EXPECT_EQ(p, SquareFree::from_square_free_decomposition(seq));
}

TEST(SquareFreeDecompositionTests, YunAlgorithm3Test) {
    auto p = (X - 4.) * (X + 3.).pow(2) * (X.pow(2) + X - 3.).pow(2) * (X.pow(2) + X).pow(4);
    auto seq = SquareFree::yun_algorithm(p);

    EXPECT_EQ(seq.size(), 4);
    EXPECT_EQ(p, SquareFree::from_square_free_decomposition(seq));
}

TEST(SquareFreeDecompositionTests, YunAlgorithm4Test) {
    auto p = (X + 2).pow(2) * (X + 1).pow(3) * X.pow(4) * (X - 1).pow(5) * (X - 2).pow(6);
    auto seq = SquareFree::yun_algorithm(p);

    EXPECT_EQ(seq.size(), 6);

    EXPECT_EQ(seq[0], one);
    EXPECT_TRUE(seq[0].is_integer());

    EXPECT_EQ(seq[1], X + 2);
    EXPECT_TRUE(seq[1].is_integer());

    EXPECT_EQ(seq[2], X + 1);
    EXPECT_TRUE(seq[2].is_integer());

    EXPECT_EQ(seq[3], X);
    EXPECT_TRUE(seq[3].is_integer());

    EXPECT_EQ(seq[4], X - 1);
    EXPECT_TRUE(seq[4].is_integer());

    EXPECT_EQ(seq[5], X - 2);
    EXPECT_TRUE(seq[5].is_integer());

    EXPECT_EQ(p, SquareFree::from_square_free_decomposition(seq));
}