/// @file square_free_decomposition_tests.cpp
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
