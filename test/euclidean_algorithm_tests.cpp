/// @file euclidean_algorithm_tests.cpp
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
#include "euclidean_algorithm.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using EuclideanAlgorithm = euclidean_algorithm<double>;
    const auto X = Polynomial::monomial(1, 1.0);
}

TEST(EuclideanAlgorithmTests, GreatestCommonDivisorTest) {
    auto p = X.pow(4) - 2 * X.pow(3) - 6 * X.pow(2) + 12 * X + 15;
    auto q = X.pow(3) + X.pow(2) - 4 * X - 4;

    EXPECT_EQ(EuclideanAlgorithm::euclidean(p, q), X + 1);
}

TEST(EuclideanAlgorithmTests, ExtendedEuclideanTest) {
    auto p = X.pow(4) - 2 * X.pow(3) - 6 * X.pow(2) + 12 * X + 15;
    auto q = X.pow(3) + X.pow(2) - 4 * X - 4;
    auto [s, t, g] = EuclideanAlgorithm::extended_euclidean(p, q); // s, t, g such that g = gcd(p, q) = s*p + t*q

    EXPECT_EQ(g, 5 * (X + 1));
    EXPECT_EQ(s, -1 * X + 3);
    EXPECT_EQ(t, X.pow(2) - 6 * X + 10);
    EXPECT_EQ(s * p + t * q, g);
}