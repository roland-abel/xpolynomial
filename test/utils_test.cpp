/// @file utils_tests.cpp
/// @brief Tests for `utils` functions.
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
#include "utils.h"

using namespace xmath;

namespace {
    constexpr auto epsilon = 1e-5;
}

// Test cases for nearly_zero function
TEST(UtilsTest, NearlyZeroWithPositiveValue) {
    EXPECT_TRUE(nearly_zero(0.00000001, epsilon));
    EXPECT_FALSE(nearly_zero(0.1, epsilon));
}

TEST(UtilsTest, NearlyZeroWithNegativeValue) {
    EXPECT_TRUE(nearly_zero(-0.00000001, epsilon));
    EXPECT_FALSE(nearly_zero(-0.1, epsilon));
}

TEST(UtilsTest, NearlyZeroWithZeroValue) {
    EXPECT_TRUE(nearly_zero(0.0, epsilon));
}

// Test cases for nearly_equal function
TEST(UtilsTest, NearlyZeroWithEqualValues) {
    EXPECT_TRUE(nearly_equal(0.1, 0.1, epsilon));
    EXPECT_TRUE(nearly_equal(-0.1, -0.1, epsilon));
}

TEST(UtilsTest, NearlyZeroWithNotEqualValues) {
    EXPECT_FALSE(nearly_equal(0.1, 0.2, epsilon));
    EXPECT_FALSE(nearly_equal(-0.1, 0.1, epsilon));
}

TEST(UtilsTest, NearlyZeroWithOppositeSigns) {
    EXPECT_FALSE(nearly_equal(0.1, -0.1, epsilon));
    EXPECT_FALSE(nearly_equal(-0.1, 0.1, epsilon));
}

// Test cases for greater_than function
TEST(UtilsTest, GreaterThanWithGreaterValues) {
    EXPECT_TRUE(greater_than(1.0, 0.9, epsilon));
    EXPECT_TRUE(greater_than(-0.1, -0.2, epsilon));
}

TEST(UtilsTest, GreaterThanWithNotGreaterValues) {
    EXPECT_FALSE(greater_than(0.9, 1.0, epsilon));
    EXPECT_FALSE(greater_than(-0.2, -0.1, epsilon));
}

TEST(UtilsTest, GreaterThanWithGreaterThanEqualValues) {
    EXPECT_FALSE(greater_than(1.0, 1.0, epsilon));
    EXPECT_FALSE(greater_than(-0.1, -0.1, epsilon));
}

// Test cases for greater_than_or_equal function
TEST(UtilsTest, GreaterThanOrEqualWithGreaterValues) {
    EXPECT_TRUE(greater_than_or_equal(1.0, 0.9, epsilon));
    EXPECT_TRUE(greater_than_or_equal(-0.1, -0.2, epsilon));
}

TEST(UtilsTest, GreaterThanOrEqualWithEqualValues) {
    EXPECT_TRUE(greater_than_or_equal(1.0, 1.0, epsilon));
    EXPECT_TRUE(greater_than_or_equal(-0.1, -0.1, epsilon));
}

TEST(UtilsTest, GreaterThanOrEqualWithNotGreaterValues) {
    EXPECT_FALSE(greater_than_or_equal(0.9, 1.0, epsilon));
    EXPECT_FALSE(greater_than_or_equal(-0.2, -0.1, epsilon));
}

// Test cases for less_than function
TEST(UtilsTest, LessThanWithLesserValues) {
    EXPECT_TRUE(less_than(0.9, 1.0, epsilon));
    EXPECT_TRUE(less_than(-0.2, -0.1, epsilon));
}

TEST(UtilsTest, LessThanWithNotLesserValues) {
    EXPECT_FALSE(less_than(1.0, 0.9, epsilon));
    EXPECT_FALSE(less_than(-0.1, -0.2, epsilon));
}

TEST(UtilsTest, LessThanWithEqualValues) {
    EXPECT_FALSE(less_than(1.0, 1.0, epsilon));
    EXPECT_FALSE(less_than(-0.1, -0.1, epsilon));
}

// Test cases for is_even function
TEST(UtilsTest, IsEvenWithEvenValues) {
    EXPECT_TRUE(is_even(2));
    EXPECT_TRUE(is_even(0));
}

TEST(UtilsTest, IsEvenWithOddValues) {
    EXPECT_FALSE(is_even(1));
    EXPECT_FALSE(is_even(-3));
}

TEST(UtilsTest, IsEvenWithZeroValue) {
    EXPECT_TRUE(is_even(0));
}

// Test cases for is_odd function
TEST(UtilsTest, IsOddWithEvenValues) {
    EXPECT_FALSE(is_odd(2));
    EXPECT_FALSE(is_odd(0));
}

TEST(UtilsTest, IsOddWithOddValues) {
    EXPECT_TRUE(is_odd(1));
    EXPECT_TRUE(is_odd(-3));
}

TEST(UtilsTest, IsOddWithZeroValue) {
    EXPECT_FALSE(is_odd(0));
}
