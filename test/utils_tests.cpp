/// @file utils_tests.cpp
/// @brief Tests for `utils` functions.
///
/// @author Roland Abel
/// @date October 8, 2023
///
/// Copyright (c) 2026 Roland Abel

#include <gtest/gtest.h>
#include <xpolynomial.h>

using namespace xmath;

namespace {
    constexpr auto epsilon = 1e-5;
}

namespace {
    constexpr auto compile_time_check = [] {
        return nearly_zero(0.00000001, epsilon) && !nearly_zero(0.1, epsilon)
            && nearly_equal(0.1, 0.1, epsilon) && !nearly_equal(0.1, 0.2, epsilon)
            && greater_than(1.0, 0.9, epsilon) && !greater_than(0.9, 1.0, epsilon)
            && greater_than_or_equal(1.0, 0.9, epsilon) && greater_than_or_equal(1.0, 1.0, epsilon)
            && less_than(0.9, 1.0, epsilon) && !less_than(1.0, 0.9, epsilon)
            && less_than_or_equal(0.9, 1.0, epsilon) && less_than_or_equal(1.0, 1.0, epsilon)
            && is_even(2) && is_even(0) && !is_even(1)
            && is_odd(1) && !is_odd(2) && is_odd(-3);
    }();
    static_assert(compile_time_check, "utils helpers are not constexpr");
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
