/// @file test_utilities.h
/// @brief Helper functions for testing.
///
/// @author Roland Abel
/// @date October 24, 2023
///
/// Copyright (c) 2026 Roland Abel

#ifndef TEST_UTILITIES_H
#define TEST_UTILITIES_H

#include <vector>
#include <algorithm>

template<typename T>
bool unique(const std::vector<T> &vec, T epsilon) {
    std::vector<T> sorted = vec;
    std::sort(sorted.begin(), sorted.end());

    auto it = std::unique(sorted.begin(), sorted.end(), [epsilon](const T &a, const T &b) {
        return std::abs(a - b) <= epsilon;
    });

    return it == sorted.end();
}

#define EXPECT_COMPLEX_NEAR(a, b, epsilon)        \
  EXPECT_NEAR((a).real(), (b).real(), epsilon);   \
  EXPECT_NEAR((a).imag(), (b).imag(), epsilon)

#define ASSERT_COMPLEX_NEAR(a, b, epsilon)        \
  ASSERT_NEAR((a).real(), (b).real(), epsilon);   \
  ASSERT_NEAR((a).imag(), (b).imag(), epsilon)

#define EXPECT_UNIQUE(vec, epsilon)        \
    EXPECT_TRUE(unique(vec, epsilon))      \

#define ASSERT_UNIQUE(vec, epsilon)        \
  ASSERT_TRUE(unique(vec, epsilon));       \

#endif // TEST_UTILITIES_H
