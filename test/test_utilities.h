/// @file test_utilities.h
/// @brief
///
/// @author Roland Abel
/// @date October 24, 2023
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

#ifndef XPOLYNOMIAL_TEST_UTILITIES_H
#define XPOLYNOMIAL_TEST_UTILITIES_H

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

#endif
