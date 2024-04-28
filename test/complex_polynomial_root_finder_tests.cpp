/// @file complex_polynomial_root_finder_tests.cpp
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

#include <gtest/gtest.h>
#include "test_utilities.h"
#include "complex_polynomial_root_finder.h"

using namespace xmath;

namespace {
    using ComplexPolynomial = complex_polynomial<double>;
    using RealPolynomial = polynomial<double, polynomial_specification<double>>;
    using RootFinder = complex_polynomial_root_finder<double>;

    constexpr auto i = std::complex(0., 1.);
    constexpr auto I = ComplexPolynomial::value_type(i);

    constexpr double epsilon = ComplexPolynomial::epsilon;
    const auto zero = ComplexPolynomial::zero();
    const auto one = ComplexPolynomial::one();

    const auto X = RealPolynomial::monomial(1, 1.0);
    const auto Y = RealPolynomial::monomial(1, 1.0);
    const auto Z = ComplexPolynomial::monomial(1, 1.0);
}

TEST(ComplexPolynomialRootFinder, QuadraticRootsOfUnityTest) {
    auto roots = RootFinder::nth_roots_of_unity(2);

    EXPECT_EQ(roots.size(), 2);
    EXPECT_COMPLEX_NEAR(roots[0], std::complex(1., 0.), epsilon);
    EXPECT_COMPLEX_NEAR(roots[1], std::complex(-1., 0.), epsilon);
}

TEST(ComplexPolynomialRootFinder, CubicRootsOfUnityTest) {
    auto roots = RootFinder::nth_roots_of_unity(3);

    EXPECT_EQ(roots.size(), 3);
    EXPECT_COMPLEX_NEAR(roots[0], std::complex(1., 0.), epsilon);
    EXPECT_COMPLEX_NEAR(roots[1], std::complex(-1. / 2., std::sqrt(3.) / 2.), epsilon);
    EXPECT_COMPLEX_NEAR(roots[2], std::complex(-1. / 2., -std::sqrt(3.) / 2.), epsilon);
}

TEST(ComplexPolynomialRootFinder, HasRootsOfUnityTest) {
    auto p = Z.pow(7) - 1.;
    auto roots_of_unity = RootFinder::nth_roots_of_unity(7);

    EXPECT_TRUE(p.has_roots(roots_of_unity));
}

TEST(ComplexPolynomialRootFinder, DurandKernerMethodTest) {
    auto p = (2.6 + i) * Z.pow(7) - 10.5;
    auto initial_points = RootFinder::nth_roots_of_unity(7);
    auto roots = RootFinder::durand_kerner_method(p, initial_points);

    EXPECT_EQ(roots.size(), 7);
    EXPECT_TRUE(p.has_roots(roots));
}

TEST(ComplexPolynomialRootFinder, AberthEhrlichMethodTest) {
    auto p = (1.2 + 3. * i) * Z.pow(8) + 23.5;
    auto initial_points = RootFinder::nth_roots_of_unity(8);
    auto roots = RootFinder::aberth_ehrlich_method(p, initial_points);

    EXPECT_EQ(roots.size(), 8);
    EXPECT_TRUE(p.has_roots(roots));
}
