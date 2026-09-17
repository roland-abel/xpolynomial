/// @file complex_polynomial_root_finder_tests.cpp
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#include <gtest/gtest.h>
#include "test_utilities.h"
#include <complex_polynomial_root_finder.h>

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
    for (const auto &z: roots) {
        EXPECT_COMPLEX_NEAR(p(z), std::complex(0., 0.), 1e-8);
    }
}

TEST(ComplexPolynomialRootFinder, AberthEhrlichMethodTest) {
    auto p = (1.2 + 3. * i) * Z.pow(8) + 23.5;
    auto initial_points = RootFinder::nth_roots_of_unity(8);
    auto roots = RootFinder::aberth_ehrlich_method(p, initial_points);

    EXPECT_EQ(roots.size(), 8);
    EXPECT_TRUE(p.has_roots(roots));
    for (const auto &z: roots) {
        EXPECT_COMPLEX_NEAR(p(z), std::complex(0., 0.), 1e-8);
    }
}

TEST(ComplexPolynomialRootFinder, AberthEhrlichMethodWithMultipleRootStaysFiniteTest) {
    // The Aberth correction has a vanishing denominator at the double root z = 1.
    auto p = (Z - 1.).pow(2);
    auto initial_points = std::vector{std::complex(1., 0.), std::complex(2., 0.)};
    auto roots = RootFinder::aberth_ehrlich_method(p, initial_points);

    ASSERT_EQ(roots.size(), 2);
    for (const auto &z: roots) {
        EXPECT_TRUE(std::isfinite(z.real()));
        EXPECT_TRUE(std::isfinite(z.imag()));
        EXPECT_COMPLEX_NEAR(z, std::complex(1., 0.), epsilon);
    }
}
