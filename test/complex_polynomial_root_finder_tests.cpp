/// @file complex_polynomial_root_finder_tests.h
///
/// @author Roland Abel
/// @date 24.10.2023

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
    auto zero = ComplexPolynomial::zero();
    auto one = ComplexPolynomial::one();

    auto X = RealPolynomial::monomial(1, 1.0);
    auto Y = RealPolynomial::monomial(1, 1.0);
    auto Z = ComplexPolynomial::monomial(1, 1.0);
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

