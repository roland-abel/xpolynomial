/// @file euclidean_algorithm_tests.h
///
/// @author Roland Abel
/// @date 08.10.2023

#include <gtest/gtest.h>
#include <cmath>
#include "euclidean_algorithm.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using EuclideanAlgorithm = euclidean_algorithm<double>;

    constexpr auto tolerance = Polynomial::tolerance;
    auto zero = Polynomial::zero();
    auto one = Polynomial::one();
    auto X = Polynomial::monomial(1, 1.0);
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
    EXPECT_EQ(s, (-1) * X + 3);
    EXPECT_EQ(t, X.pow(2) - 6 * X + 10);
    EXPECT_EQ(s * p + t * q, g);
}