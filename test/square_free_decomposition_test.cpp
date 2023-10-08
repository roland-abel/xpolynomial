/// @file square_free_decomposition_tests.h
///
/// @author Roland Abel
/// @date 08.10.2023

#include <gtest/gtest.h>
#include <cmath>
#include "polynomial.h"
#include "square_free_decomposition.h"

using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using SquareFree = square_free_decomposition<double>;
    using value_type = Polynomial::value_type;

    constexpr auto tolerance = Polynomial::tolerance;
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
