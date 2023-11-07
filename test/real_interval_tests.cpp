/// @file real_interval_tests.cpp
///
/// @author Roland Abel
/// @date 07.11.2023

#include <gtest/gtest.h>
#include "real_interval.h"

using namespace xmath;

namespace {
    using Interval = real_interval<double>;
    constexpr double epsilon = 1e-9;
}

TEST(RealIntervalTests, Constructor) {
    auto I = Interval(-1., 1.);
    EXPECT_NEAR(I.start(), -1.0, epsilon);
    EXPECT_NEAR(I.end(), 1.0, epsilon);
}