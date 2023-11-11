/// @file real_interval_tests.cpp
///
/// @author Roland Abel
/// @date 07.11.2023

#include <gtest/gtest.h>
#include "real_interval.h"

using namespace xmath;

namespace {
    using Interval = real_interval<double>;
    auto opened = real_interval<double>::interval_bounds::opened;
    auto closed = real_interval<double>::interval_bounds::closed;
    constexpr double epsilon = 1e-9;
}

TEST(RealIntervalTests, DefaultConstructorTest) {
    auto I = Interval();
    EXPECT_NEAR(I.lower(), 0.0, epsilon);
    EXPECT_NEAR(I.upper(), 1.0, epsilon);
    EXPECT_FALSE(I.is_empty());
    EXPECT_TRUE(I.is_closed());
}

TEST(RealIntervalTests, ConstructorTest) {
    auto I = Interval(-1., 1.);
    EXPECT_NEAR(I.lower(), -1.0, epsilon);
    EXPECT_NEAR(I.upper(), 1.0, epsilon);
}

TEST(RealIntervalTests, IsOpenedTest) {
    EXPECT_FALSE(Interval().is_opened());
    EXPECT_FALSE(Interval(0., 1., closed, closed).is_opened());
    EXPECT_FALSE(Interval(0., 1., opened, closed).is_opened());
    EXPECT_FALSE(Interval(0., 1., closed, opened).is_opened());
    EXPECT_TRUE(Interval(0., 1., opened, opened).is_opened());
}

TEST(RealIntervalTests, IsClosedTest) {
    EXPECT_TRUE(Interval().is_closed());
    EXPECT_FALSE(Interval(0., 1., opened, opened).is_closed());
    EXPECT_FALSE(Interval(0., 1., opened, closed).is_closed());
    EXPECT_FALSE(Interval(0., 1., closed, opened).is_closed());
    EXPECT_TRUE(Interval(0., 1., closed, closed).is_closed());
}

TEST(RealIntervalTests, IsDegenerateTest) {
    EXPECT_FALSE(Interval(1., -1.).is_degenerate());
    EXPECT_TRUE(Interval(1., 1., opened, opened).is_degenerate());
    EXPECT_TRUE(Interval(2., 2., closed, closed).is_degenerate());
}

TEST(RealIntervalTests, IsEmptyTest) {
    EXPECT_TRUE(Interval(1., -1.).is_empty());
    EXPECT_TRUE(Interval(2., 2., opened, opened).is_empty());
    EXPECT_TRUE(Interval(2., 2., closed, opened).is_empty());
    EXPECT_TRUE(Interval(2., 2., opened, closed).is_empty());
    EXPECT_FALSE(Interval(2., 2., closed, closed).is_empty());
}