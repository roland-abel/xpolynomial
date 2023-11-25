/// @file real_interval_tests.cpp
///
/// @author Roland Abel
/// @date 07.11.2023

#include <gtest/gtest.h>
#include "real_interval.h"

using namespace xmath;

namespace {
    using std::numbers::pi;
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

TEST(RealIntervalTests, IsLowerClosed) {
    EXPECT_TRUE(Interval().is_lower_closed());
    EXPECT_TRUE(Interval(0., 1., closed, closed).is_lower_closed());
    EXPECT_FALSE(Interval(0., 1., opened, closed).is_lower_closed());
    EXPECT_TRUE(Interval(0., 1., closed, opened).is_lower_closed());
    EXPECT_FALSE(Interval(0., 1., opened, opened).is_lower_closed());
}

TEST(RealIntervalTests, IsUpperClosed) {
    EXPECT_TRUE(Interval().is_lower_closed());
    EXPECT_TRUE(Interval(0., 1., closed, closed).is_upper_closed());
    EXPECT_TRUE(Interval(0., 1., opened, closed).is_upper_closed());
    EXPECT_FALSE(Interval(0., 1., closed, opened).is_upper_closed());
    EXPECT_FALSE(Interval(0., 1., opened, opened).is_upper_closed());
}

TEST(RealIntervalTests, IsLowerOpen) {
    EXPECT_FALSE(Interval().is_lower_open());
    EXPECT_FALSE(Interval(0., 1., closed, closed).is_lower_open());
    EXPECT_TRUE(Interval(0., 1., opened, closed).is_lower_open());
    EXPECT_FALSE(Interval(0., 1., closed, opened).is_lower_open());
    EXPECT_TRUE(Interval(0., 1., opened, opened).is_lower_open());
}

TEST(RealIntervalTests, IsUpperOpen) {
    EXPECT_FALSE(Interval().is_upper_open());
    EXPECT_FALSE(Interval(0., 1., closed, closed).is_upper_open());
    EXPECT_FALSE(Interval(0., 1., opened, closed).is_upper_open());
    EXPECT_TRUE(Interval(0., 1., closed, opened).is_upper_open());
    EXPECT_TRUE(Interval(0., 1., opened, opened).is_upper_open());
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

TEST(RealIntervalTests, IsHalfOpen) {
    EXPECT_FALSE(Interval(1., -1.).is_half_open());
    EXPECT_FALSE(Interval(2., 2., opened, opened).is_half_open());
    EXPECT_TRUE(Interval(2., 2., closed, opened).is_half_open());
    EXPECT_TRUE(Interval(2., 2., opened, closed).is_half_open());
    EXPECT_FALSE(Interval(2., 2., closed, closed).is_half_open());
}

TEST(RealIntervalTests, LinearTransform1) {
    const auto I = Interval(-1., 1.);
    const auto J = Interval(2., 5.);
    const auto [linear_mapping, m, c] = I.linear_transform(J);

    EXPECT_NEAR(linear_mapping(-1.), 2., epsilon);
    EXPECT_NEAR(linear_mapping(1.), 5., epsilon);
    EXPECT_NEAR(linear_mapping(0.), 3.5, epsilon);
}

TEST(RealIntervalTests, LinearTransform2) {
    const auto I = Interval(-1., 1.);
    const auto J = Interval(0., 2 * pi);
    const auto linear_mapping = I.linear_transform(J);

//    EXPECT_NEAR(linear_mapping(-1.), 0., epsilon);
//    EXPECT_NEAR(linear_mapping(-1. / 2.), pi / 2., epsilon);
//    EXPECT_NEAR(linear_mapping(0.), pi, epsilon);
//    EXPECT_NEAR(linear_mapping(1. / 2.), 3. / 2. * pi, epsilon);
//    EXPECT_NEAR(linear_mapping(1.), 2 * pi, epsilon);
}