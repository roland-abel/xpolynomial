/// @file interval_tests.cpp
/// @brief
///
/// @author Roland Abel
/// @date February 27, 2024
///
/// Copyright (c) 2026 Roland Abel

#include <numbers>
#include <gtest/gtest.h>
#include <xpolynomial.h>

using namespace xmath;

namespace {
    using std::numbers::pi;
    using Interval = interval<double>;

    const auto opened = interval<double>::interval_bounds::opened;
    const auto closed = interval<double>::interval_bounds::closed;
    constexpr double epsilon = 1e-9;
}

namespace {
    constexpr auto compile_time_check = [] {
        const auto closed_interval = Interval(2., 2., closed, closed);
        const auto opened_interval = Interval(2., 2., opened, opened);
        const auto half_open_interval = Interval(0., 1., opened, closed);

        return closed_interval.is_closed() && !closed_interval.is_empty()
            && opened_interval.is_opened() && opened_interval.is_empty() && opened_interval.is_degenerate()
            && half_open_interval.is_lower_open() && half_open_interval.is_upper_closed()
            && half_open_interval.is_half_open() && !half_open_interval.is_empty()
            && Interval(-1., 1.).length() == 2.0 && Interval().lower() == 0.0 && Interval().upper() == 1.0;
    }();
    static_assert(compile_time_check, "interval methods are not constexpr");
}

TEST(IntervalTests, DefaultCtorTest) {
    const auto I = Interval{};
    EXPECT_NEAR(I.lower(), 0.0, epsilon);
    EXPECT_NEAR(I.upper(), 1.0, epsilon);
    EXPECT_FALSE(I.is_empty());
    EXPECT_TRUE(I.is_lower_open());
    EXPECT_TRUE(I.is_upper_closed());
    EXPECT_TRUE(I.is_half_open());
}

TEST(IntervalTests, CtorTest) {
    const auto I = Interval{-1., 1.};
    EXPECT_NEAR(I.lower(), -1.0, epsilon);
    EXPECT_NEAR(I.upper(), 1.0, epsilon);

    EXPECT_TRUE(Interval().is_lower_open());
    EXPECT_TRUE(Interval().is_upper_closed());
}

TEST(IntervalTests, IsOpenedTest) {
    EXPECT_FALSE(Interval(0., 1., closed, closed).is_opened());
    EXPECT_FALSE(Interval(0., 1., opened, closed).is_opened());
    EXPECT_FALSE(Interval(0., 1., closed, opened).is_opened());
    EXPECT_TRUE(Interval(0., 1., opened, opened).is_opened());
}

TEST(IntervalTests, IsLowerClosedTest) {
    EXPECT_TRUE(Interval(0., 1., closed, closed).is_lower_closed());
    EXPECT_FALSE(Interval(0., 1., opened, closed).is_lower_closed());
    EXPECT_TRUE(Interval(0., 1., closed, opened).is_lower_closed());
    EXPECT_FALSE(Interval(0., 1., opened, opened).is_lower_closed());
}

TEST(IntervalTests, IsUpperClosedTest) {
    EXPECT_TRUE(Interval(0., 1., closed, closed).is_upper_closed());
    EXPECT_TRUE(Interval(0., 1., opened, closed).is_upper_closed());
    EXPECT_FALSE(Interval(0., 1., closed, opened).is_upper_closed());
    EXPECT_FALSE(Interval(0., 1., opened, opened).is_upper_closed());
}

TEST(IntervalTests, IsLowerOpenTest) {
    EXPECT_FALSE(Interval(0., 1., closed, closed).is_lower_open());
    EXPECT_TRUE(Interval(0., 1., opened, closed).is_lower_open());
    EXPECT_FALSE(Interval(0., 1., closed, opened).is_lower_open());
    EXPECT_TRUE(Interval(0., 1., opened, opened).is_lower_open());
}

TEST(IntervalTests, IsUpperOpenTest) {
    EXPECT_FALSE(Interval(0., 1., closed, closed).is_upper_open());
    EXPECT_FALSE(Interval(0., 1., opened, closed).is_upper_open());
    EXPECT_TRUE(Interval(0., 1., closed, opened).is_upper_open());
    EXPECT_TRUE(Interval(0., 1., opened, opened).is_upper_open());
}

TEST(IntervalTests, IsClosedTest) {
    EXPECT_FALSE(Interval(0., 1., opened, opened).is_closed());
    EXPECT_FALSE(Interval(0., 1., opened, closed).is_closed());
    EXPECT_FALSE(Interval(0., 1., closed, opened).is_closed());
    EXPECT_TRUE(Interval(0., 1., closed, closed).is_closed());
}

TEST(IntervalTests, IsDegenerateTest) {
    EXPECT_FALSE(Interval(1., -1.).is_degenerate());
    EXPECT_TRUE(Interval(1., 1., opened, opened).is_degenerate());
    EXPECT_TRUE(Interval(2., 2., closed, closed).is_degenerate());
}

TEST(IntervalTests, IsEmptyTest) {
    EXPECT_TRUE(Interval(1., -1.).is_empty());
    EXPECT_TRUE(Interval(2., 2., opened, opened).is_empty());
    EXPECT_TRUE(Interval(2., 2., closed, opened).is_empty());
    EXPECT_TRUE(Interval(2., 2., opened, closed).is_empty());
    EXPECT_FALSE(Interval(2., 2., closed, closed).is_empty());
}

TEST(IntervalTests, IsHalfOpenTest) {
    EXPECT_TRUE(Interval(1., -1.).is_half_open());
    EXPECT_FALSE(Interval(2., 2., opened, opened).is_half_open());
    EXPECT_TRUE(Interval(2., 2., closed, opened).is_half_open());
    EXPECT_TRUE(Interval(2., 2., opened, closed).is_half_open());
    EXPECT_FALSE(Interval(2., 2., closed, closed).is_half_open());
}

TEST(IntervalTests, LinearTransform1Test) {
    const auto I = Interval{-1., 1.};
    const auto J = Interval{2., 5.};
    const auto map = I.linear_transform(J);

    EXPECT_NEAR(map(-1.), 2., epsilon);
    EXPECT_NEAR(map(1.), 5., epsilon);
    EXPECT_NEAR(map(0.), 3.5, epsilon);
}

TEST(IntervalTests, LinearTransform2Test) {
    const auto I = Interval{-1., 1.};
    const auto J = Interval{0., 2 * pi};
    const auto map = I.linear_transform(J);

    EXPECT_NEAR(map(-1.), 0., epsilon);
    EXPECT_NEAR(map(-1. / 2.), pi / 2., epsilon);
    EXPECT_NEAR(map(0.), pi, epsilon);
    EXPECT_NEAR(map(1. / 2.), 3. / 2. * pi, epsilon);
    EXPECT_NEAR(map(1.), 2 * pi, epsilon);
}

TEST(IntervalTests, BisectTest) {
    const auto I = Interval{-2., 1.};
    const auto [I1, I2] = I.bisect();

    EXPECT_NEAR(I1.lower(), -2., epsilon);
    EXPECT_NEAR(I1.upper(), -0.5, epsilon);
    EXPECT_NEAR(I2.lower(), -0.5, epsilon);
    EXPECT_NEAR(I2.upper(), 1., epsilon);
}

TEST(IntervalTests, BisectHalfOpenIntervalsTest) {
    const auto I = Interval{-2., 1.};
    const auto [I1, I2] = I.bisect(opened, closed);

    EXPECT_TRUE(I1.is_lower_open());
    EXPECT_TRUE(I1.is_upper_closed());
    EXPECT_TRUE(I2.is_lower_open());
    EXPECT_TRUE(I2.is_upper_closed());
}
