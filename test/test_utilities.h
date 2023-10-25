/// @file test_utilties.h
///
/// @author Roland Abel
/// @date 24.10.2023

#ifndef XPOLYNOMIAL_TEST_UTILITIES_H
#define XPOLYNOMIAL_TEST_UTILITIES_H

#define EXPECT_COMPLEX_NEAR(a, b, abs_error)             \
  EXPECT_NEAR(a.real(), b.real(), abs_error);            \
  EXPECT_NEAR(a.imag(), b.imag(), abs_error)

#endif
