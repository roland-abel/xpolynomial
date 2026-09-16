# @file CMakefile
#
# @author Roland Abel
# @date July 02, 2024
#
# Copyright (c) 2026 Roland Abel

set(UNIT_TEST_NAME unit-tests)

find_package(GTest CONFIG REQUIRED)
enable_testing()

file(GLOB TESTS "**_tests.cpp")
add_executable(${UNIT_TEST_NAME}
        ${TESTS}
        test_main.cpp)

target_link_libraries(${UNIT_TEST_NAME}
        xpolynomial::xpolynomial
        GTest::gtest
        GTest::gtest_main
        GTest::gmock
        GTest::gmock_main
)
