# @file CMakefile
#
# @author Roland Abel
# @date July 02, 2024
#
# Copyright (c) 2026 Roland Abel

set(UNIT_TEST_NAME unit-tests)

find_package(GTest CONFIG REQUIRED)

file(GLOB TESTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*_tests.cpp")
add_executable(${UNIT_TEST_NAME}
        ${TESTS}
        test_main.cpp)

target_link_libraries(${UNIT_TEST_NAME}
        xpolynomial::xpolynomial
        GTest::gtest
)

add_test(NAME ${UNIT_TEST_NAME} COMMAND ${UNIT_TEST_NAME})
