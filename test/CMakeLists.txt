# @file CMakefile
#
# @author Roland Abel
# @date September 17, 2026
#
# Copyright (c) 2026 Roland Abel

set(UNIT_TEST_NAME unit-tests)

find_package(GTest CONFIG REQUIRED)

# Compile each public header on its own to ensure that every header is self-contained.
file(GLOB XPOLYNOMIAL_PUBLIC_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/../include/*.h")
foreach (HEADER ${XPOLYNOMIAL_PUBLIC_HEADERS})
    get_filename_component(HEADER_NAME ${HEADER} NAME_WE)
    set(HEADER_TU "${CMAKE_CURRENT_BINARY_DIR}/self_contained_${HEADER_NAME}.cpp")
    file(WRITE "${HEADER_TU}" "#include <${HEADER_NAME}.h>\n")
    add_library(${HEADER_NAME}_self_contained OBJECT "${HEADER_TU}")
    target_link_libraries(${HEADER_NAME}_self_contained PRIVATE xpolynomial::xpolynomial)
endforeach()

file(GLOB TESTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*_tests.cpp")
add_executable(${UNIT_TEST_NAME}
        ${TESTS}
        test_main.cpp)

target_link_libraries(${UNIT_TEST_NAME}
        xpolynomial::xpolynomial
        GTest::gtest
)

add_test(NAME ${UNIT_TEST_NAME} COMMAND ${UNIT_TEST_NAME})
