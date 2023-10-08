
#include_directories(
#        ../lib/google-test/googletest/include
#        ../lib/google-test/googletest/src
#)

enable_testing()

add_executable(${UNIT_TEST_NAME}
        legendre_polynomial_tests.cpp
        polynomial_tests.cpp
        matrix_tests.cpp
        euclidean_algorithm_tests.cpp
        root_finder_test.cpp
        real_polynomial_root_finder_test.cpp
        square_free_decomposition_test.cpp
)

target_link_libraries(${UNIT_TEST_NAME}
        ${LIBRARY_NAME}
        gtest
        gtest_main
)
