enable_testing()

add_executable(${UNIT_TEST_NAME}
        legendre_polynomial_tests.cpp
        polynomial_tests.cpp
        matrix_tests.cpp
        euclidean_algorithm_tests.cpp
        root_finder_test.cpp
        real_polynomial_root_finder_test.cpp
        square_free_decomposition_test.cpp
        chebyshev_polynomial_tests.cpp
        complex_polynomial_tests.cpp
        complex_polynomial_root_finder_tests.cpp
        real_interval_tests.cpp
        polynomial_interpolation_tests.cpp
)

target_include_directories(${UNIT_TEST_NAME} PUBLIC ../include)
target_link_libraries(${UNIT_TEST_NAME}
        gtest
        gtest_main
)
