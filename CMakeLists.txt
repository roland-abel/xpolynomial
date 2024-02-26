cmake_minimum_required(VERSION 3.26)
project(xpolynomial
        VERSION 1.0.0
        DESCRIPTION "Math C++ template library for polynomial calculation."
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
FetchContent_Declare(
        google-test
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

FetchContent_MakeAvailable(google-test)

set(PROJECT_NAME xpolynomial)
set(UNIT_TEST_NAME unit-tests)

add_subdirectory(samples)
add_subdirectory(tests)
