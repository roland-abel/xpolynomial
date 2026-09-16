# @file CMakefile
#
# @author Roland Abel
# @date July 02, 2024
#
# Copyright (c) 2026 Roland Abel

cmake_minimum_required(VERSION 3.22.1)
project(xpolynomial
        VERSION 1.0.0
        DESCRIPTION "Math C++ template library for polynomial calculation."
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(PROJECT_NAME xpolynomial)

add_subdirectory(examples)
add_subdirectory(test)
