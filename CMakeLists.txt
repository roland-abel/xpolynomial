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

file(GLOB XPOLYNOMIAL_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
file(GLOB XPOLYNOMIAL_TPL CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.tpp")

add_library(xpolynomial INTERFACE)
add_library(xpolynomial::xpolynomial ALIAS xpolynomial)
set_target_properties(xpolynomial PROPERTIES EXPORT_NAME xpolynomial)

target_include_directories(xpolynomial INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/xpolynomial>)

target_compile_features(xpolynomial INTERFACE cxx_std_23)

install(TARGETS xpolynomial EXPORT xpolynomialTargets)
install(FILES ${XPOLYNOMIAL_HEADERS} ${XPOLYNOMIAL_TPL} DESTINATION include/xpolynomial)
install(EXPORT xpolynomialTargets NAMESPACE xpolynomial:: DESTINATION lib/cmake/xpolynomial)

include(CMakePackageConfigHelpers)
configure_package_config_file(
        cmake/xpolynomialConfig.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/xpolynomialConfig.cmake"
        INSTALL_DESTINATION lib/cmake/xpolynomial)

write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/xpolynomialConfigVersion.cmake"
        ARCH_INDEPENDENT
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)
        
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/xpolynomialConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/xpolynomialConfigVersion.cmake"
        DESTINATION lib/cmake/xpolynomial)

add_subdirectory(examples)
add_subdirectory(test)
