/// @file complex_polynomial_root_finder.h
/// @brief
///
/// @author Roland Abel
/// @date October 20, 2023
///
/// Copyright (c) 2023 Roland Abel
///
/// This software is released under the MIT License.
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
/// THE SOFTWARE.

#include <vector>
#include <complex>
#include "complex_polynomial.h"

#ifndef COMPLEX_POLYNOMIAL_ROOT_FINDER_H_
#define COMPLEX_POLYNOMIAL_ROOT_FINDER_H_

namespace xmath {

    /// @brief
    /// @tparam T
    template<typename T>
    class complex_polynomial_root_finder {
    public:
        using complex_numbers = std::vector<std::complex<T>>;

        /// @brief Gets a vector of th n-th roots of unity.
        /// @param n Positive integer.
        /// @return The n-th roots of unity.
        static std::vector<std::complex<T>> nth_roots_of_unity(int n);

        /// @brief
        /// @param p
        /// @param initial_points
        /// @return
        static std::vector<std::complex<T>> durand_kerner_method(
                const complex_polynomial<T> &p,
                const std::vector<std::complex<T>> &initial_points,
                size_t max_iterations = 100);

        static std::vector<std::complex<T>> aberth_ehrlich_method(
                const complex_polynomial<T> &p,
                const std::vector<std::complex<T>> &initial_points,
                size_t max_iterations = 100);
    };
}

#include "complex_polynomial_root_finder.tpp"

#endif // COMPLEX_POLYNOMIAL_ROOT_FINDER_H_
