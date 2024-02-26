/// @file complex_polynomial_root_finder.h
/// @brief Defines template classes for polynomials with complex coefficients.
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

#ifndef COMPLEX_POLYNOMIAL_H_
#define COMPLEX_POLYNOMIAL_H_

#include <complex>
#include "polynomial.h"

namespace xmath {

    using std::ranges::views::transform;

    /// @brief Specialization of polynomial_specification for complex numbers.
    /// @tparam T The data type of the coefficients in the complex polynomial.
    template<typename T>
    struct polynomial_specification<std::complex<T>> {

        static_assert(std::is_floating_point_v<T>, "The type parameter must be a floating point type.");

        using value_type = std::complex<T>;
        using size_type = size_t;
        using floating_point_type = T;
        static constexpr floating_point_type epsilon = 1e-5;
        static constexpr value_type one = std::complex<T>(1, 0);
        static constexpr value_type zero = std::complex<T>(0, 0);
    };

    /// @brief Output stream operator for complex polynomials.
    /// @param os The output stream.
    /// @param p The complex polynomial to be output.
    /// @return The output stream.
    template<typename T>
    std::ostream &operator<<(std::ostream &os, const polynomial<std::complex<T>> &p) {
        for (auto coeff: p.coefficients()) {
            os << coeff;
        }
        return os;
    }

    template<typename T>
    using complex_polynomial = polynomial<std::complex<T>, polynomial_specification<std::complex<T>>>;

    template<typename T>
    using real_polynomial = polynomial<T, polynomial_specification<T>>;

    template<typename T>
    using complex_type = std::complex<T>;

    /// @brief Multiplication operator for a real polynomial and a complex number.
    /// @param p The real polynomial.
    /// @param z The complex number.
    /// @return The resulting complex polynomial.
    template<typename T>
    complex_polynomial<T> operator*(const real_polynomial<T> &p, const complex_type<T> &z) {
        return complex_polynomial<T>(p.coefficients() | transform([&](const T &coeff) {
            return coeff * z;
        }));
    }

    template<typename T>
    complex_polynomial<T> operator*(const complex_type<T> &z, const real_polynomial<T> &p) {
        return p * z;
    }

    template<typename T>
    complex_polynomial<T> operator+(const complex_polynomial<T> &p, const real_polynomial<T> &q) {
        return p + std::complex<T>(1, 0) * q;
    }

    template<typename T>
    complex_polynomial<T> operator+(const real_polynomial<T> &q, const complex_polynomial<T> &p) {
        return p + q;
    }

    /// @brief Separates a complex polynomial into its real and imaginary parts.
    /// @param p The complex polynomial to be separated.
    /// @return A pair of polynomials representing the real and imaginary parts, respectively.
    template<typename T>
    auto separate(const complex_polynomial<T> &p) -> decltype(auto);
}

#include "complex_polynomial.tpp"

#endif // COMPLEX_POLYNOMIAL_H_
