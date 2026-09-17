/// @file complex_polynomial.h
/// @brief Defines template classes for polynomials with complex coefficients.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

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

    /// @brief A polynomial with complex coefficients.
    /// @tparam T The data type of the real and imaginary parts.
    template<typename T>
    using complex_polynomial = polynomial<std::complex<T>, polynomial_specification<std::complex<T>>>;

    /// @brief A polynomial with real coefficients.
    /// @tparam T The data type of the coefficients.
    template<typename T>
    using real_polynomial = polynomial<T, polynomial_specification<T>>;

    /// @brief A complex number type.
    /// @tparam T The data type of the real and imaginary parts.
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

    /// @brief Multiplication operator for a complex number and a real polynomial.
    /// @param z The complex number.
    /// @param p The real polynomial.
    /// @return The resulting complex polynomial.
    template<typename T>
    complex_polynomial<T> operator*(const complex_type<T> &z, const real_polynomial<T> &p) {
        return p * z;
    }

    /// @brief Addition operator for a complex polynomial and a real polynomial.
    /// @param p The complex polynomial.
    /// @param q The real polynomial.
    /// @return The resulting complex polynomial.
    template<typename T>
    complex_polynomial<T> operator+(const complex_polynomial<T> &p, const real_polynomial<T> &q) {
        return p + std::complex<T>(1, 0) * q;
    }

    /// @brief Addition operator for a real polynomial and a complex polynomial.
    /// @param q The real polynomial.
    /// @param p The complex polynomial.
    /// @return The resulting complex polynomial.
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
