/// @file complex_polynomial_root_finder.h
///
/// @author Roland Abel
/// @date 20.10.2023

#ifndef COMPLEX_POLYNOMIAL_H_
#define COMPLEX_POLYNOMIAL_H_

#include <complex>
#include "polynomial.h"

namespace xmath {

    namespace {
        using std::ranges::views::transform;
    }

    template<typename T>
    struct polynomial_specification<std::complex<T>> {
        using value_type = std::complex<T>;
        using size_type = size_t;
        using floating_point_type = T;
        static constexpr floating_point_type epsilon = 1e-5;
        static constexpr value_type one = std::complex<T>(1, 0);
        static constexpr value_type zero = std::complex<T>(0, 0);
    };

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

    template<typename T>
    complex_polynomial<T> operator*(const real_polynomial<T> &p, const complex_type<T> &z) {
        return complex_polynomial<T>(p.coefficients() | transform([&](const T &coeff) {
            return coeff * z;
        }));
    }

    template<typename T>
    inline complex_polynomial<T> operator*(const complex_type<T> &z, const real_polynomial<T> &p) {
        return p * z;
    }

    template<typename T>
    inline complex_polynomial<T> operator+(const complex_polynomial<T> &p, const real_polynomial<T> &q) {
        return p + std::complex<T>(1, 0) * q;
    }

    template<typename T>
    inline complex_polynomial<T> operator+(const real_polynomial<T> &q, const complex_polynomial<T> &p) {
        return p + q;
    }

    /// @brief
    /// @tparam T
    /// @param polynomial
    /// @return
    template<typename T>
    auto separate(const complex_polynomial<T> &polynomial) -> decltype(auto);
}

#include "complex_polynomial.tpp"

#endif // COMPLEX_POLYNOMIAL_H_
