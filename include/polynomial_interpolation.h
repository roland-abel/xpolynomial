/// @file polynomial_interpolation.h
/// @brief Lagrange polynomial interpolation.
///
/// @author Roland Abel
/// @date September 17, 2026
///
/// Copyright (c) 2026 Roland Abel

#ifndef POLYNOMIAL_INTERPOLATION_H_
#define POLYNOMIAL_INTERPOLATION_H_

#include <vector>
#include <optional>
#include <type_traits>
#include "polynomial.h"

namespace xmath {

    /// @brief Provides methods for polynomial interpolation using Lagrange basis.
    /// @tparam T The data type of the coefficients in the polynomials.
    template<typename T>
    class polynomial_interpolation {
    public:
        using polynomial_sequence = std::vector<polynomial<T>>;
        using value_type = typename polynomial<T>::value_type;
        using values_type = typename polynomial<T>::values_type;

        static_assert(std::is_floating_point_v<T>, "The type parameter must be a floating point type.");

        /// @brief Gets the Lagrange basis for polynomials of degree k = xs.size().
        /// @param xs The distinct x points x_0, ..., x_k.
        /// @return The Lagrange basis for polynomials. For each polynomial l_j is l_j(x_m) = 0 if j != m and l_j(x_j) = 1.
        static polynomial_sequence lagrange_basis(const std::vector<T> &xs);

        /// @brief Performs Lagrange interpolation to compute the polynomial that passes through a set of given points.
        /// @param x_values The vector containing the x coordinates of the interpolation points.
        /// @param y_values The vector containing the corresponding y coordinates of the interpolation points.
        /// @return A polynomial representing the Lagrange interpolation polynomial.
        static std::optional<polynomial<T>> lagrange_interpolation(const values_type &x_values, const values_type &y_values);
    };
}

#include "polynomial_interpolation.tpp"

#endif // POLYNOMIAL_INTERPOLATION_H_