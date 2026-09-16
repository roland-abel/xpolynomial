/// @file legendre_polynomial.h
/// @brief Creates Legendre polynomials.
///
/// @author Roland Abel
/// @date October 8, 2023
///
/// Copyright (c) 2026 Roland Abel

#ifndef LEGENDRE_POLYNOMIAL_H_
#define LEGENDRE_POLYNOMIAL_H_

#include <vector>
#include "polynomial.h"

namespace xmath {

    /// @brief Provides methods useful for dealing with Legendre polynomials.
    /// @tparam T The data type of the coefficients.
    template<typename T>
    class legendre_polynomial {
    public:
        using polynomial_sequence = std::vector<polynomial<T>>;

        /// @brief Creates the Legendre polynomial P_n for the given order.
        /// @param order The order of the Legendre polynomial.
        /// @return The Legendre polynomial P_n for the given order.
        static polynomial<T> create(size_t order);

    private:
        /// Cache of the calculated Legendre polynomials.
        static polynomial_sequence legendre_polynomial_;
    };
}

#include "legendre_polynomial.tpp"

#endif // LEGENDRE_POLYNOMIAL_H_