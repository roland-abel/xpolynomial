/// @file legendre_polynomial.h
///
/// @author Roland Abel
/// @date 08.10.2023

#ifndef LEGENDRE_POLYNOMIAL_H_
#define LEGENDRE_POLYNOMIAL_H_

#include <vector>
#include "polynomial.h"

namespace xmath {

    /// @brief Provides methods useful for dealing with Legendre polynomials.
    template<typename T>
    class legendre_polynomial {
    public:
        using polynomial_sequence = std::vector<polynomial<T>>;

        /// @brief Creates the Legendre polynomial P_n for the given order.
        /// @param order The order of the polynomial.
        /// @return The Legendre polynomial P_n for the given order.
        static polynomial<T> create(size_t order);

    private:
        static polynomial_sequence legendre_polynomial_;
    };
}

#include "legendre_polynomial.tpp"

#endif // LEGENDRE_POLYNOMIAL_H_