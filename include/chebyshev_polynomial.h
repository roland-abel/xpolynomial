/// @file chebyshev_polynomial.h
///
/// @author Roland Abel
/// @date 08.10.2023

#ifndef CHEBYSHEV_POLYNOMIAL_H_
#define CHEBYSHEV_POLYNOMIAL_H_

#include <vector>
#include "polynomial.h"

namespace xmath {

    /// @brief Provides methods useful for dealing with Chebyshev polynomials of the first and second kind.
    template<typename T>
    class chebyshev_polynomial {
    public:
        using polynomial_sequence = std::vector<polynomial<T>>;
        using value_type = polynomial<T>::value_type;
        using values_type = polynomial<T>::values_type;

        /// @brief Creates first kind Chebyshev polynomial T_n for the given order.
        /// @param order The order of the polynomial.
        /// @return The Chebyshev polynomial T_n of the first kind.
        static polynomial<T> create_1st_kind(size_t order);

        /// @brief
        /// @param order
        /// @return
        static values_type roots_1st_kind(size_t order);

    private:
        static polynomial_sequence chebyshev_1st_kind_polynomials_;
    };
}

#include "chebyshev_polynomial.tpp"

#endif // CHEBYSHEV_POLYNOMIAL_H_