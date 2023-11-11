/// @file chebyshev_polynomial.h
///
/// @author Roland Abel
/// @date 08.10.2023

#ifndef CHEBYSHEV_POLYNOMIAL_H_
#define CHEBYSHEV_POLYNOMIAL_H_

#include <vector>
#include "real_interval.h"
#include "polynomial.h"

namespace xmath {

    /// @brief Provides methods useful for dealing with Chebyshev polynomials of the first and second kind.
    /// @tparam T The data type of the coefficients.
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

        /// @brief Computes the Chebyshev nodes for a given order within a specified interval.
        /// @param order  The order of the Chebyshev nodes.
        /// @return A vector containing the Chebyshev nodes for the specified order within the interval.
        static values_type chebyshev_nodes(size_t order, const real_interval<T>& I);

        /// Calculates the Chebyshev series of order n for points in `xs` and the given a list of coefficients `alphas`.
        /// @param alphas A vector containing the coefficients for the Chebyshev series polynomial.
        /// @return The Chebyshev series of order n for the given coefficients.
        static polynomial<T> chebyshev_series(const values_type &alphas);

        /// Evaluates the Chebyshev series for the given coefficients and input value by using
        /// Clenshaw algorithm.
        /// @param alphas A vector containing the coefficients for the Chebyshev series polynomial.
        /// @param xs The input x points.
        /// @return The Chebyshev series of order n for the given coefficients for the `xs`.
        static value_type clenshaw(const values_type &alphas, const value_type &x);

    private:
        static polynomial_sequence chebyshev_1st_kind_polynomials_;
    };
}

#include "chebyshev_polynomial.tpp"

#endif // CHEBYSHEV_POLYNOMIAL_H_