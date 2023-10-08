//
// Created by abel on 06.10.2023.
//

#ifndef SQUARE_FREE_DECOMPOSITION_H_
#define SQUARE_FREE_DECOMPOSITION_H_

#include <vector>
#include "polynomial.h"

namespace xmath {

    /// @brief A square-free decomposition or square-free factorization of a polynomial is a factorization
    /// into powers of square-free polynomials p = q_1 q_2^2 q_3^3 ... qn^n where those of the q_i are
    /// non-constant and pairwise coprime square-free polynomials.
    /// See: https://en.wikipedia.org/wiki/Square-free_polynomial
    /// @tparam T The coefficient's type.
    template<typename T>
    class square_free_decomposition {
    public:
        using value_type = polynomial<T>::value_type;
        using polynomial_sequence = std::vector<polynomial<T>>;

        /// @brief Checks if the given polynomial is square-free. A polynomial p is square-free if and only
        /// if the greatest common division of the polynomial p and its derivative p' is constant.
        /// @return True if the polynomial p is square-free; otherwise false.
        static bool is_square_free(const polynomial<T> &p);

        /// @brief Gets the square-free decomposition of the polynomial p.
        /// @param p The polynomial.
        /// @return The square-free decomposition of p.
        static polynomial_sequence yun_algorithm(const polynomial<T> &p);

        /// @brief Gets the polynomial from the square-free decomposition.
        /// @param decomposition The square-free decomposition.
        /// @return The polynomial which square-free decomposition is given.
        static polynomial<T> from_square_free_decomposition(const polynomial_sequence& decomposition);
    };
}

#include "square_free_decomposition.tpp"

#endif // SQUARE_FREE_DECOMPOSITION_H_
