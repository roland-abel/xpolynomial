/// @file chebyshev_polynomial.tpp
///
/// @author Roland Abel
/// @date 08.10.2023

#include "chebyshev_polynomial.h"

namespace xmath {

    template<typename T>
    chebyshev_polynomial<T>::polynomial_sequence create_1st_kind_chebyshev_polynomials() {
        const auto &one = polynomial<T>::one();
        const auto &X = polynomial<T>::monomial(1);

        return std::vector<polynomial<T>>(
                {
                        one,
                        X,
                        2 * X.pow(2) - 1,
                        4 * X.pow(3) - 3 * X,
                        8 * X.pow(4) - 8 * X.pow(2) + 1,
                });
    }

    template<typename T> chebyshev_polynomial<T>::polynomial_sequence
            chebyshev_polynomial<T>::chebyshev_1st_kind_polynomials_ = create_1st_kind_chebyshev_polynomials<T>();

    template<typename T>
    polynomial<T> chebyshev_polynomial<T>::create_1st_kind(size_t n) {
        static auto zero = polynomial<T>::zero();
        static auto X = polynomial<T>::monomial(1);

        auto &chebyshev = chebyshev_polynomial<T>::chebyshev_1st_kind_polynomials_;

        if (n < chebyshev.size()) {
            return chebyshev[n];
        }

        auto T_n = chebyshev[n - 2]; // T_n
        auto T_n1 = chebyshev[n - 1]; // T_{n+1}

        auto T_n2 = zero;   // T_{n+2}
        for (size_t i = chebyshev.size(); i <= n; ++i) {

            T_n2 = 2 * X * T_n1 - T_n;
            chebyshev.push_back(T_n2);

            T_n1 = T_n2;
            T_n = T_n1;
        }
        return T_n2;
    }
}