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
                        16 * X.pow(5) - 20 * X.pow(3) + 5 * X,
                        32 * X.pow(6) - 48 * X.pow(4) + 18 * X.pow(2) - 1,
                        64 * X.pow(7) - 112 * X.pow(5) + 56 * X.pow(3) - 7 * X,
                        128 * X.pow(8) - 256 * X.pow(6) + 160 * X.pow(4) - 32 * X.pow(2) + 1,
                        256 * X.pow(9) - 576 * X.pow(7) + 432 * X.pow(5) - 120 * X.pow(3) + 9 * X,
                        512 * X.pow(10) - 1280 * X.pow(8) + 1120 * X.pow(6) - 400 * X.pow(4) + 50 * X.pow(2) - 1,
                });
    }

    template<typename T>
    chebyshev_polynomial<T>::polynomial_sequence
            chebyshev_polynomial<T>::chebyshev_1st_kind_polynomials_ = create_1st_kind_chebyshev_polynomials<T>();

    template<typename T>
    polynomial<T> chebyshev_polynomial<T>::create_1st_kind(size_t order) {
        static auto zero = polynomial<T>::zero();
        static auto X = polynomial<T>::monomial(1);

        auto &chebyshev = chebyshev_polynomial<T>::chebyshev_1st_kind_polynomials_;

        if (order < chebyshev.size()) {
            return chebyshev[order];
        }

        auto T_nm1 = chebyshev[order - 1]; // T_{n-1}
        auto T_nm2 = chebyshev[order - 2]; // T_{n+1}

        auto T_n = zero;   // T_n
        for (size_t i = chebyshev.size(); i <= order; ++i) {

            T_n = 2 * X * T_nm1 - T_nm2;
            chebyshev.push_back(T_n);

            T_nm1 = T_n;
            T_nm2 = T_nm1;
        }
        return T_n;
    }
}