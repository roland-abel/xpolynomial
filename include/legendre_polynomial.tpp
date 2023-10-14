/// @file legendre_polynomial.tpp
///
/// @author Roland Abel
/// @date 08.10.2023

#include "utils.h"
#include "legendre_polynomial.h"

namespace xmath {

    template<typename T>
    std::vector<polynomial<T>> create_legendre_polynomials() {
        const auto &one = polynomial<T>::one();
        const auto &X = polynomial<T>::monomial(1);

        return std::vector<polynomial<T>>(
                {
                        one,
                        X,
                        .5 * (3 * X.pow(2) - 1),
                        .5 * (5 * X.pow(3) - 3 * X),
                        (1. / 8) * (35. * X.pow(4) - 30 * X.pow(2) + 3),
                        (1. / 8) * (63 * X.pow(5) - (70 * X.pow(3)) + (15 * X)),
                        (1. / 16) * (231 * X.pow(6) - 315 * X.pow(4) + 105 * X.pow(2) - 5),
                        (1. / 16) * (429 * X.pow(7) - 693 * X.pow(5) + 315 * X.pow(3) - 35 * X),
                        (1. / 128) * (6435 * X.pow(8) - 12012 * X.pow(6) + 6930 * X.pow(4) - 1260 * X.pow(2) + 35)
                });
    }

    template<typename T>
    legendre_polynomial<T>::polynomial_sequence
            legendre_polynomial<T>::legendre_polynomial_ = create_legendre_polynomials<T>();

    template<typename T>
    polynomial<T> legendre_polynomial<T>::create(size_t order) {
        static auto zero = polynomial<T>::zero();
        static auto X = polynomial<T>::monomial(1);

        auto &legendre = legendre_polynomial<T>::legendre_polynomial_;

        if (order < legendre.size()) {
            return legendre[order];
        }

        auto P_nm1 = legendre[order - 1]; // P_{n-1}
        auto P_nm2 = legendre[order - 2]; // P_{n-2}

        auto P_n = zero;   // P_n
        for (size_t i = legendre.size(); i <= order; ++i) {
            const auto n = static_cast<double>(order);

            P_n = ((2. * n - 1.) / n) * X * P_nm1 - ((n - 1.) / n) * P_nm2;
            legendre.push_back(P_n);

            P_nm1 = P_n;
            P_nm2 = P_nm1;
        }

        return P_n;
    }
}
