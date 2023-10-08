/// @file legendre_polynomial.tpp
///
/// @author Roland Abel
/// @date 08.10.2023

#include "utils.h"
#include "legendre_polynomial.h"

namespace xmath {

    template<typename T>
    std::vector<polynomial<T>> create_legendre_polynomials() {
        auto one = polynomial<T>::one();
        auto X = polynomial<T>::monomial(1);
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
    auto legendre_polynomials = create_legendre_polynomials<T>();

    template<typename T>
    legendre_polynomial<T>::legendre_polynomial(size_t order)
            : polynomial<T>(create(order)) {
    }

    template<typename T>
    const std::vector<typename polynomial<T>::value_type> &legendre_polynomial<T>::weights() const {
        return weight_;
    }

    template<typename T>
    const std::vector<typename polynomial<T>::value_type> &legendre_polynomial<T>::roots() const {
        return roots_;
    }

    template<typename T>
    polynomial<T> legendre_polynomial<T>::create(size_t order) {
        static auto zero = polynomial<T>::zero();
        static auto X = polynomial<T>::monomial(1);

        if (order < legendre_polynomials<T>.size()) {
            return legendre_polynomials<T>[order];
        }

        auto pnm2 = legendre_polynomials<T>[order - 2]; // P_{n-2}
        auto pnm1 = legendre_polynomials<T>[order - 1]; // P_{n-1}

        auto pn = zero;
        for (size_t i = legendre_polynomials<T>.size(); i <= order; ++i) {
            auto n = static_cast<double >(order);

            pn = ((2. * n - 1.) / n) * X * pnm1 - ((n - 1.) / n) * pnm2;
            legendre_polynomials<T>.push_back(pn);

            pnm2 = pnm1;
            pnm1 = pn;
        }

        return pn;
    }
}
