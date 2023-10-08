/// @file euclidean_algorithm.tpp
///
/// @author Roland Abel
/// @date 08.10.2023

#include "euclidean_algorithm.h"

namespace xmath {

    template<typename T>
    polynomial<T> euclidean_algorithm<T>::euclidean(const polynomial<T> &p, const polynomial<T> &q) {
        auto a = p;
        auto b = q;

        while (!b.is_zero()) {
            auto r = a % b; // remainder
            a = b;
            b = r;
        }
        return a.normalize();
    }

    template<typename T>
    std::tuple<polynomial<T>, polynomial<T>, polynomial<T>>
    euclidean_algorithm<T>::extended_euclidean(const polynomial<T> &p, const polynomial<T> &q) {
        auto zero = polynomial<T>::zero();
        auto one = polynomial<T>::one();

        auto a = p;
        auto b = q;

        auto a1 = one;
        auto a2 = zero;

        auto b1 = zero;
        auto b2 = one;

        while (!b.is_zero()) {
            auto [quotient, remainder] = a.divide(b);

            a = b;
            b = remainder;

            auto r1 = a1 - quotient * b1;
            auto r2 = a2 - quotient * b2;

            a1 = b1;
            a2 = b2;

            b1 = r1;
            b2 = r2;
        }
        return std::make_tuple(a1, a2, a);
    }
}