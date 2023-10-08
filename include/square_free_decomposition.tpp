/// @file square_free_decomposition.tpp
///
/// @author Roland Abel
/// @date 08.10.2023

#include "euclidean_algorithm.h"
#include "square_free_decomposition.h"

namespace xmath {

    template<typename T>
    bool square_free_decomposition<T>::is_square_free(const polynomial<T> &p) {
        return euclidean_algorithm<T>::euclidean(p, p.derive()).is_constant();
    }

    template<typename T>
    polynomial<T> square_free_decomposition<T>::from_square_free_decomposition(
            const square_free_decomposition::polynomial_sequence &decomposition) {

        auto q = polynomial<T>::one();
        for (int k = 0; k < decomposition.size(); ++k) {
            q *= decomposition[k].pow(k + 1);
        }
        return q;
    }

    template<typename T>
    square_free_decomposition<T>::polynomial_sequence
    square_free_decomposition<T>::yun_algorithm(const polynomial<T> &p) {
        auto decomposition = xmath::square_free_decomposition<T>::polynomial_sequence();

        auto p_prim = p.derive();
        auto q = euclidean_algorithm<T>::euclidean(p, p_prim);
        auto b = p / q;
        auto c = p_prim / q;
        auto d = c - b.derive();

        q = euclidean_algorithm<T>::euclidean(b, d);
        decomposition.push_back(q);

        while (!d.is_zero()) {
            b = b / q;
            c = d / q;
            d = c - b.derive();

            q = euclidean_algorithm<T>::euclidean(b, d);
            decomposition.push_back(q);
        }
        return decomposition;
    }
}

