/// @file chebyshev_polynomial.tpp
///
/// @author Roland Abel
/// @date 08.10.2023

#include <ranges>
#include "chebyshev_polynomial.h"

namespace xmath {

    namespace {
        using std::numbers::pi;
        using std::views::iota;
        using std::views::transform;
    }

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
        return create_1st_kind(order, chebyshev_polynomial<T>::chebyshev_1st_kind_polynomials_);
    }

    template<typename T>
    polynomial<T> chebyshev_polynomial<T>::create_1st_kind(
            size_t order,
            chebyshev_polynomial<T>::polynomial_sequence &chebyshev_cache) {
        static auto zero = polynomial<T>::zero();
        static auto X = polynomial<T>::monomial(1);

        if (order < chebyshev_cache.size()) {
            return chebyshev_cache[order];
        } else if (order == 0) {
            chebyshev_cache.push_back(polynomial<T>::one());
            return chebyshev_cache[0];
        } else if (order == 1) {
            chebyshev_cache.push_back(X);
            return chebyshev_cache[1];
        } else {
            auto T_nm1 = create_1st_kind(order - 1, chebyshev_cache); // T_{n-1}
            auto T_nm2 = create_1st_kind(order - 2, chebyshev_cache); // T_{n-2}

            auto T_n = 2 * X * T_nm1 - T_nm2;
            chebyshev_cache.push_back(T_n);

            return T_n;
        }
    }

    template<typename T>
    chebyshev_polynomial<T>::values_type chebyshev_polynomial<T>::chebyshev_nodes(
            size_t order,
            const real_interval<T> &I) {

        if (order == 0) {
            xmath::chebyshev_polynomial<T>::values_type();
        }

        const auto pi_half = std::numbers::pi / 2.;
        auto kth_node = [pi_half, order](const auto &k) {
            return std::cos(pi_half * ((2. * k + 1.) / order));
        };

        const auto A = .5 * (I.lower() + I.upper());
        const auto B = .5 * (I.upper() - I.lower());

        auto z_node = [A, B, &I](const auto &x) {
            return A + B * x;
        };
        return to_vector(iota(0, static_cast<int>(order)) | transform(kth_node) | transform(z_node));
    }

    template<typename T>
    chebyshev_polynomial<T>::value_type chebyshev_polynomial<T>::clenshaw(
            const values_type &alphas,
            const value_type &x) {
        auto beta1 = polynomial<T>::spec::zero;
        auto beta2 = polynomial<T>::spec::zero;

        for (int k = alphas.size() - 1; k > 0; k--) {
            auto beta = alphas[k] + 2. * x * beta1 - beta2;

            beta2 = beta1;
            beta1 = beta;
        }
        return alphas[0] + x * beta1 - beta2;
    }

    template<typename T>
    polynomial<T> chebyshev_polynomial<T>::chebyshev_series(const values_type &alphas) {
        static auto X = polynomial<T>::monomial(1);

        auto beta1 = polynomial<T>::zero();
        auto beta2 = polynomial<T>::zero();

        for (int k = alphas.size() - 1; k > 0; k--) {
            polynomial<T> beta = alphas[k] + 2. * X * beta1 - beta2;

            beta2 = beta1;
            beta1 = beta;
        }
        return alphas[0] + X * beta1 - beta2;
    }
}