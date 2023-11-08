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

    template<typename T>
    chebyshev_polynomial<T>::values_type chebyshev_polynomial<T>::roots_1st_kind(size_t order) {
        auto kth_root = [order](const auto &k) {
            return std::cos((2. * k + 1.) * std::numbers::pi / (2 * order));
        };

        return order == 0
               ? xmath::chebyshev_polynomial<T>::values_type()
               : to_vector(iota(0, static_cast<int>(order)) | transform(kth_root));
    }

    template<typename T>
    chebyshev_polynomial<T>::value_type chebyshev_polynomial<T>::clenshaw(
            const values_type &alphas,
            const value_type &x) {

        const auto n = alphas.size();
        std::function<T(const int k)> beta = [&](const int k) {
            if (k == n) {
                return alphas[n];
            } else if (k == n - 1) {
                return alphas[n - 1] + 2 * x * beta(k + 1);
            } else {
                return alphas[k] + 2 * x * beta(k + 1) - beta(k + 2);
            }
        };

        auto result = alphas[0] + x * beta(1) - beta(2);
        return result;
    }
    
    template<typename T>
    polynomial<T> chebyshev_polynomial<T>::chebyshev_series(const values_type &alphas) {
        auto series = polynomial<T>::zero();
        for (int k = 0; k < alphas.size(); ++k) {
            series += alphas[k] * chebyshev_polynomial<T>::create_1st_kind(k);
        }
        return series;
    }
}