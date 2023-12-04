/// @file chebyshev_polynomial.tpp
/// @brief
///
/// @author Roland Abel
/// @date October 8, 2023
///
/// Copyright (c) 2023 Roland Abel
///
/// This software is released under the MIT License.
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
/// THE SOFTWARE.

#include <ranges>
#include <numeric>
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
    polynomial<T> chebyshev_polynomial<T>::create_1st_kind(size_t N) {
        return create_1st_kind(N, chebyshev_polynomial<T>::chebyshev_1st_kind_polynomials_);
    }

    template<typename T>
    polynomial<T> chebyshev_polynomial<T>::create_1st_kind(
            size_t N,
            chebyshev_polynomial<T>::polynomial_sequence &chebyshev_cache) {
        static auto zero = polynomial<T>::zero();
        static auto X = polynomial<T>::monomial(1);

        if (N < chebyshev_cache.size()) {
            return chebyshev_cache[N];
        } else if (N == 0) {
            chebyshev_cache.push_back(polynomial<T>::one());
            return chebyshev_cache[0];
        } else if (N == 1) {
            chebyshev_cache.push_back(X);
            return chebyshev_cache[1];
        } else {
            auto T_nm1 = create_1st_kind(N - 1, chebyshev_cache); // T_{n-1}
            auto T_nm2 = create_1st_kind(N - 2, chebyshev_cache); // T_{n-2}

            auto T_n = 2 * X * T_nm1 - T_nm2;
            chebyshev_cache.push_back(T_n);

            return T_n;
        }
    }

    template<typename T>
    chebyshev_polynomial<T>::values_type chebyshev_polynomial<T>::chebyshev_nodes(
            size_t N,
            const real_interval<T> &I) {

        if (N == 0) {
            xmath::chebyshev_polynomial<T>::values_type();
        }

        const auto pi_half = std::numbers::pi / 2.;
        auto kth_node = [pi_half, N](const auto &k) {
            return std::cos(pi_half * ((2. * k - 1.) / N));
        };

        const auto A = .5 * (I.lower() + I.upper());
        const auto B = .5 * (I.upper() - I.lower());

        auto z_node = [A, B](const auto &x) {
            return A + B * x;
        };
        return to_vector(iota(1, static_cast<int>(N + 1))
                         | transform(kth_node)
                         | transform(z_node));
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

    template<typename T>
    chebyshev_polynomial<T>::value_type chebyshev_polynomial<T>::chebyshev_quadrature(
            std::function<value_type(value_type)> func,
            uint32_t N) {

        auto value = polynomial<T>::spec::zero;
        for (int i = 1; i <= N; ++i) {
            value += func(std::cos((2. * i - 1.) * pi / (2. * N)));
        }
        return (pi / N) * value;
    }
}
