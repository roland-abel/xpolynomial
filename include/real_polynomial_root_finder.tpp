/// @file real_polynomial_root_finder.tpp
/// @brief Root finder class of polynomials with real coefficients.
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

#pragma once

#include <valarray>
#include <numeric>
#include "utils.h"
#include "root_finder.h"
#include "square_free_decomposition.h"
#include "real_polynomial_root_finder.h"

namespace xmath {
    template<typename T, typename FP=double_t>
    size_t sign_changes(const std::vector<T> &sequence, FP epsilon = 1e-5) {
        auto changes = 0;
        auto size = sequence.size();

        if (size <= 1) {
            return 0; // A sequence with only one value hasn't a sign change.
        }

        int prev_sign = sequence[0] >= 0 ? 1 : -1;

        for (int i = 1; i < size; ++i) {
            if (nearly_zero<T>(sequence[i], epsilon)) {
                continue;
            }

            int current_sign = sequence[i] >= 0 ? 1 : -1;
            if (current_sign != prev_sign) {
                changes++;
                prev_sign = current_sign;
            }
        }
        return changes;
    }

    template<typename T>
    std::optional<std::tuple<T, T> >
    real_polynomial_root_finder<T>::quadratic_roots(const polynomial<T> &p) {
        if (!p.is_quadratic()) {
            return {};
        }

        const auto a = p[2];
        const auto b = p[1];
        const auto c = p[0];

        auto d = b * b - 4 * a * c;

        if (d < 0) {
            return {};
        }

        return std::make_tuple(
            (-b + std::sqrt(d)) / (2. * a),
            (-b - std::sqrt(d)) / (2. * a));
    }

    template<typename T>
    bool real_polynomial_root_finder<T>::has_cubic_normal_form(const polynomial<T> &p) {
        return p.is_cubic()
               && nearly_equal<T>(p[3], 1)
               && nearly_zero<T>(p[2], 1);
    }

    template<typename T>
    std::vector<T> real_polynomial_root_finder<T>::cubic_roots(const polynomial<T> &p) {
        if (!p.is_cubic()) {
            return {};
        }

        // Transform into the cubic normal form polynomial Y^3 + aY + b by the substitution X = Y - p[2]/3.
        const auto q = p.normalize();
        const auto a = (1. / 3.) * (3. * q[1] - std::pow(q[2], 2));
        const auto b = (1. / 27.) * (2. * std::pow(q[2], 3) - 9. * q[1] * q[2] + 27. * q[0]);

        auto p_normal = polynomial<T>({b, a, 0, 1});

        // Gets the real roots of the normal form polynomial
        auto roots = cubic_normal_form_roots(p_normal);

        // Re-substitution the roots
        return to_vector(roots | std::views::transform([&q](T r) { return r - q[2] / 3.; }));
    }

    template<typename T>
    std::vector<T> real_polynomial_root_finder<T>::cubic_normal_form_roots(const polynomial<T> &p) {
        if (!has_cubic_normal_form(p)) {
            return {};
        }

        // The polynomial p is given by the normal form p = X^3 + aX + b
        const auto a = p[1];
        const auto b = p[0];

        const auto w = (std::pow(b, 2) / 4.) + (std::pow(a, 3) / 27.);
        if (nearly_zero(w)) {
            if (nearly_zero(b)) {
                return {0., 0., 0.};
            }

            auto r1 = std::sqrt(-a / 3.);
            auto r2 = -2 * r1;

            if (b > 0) {
                return {r2, r1, r1};
            }
            if (b < 0) {
                return {r2, -r1, -r1};
            }
        }

        if (w > 0) {
            // There are one real root and two conjugate imaginary roots
            const auto A = std::pow(-b / 2. + std::sqrt(w), 1. / 3.);
            const auto B = std::pow(-b / 2. - std::sqrt(w), 1. / 3.);

            return {A + B};
        }

        // There are three unequal real roots (w < 0)
        auto t = std::sqrt(-(std::pow(b, 2) / 4.) / (std::pow(a, 3) / 27.));
        auto phi = std::acos(b > 0 ? -t : t);

        auto root = [a, phi](const int k) {
            return 2 * std::sqrt(-a / 3.) * std::cos((phi + 2 * k * std::numbers::pi) / 3.);
        };

        return {root(0), root(1), root(2)};
    }

    template<typename T>
    std::optional<T> real_polynomial_root_finder<T>::newton_raphson(
        const polynomial<T> &p,
        value_type initial,
        int max_iterations,
        value_type tolerance) {
        return root_finder<T>::newton_raphson(
            p,
            p.derive(),
            initial,
            max_iterations,
            tolerance);
    }

    template<typename T>
    size_t real_polynomial_root_finder<T>::sign_changes(const polynomial<T> &p) {
        return xmath::sign_changes(p.coefficients(), polynomial<T>::epsilon);
    }

    template<typename T>
    std::optional<T> real_polynomial_root_finder<T>::cauchy_bounds(const polynomial<T> &p) {
        if (p.is_zero()) {
            return {};
        }

        const auto zero = polynomial<T>::spec::zero;
        const auto one = polynomial<T>::spec::one;

        auto num_coefficients = p.coefficients().size();
        auto max = zero;

        for (auto i = 0; i < num_coefficients - 1; ++i) {
            max = std::max(max, std::abs(p[i]));
        }
        return one + max / std::abs(p.leading_coefficient());
    }

    template<typename T>
    std::optional<T> real_polynomial_root_finder<T>::lagrange_bounds(const polynomial<T> &p) {
        if (p.is_zero()) {
            return {};
        }
        const auto zero = polynomial<T>::spec::zero;
        const auto one = polynomial<T>::spec::one;

        const auto lc = p.leading_coefficient();
        const auto sum = std::accumulate(
            p.coefficients().begin(),
            p.coefficients().end() - 1, zero, [=](auto s, auto c) {
                return s + std::abs(c / lc);
            });
        return std::max(one, sum);
    }

    template<typename T>
    typename real_polynomial_root_finder<T>::polynomial_sequence
    real_polynomial_root_finder<T>::sturm_sequence(const polynomial<T> &p) {
        auto seq = std::vector<polynomial<T> >();

        seq.push_back(p);
        seq.push_back(p.derive());

        size_t i = 1;
        while (!seq[i].is_constant()) {
            seq.push_back(-(seq[i - 1] % seq[i]));
            i++;
        }

        return seq;
    }

    template<typename T>
    std::vector<short> real_polynomial_root_finder<T>::sign_variations(
        const polynomial_sequence &seq,
        const value_type &x) {
        auto variations = std::vector<short>();
        for (auto p: seq) {
            auto y = p(x);
            if (!nearly_zero(y)) {
                variations.push_back(y < 0 ? -1 : 1);
            }
        }
        return variations;
    }

    template<typename T>
    std::optional<int> real_polynomial_root_finder<T>::number_distinct_roots(const polynomial<T> &p, const real_interval<T> &I) {
        if (!I.is_lower_open() || !I.is_upper_closed()) {
            return {};
        }

        if (!square_free_decomposition<T>::is_square_free(p)) {
            return {};
        }

        auto seq = sturm_sequence(p);
        return xmath::sign_changes(sign_variations(seq, I.lower()))
               - xmath::sign_changes(sign_variations(seq, I.upper()));
    }

    template<typename T>
    std::optional<int> real_polynomial_root_finder<T>::number_distinct_roots(const polynomial<T> &p) {
        return cauchy_bounds(p).and_then([&p](T bound) {
            return number_distinct_roots(p, real_interval<T>(
                                             -bound,
                                             bound,
                                             real_interval<T>::interval_bounds::opened,
                                             real_interval<T>::interval_bounds::closed));
        });
    }

    template<typename T>
    std::vector<real_interval<T> > real_polynomial_root_finder<T>::root_isolation(const polynomial<T> &p) {
        auto intervals = std::vector<real_interval<T> >();
        if (p.is_constant() || !square_free_decomposition<T>::is_square_free(p)) {
            return intervals;
        }

        auto seq = sturm_sequence(p); // The Sturm's polynomial sequence.

        auto number_of_roots = [&seq](const real_interval<T> &I) {
            // Gets the number of roots of the polynomial p within the real_interval I.
            return xmath::sign_changes(sign_variations(seq, I.lower())) - xmath::sign_changes(sign_variations(seq, I.upper()));
        };

        std::function<void(const real_interval<T> &)> root_isolation = [&, p, number_of_roots](const real_interval<T> &I) {
            const auto epsilon = polynomial<T>::epsilon;
            auto r = number_of_roots(I);
            if (r == 0) {
                // nothing to do
            } else if (r == 1) {
                // Interval contains exactly one root
                intervals.push_back(I);
            } else {
                const auto c = (I.lower() + I.upper()) / static_cast<T>(2.);
                const auto t = nearly_zero(p(c), epsilon) ? epsilon : 0;

                const auto J_left = real_interval(I.lower(), c + t);
                const auto J_right = real_interval(c + t, I.upper());

                root_isolation(J_left);
                root_isolation(J_right);
            }
        };

        if (auto bounds = cauchy_bounds(p); bounds.has_value()) {
            const auto a = bounds.value();
            root_isolation(real_interval<T>(-a, a));
        }

        return intervals;
    }

    template<typename T>
    std::tuple<typename real_polynomial_root_finder<T>::roots_type, typename real_polynomial_root_finder<T>::multiplicities_type>
    real_polynomial_root_finder<T>::find_roots(const polynomial<T> &p, T epsilon) {
        auto roots = std::vector<value_type>();
        auto multiplicities = std::vector<unsigned short>();

        auto square_free_seq = square_free_decomposition<T>::yun_algorithm(p).value();
        for (int k = 0; k < square_free_seq.size(); ++k) {
            auto q = square_free_seq[k];
            auto intervals = root_isolation(q);

            for (auto I: intervals) {
                multiplicities.push_back(k + 1);

                auto root = root_finder<T>::bisection(q, I, epsilon);
                if (root.has_value()) {
                    roots.push_back(root.value());
                }
            }
        }
        return std::make_tuple(roots, multiplicities);
    }
}
