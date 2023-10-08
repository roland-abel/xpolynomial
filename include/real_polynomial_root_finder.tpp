//
// Created by abel on 26.08.2023.
//

#include <limits>
#include <valarray>
#include <numeric>
#include "utils.h"
#include "root_finder.h"
#include "square_free_decomposition.h"
#include "real_polynomial_root_finder.h"

namespace xmath {

    template<typename T>
    std::tuple<typename polynomial<T>::value_type, typename polynomial<T>::value_type>
    real_polynomial_root_finder<T>::quadratic_roots(const polynomial<T> &p) {
        const auto NaN = std::numeric_limits<T>::quiet_NaN();

        if (!p.is_quadratic()) {
            return std::make_pair(NaN, NaN);
        }

        const auto a = p[2];
        const auto b = p[1];
        const auto c = p[0];

        auto d = b * b - 4 * a * c;

        if (d < 0) {
            return std::make_pair(NaN, NaN);
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
    std::tuple<typename polynomial<T>::value_type, typename polynomial<T>::value_type, typename polynomial<T>::value_type>
    real_polynomial_root_finder<T>::cubic_roots(const polynomial<T> &p) {
        const auto NaN = std::numeric_limits<T>::quiet_NaN();

        if (!p.is_cubic()) {
            return std::make_tuple(NaN, NaN, NaN);
        }

        // Transform into the cubic normal form polynomial Y^3 + aY + b by the substitution X = Y - p[2]/3.
        const auto q = p.normalize();
        const auto a = (1. / 3.) * (3. * q[1] - std::pow(q[2], 2));
        const auto b = (1. / 27.) * (2. * std::pow(q[2], 3) - 9. * q[1] * q[2] + 27. * q[0]);

        auto p_normal = polynomial<T>({b, a, 0, 1});

        // Gets the real roots of the normal form polynomial
        auto roots = cubic_normal_form_roots(p_normal);

        // Re-substitution the roots
        auto rsub = [&q](T r) { return r - q[2] / 3.; };
        return std::make_tuple(
                rsub(std::get<0>(roots)),
                rsub(std::get<1>(roots)),
                rsub(std::get<2>(roots)));
    }

    template<typename T>
    std::tuple<typename polynomial<T>::value_type, typename polynomial<T>::value_type, typename polynomial<T>::value_type>
    real_polynomial_root_finder<T>::cubic_normal_form_roots(const polynomial<T> &p) {
        const auto NaN = std::numeric_limits<T>::quiet_NaN();

        if (!has_cubic_normal_form(p)) {
            return std::make_tuple(NaN, NaN, NaN);
        }

        // The polynomial p is given by the normal form p = X^3 + aX + b
        const auto a = p[1];
        const auto b = p[0];

        const auto w = (std::pow(b, 2) / 4.) + (std::pow(a, 3) / 27.);
        if (nearly_zero(w)) {
            // There are three real roots of which at least two are equal
            if (nearly_zero(b)) {
                return std::make_tuple(0., 0., 0.);
            }

            auto r1 = std::sqrt(-a / 3.);
            auto r2 = -2 * r1;

            if (b > 0) {
                return std::make_tuple(r2, r1, r1);
            }
            if (b < 0) {
                return std::make_tuple(r2, -r1, -r1);
            }
        }

        if (w > 0) {
            // There are one real root and two conjugate imaginary roots
            const auto A = std::pow(-b / 2. + std::sqrt(w), 1. / 3.);
            const auto B = std::pow(-b / 2. - std::sqrt(w), 1. / 3.);

            return std::make_tuple(A + B, NaN, NaN);
        }

        // There are three unequal real roots (w < 0)
        auto t = std::sqrt(-(std::pow(b, 2) / 4.) / (std::pow(a, 3) / 27.));
        auto phi = std::acos(b > 0 ? -t : t);

        auto root = [a, phi](int k) {
            return 2 * std::sqrt(-a / 3.) * std::cos((phi + 2 * k * std::numbers::pi) / 3.);
        };

        return std::make_tuple(root(0), root(1), root(2));
    }

    template<typename T>
    polynomial<T>::value_type real_polynomial_root_finder<T>::newton_raphson(
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
        return xmath::sign_changes(p.coefficients(), polynomial<T>::tolerance);
    }

    template<typename T>
    polynomial<T>::value_type real_polynomial_root_finder<T>::cauchy_bound(const polynomial<T> &p) {
        const auto NaN = std::numeric_limits<T>::quiet_NaN();

        if (p.is_zero()) {
            return NaN;
        }

        auto num_coefficients = p.coefficients().size();
        auto max = .0;

        for (auto i = 0; i < num_coefficients - 1; ++i) {
            max = std::max((double) max, (double) std::abs(p[i]));
        }

        return 1 + max / std::abs(p.leading_coefficient());
    }

    template<typename T>
    real_polynomial_root_finder<T>::polynomial_sequence
    real_polynomial_root_finder<T>::sturm_sequence(const polynomial<T> &p) {
        auto seq = std::vector<polynomial<T>>();

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
    int real_polynomial_root_finder<T>::number_distinct_roots(const polynomial<T> &p, const interval<T> &I) {
        auto NaN = std::numeric_limits<int>::quiet_NaN();

        if (!square_free_decomposition<T>::is_square_free(p)) {
            return NaN;
        }

        auto seq = sturm_sequence(p);
        return xmath::sign_changes(sign_variations(seq, I.start()))
               - xmath::sign_changes(sign_variations(seq, I.end()));
    }

    template<typename T>
    int real_polynomial_root_finder<T>::number_distinct_roots(const polynomial<T> &p) {
        auto bound = cauchy_bound(p);
        return number_distinct_roots(p, interval<T>(-bound, bound));
    }

    template<typename T>
    std::vector<interval<T>> real_polynomial_root_finder<T>::root_isolation(const polynomial<T> &p) {

        auto intervals = std::vector<interval<T>>();
        if (p.is_constant()) {
            return intervals;
        }

        auto seq = sturm_sequence(p); // The Sturm's polynomial sequence.

        // Gets the number of roots of the polynomial p within the interval I.
        auto number_of_roots = [&seq](const interval<T> &I) {
            return xmath::sign_changes(sign_variations(seq, I.start()))
                   - xmath::sign_changes(sign_variations(seq, I.end()));
        };

        std::function<void(const interval<T> &)> root_isolation = [&, number_of_roots](const interval<T> &I) {
            auto r = number_of_roots(I);
            if (r == 0) {
                // nothing to do
            } else if (r == 1) {
                // Interval contains exactly one root
                intervals.push_back(I);
            } else {
                const auto J = I.bisect();

                root_isolation(J.first);
                root_isolation(J.second);
            }
        };

        auto a = cauchy_bound(p);
        root_isolation(interval<T>(-a, a));

        return intervals;
    }

    template<typename T>
    std::vector<typename polynomial<T>::value_type> real_polynomial_root_finder<T>::find_roots(
            const polynomial<T> &p,
            value_type tolerance) {
        auto roots = std::vector<typename polynomial<T>::value_type>();
        auto intervals = root_isolation(p);

        for (auto I: intervals) {
            roots.push_back(root_finder<T>::bisection(p, I, 1e-9));
        }
        return roots;
    }
}
