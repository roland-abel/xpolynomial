/// @file polynomial.tpp
/// @brief This file contains the implementation of the polynomial class.
///
/// @author Roland Abel
/// @date August 19, 2023
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

#include <utility>
#include <limits>
#include <ranges>
#include <algorithm>
#include <iostream>
#include "utils.h"
#include "polynomial.h"

namespace xmath {

    namespace {
        using std::views::transform;
        using std::ranges::for_each;
    }

    template<typename T>
    polynomial<T>::polynomial() {
        coeffs_.push_back(0);
    }

    template<typename T>
    polynomial<T>::polynomial(polynomial<T> &&p) noexcept
            : coeffs_(std::move(p.coeffs_)) {
    }

    template<typename T>
    polynomial<T>::polynomial(std::initializer_list<value_type> coeffs)
            : coeffs_(coeffs) {
        trim_coefficients();
    }

    template<typename T>
    polynomial<T>::polynomial(const values_type &coeffs)
            : coeffs_(std::move(coeffs)) {
        trim_coefficients();
    }

    template<typename T>
    polynomial<T>::polynomial(size_type degree)
            : coeffs_(std::max(static_cast<size_type>(1), degree + 1)) {
    }

    template<typename T>
    polynomial<T>::polynomial(const std::ranges::range auto &range)
            : polynomial<T>(std::vector<T>{range.begin(), range.end()}) {
    }

    template<typename T>
    polynomial<T> polynomial<T>::zero() {
        return {spec::zero};
    }

    template<typename T>
    polynomial<T> polynomial<T>::one() {
        return {spec::one};
    }

    template<typename T>
    bool polynomial<T>::nearly_equal(value_type a, value_type b) {
        return xmath::nearly_equal<value_type, floating_point_type>(a, b, epsilon);
    }

    template<typename T>
    bool polynomial<T>::nearly_zero(value_type a) {
        return xmath::nearly_zero<value_type, floating_point_type>(a, epsilon);
    }

    template<typename T>
    polynomial<T> polynomial<T>::monomial(size_type degree, value_type coeff) {
        values_type coeffs(degree + 1);
        coeffs[degree] = coeff;

        return polynomial<T>(coeffs);
    }

    template<typename T>
    polynomial<T>::size_type polynomial<T>::degree() const {
        return coeffs_.size() - 1;
    }

    template<typename T>
    polynomial<T> &polynomial<T>::trim_coefficients() {
        auto coeff = coeffs_.rbegin();
        while (coeffs_.size() > 1 && nearly_zero(*coeff)) {
            coeffs_.pop_back();
            ++coeff;
        }
        return *this;
    }

    template<typename T>
    bool polynomial<T>::is_zero() const {
        return degree() == 0 && nearly_zero(leading_coefficient());
    }

    template<typename T>
    bool polynomial<T>::is_one() const {
        return degree() == 0 && nearly_equal(leading_coefficient(), spec::one);
    }

    template<typename T>
    bool polynomial<T>::is_constant() const {
        return degree() == 0;
    }

    template<typename T>
    bool polynomial<T>::is_linear() const {
        return degree() <= 1;
    }

    template<typename T>
    bool polynomial<T>::is_quadratic() const {
        return degree() == 2;
    }

    template<typename T>
    bool polynomial<T>::is_cubic() const {
        return degree() == 3;
    }

    template<typename T>
    bool polynomial<T>::is_normalized() const {
        return nearly_equal(leading_coefficient(), spec::one);
    }

    template<typename T>
    bool polynomial<T>::is_integer() const {
        return std::all_of(coefficients().begin(), coefficients().end(), [](double c) {
            return nearly_equal(c, std::round(c));
        });
    }

    template<typename T>
    polynomial<T> polynomial<T>::to_integer() const {
        auto q = polynomial<T>(degree());
        for (auto i = 0; i < degree() + 1; ++i) {
            q[i] = std::round(at(i));
        }
        return q.trim_coefficients();
    }

    template<typename T>
    polynomial<T>::value_type polynomial<T>::leading_coefficient() const {
        return coeffs_.back();
    }

    template<typename T>
    inline const polynomial<T>::values_type &polynomial<T>::coefficients() const {
        return coeffs_;
    }

    template<typename T>
    inline polynomial<T>::value_type polynomial<T>::at(size_type index) const {
        return index > degree() ? spec::zero : coeffs_[index];
    }

    template<typename T>
    inline polynomial<T>::value_type &polynomial<T>::at(size_type index) {
        return coeffs_[index];
    }

    template<typename T>
    inline polynomial<T>::value_type polynomial<T>::operator[](size_type index) const {
        return at(index);
    }

    template<typename T>
    inline polynomial<T>::value_type &polynomial<T>::operator[](size_type index) {
        return at(index);
    }

    template<typename T>
    bool polynomial<T>::operator==(const polynomial<T> &p) const {
        auto is_equal = [](value_type a, value_type b) { return nearly_equal(a, b); };
        return (p.degree() == degree())
               && std::equal(coefficients().cbegin(), coefficients().cend(),
                             p.coefficients().cbegin(), is_equal);
    }

    template<typename T>
    polynomial<T>::value_type polynomial<T>::operator()(value_type x) const {
        return evaluate(x);
    }

    template<typename T>
    bool polynomial<T>::operator!=(const polynomial<T> &p) const {
        return !(*this == p);
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator+(value_type scalar) const {
        auto coeffs = coefficients();
        coeffs[0] += scalar;

        return polynomial<T>(coeffs);
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator+() const {
        return *this;
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator+=(value_type scalar) {
        at(0) += scalar;
        return *this;
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator-(value_type scalar) const {
        return operator+(-scalar);
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator-() const {
        return (-1) * *this;
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator-=(value_type scalar) {
        at(0) -= scalar;
        return *this;
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator*(value_type scalar) const {
        return polynomial<T>(coefficients() | transform([&](const T &c) {
            return c * scalar;
        }));
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator*=(value_type scalar) {
        for_each(coeffs_, [&](value_type &c) {
            c *= scalar;
        });
        return *this;
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator/(value_type scalar) const {
        return polynomial<T>(coefficients() | transform([&](const T &c) {
            return c / scalar;
        }));
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator/=(value_type scalar) {
        for_each(coeffs_, [&](value_type &c) {
            c /= scalar;
        });
        return *this;
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator+(const polynomial<T> &p) const {
        auto sum = polynomial<T>(std::max(p.degree(), degree()));
        for (auto i = 0; i < sum.degree() + 1; ++i) {
            sum[i] = at(i) + p[i];
        }
        return sum.trim_coefficients();
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator-(const polynomial<T> &p) const {
        return *this + (-1) * p;
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator*(const polynomial<T> &p) const {
        auto product = polynomial<T>(p.degree() + degree());
        for (auto i = 0; i < p.degree() + 1; ++i) {
            for (auto j = 0; j < degree() + 1; ++j) {
                product[i + j] += p[i] * at(j);
            }
        }
        return product.trim_coefficients();
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator+=(const polynomial<T> &p) {
        coeffs_.resize(std::max(p.degree(), degree()) + 1);
        for (auto i = 0; i < coeffs_.size() + 1; ++i) {
            at(i) += p[i];
        }
        return trim_coefficients();
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator-=(const polynomial<T> &p) {
        coeffs_.resize(std::max(p.degree(), degree()) + 1);
        for (auto i = 0; i < coeffs_.size() + 1; ++i) {
            at(i) -= p[i];
        }
        return trim_coefficients();
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator*=(const polynomial<T> &p) {
        coeffs_ = std::move(((*this) * p).coeffs_);
        return *this;
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator/=(const polynomial<T> &p) {
        coeffs_ = std::move(((*this) / p).coeffs_);
        return *this;
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator/(const polynomial<T> &p) const {
        return get<0>(divide(p));
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator%(const polynomial<T> &p) const {
        return get<1>(divide(p));
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator%=(const polynomial<T> &p) {
        coeffs_ = std::move(((*this) % p).coeffs_);
        return *this;
    }

    template<typename T>
    polynomial<T> polynomial<T>::compose(const polynomial<T> &q) const {
        auto composition = zero();
        for (size_t i = 0; i < degree() + 1; ++i) {
            composition = composition + at(i) * q.pow(i);
        }
        return composition;
    }

    template<typename T>
    polynomial<T> polynomial<T>::pow(unsigned int exponent) const {
        if (exponent == 0) {
            return one();
        } else if (exponent == 1) {
            return *this;
        } else if (exponent == 2) {
            return *this * *this;
        }

        auto power = one();
        auto base = *this;

        while (exponent > 0) {
            if (exponent % 2 == 1) {
                power = power * base;
            }
            base = base * base;
            exponent /= 2;
        }
        return power;
    }

    template<typename T>
    polynomial<T> polynomial<T>::derive() const {
        if (is_constant()) {
            return zero();
        }

        auto derivative = polynomial<T>(degree() - 1);
        auto index_coeff_pairs =
                std::views::iota((size_t) 1, degree() + 1) | transform([&](int index) {
                    return std::make_pair(index, at(index));
                });

        auto derive = [&](const auto &index_value) {
            auto [exponent, coeff] = index_value;
            derivative[exponent - 1] = T(exponent) * coeff;
        };

        for_each(index_coeff_pairs, derive);
        return derivative.trim_coefficients();
    }

    template<typename T>
    polynomial<T> polynomial<T>::integrate() const {
        if (is_constant()) {
            return {0, at(0)};
        }

        auto primitive = polynomial<T>(degree() + 1);
        auto index_coeff_pairs =
                std::views::iota((size_t) 1, primitive.degree() + 1) | transform([&](int index) {
                    return std::make_pair(index, at(index - 1));
                });

        auto anti_derive = [&](const auto &index_value) {
            auto [exponent, coeff] = index_value;
            primitive[exponent] = (value_type)(1. / exponent) * at(exponent - 1);
        };

        for_each(index_coeff_pairs, anti_derive);
        return primitive;
    }

    template<typename T>
    polynomial<T>::value_type polynomial<T>::evaluate(value_type x) const {
        auto value = leading_coefficient();
        for (long i = static_cast<long>(degree()) - 1; i >= 0; --i) {
            value = value * x + at(i);
        }
        return value;
    }

    template<typename T>
    polynomial<T> operator+(typename polynomial<T>::value_type scalar, const polynomial<T> &polynomial) {
        return polynomial.operator+(scalar);
    }

    template<typename T>
    polynomial<T> operator*(typename polynomial<T>::value_type scalar, const polynomial <T> &polynomial) {
        return polynomial.operator*(scalar);
    }

    template<typename T>
    polynomial<T> polynomial<T>::from_roots(const values_type &roots) {
        const auto num_roots = roots.size();
        if (num_roots == 0) {
            return one();
        } else if (num_roots == 1) {
            return polynomial<T>({-roots.front(), spec::one});
        }

        auto polynomial = one();
        for (auto root: roots) {
            polynomial = polynomial * (monomial(1) - root);
        }
        return polynomial.normalize();
    }

    template<typename T>
    polynomial<T> polynomial<T>::normalize() const {
        return *this / leading_coefficient();
    }

    template<typename T>
    bool polynomial<T>::has_root(const value_type &value) const {
        return nearly_zero(evaluate(value));
    }

    template<typename T>
    bool polynomial<T>::has_roots(const values_type &values) const {
        return std::all_of(values.begin(), values.end(), [&](const auto &value) {
            return has_root(value);
        });
    }

    template<typename T>
    std::tuple<polynomial<T>, polynomial<T>> polynomial<T>::divide(const polynomial<T> &divisor) const {
        auto q = zero();
        auto r = *this;

        long delta = static_cast<long>(r.degree() - divisor.degree());
        while (r != zero() && delta >= 0) {
            const auto t = (r.leading_coefficient() / divisor.leading_coefficient()) * monomial(delta, spec::one);

            q = q + t;
            r = r - divisor * t;

            delta = static_cast<long>(r.degree() - divisor.degree());
        }
        return std::make_tuple(q, r);
    }
}
