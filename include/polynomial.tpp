/// @file polynomial.tpp
/// @brief This file contains the implementation of the polynomial class.
///
/// @author Roland Abel
/// @date 19.08.2023

#include <utility>
#include <limits>
#include <algorithm>
#include <iostream>
#include "utils.h"
#include "polynomial.h"

namespace xmath {

    template<typename T>
    polynomial<T>::polynomial() {
        coeffs_.push_back(0);
    }

    template<typename T>
    polynomial<T>::polynomial(std::initializer_list<coeff_type> coeffs)
            : coeffs_(coeffs) {
        trim_coefficients();
    }

    template<typename T>
    polynomial<T>::polynomial(coeffs_type coeffs)
            : coeffs_(std::move(coeffs)) {
        trim_coefficients();
    }

    template<typename T>
    polynomial<T>::polynomial(size_type degree)
            : coeffs_(std::max(static_cast<size_type>(1), degree + 1)) {
    }

    template<typename T>
    polynomial<T> polynomial<T>::zero() {
        return {0};
    }

    template<typename T>
    polynomial<T> polynomial<T>::one() {
        return {1};
    }

    template<typename T>
    polynomial<T> polynomial<T>::monomial(size_type degree, coeff_type coeff) {
        coeffs_type coeffs(degree + 1);
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
        while (coeffs_.size() > 1 && nearly_zero<value_type>(*coeff, tolerance)) {
            coeffs_.pop_back();
            ++coeff;
        }
        return *this;
    }

    template<typename T>
    bool polynomial<T>::is_zero() const {
        return degree() == 0 && nearly_zero<value_type>(leading_coefficient(), tolerance);
    }

    template<typename T>
    bool polynomial<T>::is_one() const {
        return degree() == 0 && nearly_equal(leading_coefficient(), 1., tolerance);
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
    polynomial<T>::coeff_type polynomial<T>::leading_coefficient() const {
        return coeffs_.back();
    }

    template<typename T>
    const polynomial<T>::coeffs_type &polynomial<T>::coefficients() const {
        return coeffs_;
    }

    template<typename T>
    polynomial<T>::coeff_type polynomial<T>::operator[](size_type index) const {
        return index > degree() ? (coeff_type) 0.0 : coeffs_[index];
    }

    template<typename T>
    polynomial<T>::coeff_type &polynomial<T>::operator[](size_type index) {
        return coeffs_[index];
    }

    template<typename T>
    bool polynomial<T>::operator==(const polynomial<T> &p) const {
        return (p.degree() == degree())
               && std::equal(
                coeffs_.cbegin(),
                coeffs_.cend(),
                p.coeffs_.cbegin(),
                [](coeff_type a, coeff_type b) {
                    return nearly_equal<value_type>(a, b, tolerance);
                });
    }

    template<typename T>
    polynomial<T>::coeff_type polynomial<T>::operator()(coeff_type x) const {
        return evaluate(x);
    }

    template<typename T>
    bool polynomial<T>::operator!=(const polynomial<T> &p) const {
        return !(*this == p);
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator+(coeff_type scalar) const {
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
        coeffs_[0] += scalar;
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
        coeffs_[0] -= scalar;
        return *this;
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator*(value_type scalar) const {
        auto result = polynomial<T>(coeffs_);
        for (auto &coefficient: result.coeffs_) {
            coefficient *= scalar;
        }
        return result;
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator*=(value_type scalar) {
        for (auto &coefficient: coeffs_) {
            coefficient *= scalar;
        }
        return *this;
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator/(value_type scalar) const {
        auto result = polynomial<T>(coeffs_);
        for (auto &coefficient: result.coeffs_) {
            coefficient /= scalar;
        }
        return result;
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator/=(value_type scalar) {
        for (auto &coefficient: coeffs_) {
            coefficient /= scalar;
        }
        return *this;
    }

    template<typename T>
    polynomial<T> polynomial<T>::operator+(const polynomial<T> &p) const {
        auto sum = polynomial<T>(std::max(p.degree(), degree()));
        for (auto i = 0; i < sum.degree() + 1; ++i) {
            sum[i] = operator[](i) + p[i];
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
                product[i + j] += p[i] * coeffs_[j];
            }
        }
        return product.trim_coefficients();
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator+=(const polynomial<T> &p) {
        coeffs_.resize(std::max(p.degree(), degree()) + 1);
        for (auto i = 0; i < coeffs_.size() + 1; ++i) {
            coeffs_[i] += p[i];
        }
        return trim_coefficients();
    }

    template<typename T>
    polynomial<T> &polynomial<T>::operator-=(const polynomial<T> &p) {
        coeffs_.resize(std::max(p.degree(), degree()) + 1);
        for (auto i = 0; i < coeffs_.size() + 1; ++i) {
            coeffs_[i] -= p[i];
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
            composition = composition + coeffs_[i] * q.pow(i);
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
        for (auto i = 0; i < derivative.degree() + 1; ++i) {
            derivative[i] = (coeff_type)(i + 1) * coeffs_[i + 1];
        }
        return derivative.trim_coefficients();
    }

    template<typename T>
    polynomial<T> polynomial<T>::integrate() const {
        if (is_constant()) {
            return zero();
        }

        auto integral = polynomial<T>(degree() + 1);
        for (auto i = 1; i < integral.degree() + 1; ++i) {
            integral[i] = (coeff_type)(1. / i) * coeffs_[i - 1];
        }
        return integral;
    }

    template<typename T>
    polynomial<T>::value_type polynomial<T>::evaluate(value_type x) const {
        auto value = leading_coefficient();
        for (long i = static_cast<long>(degree()) - 1; i >= 0; --i) {
            value = value * x + coeffs_[i];
        }
        return value;
    }

    template<typename T>
    polynomial<T> operator+(typename polynomial<T>::value_type scalar, const polynomial <T> &polynomial) {
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
            return polynomial<T>({-roots.front(), 1});
        }

        auto polynomial = one();
        for (auto root: roots) {
            polynomial = polynomial * (monomial(1) - root);
        }
        return polynomial;
    }

    template<typename T>
    polynomial<T> polynomial<T>::normalize() const {
        return *this / leading_coefficient();
    }

    template<typename T>
    bool polynomial<T>::has_root(const value_type &value) const {
        return nearly_zero<value_type>(evaluate(value), tolerance);
    }

    template<typename T>
    bool polynomial<T>::has_roots(const values_type &values) const {
        return std::all_of(values.begin(), values.end(), [&](const auto &value) {
            return has_root(value);
        });
    }

    template<typename T>
    std::tuple<polynomial<T>, polynomial<T>> polynomial<T>::divide(const polynomial<T> &divisor) const {
        auto quotient = zero();
        auto remainder = *this;

        long delta = static_cast<long>(remainder.degree() - divisor.degree());
        while (remainder != zero() && delta >= 0) {
            auto divider = (remainder.leading_coefficient() / divisor.leading_coefficient()) * monomial(delta, 1.0);
            quotient = quotient + divider;
            remainder = remainder - divisor * divider;

            delta = static_cast<long>(remainder.degree() - divisor.degree());
        }
        return std::make_tuple(quotient, remainder);
    }
}
