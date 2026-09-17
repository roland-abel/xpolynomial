/// @file root_finder.tpp
/// @brief Root finder class using various numerical methods.
///
/// @author Roland Abel
/// @date October 8, 2023
///
/// Copyright (c) 2026 Roland Abel

#ifndef ROOT_FINDER_TPP_H_
#define ROOT_FINDER_TPP_H_

#include "utils.h"
#include "root_finder.h"

namespace xmath {

    template<typename T>
    std::optional<T> root_finder<T>::bisection(
            const std::function<value_type(value_type)> &func,
            const interval<value_type> &I,
            value_type epsilon) {

        if (greater_than_or_equal(func(I.lower()) * func(I.upper()), 0., epsilon)) {
            return {}; // incorrect endpoints a and b.
        }

        auto bisection = [&](value_type a, value_type b) {
            auto c = a;
            while (!nearly_zero<value_type>(b - a, epsilon)) {
                // Choose the middle point as the new estimation for the root
                c = (b + a) / 2.;
                if (nearly_zero<value_type>(func(c))) {
                    return c;
                } else if (func(c) * func(a) < 0) {
                    b = c;
                } else {
                    a = c;
                }
            }
            return c;
        };
        return bisection(I.lower(), I.upper());
    }

    template<typename T>
    std::optional<T> root_finder<T>::regula_falsi(
            const std::function<value_type(value_type)> &func,
            const interval<value_type> &I,
            value_type epsilon) {

        if (greater_than_or_equal(func(I.lower()) * func(I.upper()), 0., epsilon)) {
            return {}; // incorrect endpoints a and b.
        }

        auto regula_falsi = [&](value_type a, value_type b) {
            auto c = a;
            while ((b - a) >= epsilon) {
                c = (a * func(b) - b * func(a)) / (func(b) - func(a));

                if (nearly_zero<value_type>(func(c), epsilon)) {
                    return c;
                } else if (func(c) * func(a) < 0) {
                    b = c;
                } else {
                    a = c;
                }
            }
            return c;
        };
        return regula_falsi(I.lower(), I.upper());
    }

    template<typename T>
    std::optional<T> root_finder<T>::newton_raphson(
            const std::function<value_type(const value_type &)> &func,
            const std::function<value_type(const value_type &)> &derive,
            value_type initial,
            int max_iterations,
            value_type epsilon) {

        auto num_itr = 1;
        auto x = initial;

        auto y = func(x);
        auto dfdx = derive(x);

        while (std::abs(y) >= epsilon && num_itr < max_iterations) {
            if (nearly_zero(dfdx, epsilon)) {
                return {};
            }

            x = x - (y / dfdx);
            y = func(x);
            dfdx = derive(x);

            ++num_itr;
        }
        return x;
    }
}



#endif // ROOT_FINDER_TPP_H_
