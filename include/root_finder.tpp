//
// Created by abel on 26.08.2023.
//

#include <limits>
#include <valarray>
#include <numeric>
#include "utils.h"
#include "root_finder.h"

namespace xmath {

    template<typename T>
    root_finder<T>::value_type root_finder<T>::bisection(
            const std::function<value_type(value_type)> &func,
            const interval<value_type> &I,
            value_type tolerance) {

        const auto NaN = std::numeric_limits<T>::quiet_NaN();

        if (func(I.start()) * func(I.end()) >= tolerance) {
            return NaN; // incorrect endpoints a and b.
        }

        auto bisection = [&](value_type a, value_type b) {
            auto c = a;
            while ((b - a) >= tolerance) {
                // Choose the middle point as the new estimation for the root
                c = (b + a) / 2.;

                if (nearly_zero<value_type>(func(c), tolerance)) {
                    return c;
                } else if (func(c) * func(a) < 0) {
                    b = c;
                } else {
                    a = c;
                }
            }
            return c;
        };
        return bisection(I.start(), I.end());
    }

    template<typename T>
    root_finder<T>::value_type root_finder<T>::regula_falsi(
            const std::function<value_type(value_type)> &func,
            const interval<value_type> &I,
            value_type tolerance) {

        const auto NaN = std::numeric_limits<T>::quiet_NaN();

        if (func(I.start()) * func(I.end()) >= tolerance) {
            return NaN; // incorrect endpoints a and b.
        }

        auto regula_falsi = [&](value_type a, value_type b) {
            auto c = a;
            while ((b - a) >= tolerance) {
                c = (a * func(b) - b * func(a)) / (func(b) - func(a));

                if (nearly_zero<value_type>(func(c), tolerance)) {
                    return c;
                } else if (func(c) * func(a) < 0) {
                    b = c;
                } else {
                    a = c;
                }
            }
            return c;
        };
        return regula_falsi(I.start(), I.end());
    }

    template<typename T>
    root_finder<T>::value_type root_finder<T>::newton_raphson(
            const std::function<value_type(value_type)> &func,
            const std::function<value_type(value_type)> &derive,
            value_type initial,
            int max_iterations,
            value_type tolerance) {

        const auto NaN = std::numeric_limits<T>::quiet_NaN();

        auto num_itr = 1;
        auto x = initial;

        auto y = func(x);
        auto dfdx = derive(x);

        while (std::abs(y) >= tolerance && num_itr < max_iterations) {
            if (nearly_zero(dfdx, tolerance)) {
                return NaN;
            }

            x = x - (y / dfdx);
            y = func(x);
            dfdx = derive(x);

            ++num_itr;
        }
        return x;
    }
}


