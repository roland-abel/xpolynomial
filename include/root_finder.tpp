/// @file root_finder.tpp
/// @brief Root finder class using various numerical methods.
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

#include <limits>
#include <valarray>
#include <numeric>
#include "utils.h"
#include "root_finder.h"

namespace xmath {

    template<typename T>
    std::optional<T> root_finder<T>::bisection(
            const std::function<value_type(value_type)> &func,
            const real_interval<value_type> &I,
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
            const real_interval<value_type> &I,
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


