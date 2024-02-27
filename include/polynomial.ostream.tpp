/// @file polynomial.tpp
/// @brief Implementation of operator<< for the polynomial class.
///
/// @author Roland Abel
/// @date October 11, 2023
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

#include <sstream>

namespace xmath {
    template<typename T>
    std::string polynomial<T>::to_string() const {
        std::stringstream os;
        os << *this;
        return os.str();
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &os, const polynomial<T> &p) {
        const auto epsilon = polynomial<T>::epsilon;
        static const auto X = "x";
        static const auto Power = "^";
        static const auto Space = " ";
        static const auto Plus = "+";
        static const auto Minus = "-";

        auto to_string = [epsilon](const typename polynomial<T>::value_type coeff) {
            return nearly_zero(coeff, epsilon) ? 0.0 : coeff;
        };

        if (p.is_constant()) {
            os << to_string(p.leading_coefficient());
        } else {
            auto is_first_term = true;
            for (int k = p.degree(); k >= 0; --k) {
                auto coeff = to_string(p[k]);
                if (nearly_zero<T>(coeff, epsilon)) {
                    continue;
                }

                auto sign = p[k] >= 0 ? Plus : Minus;
                if (is_first_term) {
                    is_first_term = false;
                    if (p[k] < 0) {
                        os << sign;
                    }
                } else {
                    os << Space << sign << Space;
                }

                auto abs_coeff = std::fabs(coeff);
                auto is_one = nearly_equal<T>(abs_coeff, polynomial<T>::spec::one, epsilon);

                if (k == 0 || !is_one) {
                    os << abs_coeff;
                }

                if (k > 0) {
                    os << X;
                    if (k > 1) {
                        os << Power << k;
                    }
                }
            }
        }
        return os;
    }
}
