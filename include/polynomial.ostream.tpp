/// @file polynomial.tpp
/// @brief Implementation of operator<< for the polynomial class.
///
/// @author Roland Abel
/// @date 11.10.2023.

#include <sstream>
#include "polynomial.h"

namespace xmath {

    template<typename T>
    std::string polynomial<T>::to_string() const {
        std::stringstream os;
        os << *this;
        return os.str();
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &os, const polynomial<T> &p) {
        const auto eps = polynomial<T>::tolerance;
        const auto X = "x";
        const auto Power = "^";
        const auto Space = " ";
        const auto Plus = "+";
        const auto Minus = "-";

        auto to_string = [eps](const polynomial<T>::value_type coeff) {
            return nearly_zero(coeff, eps) ? 0.0 : coeff;
        };

        if (p.is_constant()) {
            os << to_string(p.leading_coefficient());
        } else {
            auto is_first_term = true;
            for (int k = p.degree(); k >= 0; --k) {
                auto coeff = to_string(p[k]);
                if (nearly_zero(coeff, eps)) {
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
                auto is_one = nearly_equal(abs_coeff, 1.0, eps);

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

