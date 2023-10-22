/// @file complex_polynomial_root_finder.tpp
///
/// @author Roland Abel
/// @date 20.10.2023

#include "complex_polynomial.h"

namespace xmath {

    template<typename T>
    using real_values_type = polynomial<T, polynomial_specification<T>>::values_type;

    template<typename T>
    auto separate(const complex_polynomial<T> &polynomial) -> decltype(auto) {
        auto real_coeffs = real_values_type<T>();
        auto imag_coeffs = real_values_type<T>();

        for (auto z: polynomial.coefficients()) {
            real_coeffs.push_back(z.real());
            imag_coeffs.push_back(z.imag());
        }
        return std::pair<real_polynomial < T>, real_polynomial < T >> (real_coeffs, imag_coeffs);
    }
}