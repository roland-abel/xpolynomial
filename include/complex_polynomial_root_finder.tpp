/// @file complex_polynomial_root_finder.tpp
///
/// @author Roland Abel
/// @date 20.10.2023

#include <ranges>
#include <algorithm>
#include <vector>
#include "complex_polynomial_root_finder.h"

namespace xmath {

    template<typename T>
    std::vector<std::complex<T>> complex_polynomial_root_finder<T>::nth_roots_of_unity(int n) {
        auto I = std::complex<T>(0., 1.);
        auto get_root = [=](int k) {
            return std::complex<T>(
                    std::cos((2. / (T) n) * k * std::numbers::pi),
                    std::sin((2. / (T) n) * k * std::numbers::pi));
        };

        std::vector<std::complex<T>> roots{};
        for (auto r: std::ranges::iota_view{0, n} | std::views::transform(get_root)) {
            roots.push_back(r);
        }
        return roots;
    }
}