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

    template<typename T>
    std::vector<std::complex<T>> complex_polynomial_root_finder<T>::durand_kerner_method(
            const complex_polynomial<T> &p,
            const std::vector<std::complex<T>> &initial_points,
            size_t max_iterations) {

        auto are_almost_roots = [&p](const std::vector<std::complex<T>> z_pts) {
            return std::ranges::all_of(z_pts.cbegin(), z_pts.cend(), [&p](const auto &z) {
                return nearly_zero(std::abs(p(z)));
            });
        };

        if (initial_points.size() != p.degree()) {
            return {};
        }

        size_t i = 0;
        auto p_norm = p.normalize();
        auto roots = initial_points;

        while (++i < max_iterations && !are_almost_roots(roots)) {
            auto g = complex_polynomial<T>::from_roots(roots);
            auto q = g.derive();

            // Weierstrass’ correction
            std::for_each(roots.begin(), roots.end(), [&p_norm, &q](auto &z) {
                z = z - (p_norm(z) / q(z));
            });
        }
        return roots;
    }
}