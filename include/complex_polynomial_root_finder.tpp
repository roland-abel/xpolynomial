/// @file complex_polynomial_root_finder.tpp
///
/// @author Roland Abel
/// @date 20.10.2023

#include <ranges>
#include <numeric>
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

        if (initial_points.size() != p.degree()) {
            return {};
        }

        auto p_norm = p.normalize();
        auto approx_roots = initial_points;

        auto are_almost_roots = [&p](const std::vector<std::complex<T>> z_pts) {
            return std::ranges::all_of(z_pts.cbegin(), z_pts.cend(), [&p](const auto &z) {
                return nearly_zero(std::abs(p(z)));
            });
        };

        size_t iteration = 0;
        while (++iteration < max_iterations && !are_almost_roots(approx_roots)) {
            auto g = complex_polynomial<T>::from_roots(approx_roots);
            auto q = g.derive();

            approx_roots = to_vector(approx_roots | std::views::transform([&p_norm, &q](auto &z) {
                // Weierstrass’ correction
                return z - (p_norm(z) / q(z));
            }));
        }
        return approx_roots;
    }

    template<typename T>
    std::vector<std::complex<T>>
    complex_polynomial_root_finder<T>::aberth_ehrlich_method(
            const complex_polynomial<T> &p,
            const std::vector<std::complex<T>> &initial_points,
            size_t max_iterations) {

        if (initial_points.size() != p.degree()) {
            return {};
        }

        auto p_norm = p.normalize();
        auto p_prim = p_norm.derive();
        auto approx_roots = initial_points;

        auto are_almost_roots = [&p](const std::vector<std::complex<T>> z_pts) {
            return std::ranges::all_of(z_pts.cbegin(), z_pts.cend(), [&p](const auto &z) {
                return nearly_zero(std::abs(p(z)));
            });
        };

        size_t iteration = 0;
        while (++iteration < max_iterations && !are_almost_roots(approx_roots)) {

            approx_roots = to_vector(approx_roots | std::views::transform([&](const auto &z) {
                auto S = [&](const auto &z) {
                    auto r = approx_roots
                             | std::views::filter([&z](const auto &w) { return z != w; })
                             | std::views::transform([&z](const auto &w) { return 1. / (z - w); });

                    return std::accumulate(r.begin(), r.end(), std::complex<T>());
                };
                return z - p_norm(z) * (1. / (p_prim(z) - p_norm(z) * S(z)));
            }));
        }
        return approx_roots;
    }
}