/// @file complex_polynomial_root_finder.tpp
/// @brief Root finder class for polynomials with complex coefficients.
///
/// @author Roland Abel
/// @date October 20, 2023
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

#include <ranges>
#include <numeric>
#include <algorithm>
#include <vector>
#include "complex_polynomial_root_finder.h"

namespace xmath {
    using std::views::filter;
    using std::views::transform;
    using std::ranges::all_of;

    template<typename T>
    std::vector<std::complex<T> > complex_polynomial_root_finder<T>::nth_roots_of_unity(int n) {
        auto I = std::complex<T>(0., 1.);
        auto get_root = [=](int k) {
            return std::complex<T>(
                std::cos(2. / static_cast<T>(n) * k * std::numbers::pi),
                std::sin(2. / static_cast<T>(n) * k * std::numbers::pi));
        };

        std::vector<std::complex<T> > roots{};
        for (auto r: std::ranges::iota_view{0, n} | transform(get_root)) {
            roots.push_back(r);
        }
        return roots;
    }

    template<typename T>
    std::vector<std::complex<T> > complex_polynomial_root_finder<T>::durand_kerner_method(
        const complex_polynomial<T> &p,
        const std::vector<std::complex<T> > &initial_points,
        size_t max_iterations) {
        if (initial_points.size() != p.degree()) {
            return {};
        }

        auto p_norm = p.normalize();
        auto approx_roots = initial_points;

        auto are_almost_roots = [&p](const std::vector<std::complex<T> > z_pts) {
            return all_of(z_pts.cbegin(), z_pts.cend(), [&p](const auto &z) {
                return nearly_zero(std::abs(p(z)));
            });
        };

        size_t iteration = 0;
        while (++iteration < max_iterations && !are_almost_roots(approx_roots)) {
            auto g = complex_polynomial<T>::from_roots(approx_roots);
            auto q = g.derive();

            approx_roots = to_vector(approx_roots | transform([&p_norm, &q](auto &z) {
                // Weierstrass’ correction
                return z - (p_norm(z) / q(z));
            }));
        }
        return approx_roots;
    }

    template<typename T>
    std::vector<std::complex<T> >
    complex_polynomial_root_finder<T>::aberth_ehrlich_method(
        const complex_polynomial<T> &p,
        const std::vector<std::complex<T> > &initial_points,
        size_t max_iterations) {
        if (initial_points.size() != p.degree()) {
            return {};
        }

        auto p_norm = p.normalize();
        auto p_prim = p_norm.derive();
        auto approx_roots = initial_points;

        auto are_almost_roots = [&p](const std::vector<std::complex<T> > z_pts) {
            return all_of(z_pts.cbegin(), z_pts.cend(), [&p](const auto &z) {
                return nearly_zero(std::abs(p(z)));
            });
        };

        size_t iteration = 0;
        while (++iteration < max_iterations && !are_almost_roots(approx_roots)) {
            approx_roots = to_vector(approx_roots | transform([&](const auto &z) {
                auto S = [&](const auto &s) {
                    auto r = approx_roots
                             | filter([&s](const auto &w) { return s != w; })
                             | transform([&s](const auto &w) { return 1. / (s - w); });

                    return std::accumulate(r.begin(), r.end(), std::complex<T>());
                };
                return z - p_norm(z) * (1. / (p_prim(z) - p_norm(z) * S(z)));
            }));
        }
        return approx_roots;
    }
}
