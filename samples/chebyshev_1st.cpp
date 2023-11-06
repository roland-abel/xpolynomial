/// @file chebyshev_1st.cpp
/// @brief Example for Chebyshev polynomials of the 1st kind.
///
/// @author Roland Abel
/// @date 14.10.2023

#include "chebyshev_polynomial.h"

using namespace std;
using namespace xmath;

namespace {
    using ChebyshevPolynomial = chebyshev_polynomial<float>;
    auto X = polynomial<double>::monomial(1, 1.0);
}

int main() {
    const auto max_order = 10;

    cout << "Chebyshev Polynomials of 1st kind: " << endl;
    for (int n = 0; n <= max_order; ++n) {
        auto T_n = ChebyshevPolynomial::create_1st_kind(n);
        cout << "T_" << n << ": " << T_n << endl;
    }
    return 0;
}