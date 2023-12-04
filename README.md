# xpolynomial

C++ template project provides a flexible way to work with polynomials and calculate their real roots.

## Overview

This C++ template project provides a flexible way of working with polynomials and the 
calculation of their real roots.

## Features

### Polynomial

The `polynomial<T>` class allows you to create and manipulate polynomials in a single variable. 
This template class is designed to handle polynomials with coefficients of various numeric 
data types and provides a wide range of functions for working with polynomials, including 
basic operations, evaluation, normalization, and more.

```c++
#include <numeric>
#include <vector>
#include <iostream>
#include "polynomial.h"

using namespace std;
using namespace xmath;

namespace {
    auto X = polynomial<double>::monomial(1, 1.0);
}

int main() {
    // Create a 4th degree polynomial
    auto p = 3 * X.pow(4) - 2.5 * X.pow(3) + X.pow(2) - X + 1;

    // Evaluate the polynomial for a range of values
    vector<double> values(10);
    iota(values.begin(), values.end(), 1);

    for (auto x: values) {
        cout << "p(" << x << ") = " << p(x) << endl;
    }
    return 0;
}
```

### Euclidean Algorithm

The primary functions of the `euclidean_algorithm<>` template are to compute the greatest common divisor (gcd)
of two polynomials and the extended Euclidean algorithm for polynomials.

```c++
#include <iostream>
#include "euclidean_algorithm.h"

using namespace std;
using namespace xmath;

namespace {
    auto X = polynomial<double>::monomial(1, 1.0);
    using Euclidean = euclidean_algorithm<double>;
}

int main() {
    auto p = X.pow(4) - 2 * X.pow(3) - 6 * X.pow(2) + 12 * X + 15;
    auto q = X.pow(3) + X.pow(2) - 4 * X - 4;
    
    // s, t, g such that g = gcd(p, q) = s*p + t*q
    auto [s, t, g] = Euclidean::extended_euclidean(p, q); 

    cout << "p = " << p.to_string() << endl
         << "q = " << q.to_string() << endl << endl
         << "g = gcd(p, q) = " << g << endl
         << "s = " << s << endl
         << "t = " << t << endl
         << "g = s*p + t*q = " << s * p + t * q << endl;

    return 0;
}
```

### Chebyshev Polynomial

The `chebyshev_polynomial<>` template class provides a set of methods for working with Chebyshev polynomials
of the first kind. These polynomials are an integral part of several mathematical applications, and the 
class facilitates their generation, evaluation, and interpolation.

```c++
#include "chebyshev_polynomial.h"
#include "real_polynomial_root_finder.h"

using namespace std;
using namespace xmath;

namespace {
    using ChebyshevPolynomial = chebyshev_polynomial<double>;
    using RootFinder = real_polynomial_root_finder<double>;
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
```

### Root Finding Algorithm

The `real_polynomial_root_finder<>` and `complex_polynomial_root_finder<>` classes are utility classes specifically
designed for finding roots of polynomials with real and complex coefficients, respectively.
These classes provide various numerical methods to obtain accurate root approximations.

```c++
#include "real_polynomial_root_finder.h"

using namespace std;
using namespace xmath;

namespace {
    using RootFinder = real_polynomial_root_finder<double>;
    auto X = polynomial<double>::monomial(1, 1.0);
}

int main() {
    auto p = (X + 3.).pow(3) * (X - 1.) * (X.pow(2) + X - 2).pow(3);
    auto [roots, multiplicities] = RootFinder::find_roots(p);

    cout << "Polynomial: " << p << endl << endl;
    for (int k = 0; k < roots.size(); ++k) {
        cout << "Root: r[" << k << "] = " << roots[k] << ", Multiplicity: "
             << multiplicities[k] << endl;
    }
    return 0;
}
```

### Square-free Decomposition

The `square_free_decomposition<>` class determines whether a polynomial is square-free and uses Yun's algorithm
for square-free decomposition. In addition, this class provides a method for constructing the original polynomial
from its square-free decomposition.

```c++
#include <iostream>
#include "square_free_decomposition.h"

using namespace std;
using namespace xmath;

namespace {
    auto X = polynomial<double>::monomial(1, 1.0);
    using SquareFree = square_free_decomposition<double>;
}

int main() {
    auto p = X * (X - 4) * (X + 3.).pow(2) * (X.pow(2) + X).pow(3) * (X.pow(2) + X).pow(4);
    cout << "p = " << p << endl;

    auto square_free_seq = SquareFree::yun_algorithm(p).value();
    for (int k = 0; k < square_free_seq.size(); ++k) {
        cout << "q" << k << " = " << square_free_seq[k] << endl;
    }
    return 0;
}
```

### Polynomial Interpolation

The `polynomial_interpolation<>` class provides functionality for Lagrange polynomial interpolation,
a method to construct a polynomial that passes through a given set of data points.

```c++
#include <iostream>
#include "polynomial_interpolation.h"

using namespace std;
using namespace xmath;

namespace {
    using Polynomial = polynomial<double>;
    using Interpolation = polynomial_interpolation<double>;
}

int main() {
    const auto x_values = {-3., -2., -1., 0., 1., 2., 3.};
    const auto y_values = {-2.4, -1.5, 1.1, 2.5, -3.6, -1.25, -2.1};

    const auto p = Interpolation::lagrange_interpolation(x_values, y_values).value();
    cout << "p(x) = " << p << endl;

    for (auto x: x_values) {
        cout << "p(" << x << ") = " << p(x) << endl;
    }
    return 0;
}
```

## Author

Roland Abel

## License

This software is available under the following licenses:

See [MIT License (MIT)](LICENSE)
