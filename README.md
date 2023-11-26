# xpolynomial

C++ template project provides a flexible way to work with polynomials and calculate their real roots.

## Overview

This C++ template library provides functionality for working with polynomials.
It includes basic operations such as addition, subtraction, multiplication,
and division, as well as polynomial evaluation.

## Features

### Polynomial

The `polynomial<T>` class allows you to create and manipulate polynomials in a single variable. This template class
designed to handle polynomials with coefficients of various numerical data types and provides a wide range of
functionalities for working with polynomials, including basic operations, evaluation, normalization, and more.

### Euclidean Algorithm

The primary functionalities of the `euclidean_algorithm<>` template is computing the greatest common divisor (gcd)
of two polynomials and the extended Euclidean algorithm for polynomials.

### Chebyshev Polynomial

The `chebyshev_polynomial<>` template class provides a set of methods for working with Chebyshev polynomials
of the first kind. These polynomials are integral in various mathematical applications, and the class facilitates
their generation, evaluation, and interpolation.

### Root Finding Algorithm

The `real_polynomial_root_finder<>` and `complex_polynomial_root_finder<>` classes are utility classes specifically
designed for finding roots of polynomials with real and complex coefficients, respectively.
These classes provide various numerical methods to achieve accurate root approximations.

Sturm's method plays a crucial role in identifying, counting, and isolating roots within intervals.
Its effectiveness lies in constructing a sequence of polynomials whose sign variations provide information
about the distribution of roots. The algorithm's complexity is polynomial in nature, making it suitable for a wide
range of polynomial degrees.

### Square-free Decomposition

The `square_free_decomposition<>` class identifies whether a polynomial is square-free and employs Yun's algorithm
for square-free decomposition. Additionally, it provides a method to reconstruct the original polynomial
from its square-free decomposition.

### Polynomial Interpolation

The `polynomial_interpolation<>` class provides functionality for Lagrange polynomial interpolation,
a method to construct a polynomial that passes through a given set of data points.


## Usage

To use this project, you can include the `polynomial<double>` and `real_polynomial_root_finder<double>` classes in your
C++ code and use them as follows:

```cpp
#include "real_polynomial_root_finder.h"

using namespace std;
using namespace xmath;

namespace {
    using RootFinder = real_polynomial_root_finder<double>;
    auto X = polynomial<double>::monomial(1, 1.0);
}

int main() {
    // Create a polynomial
    auto p = (X + 1 / 3.).pow(3) * (X - 1.5) * (X.pow(2) + X - 2).pow(2);
    auto [roots, multiplicities] = RootFinder::find_roots(p);

    cout << "Polynomial: " << p << endl << endl;
    for (int k = 0; k < roots.size(); ++k) {
        cout << "Root: r[" << k << "] = " << roots[k]
             << ", Multiplicity: " << multiplicities[k] << endl;
    }
    return 0;
}
```

## Usage

todo:

## Author

Roland Abel

## License

This software is available under the following licenses:

See [MIT License (MIT)](LICENSE)
