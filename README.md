# xpolynomial

xpolynomial is a header-only C++23 template library for creating and manipulating polynomials
and calculating their roots.

## Overview

This C++ template library provides a flexible way of working with polynomials in a single
variable. It supports coefficients of arbitrary numeric types and offers algorithms built
around polynomials:

- **Polynomials** the `polynomial<T>` and `complex_polynomial<T>` template classes for real
  and complex-valued polynomials.
- **Root finding** for real and complex polynomials, including root multiplicities.
- **Further algorithms** such as square-free decomposition (Yun), polynomial interpolation,
  Chebyshev and Legendre polynomials, and extended Euclidean algorithms for polynomials.
- **A parser** to construct polynomials from string expressions such as `"2*X^3 - 3*X + 1"`.

The library is:

- **Header-only** - nothing has to be compiled or linked; installing it just copies headers.
- **Dependency-free** - it requires only the C++ standard library.
- **C++23** - a C++23-capable compiler is required.

The complete API is available through a single umbrella header:

```c++
#include <xpolynomial.h>
```

Individual headers such as `#include <polynomial.h>` can also be included.

### API Overview

| Header | Provides |
|---|---|
| `<polynomial.h>` | `polynomial<T>` with arithmetic, evaluation, calculus and polynomial division |
| `<complex_polynomial.h>` | `complex_polynomial<T>`, `real_polynomial<T>`, utils like `separate()` |
| `<interval.h>` | `interval<T>` with open/closed boundary handling and bisection |
| `<real_polynomial_root_finder.h>` | real root finding (quadratic/cubic formulas, Newton, Sturm isolation) |
| `<complex_polynomial_root_finder.h>` | Durand-Kerner, Aberth-Ehrlich and roots of unity |
| `<root_finder.h>` | generic bisection / regula falsi / Newton-Raphson on arbitrary callables |
| `<square_free_decomposition.h>` | Yun's square-free decomposition |
| `<polynomial_interpolation.h>` | Lagrange interpolation |
| `<chebyshev_polynomial.h>` | Chebyshev polynomials, nodes, Clenshaw, Gauss quadrature |
| `<legendre_polynomial.h>` | Legendre polynomials |
| `<euclidean_algorithm.h>` | (extended) Euclidean algorithm for polynomials |
| `<polynomial_parser.h>` | `parser::parse_polynomial` to build polynomials from strings |
| `<utils.h>` | numeric helpers (`nearly_*`, comparisons) |

A note on error handling: functions returning `std::optional` or `std::expected` report
failures without throwing (e.g. `parse_polynomial`, `find_roots`, root finders). Always
check the result before calling `.value()`, which throws on error.

## Integration

The library is built with CMake and is consumed either via `FetchContent` or by installing it
and using `find_package`.

### Requirements

- CMake 3.22 or newer
- A compiler with C++23 support (e.g. GCC 14 or newer, Clang 17 or newer, MSVC 19.3x)

### FetchContent

```cmake
include(FetchContent)

FetchContent_Declare(xpolynomial
    GIT_REPOSITORY https://github.com/roland-abel/xpolynomial.git
    GIT_TAG main)
FetchContent_MakeAvailable(xpolynomial)

add_executable(demo demo.cpp)
target_link_libraries(demo PRIVATE xpolynomial::xpolynomial)
```

### Install and find_package

Build and install the library:

```sh
cmake -S xpolynomial -B build
cmake --build build
cmake --install build --prefix /path/to/install
```

Use it in your project:

```cmake
find_package(xpolynomial CONFIG REQUIRED)

add_executable(demo demo.cpp)
target_link_libraries(demo PRIVATE xpolynomial::xpolynomial)
```

The `xpolynomial::xpolynomial` target provides all headers, the C++23 language standard and
the umbrella header `<xpolynomial.h>`.

### Building the tests

The project requires GoogleTest for its test suite. Point `CMAKE_PREFIX_PATH` to an
installation of GoogleTest and run the tests afterwards. Tests are enabled by default and can
be turned off with `-DBUILD_TESTING=OFF`; examples with `-DXPOLYNOMIAL_BUILD_EXAMPLES=OFF`.

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/gtest
cmake --build build
ctest --test-dir build --output-on-failure
```

### Example

The following program parses a polynomial from a string and computes its real roots via the
umbrella header:

```c++
#include <xpolynomial.h>
#include <iostream>

using namespace std;
using namespace xmath;

int main() {
    auto poly = parser::parse_polynomial("X^3 - 3*X + 1").value();
    cout << "p(x) = " << poly << endl;

    auto [roots, multiplicities] = real_polynomial_root_finder<double>::find_roots(poly);
    for (std::size_t k = 0; k < roots.size(); ++k) {
        cout << "root r[" << k << "] = " << roots[k]
             << ", multiplicity = " << multiplicities[k] << endl;
    }
    return 0;
}
```

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
#include <xpolynomial.h>

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

#### Compile-time evaluation

The methods of the real `polynomial<T>` class (construction, arithmetic, evaluation, calculus,
division and comparison) are `constexpr` and can be evaluated at compile time. A polynomial
object itself cannot have static storage duration, because the coefficients live in a
`std::vector` whose allocation would have to persist; constant evaluation is therefore only
possible within a constant expression, for example inside a `constexpr` function or lambda.

```c++
constexpr bool check() {
    const auto X = polynomial<double>::monomial(1, 1.0);
    const auto p = (X - 1.0).pow(2);          // (X - 1)^2
    return p(3.0) == 4.0                      // evaluation
        && p.derive() == 2.0 * X - 2.0        // calculus
        && p.is_root(1.0);                    // root check
}

static_assert(check());
```

`complex_polynomial<T>` is not `constexpr`, because `std::abs(std::complex<T>)` is not a
constant expression.

### Interval

The `interval<T>` class represents a real interval with configurable open/closed boundary
types. It supports the usual queries (`is_open`, `is_degenerate`, `is_empty`, ...), bisection
and a linear transform to another interval.

```c++
#include <iostream>
#include <xpolynomial.h>

using namespace std;
using namespace xmath;

int main() {
    auto I = interval(0., 1.);
    cout << "I = [" << I.lower() << ", " << I.upper() << "]" << endl;
    cout << "length = " << I.length() << endl;
    cout << "degenerate? " << boolalpha << I.is_degenerate() << endl;

    auto [left, right] = I.bisect();
    cout << "[" << left.lower() << ", " << left.upper() << "] "
         << "[" << right.lower() << ", " << right.upper() << "]" << endl;
    return 0;
}
```

### Euclidean Algorithm

The primary functions of the `euclidean_algorithm<>` template are to compute the greatest common divisor (gcd)
of two polynomials and the extended Euclidean algorithm for polynomials.

```c++
#include <iostream>
#include <xpolynomial.h>

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
#include <xpolynomial.h>

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

### Legendre Polynomial

The `legendre_polynomial<>` class computes the Legendre polynomials `P_n(x)` via a recurrence
relation. Results are cached, so repeated calls are cheap.

```c++
#include <iostream>
#include <xpolynomial.h>

using namespace std;
using namespace xmath;

int main() {
    for (int n = 0; n <= 4; ++n) {
        cout << "P_" << n << ": " << legendre_polynomial<double>::create(n) << endl;
    }
    return 0;
}
```

### Root Finding Algorithm

The `real_polynomial_root_finder<>` and `complex_polynomial_root_finder<>` classes are utility classes specifically
designed for finding roots of polynomials with real and complex coefficients, respectively.
These classes provide various numerical methods to obtain root approximations.

```c++
#include <xpolynomial.h>

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

### Complex Polynomial

The `complex_polynomial<T>` class is a `polynomial<std::complex<T>>` with specialized
coefficient handling. It interoperates with real polynomials, and the `separate()` function
splits a complex polynomial into its real and imaginary parts.

```c++
#include <complex>
#include <iostream>
#include <xpolynomial.h>

using namespace std;
using namespace xmath;

namespace {
    using ComplexPolynomial = complex_polynomial<double>;
    using RealPolynomial = polynomial<double>;

    constexpr auto i = std::complex(0., 1.);
    auto Z = ComplexPolynomial::monomial(1, 1.0);
}

int main() {
    auto p = ComplexPolynomial({1. - i, 2. + 3. * i, -1.}); // (2+3i)*x^2 + (1-i)*x - 1
    cout << "p = " << p << endl;
    cout << "p(1) = " << p(1.) << endl;

    auto [real, imag] = separate(p);
    cout << "Re(p) = " << real << ", Im(p) = " << imag << endl;
    return 0;
}
```

### Complex Polynomial Root Finding

The `complex_polynomial_root_finder<>` class finds the complex roots of a polynomial with the
Durand-Kerner and Aberth-Ehrlich methods, and provides the computed roots of unity.

```c++
#include <complex>
#include <iostream>
#include <xpolynomial.h>

using namespace std;
using namespace xmath;

int main() {
    using RootFinder = complex_polynomial_root_finder<double>;

    auto roots = RootFinder::nth_roots_of_unity(5);
    cout << "5-th roots of unity:" << endl;
    for (auto z: roots) {
        cout << "(" << z.real() << ", " << z.imag() << ")" << endl;
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
#include <xpolynomial.h>

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
#include <xpolynomial.h>

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

### Polynomial Parser

The `parser::parse_polynomial` function constructs a polynomial from a string expression. It
supports `+`, `-`, `*`, `/`, `^`, parentheses, the variable `X` and unary signs. On failure it
returns a `std::expected` with a descriptive `error_t`, so no exception is thrown.

```c++
#include <iostream>
#include <xpolynomial.h>

using namespace std;
using namespace xmath;

int main() {
    auto result = parser::parse_polynomial("2*X^3 - 3*X + 1");
    if (result.has_value()) {
        cout << "p(x) = " << result.value() << endl;
    } else {
        cerr << "parse error" << endl;
    }

    auto invalid = parser::parse_polynomial("2*X^^3");
    if (!invalid.has_value()) {
        cerr << "expected error for invalid expression" << endl;
    }
    return 0;
}
```

## Author

Roland Abel

## Changelog

- **2026** - Parser, complex/real root finders, Chebyshev, Legendre and interval support;
  header-only C++23 rewrite.
- **2024** - Initial release with polynomial arithmetic, Euclidean algorithm and root finding.

## License

This software is available under the following licenses:

See [MIT License (MIT)](LICENSE)
