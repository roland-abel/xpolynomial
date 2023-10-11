# xpolynomial
This C++ template project provides a flexible way to work with polynomials and calculate their real roots.

## Features

### Polynomial Class
The `polynomial<T>` class allows you to create and manipulate polynomials in a single variable. Here are some of the available functions:

- `polynomial<T>()`: Constructor the zero polynomial.
- `polynomial<T>(std::initializer_list<coeff_type> coeffs)`: Constructor that creates a polynomial with the given coefficients.
- `value_type evaluate(value_type x) const`: Computes the value of the polynomial for a given variable `x`.
- `polynomial<T> derive() const`: Computes the derivative of the polynomial.
- `polynomial<T> operator+(const polynomial<T> &p) const`: Adds two polynomials together.
- `polynomial<T> operator-(const polynomial<T> &p) const`: Subtracts one polynomial from another.
- `polynomial<T> operator*(const polynomial<T> &p) const`: Multiplies two polynomials.
- `std::string to_string() const`: Returns the polynomial as a string.

### Root Finder Class
The `real_polynomial_root_finder<T>` class finds all real roots of a given polynomial using the Bisection method. Here are some of the available functions:

- `int number_distinct_roots(const polynomial<T> &p, const interval<T> &I)`: Calculates the number of distinct real roots for a given Interval.
- `std::tuple<roots_type, multiplicities_type> find_roots(const polynomial<T> &p)`: Find all the real roots with the corresponding multiplicities.

### Sturm Method for Counting Roots

The number of real roots of a polynomial can be determined using the Sturm method. This method involves creating a sequence of polynomials, known as Sturm sequences, and counting the sign changes in this sequence. The number of sign changes corresponds to the number of distinct real roots of the polynomial.


## Usage

To use this project, you can include the `polynomial<double>` and `real_polynomial_root_finder<double>` classes in your C++ code and use them as follows:

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

## License

This software is available under the following licenses:

See [MIT License (MIT)](LICENSE)
