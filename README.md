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


| Method                                              | Description                                              | Time Complexity           |
|-----------------------------------------------------|----------------------------------------------------------|---------------------------|
| `polynomial()`                                      | Default constructor creating the zero polynomial.        | O(1)                      |
| `polynomial(polynomial&& p)`                        | Move constructor.                                       | O(1)                      |
| `polynomial(std::initializer_list coeffs)`          | Constructor from an initializer list of coefficients.   | O(n)                      |
| `polynomial(const values_type& coeffs)`             | Constructor from a container of coefficients.           | O(n)                      |
| `polynomial(size_type degree)`                      | Constructor to create a polynomial of given degree.     | O(degree)                 |
| `polynomial(const std::ranges::range auto& range)`  | Constructor from a range of coefficients.               | O(n)                      |
| `zero()`                                            | Static method to create the zero polynomial.            | O(1)                      |
| `one()`                                             | Static method to create the polynomial 1.               | O(1)                      |
| `nearly_equal(value_type a, value_type b)`          | Checks if two values are nearly equal.                  | O(1)                      |
| `nearly_zero(value_type a)`                         | Checks if a value is nearly zero.                       | O(1)                      |
| `monomial(size_type degree, value_type coeff)`      | Creates a monomial of given degree and coefficient.     | O(degree)                 |
| `degree()`                                          | Returns the degree of the polynomial.                   | O(1)                      |
| `trim_coefficients()`                               | Removes leading zero coefficients.                      | O(degree)                 |
| `is_zero()`                                         | Checks if the polynomial is the zero polynomial.        | O(degree)                 |
| `is_one()`                                          | Checks if the polynomial is the constant 1.             | O(degree)                 |
| `is_constant()`                                     | Checks if the polynomial is a constant.                 | O(1)                      |
| `is_linear()`                                       | Checks if the polynomial is linear.                     | O(1)                      |
| `is_quadratic()`                                    | Checks if the polynomial is quadratic.                  | O(1)                      |
| `is_cubic()`                                        | Checks if the polynomial is cubic.                      | O(1)                      |
| `is_normalized()`                                   | Checks if the polynomial is normalized.                | O(1)                      |
| `leading_coefficient()`                             | Returns the leading coefficient.                        | O(1)                      |
| `coefficients()`                                    | Returns the coefficients of the polynomial.             | O(1)                      |
| `at(size_type index) const`                         | Returns the coefficient at the given index.             | O(1)                      |
| `at(size_type index)`                               | Returns a reference to the coefficient at the index.    | O(1)                      |
| `operator[](size_type index) const`                 | Returns the coefficient at the given index.             | O(1)                      |
| `operator[](size_type index)`                       | Returns a reference to the coefficient at the index.    |
| `operator==(const polynomial& p) const`             | Checks if two polynomials are equal.                    | O(n)                      |
| `operator()(value_type x) const`                    | Evaluates the polynomial at a given point.              | O(degree)                 |
| `operator!=(const polynomial& p) const`             | Checks if two polynomials are not equal.                | O(n)                      |
| `operator+(value_type scalar) const`                | Adds a scalar to the polynomial.                         | O(n)                      |
| `operator+()`                                       | Unary plus operator.                                    | O(1)                      |
| `operator+=(value_type scalar)`                     | Adds a scalar to the polynomial in place.               | O(n)                      |
| `operator-(value_type scalar) const`                | Subtracts a scalar from the polynomial.                 | O(n)                      |
| `operator-()`                                       | Unary minus operator.                                   | O(n)                      |
| `operator-=(value_type scalar)`                     | Subtracts a scalar from the polynomial in place.        | O(n)                      |
| `operator*(value_type scalar) const`                | Multiplies the polynomial by a scalar.                  | O(n)                      |
| `operator*=(value_type scalar)`                     | Multiplies the polynomial by a scalar in place.         | O(n)                      |
| `operator/(value_type scalar) const`                | Divides the polynomial by a scalar.                     | O(n)                      |
| `operator/=(value_type scalar)`                     | Divides the polynomial by a scalar in place.            | O(n)                      |
| `operator+(const polynomial& p) const`              | Adds another polynomial.                                | O(max(n, m))              |
| `operator-(const polynomial& p) const`              | Subtracts another polynomial.                           | O(max(n, m))              |
| `operator*(const polynomial& p) const`              | Multiplies by another polynomial.                       | O(n * m)                  |
| `operator+=(const polynomial& p)`                   | Adds another polynomial in place.                       | O(max(n, m))              |
| `operator-=(const polynomial& p)`                   | Subtracts another polynomial in place.                  | O(max(n, m))              |
| `operator*=(const polynomial& p)`                   | Multiplies by another polynomial in place.              | O(n * m)                  |
| `operator/=(const polynomial& p)`                   | Divides by another polynomial in place.                 | O((n - m) * m)            |
| `operator/(const polynomial& p) const`              | Divides by another polynomial.                          | O((n - m) * m)            |
| `operator%(const polynomial& p) const`              | Computes the remainder when divided by another polynomial.| O((n - m) * m)            |
| `operator%=(const polynomial& p)`                   | Computes the remainder when divided in place by another polynomial.| O((n - m) * m)            |
| `compose(const polynomial& q) const`                | Composes the polynomial with another.                   | O(degree^2)               |
| `pow(unsigned int exponent) const`                  | Raises the polynomial to a given power.                 | O(log(exponent))          |
| `derive() const`                                    | Computes the derivative of the polynomial.              | O(degree)                 |
| `integrate() const`                                 | Computes the indefinite integral of the polynomial.     | O(degree)                 |
| `evaluate(value_type x) const`                      | Evaluates the polynomial at a given point.              | O(degree)                 |
| `operator+(value_type scalar, const polynomial& p)` | Adds a scalar to the polynomial.                         | O(n)                      |
| `operator*(value_type scalar, const polynomial& p)` | Multiplies the polynomial by a scalar.                  | O(n)                      |
| `from_roots(const values_type& roots)`              | Creates a polynomial from its roots.                    | O(roots.size()^2)         |
| `normalize() const`                                 | Normalizes the coefficients of the polynomial.         | O(n)                      |
| `has_root(const value_type& value) const`           | Checks if a given value is a root of the polynomial.   | O(degree)                 |
| `has_roots(const values_type& values) const`        | Checks if a set of values are roots of the polynomial. | O(degree * values.size()) |
| `divide(const polynomial& divisor) const`           | Divides the polynomial by another polynomial.          | O((n - m) * m)            |


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

| Method                                                                                                 | Description                                                                                 | Complexity                                 |
|--------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------|--------------------------------------------|
| `sign_changes(const std::vector<T>& sequence, T epsilon = 1e-5)`                                       | Counts the number of sign changes in a sequence.                                          | O(n) (where n is the size of the sequence) |
| `quadratic_roots(const polynomial<T>& p)`                                                              | Finds the roots of a quadratic polynomial.                                                 | O(1)                                       |
| `has_cubic_normal_form(const polynomial<T>& p)`                                                        | Checks if a cubic polynomial is in normal form.                                            | O(1)                                       |
| `cubic_roots(const polynomial<T>& p)`                                                                  | Finds the roots of a cubic polynomial.                                                     | O(1)                                       |
| `cubic_normal_form_roots(const polynomial<T>& p)`                                                      | Finds the roots of a cubic polynomial in normal form.                                      | O(1)                                       |
| `newton_raphson(const polynomial<T>& p, value_type initial, int max_iterations, value_type tolerance)` | Applies the Newton-Raphson method to find a root of a polynomial.                           | O(max_iterations)                          |
| `sign_changes(const polynomial<T>& p)`                                                                 | Counts the number of sign changes in the Sturm sequence of a polynomial.                     | O(n^2)                                     |
| `cauchy_bounds(const polynomial<T>& p)`                                                                | Computes the Cauchy bounds for a polynomial.                                               | O(n)                                       |
| `lagrange_bounds(const polynomial<T>& p)`                                                              | Computes the Lagrange bounds for a polynomial.                                             | O(n)                                       |
| `sturm_sequence(const polynomial<T>& p)`                                                               | Generates the Sturm sequence of a polynomial.                                              | O(n^2)                                     |
| `sign_variations(const polynomial_sequence& seq, const value_type& x)`                                 | Counts the number of sign variations in a Sturm sequence at a given point.                 | O(seq.size())                              |
| `number_distinct_roots(const polynomial<T>& p, const real_interval<T>& I)`                             | Counts the number of distinct roots of a polynomial in a given interval using Sturm's method.| O(seq.size())                              |
| `number_distinct_roots(const polynomial<T>& p)`                                                        | Counts the number of distinct roots of a polynomial using Sturm's method.                   | O(seq.size())                              |
| `root_isolation(const polynomial<T>& p)`                                                               | Performs root isolation using Sturm's method.                                              | O(seq.size() * log(a))                     |
| `find_roots(const polynomial<T>& p, const T precision)`                                                | Finds the roots of a polynomial with given precision using root isolation and bisection.   |                                            |



### Square-free Decomposition

The `square_free_decomposition<>` class identifies whether a polynomial is square-free and employs Yun's algorithm
for square-free decomposition. Additionally, it provides a method to reconstruct the original polynomial
from its square-free decomposition.

| Method                           | Description                                                                           | Time Complexity                                            |
|----------------------------------|---------------------------------------------------------------------------------------|------------------------------------------------------------|
| `is_square_free`                 | Checks if the given polynomial is square-free.                                        | $O(N^2)$, where N is the degree of the polynomial.         |
| `from_square_free_decomposition` | Reconstructs the original polynomial from its square-free decomposition.              | O(N^2), where N is the size of the decomposition sequence. |
| `yun_algorithm`                  | Computes the square-free decomposition of the given polynomial using Yun's algorithm. | O(N^3), where N is the degree of the polynomial.           |

### Polynomial Interpolation

The `polynomial_interpolation<>` class provides functionality for Lagrange polynomial interpolation,
a method to construct a polynomial that passes through a given set of data points.

| Method                                                                                                                        | Description                                                                                         | Time Complexity                                              |
|-------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------|--------------------------------------------------------------|
| `lagrange_basis(const std::vector<T> &xs)`      | Computes the Lagrange basis polynomials for the given set of interpolation points.                  | Time: O(N^2), where N is the number of interpolation points. |
| `lagrange_interpolation(const values_type &x_values, const values_type &y_values)` | Performs Lagrange interpolation to find the polynomial that passes through the given set of points. | Time: O(N^2), where N is the number of interpolation points. |


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
