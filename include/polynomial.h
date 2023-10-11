/// @file polynomial.h
/// @brief This file contains the declaration of the polynomial class and related functions.
///
/// @author Roland Abel
/// @date 19.08.2023

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <utility>
#include <iostream>
#include <vector>
#include "matrix.h"

namespace xmath {

    /// @brief Represents a polynomial with coefficients from template parameter type T.
    template<typename T>
    class polynomial {
    public:
        using coeffs_type = std::vector<T>;
        using values_type = std::vector<T>;
        using coeff_type = coeffs_type::value_type;
        using value_type = values_type::value_type;
        using size_type = coeffs_type::size_type;

        static constexpr value_type tolerance = 1e-5;

    public:
        /// @brief Default constructor. Creates a zero polynomial.
        explicit polynomial<T>();

        /// @brief Constructor that takes a list of coefficients.
        /// @param coeffs The coefficients of the polynomial in descending order.
        polynomial(std::initializer_list<coeff_type> coeffs);

        /// @brief Constructor that takes a list of coefficients.
        /// @param coeffs The coefficients of the polynomial in descending order.
        explicit polynomial(coeffs_type coeffs);

        /// @brief Copy constructor.
        /// @param p The polynomial to be copied.
        polynomial(const polynomial<T> &p) = default;

        /// @brief Destructor.
        ~polynomial() = default;

        /// @brief Checks if the polynomial is the zero polynomial.
        /// @return True if the polynomial is zero, false otherwise.
        bool is_zero() const;

        /// @brief Checks if the polynomial is the 1 polynomial.
        /// @return True if the polynomial is 1, false otherwise.
        bool is_one() const;

        /// @brief Checks if the polynomial is a constant (degree 0) polynomial.
        /// @return True if the polynomial is a constant, false otherwise.
        bool is_constant() const;

        /// @brief Checks if the polynomial is a linear (degree 1) polynomial.
        /// @return True if the polynomial is linear, false otherwise.
        bool is_linear() const;

        /// @brief Checks if the polynomial is a quadratic (degree 2) polynomial.
        /// @return True if the polynomial is quadratic, false otherwise.
        bool is_quadratic() const;

        /// @brief
        /// @return
        bool is_cubic() const;

        /// @brief Returns the degree of the polynomial.
        /// @return The degree of the polynomial.
        size_type degree() const;

        /// @brief Returns the leading coefficient of the polynomial.
        /// @return The leading coefficient.
        coeff_type leading_coefficient() const;

        /// @brief Returns the coefficients of th polynomial.
        /// @return The coefficients.
        const coeffs_type &coefficients() const;

        /// @brief Evaluates the polynomial at a given value.
        /// @param x The value at which to evaluate the polynomial.
        /// @return The result of the evaluation.
        value_type evaluate(value_type x) const;

        /// @brief Creates a zero polynomial.
        /// @return A zero polynomial.
        static polynomial<T> zero();

        /// @brief Creates a polynomial representing the constant 1.
        /// @return A polynomial representing 1.
        static polynomial<T> one();

        /// @brief Creates a monomial of the given degree.
        /// @param degree The degree of the monomial.
        /// @param coeff The coefficient.
        /// @return
        static polynomial<T> monomial(size_type degree, coeff_type coeff = (coeff_type) 1);

        /// @brief Constructs a polynomial from with the given roots.
        /// @param roots The vector of roots of the polynomial.
        /// @return The polynomial with the given roots.
        static polynomial<T> from_roots(const values_type &roots);

        /// @brief
        /// @return
        polynomial<T> normalize() const;

        /// @brief
        /// @param p
        /// @return
        polynomial<T> &operator=(const polynomial<T> &p) = default;

        /// @brief Access operator for coefficient at a given index.
        /// @param index The index of the coefficient.
        /// @return The coefficient at the specified index.
        coeff_type operator[](size_type index) const;

        /// @brief Accesses the coefficient at the specified index for modification.
        /// @param index The index of the coefficient.
        /// @return A reference to the coefficient at the specified index.
        coeff_type &operator[](size_type index);

        /// @brief Checks if two polynomials are equal.
        /// @param p The polynomial to compare with.
        /// @return True if the polynomials are equal, false otherwise.
        bool operator==(const polynomial<T> &p) const;

        /// @brief Checks if two polynomials are not equal.
        /// @param p The polynomial to compare with.
        /// @return True if the polynomials are not equal, false otherwise.
        bool operator!=(const polynomial<T> &p) const;

        /// @brief Evaluates the polynomial at a given value.
        /// @param x The value at which to evaluate the polynomial.
        /// @return The result of the evaluation.
        value_type operator()(value_type x) const;

        /// @brief Adds a scalar to the polynomial.
        /// @param scalar The scalar value to be added.
        /// @return The polynomial resulting from adding the scalar to each coefficient.
        polynomial<T> operator+(value_type scalar) const;

        /// @brief
        /// @return
        polynomial<T> operator+() const;

        /// @brief Adds a scalar to the this polynomial.
        /// @param scalar
        /// @return The modified this polynomial.
        polynomial<T> &operator+=(value_type scalar);

        /// @brief Subtracts a scalar from the polynomial.
        /// @param scalar The scalar value to be subtracted.
        /// @return The polynomial resulting from subtracting the scalar from each coefficient.
        polynomial<T> operator-(value_type scalar) const;

        /// @brief
        /// @return
        polynomial<T> operator-() const;

        /// @brief Subtracts a scalar to the this polynomial.
        /// @param scalar The scalar value to be subtracted.
        /// @return The modified this polynomial. .
        polynomial<T> &operator-=(value_type scalar);

        /// @brief Multiplies the polynomial by a scalar.
        /// @param scalar The scalar value to be multiplied with.
        /// @return The polynomial resulting from multiplying each coefficient by the scalar.
        polynomial<T> operator*(value_type scalar) const;

        /// @brief Multiplies the polynomial by a scalar.
        /// @param scalar The scalar value to be multiplied with.
        /// @return The modified this polynomial.
        polynomial<T> &operator*=(value_type scalar);

        /// @brief
        /// @param scalar
        /// @return
        polynomial<T> operator/(value_type scalar) const;

        /// @brief
        /// @param scalar
        /// @return The modified this polynomial.
        polynomial<T> &operator/=(value_type scalar);

        /// @brief Adds two polynomials.
        /// @param p The polynomial to be added.
        /// @return The polynomial resulting from adding the coefficients of this polynomial
        /// and the input polynomial.
        polynomial<T> operator+(const polynomial<T> &p) const;

        /// @brief Subtracts a polynomial from another polynomial.
        /// @param p The polynomial to be subtracted.
        /// @return The polynomial resulting from subtracting the coefficients of the input polynomial
        /// from this polynomial.
        polynomial<T> operator-(const polynomial<T> &p) const;

        /// @brief Multiplies two polynomials.
        /// @param p The polynomial to be multiplied with.
        /// @return The polynomial resulting from multiplying the coefficients of this polynomial
        /// and the input polynomial.
        polynomial<T> operator*(const polynomial<T> &p) const;

        /// @brief
        /// @param p
        /// @return
        polynomial<T> &operator+=(const polynomial<T> &p);

        /// @brief
        /// @param p
        /// @return
        polynomial<T> &operator-=(const polynomial<T> &p);

        /// @brief
        /// @param p
        /// @return
        polynomial<T> &operator*=(const polynomial<T> &p);

        /// @brief
        /// @param p
        /// @return The modified this polynomial.
        polynomial<T> operator/(const polynomial<T> &p) const;

        /// @brief
        /// @param p
        /// @return
        polynomial<T> &operator/=(const polynomial<T> &p);

        /// @brief
        /// @param p
        /// @return The modified this polynomial.
        polynomial<T> operator%(const polynomial<T> &p) const;

        /// @brief
        /// @param p
        /// @return
        polynomial<T> &operator%=(const polynomial<T> &p);

        /// @brief
        /// @param os
        /// @param obj
        /// @return
        template<typename C>
        friend std::ostream &operator<<(std::ostream &os, const polynomial<C> &p);

        /// @brief Divides two polynomials and returns the quotient and remainder (Euclidean Division).
        /// @param dividend The polynomial to be divided.
        /// @param divisor The polynomial by which to divide.
        /// @return The result (q, r) of polynomial division containing quotient and remainder
        /// such that polynomial = divisor * q + r.
        std::tuple<polynomial<T>, polynomial<T>> divide(const polynomial<T> &divisor) const;

        /// @brief Composes the current polynomial with another polynomial, resulting in the composition p(q(x)).
        /// @param q The polynomial q(x) to be composed with the current polynomial p.
        /// @return The polynomial resulting from the composition p(q(x)).
        polynomial<T> compose(const polynomial<T> &q) const;

        /// @brief Raises the polynomial to a given non-negative integer exponent.
        /// @param exponent The exponent to raise the polynomial to.
        /// @return The resulting polynomial after raising to the given exponent.
        polynomial<T> pow(unsigned int exponent) const;

        /// @brief Computes the derive of the polynomial.
        /// @return The derive polynomial.
        polynomial<T> derive() const;

        /// @brief Computes the indefinite integral of the polynomial.
        /// @return The integral polynomial.
        polynomial<T> integrate() const;

        /// @brief Checks if the polynomial has the given root.
        /// @param value The value to check.
        /// @return True if the polynomial has the given root; otherwise false.
        bool has_root(const value_type &value) const;

        /// @brief Checks if the polynomial has given roots.
        /// @param values The values to check.
        /// @return True if the polynomial has all the given roots; otherwise false.
        bool has_roots(const values_type &values) const;

    private:
        /// @brief Constructor that creates a zero polynomial of a given degree.
        explicit polynomial<T>(size_type degree);

        /// @brief Trims leading zero coefficients from the polynomial.
        polynomial<T> &trim_coefficients();

    private:
        coeffs_type coeffs_;
    };
}

#include "polynomial.tpp"
#include "polynomial.ostream.tpp"

#endif // POLYNOMIAL_H_
