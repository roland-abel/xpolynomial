/// @file matrix.h
///
/// @author Roland Abel
/// @date 08.10.2023

#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
#include <functional>
#include <ostream>

namespace xmath {

    template<typename T>
    class matrix {
    public:
        using coeffs_type = std::vector<T>;
        using coeff_type = coeffs_type::value_type;
        using size_type = coeffs_type::size_type;

        static constexpr coeff_type tolerance = 1e-5;

        class MatrixProxy {
        public:
            MatrixProxy(const matrix<T> &mat, size_type row);

            coeff_type operator[](size_type col);

        private:
            const matrix &inner_matrix_;
            size_type row_;
        };

    public:

        explicit matrix(size_type num_rows = 0, size_type num_cols = 0);

        matrix(size_type num_rows, size_type num_cols, coeffs_type coeffs);

        matrix(size_type num_rows, size_type num_cols, const coeff_type &value);

        explicit matrix(const std::vector<std::vector<coeff_type>> &coeffs);

        matrix(const matrix &mat) = default;

        ~matrix() = default;

    public:
        constexpr inline size_type rows() const { return num_rows_; }

        constexpr inline size_type cols() const { return num_cols_; }

        constexpr inline size_type index(size_type row, size_type col) const { return row * num_cols_ + col; }

    public:
        template<typename C>
        friend std::ostream &operator<<(std::ostream &os, const matrix<C> &mat);

        MatrixProxy operator[](size_type row) const;

        coeff_type operator()(size_type idx) const;

        coeff_type &operator()(size_type idx);

        coeff_type operator()(size_type row, size_type col) const;

        coeff_type &operator()(size_type row, size_type col);

        bool operator==(const matrix &mat) const;

        bool operator!=(const matrix &mat) const;

        matrix<T> operator+(const matrix &mat) const;

        matrix<T> operator-(const matrix &mat) const;

        matrix<T> operator*(const matrix &mat) const;

        matrix<T> operator*(const double &scalar) const;

        matrix<T> operator/(const double &scalar) const;

    public:
        coeff_type at(size_type row, size_type col) const;

        coeff_type &at(size_type row, size_type col);

        bool is_square() const noexcept;

        bool is_empty() const noexcept;

        bool is_symmetrical() const noexcept;

        void resize(size_type rows, size_type cols);

        matrix<T> transpose() const;

    private:
        void check_dimension();

        matrix<T> make_matrix(const coeff_type &value = static_cast<coeff_type >(0.0)) const;

        matrix<T> apply(const std::function<coeff_type(const coeff_type &, const size_type &)> &func) const;

        // The number of rows of the matrix.
        size_type num_rows_ = 0;

        // The number of columns of the matrix,
        size_type num_cols_ = 0;

        // The matrix's data container.
        coeffs_type coeffs_;
    };
}

#include "matrix.tpp"

#endif // MATRIX_H_