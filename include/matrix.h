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
    struct matrix_specification {
        using value_type = T;
        using size_type = size_t;
        using floating_point_type = value_type;
        static constexpr floating_point_type epsilon = 1e-5;
        static constexpr value_type one = 1.0;
        static constexpr value_type zero = 0.0;
    };

    template<>
    struct matrix_specification<double> {
        using value_type = double;
        using size_type = size_t;
        using floating_point_type = value_type;
        static constexpr floating_point_type epsilon = 1e-5;
        static constexpr value_type one = 1.0;
        static constexpr value_type zero = 0.0;
    };

    template<typename T, typename S = matrix_specification<T>>
    class matrix {
    };

    template<typename T>
    class matrix<T, matrix_specification<T>> {
    public:
        using spec = matrix_specification<T>;
        using value_type = spec::value_type;
        using size_type = spec::size_type;
        using floating_point_type = spec::floating_point_type;
        using values_type = std::vector<T>;
        static constexpr floating_point_type epsilon = spec::epsilon;

        class MatrixProxy {
        public:
            MatrixProxy(const matrix<T> &mat, size_type row);

            value_type operator[](size_type col);

        private:
            const matrix<T> &inner_matrix_;
            size_type row_;
        };

        explicit matrix(size_type num_rows = 0, size_type num_cols = 0);

        matrix(size_type num_rows, size_type num_cols, std::initializer_list<value_type> coeffs);

        matrix(size_type num_rows, size_type num_cols, const values_type &coeffs);

        matrix(size_type num_rows, size_type num_cols, const value_type &value);

        explicit matrix(const std::vector<values_type> &coeffs);

        explicit matrix(size_type num_rows, size_type num_cols, const std::ranges::range auto &range);

        matrix(const matrix<T> &mat) = default;

        ~matrix() = default;

        constexpr inline size_type rows() const { return num_rows_; }

        constexpr inline size_type cols() const { return num_cols_; }

        constexpr inline size_type index(size_type row, size_type col) const { return row * num_cols_ + col; }

        /// @brief Returns the coefficients of the polynomial.
        /// @return The coefficients.
        const values_type &coefficients() const;

        template<typename C>
        friend std::ostream &operator<<(std::ostream &os, const matrix<C> &mat);

        MatrixProxy operator[](size_type row) const;

        value_type operator()(size_type idx) const;

        value_type &operator()(size_type idx);

        value_type operator()(size_type row, size_type col) const;

        value_type &operator()(size_type row, size_type col);

        bool operator==(const matrix<T> &mat) const;

        bool operator!=(const matrix<T> &mat) const;

        matrix<T> operator+(const matrix<T> &M) const;

        matrix<T> operator-(const matrix<T> &mat) const;

        matrix<T> operator*(const matrix<T> &mat) const;

        matrix<T> operator*(const double &scalar) const;

        matrix<T> operator/(const double &scalar) const;

    public:
        value_type at(size_type row, size_type col) const;

        value_type &at(size_type row, size_type col);

        bool is_square() const noexcept;

        bool is_empty() const noexcept;

        bool is_zero() const noexcept;

        bool is_symmetrical() const noexcept;

        void resize(size_type rows, size_type cols);

        matrix<T> transpose() const;

    private:
        void check_dimension();

        // The number of rows of the matrix.
        size_type num_rows_ = 0;

        // The number of columns of the matrix,
        size_type num_cols_ = 0;

        // The matrix's data container.
        values_type coeffs_;
    };
}

#include "matrix.tpp"

#endif // MATRIX_H_