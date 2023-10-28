/// @file matrix.tpp
///
/// @author Roland Abel
/// @date 08.10.2023

#include <utility>
#include "matrix.h"
#include "utils.h"

namespace xmath {

    template<typename T>
    matrix<T>::MatrixProxy::MatrixProxy(const matrix &mat, size_type row)
            : inner_matrix_(mat), row_(row) {
    }

    template<typename T>
    matrix<T>::value_type matrix<T>::MatrixProxy::operator[](size_type col) {
        return inner_matrix_(inner_matrix_.index(row_, col));
    }

    template<typename T>
    matrix<T>::MatrixProxy matrix<T>::operator[](size_type row) const {
        return {*this, row};
    }

    template<typename T>
    matrix<T>::matrix(size_type num_rows, size_type num_cols)
            : matrix(num_rows, num_cols, 0.0) {
    }

    template<typename T>
    matrix<T>::matrix(size_type num_rows, size_type num_cols, values_type coeffs)
            : num_rows_(num_rows), num_cols_(num_cols), coeffs_(std::move(coeffs)) {
        check_dimension();
    }

    template<typename T>
    matrix<T>::matrix(size_type num_rows, size_type num_cols, const value_type &value)
            : num_rows_(num_rows), num_cols_(num_cols), coeffs_(num_rows * num_cols, value) {
    }

    template<typename T>
    matrix<T>::matrix(const std::vector<std::vector<value_type>> &coeffs)
            : num_rows_(coeffs.size()) {
        num_cols_ = 0;
        for (const auto &coefficient: coeffs) {
            num_cols_ = std::max(num_cols_, coefficient.size());
        }

        coeffs_.resize(num_rows_ * num_cols_, 0.0);

        for (auto row = 0; row < coeffs.size(); ++row) {
            auto inner_values = coeffs.at(row);
            for (auto col = 0; col < inner_values.size(); ++col) {
                coeffs_.at(index(row, col)) = inner_values.at(col);
            }
        }
    }

    template<typename C>
    std::ostream &operator<<(std::ostream &os, const matrix<C> &mat) {
        const char fill = ' ';

        for (typename matrix<C>::size_type i = 0; i < mat.rows(); ++i) {
            for (typename matrix<C>::size_type j = 0; j < mat.cols(); ++j) {
                os << mat(i, j) << fill;
            }
            os << std::endl;
        }
        return os;
    }

    template<typename T>
    matrix<T> matrix<T>::transpose() const {
        matrix mat(cols(), rows());
        for (size_type i = 0; i < mat.rows(); i++) {
            for (size_type j = 0; j < mat.cols(); j++) {
                mat(i, j) = at(j, i);
            }
        }
        return mat;
    }

    template<typename T>
    matrix<T> matrix<T>::operator+(const matrix &mat) const {
        return apply([&](const double &_, const size_type &idx) { return this->operator()(idx) + mat(idx); });
    }

    template<typename T>
    matrix<T> matrix<T>::operator-(const matrix &mat) const {
        return apply([&](const double &_, const size_type &idx) { return this->operator()(idx) - mat(idx); });
    }

    template<typename T>
    matrix<T> matrix<T>::operator*(const matrix &mat) const {
        auto m = matrix(rows(), mat.cols());
        for (size_type i = 0; i < rows(); ++i) {
            for (size_type k = 0; k < mat.cols(); ++k) {
                for (size_type j = 0; j < cols(); ++j)
                    m(i, k) += at(i, j) * mat(j, k);
            }
        }
        return m;
    }

    template<typename T>
    matrix<T> matrix<T>::make_matrix(const value_type &value) const {
        return {rows(), cols(), value};
    }

    template<typename T>
    matrix<T> matrix<T>::apply(const std::function<value_type(const value_type &, const size_type &)> &func) const {
        auto mat = make_matrix();
        for (auto [itr, end, idx] = std::tuple{coeffs_.cbegin(), coeffs_.cend(), 0};
             idx < coeffs_.size(); ++itr, ++idx) {
            mat(idx) = func(*itr, idx);
        }
        return mat;
    }

    template<typename T>
    matrix<T> operator*(const matrix<T> &lhs, const matrix<T> &rhs) {
        return lhs.operator*(rhs);
    }

    template<typename T>
    matrix<T> matrix<T>::operator*(const double &scalar) const {
        return apply([&](const double &_, const size_type &idx) { return scalar * this->operator()(idx); });
    }

    template<typename T>
    matrix<T> matrix<T>::operator/(const double &scalar) const {
        return apply([&](const double &_, const size_type &idx) { return this->operator()(idx) / scalar; });
    }

    template<typename T>
    void matrix<T>::check_dimension() {
        if (num_rows_ * num_cols_ != coeffs_.size()) {
            throw std::runtime_error(
                    "The dimension numbers of rows columns don't fit the number of initial values.");
        }
    }

    template<typename T>
    matrix<T> operator*(const double &value, const matrix<T> &mat) {
        return mat * value;
    }


    template<typename T>
    matrix<T> operator/(const double &value, const matrix<T> &mat) {
        return mat;
    }

    template<typename T>
    bool matrix<T>::is_symmetrical() const noexcept {
        if (!is_square()) {
            return false;
        }

        if (rows() == 1) {
            return true;
        }

        for (size_type i = 0; i < rows(); i++) {
            for (size_type j = i + 1; j < cols(); j++) {
                if (!nearly_equal<value_type>(at(i, j), at(j, i), epsilon)) {
                    return false;
                }
            }
        }
        return true;
    }

    template<typename T>
    matrix<T>::value_type matrix<T>::operator()(size_type idx) const {
        return coeffs_.at(idx);
    }

    template<typename T>
    matrix<T>::value_type &matrix<T>::operator()(size_type idx) {
        return coeffs_.at(idx);
    }

    template<typename T>
    matrix<T>::value_type matrix<T>::operator()(size_type row, size_type col) const {
        return coeffs_.at(index(row, col));
    }

    template<typename T>
    matrix<T>::value_type &matrix<T>::operator()(size_type row, size_type col) {
        return coeffs_.at(index(row, col));
    }

    template<typename T>
    matrix<T>::value_type matrix<T>::at(size_type row, size_type col) const {
        return operator()(row, col);
    }

    template<typename T>
    matrix<T>::value_type &matrix<T>::at(size_type row, size_type col) {
        return operator()(row, col);
    }

    template<typename T>
    bool matrix<T>::is_square() const noexcept {
        return num_rows_ == num_cols_;
    }

    template<typename T>
    bool matrix<T>::is_empty() const noexcept {
        return coeffs_.empty();
    }

    template<typename T>
    void matrix<T>::resize(size_type rows, size_type cols) {
        num_rows_ = rows;
        num_cols_ = cols;
    }

    template<typename T>
    bool matrix<T>::operator==(const matrix &mat) const {
        if (mat.rows() != rows() || mat.cols() != cols())
            return false;

        return std::equal(
                coeffs_.begin(),
                coeffs_.end(),
                mat.coeffs_.begin(),
                [](value_type x, value_type y) { return nearly_equal<value_type>(x, y, epsilon); });
    }

    template<typename T>
    bool matrix<T>::operator!=(const matrix &mat) const {
        return !(*this == mat);
    }
}