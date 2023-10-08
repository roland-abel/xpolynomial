/// @file legendre_polynomial.h
///
/// @author Roland Abel
/// @date 08.10.2023

#ifndef LEGENDRE_POLYNOMIAL_H_
#define LEGENDRE_POLYNOMIAL_H_

#include <vector>
#include "polynomial.h"

namespace xmath {

    template<typename T>
    class legendre_polynomial : public polynomial<T> {
    public:
        using value_t = polynomial<T>::value_type;

    public:
        explicit legendre_polynomial(size_t order);

    public:
        const std::vector<value_t> &weights() const;

        const std::vector<value_t> &roots() const;

    private:
        static polynomial<T> create(size_t order);

    private:
        std::vector<value_t> weight_;
        std::vector<value_t> roots_;
    };
}

#include "legendre_polynomial.tpp"

#endif // LEGENDRE_POLYNOMIAL_H_