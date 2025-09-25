// Copyright © 2018–2025 Dario Izzo (dario.izzo@gmail.com),
// Francesco Biscani (bluescarni@gmail.com),
// Sean Cowan (lambertarc@icloud.com)
//
// This file is part of the audi library.
//
// The audi library is free software: you can redistribute it and/or modify
// it under the terms of either:
//   - the GNU General Public License as published by the Free Software
//     Foundation, either version 3 of the License, or (at your option)
//     any later version, or
//   - the GNU Lesser General Public License as published by the Free
//     Software Foundation, either version 3 of the License, or (at your
//     option) any later version.
//
// The audi library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License and the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// and the GNU Lesser General Public License along with the audi library.
// If not, see <https://www.gnu.org/licenses/>.

#ifndef AUDI_VECTORIZED_HPP
#define AUDI_VECTORIZED_HPP

#include <algorithm>
#include <cstddef>
#include <exception>
#include <initializer_list>
#include <type_traits>
#include <vector>

#include <boost/serialization/vector.hpp>

#include <obake/math/fma3.hpp>
#include <obake/math/negate.hpp>
#include <obake/math/pow.hpp>

#include <audi/back_compatibility.hpp>
#include <audi/type_traits.hpp>

// The streaming operator will only output the first MAX_STREAMED_COMPONENTS elements of the vector
#define MAX_STREAMED_COMPONENTS 5u

namespace audi
{

// This class implements a vectorized coefficient to be used as the coefficient type in an obake polynomial
// The coefficient is, essentially a vector of doubles [a0,a1,a2,...an] on which all arithmetic operations and
// function calls operate element-wise.

template <typename T>
struct vectorized {
    // type to be vectorized
    using v_type = T;

    // enables operators between vectorized types and arithmetic types (extended to include mppp::real128)
    template <typename T1>
    using operator_enabler
        = enable_if_t<(std::is_same<T1, vectorized<T>>::value || audi::is_arithmetic<T1>::value), int>;

    // Default constructor. Constructs [0.]
    vectorized() : m_c{T(0.)} {};

    // Constructor from other, compatible, types
    template <typename T1, std::enable_if_t<std::is_constructible_v<T, T1>, int> = 0>
    explicit vectorized(T1 a) : m_c{static_cast<T>(a)}
    {
    }

    // Constructor from an std::vector
    explicit vectorized(const std::vector<T> &c) : m_c(c)
    {
        if (m_c.size() == 0) {
            throw std::invalid_argument("Cannot build an empty coefficient_v (lvalue)");
        }
    };
    // Constructor from an std::vector r value
    explicit vectorized(std::vector<T> &&c) : m_c(c)
    {
        if (m_c.size() == 0) {
            throw std::invalid_argument("Cannot build an empty coefficient_v (rvalue)");
        }
    };
    explicit vectorized(std::initializer_list<T> c1) : m_c(c1)
    {
        if (m_c.size() == 0) {
            throw std::invalid_argument("Cannot build an empty coefficient_v (initializer)");
        }
    }
    // ----------------- Juice implementation of operators. It deals with the case [b1] op [a1,a2,..an].
    vectorized<T> &operator+=(const vectorized<T> &d1)
    {
        if (d1.size() == this->size()) {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(), std::plus<T>());
            return *this;
        } else if (d1.size() == 1u) {
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(),
                           [&d1](const T &x) { return x + d1.m_c[0]; });
            return *this;
        } else if (this->size() == 1u) {
            T scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(),
                           [scalar](const T &x) { return x + scalar; });
            return *this;
        }
        throw std::invalid_argument("Coefficients of different sizes in +");
    }
    vectorized<T> &operator-=(const vectorized<T> &d1)
    {
        if (d1.size() == this->size()) {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(), std::minus<T>());
            return *this;
        } else if (d1.size() == 1u) {
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(),
                           [&d1](const T &x) { return x - d1.m_c[0]; });
            return *this;
        } else if (this->size() == 1u) {
            T scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(),
                           [scalar](const T &x) { return scalar - x; });
            return *this;
        }
        throw std::invalid_argument("Coefficients of different sizes in -");
    }
    vectorized<T> &operator*=(const vectorized<T> &d1)
    {
        if (d1.size() == this->size()) {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(), std::multiplies<T>());
            return *this;
        } else if (d1.size() == 1u) {
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(),
                           [&d1](const T &x) { return x * d1.m_c[0]; });
            return *this;
        } else if (this->size() == 1u) {
            T scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(),
                           [scalar](const T &x) { return scalar * x; });
            return *this;
        }
        throw std::invalid_argument("Coefficients of different sizes in *");
    }
    vectorized<T> &operator/=(const vectorized<T> &d1)
    {
        if (d1.size() == this->size()) {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(), std::divides<T>());
            return *this;
        } else if (d1.size() == 1u) {
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(),
                           [&d1](const T &x) { return x / d1.m_c[0]; });
            return *this;
        } else if (this->size() == 1u) {
            T scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(),
                           [scalar](const T &x) { return scalar / x; });
            return *this;
        }
        throw std::invalid_argument("Coefficients of different sizes in /");
    }
    template <typename T1>
    vectorized<T> &operator/=(const T1 &d1)
    {
        std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(), [&d1](const T &x) { return x / d1; });
        return *this;
    }
    vectorized<T> operator-() const
    {
        vectorized<T> retval(m_c);
        transform(retval.m_c.begin(), retval.m_c.end(), retval.m_c.begin(), std::negate<T>());
        return retval;
    }
    T operator[](const typename std::vector<T>::size_type idx) const
    {
        if (m_c.size() == 1u) {
            return m_c[0];
        }
        return m_c[idx];
    }
    friend std::ostream &operator<<(std::ostream &os, const vectorized<T> &d)
    {
        os << "[";
        if (d.size() <= MAX_STREAMED_COMPONENTS) {
            for (auto i = 0u; i < d.size() - 1u; ++i) {
                os << d.m_c[i] << ", ";
            }
            os << d.m_c[d.m_c.size() - 1u] << "]";
        } else {
            for (auto i = 0u; i < MAX_STREAMED_COMPONENTS; ++i) {
                os << d.m_c[i] << ", ";
            }
            os << "... ]";
        }
        return os;
    }
    typename std::vector<T>::size_type size() const
    {
        return m_c.size();
    }
    typename std::vector<T>::const_iterator begin() const
    {
        return m_c.begin();
    }
    typename std::vector<T>::const_iterator end() const
    {
        return m_c.end();
    }
    typename std::vector<T>::iterator begin()
    {
        return m_c.begin();
    }
    typename std::vector<T>::iterator end()
    {
        return m_c.end();
    }
    void resize(typename std::vector<T>::size_type new_size)
    {
        m_c.resize(new_size);
    }
    void resize(typename std::vector<T>::size_type new_size, T val)
    {
        m_c.resize(new_size, val);
    }
    void set_value(const typename std::vector<T>::size_type idx, T val)
    {
        m_c[idx] = val;
    }
    const std::vector<T> &get_v() const
    {
        return m_c;
    }
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar &m_c;
    }

    // The container
    std::vector<T> m_c;
}; // end of struct vectorized

/// Type is a vectorized
/**
 * Checks whether T is a vectorized type. Provides the member constant value which is equal to true,
 * if T is the type vectorized<U> for any U.
 *
 * \tparam T a type to check
 */

template <typename T>
struct is_vectorized : std::false_type {
};
template <typename T>
struct is_vectorized<vectorized<T>> : std::true_type {
};

// ------------------- Binary arithmetic operators implemented using the available operators +=,-=, etc.
template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && is_vectorized<T2>::value), int> = 0>
inline T1 operator+(T1 d1, const T2 &d2)
{
    d1 += d2;
    return d1;
}
template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && audi::is_arithmetic<T2>::value), int> = 0>
inline T1 operator+(T1 d1, const T2 &d2v)
{
    T1 d2(d2v);
    return d1 + d2;
}
template <typename T1, typename T2, enable_if_t<(audi::is_arithmetic<T1>::value && is_vectorized<T2>::value), int> = 0>
T2 operator+(T1 d1v, const T2 &d2)
{
    T2 d1(d1v);
    return d1 + d2;
}

template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && is_vectorized<T2>::value), int> = 0>
inline T1 operator-(T1 d1, const T2 &d2)
{
    d1 -= d2;
    return d1;
}
template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && audi::is_arithmetic<T2>::value), int> = 0>
inline T1 operator-(T1 d1, const T2 &d2v)
{
    T1 d2(d2v);
    return d1 - d2;
}
template <typename T1, typename T2, enable_if_t<(audi::is_arithmetic<T1>::value && is_vectorized<T2>::value), int> = 0>
inline T2 operator-(T1 d1v, const T2 &d2)
{
    T2 d1(d1v);
    return d1 - d2;
}

template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && is_vectorized<T2>::value), int> = 0>
inline T1 operator*(T1 d1, const T2 &d2)
{
    d1 *= d2;
    return d1;
}
template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && audi::is_arithmetic<T2>::value), int> = 0>
inline T1 operator*(T1 d1, const T2 &d2v)
{
    T1 d2(d2v);
    return d1 * d2;
}
template <typename T1, typename T2, enable_if_t<(audi::is_arithmetic<T1>::value && is_vectorized<T2>::value), int> = 0>
inline T2 operator*(T1 d1v, const T2 &d2)
{
    T2 d1(d1v);
    return d1 * d2;
}

template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && is_vectorized<T2>::value), int> = 0>
inline T1 operator/(T1 d1, const T2 &d2)
{
    d1 /= d2;
    return d1;
}
template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && audi::is_arithmetic<T2>::value), int> = 0>
inline T1 operator/(T1 d1, const T2 &d2v)
{
    T1 d2(d2v);
    return d1 / d2;
}
template <typename T1, typename T2, enable_if_t<(audi::is_arithmetic<T1>::value && is_vectorized<T2>::value), int> = 0>
inline T2 operator/(T1 d1v, const T2 &d2)
{
    T2 d1(d1v);
    return d1 / d2;
}

// Various comparison operators
template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && is_vectorized<T2>::value), int> = 0>
inline bool operator==(const T1 &d1, const T2 &d2)
{
    if (d1.size() == d2.size()) {
        return d1.m_c == d2.m_c;
    } else if (d1.size() == 1u) {
        return std::all_of(d2.begin(), d2.end(), [d1](const typename T1::v_type &x) { return x == d1[0]; });
    } else if (d2.size() == 1u) {
        return std::all_of(d1.begin(), d1.end(), [d2](const typename T1::v_type &x) { return x == d2[0]; });
    }
    return false;
}

template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && audi::is_arithmetic<T2>::value), int> = 0>
inline bool operator==(const T1 &d1, const T2 &d2v)
{
    T1 d2(d2v);
    return d1 == d2;
}
template <typename T1, typename T2, enable_if_t<(audi::is_arithmetic<T1>::value && is_vectorized<T2>::value), int> = 0>
inline bool operator==(const T1 &d1v, const T2 &d2)
{
    T2 d1(d1v);
    return d1 == d2;
}

template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value || is_vectorized<T2>::value), int> = 0>
inline bool operator!=(const T1 &d1, const T2 &d2)
{
    return !(d1 == d2);
}

template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && is_vectorized<T2>::value), int> = 0>
inline bool operator>(const T1 &d1, const T2 &d2)
{
    if (d1.size() == d2.size()) {
        return d1.m_c > d2.m_c;
    } else if (d1.size() == 1u) {
        return std::all_of(d2.begin(), d2.end(), [d1](const typename T1::v_type &x) { return d1[0] > x; });
    } else if (d2.size() == 1u) {
        return std::all_of(d1.begin(), d1.end(), [d2](const typename T1::v_type &x) { return x > d2[0]; });
    }
    return false;
}
template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && audi::is_arithmetic<T2>::value), int> = 0>
inline bool operator>(const T1 &d1, const T2 &d2v)
{
    T1 d2(d2v);
    return d1 > d2;
}
template <typename T1, typename T2, enable_if_t<(audi::is_arithmetic<T1>::value && is_vectorized<T2>::value), int> = 0>
inline bool operator>(const T1 &d1v, const T2 &d2)
{
    T2 d1(d1v);
    return d1 > d2;
}

template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && is_vectorized<T2>::value), int> = 0>
inline bool operator<(const T1 &d1, const T2 &d2)
{
    if (d1.size() == d2.size()) {
        return d1.m_c < d2.m_c;
    } else if (d1.size() == 1u) {
        return std::all_of(d2.begin(), d2.end(), [d1](const typename T1::v_type &x) { return d1[0] < x; });
    } else if (d2.size() == 1u) {
        return std::all_of(d1.begin(), d1.end(), [d2](const typename T1::v_type &x) { return x < d2[0]; });
    }
    return false;
}
template <typename T1, typename T2, enable_if_t<(is_vectorized<T1>::value && audi::is_arithmetic<T2>::value), int> = 0>
inline bool operator<(const T1 &d1, const T2 &d2v)
{
    T1 d2(d2v);
    return d1 < d2;
}
template <typename T1, typename T2, enable_if_t<(audi::is_arithmetic<T1>::value && is_vectorized<T2>::value), int> = 0>
inline bool operator<(const T1 &d1v, const T2 &d2)
{
    T2 d1(d1v);
    return d1 < d2;
}

// ---------------------   impl functions needed for vectorized_double to pass the obake::is_cf type trait
template <typename T>
inline bool is_zero(const vectorized<T> &v)
{
    return std::all_of(v.begin(), v.end(), [](const T &x) { return x == 0.; });
}

template <typename T>
inline vectorized<T> diff(const vectorized<T> &in, const std::string &)
{
    return vectorized<T>(std::vector<T>(in.size(), 0.));
}

template <typename T, typename U>
inline vectorized<T> pow(const vectorized<T> &c, const U &exp)
{
    auto retval(c);
    std::transform(retval.begin(), retval.end(), retval.begin(), [exp](const T &x) { return obake::pow(x, exp); });
    return retval;
}

template <typename T>
inline void negate(vectorized<T> &c) noexcept(noexcept(obake::negate(std::declval<T &>())))
{
    for (auto &x : c) {
        obake::negate(x);
    }
}

template <typename T>
inline std::size_t byte_size(const vectorized<T> &c) noexcept
{
    return sizeof(T) * c.size() + sizeof(c);
}

template <typename T>
inline void fma3(vectorized<T> &ret, const vectorized<T> &x, const vectorized<T> &y)
{
    const auto x_size = x.size(), y_size = y.size(), ret_size = ret.size();
    constexpr auto use_fma = obake::is_mult_addable_v<T &, const T &, const T &>;

    // Most frequent case: all dimensions are equal case: n,n,n
    if ((x_size == y_size) && (x_size == ret_size)) {
        auto ret_it = ret.begin();
        for (decltype(x.size()) i = 0; i < x_size; ++i, ++ret_it) {
            if constexpr (use_fma) {
                obake::fma3(*ret_it, x[i], y[i]);
            } else {
                *ret_it += x[i] * y[i];
            }
        }
    } else if (ret_size == 1u) { // We are in the case 1, .., ..
        if (x_size == 1u) {      // We are in the case 1, 1, n with n > 1
            ret.resize(y_size, ret[0]);
            auto ret_it = ret.begin();
            const auto x0 = x[0];
            for (decltype(y.size()) i = 0; i < y_size; ++i, ++ret_it) {
                if constexpr (use_fma) {
                    obake::fma3(*ret_it, x0, y[i]);
                } else {
                    *ret_it += x0 * y[i];
                }
            }
        } else if (y_size == 1u) { // We are in the case 1, n, 1 with n > 1
            ret.resize(x_size, ret[0]);
            auto ret_it = ret.begin();
            const auto y0 = y[0];
            for (decltype(x.size()) i = 0; i < x_size; ++i, ++ret_it) {
                if constexpr (use_fma) {
                    obake::fma3(*ret_it, x[i], y0);
                } else {
                    *ret_it += x[i] * y0;
                }
            }
        } // We are in the case 1, n, n with n > 1
        else if (y_size == x_size) {
            ret.resize(x_size, ret[0]);
            auto ret_it = ret.begin();
            for (decltype(x.size()) i = 0; i < x_size; ++i, ++ret_it) {
                if constexpr (use_fma) {
                    obake::fma3(*ret_it, x[i], y[i]);
                } else {
                    *ret_it += x[i] * y[i];
                }
            }
        } else {
            throw std::invalid_argument("Coefficients of different sizes in fma3");
        }
    } else if (x_size == 1u) { // We are in the case n, 1
        if (y_size == 1u) {    // We are in the case n, 1, 1
            auto ret_it = ret.begin();
            const auto x0 = x[0];
            const auto y0 = y[0];
            for (decltype(ret.size()) i = 0u; i < ret_size; ++i, ++ret_it) {
                if constexpr (use_fma) {
                    obake::fma3(*ret_it, x0, y0);
                } else {
                    *ret_it += x0 * y0;
                }
            }
        } else if (y_size == ret_size) { // We are in the case n, 1, n
            auto ret_it = ret.begin();
            const auto x0 = x[0];
            for (decltype(y.size()) i = 0; i < y_size; ++i, ++ret_it) {
                if constexpr (use_fma) {
                    obake::fma3(*ret_it, x0, y[i]);
                } else {
                    *ret_it += x0 * y[i];
                }
            }
        } else {
            throw std::invalid_argument("Coefficients of different sizes in fma3");
        }

    } else if (y_size == 1u) {    // We are in the case n, n, 1
        if (ret_size == x_size) { // We are in the case n, n, 1
            auto ret_it = ret.begin();
            const auto y0 = y[0];
            for (decltype(x.size()) i = 0; i < x_size; ++i, ++ret_it) {
                if constexpr (use_fma) {
                    obake::fma3(*ret_it, x[i], y0);
                } else {
                    *ret_it += x[i] * y0;
                }
            }
        } else {
            throw std::invalid_argument("Coefficients of different sizes in fma3");
        }
    }
}

} // namespace audi

#undef MAX_STREAMED_COMPONENTS
#endif
