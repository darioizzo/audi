#ifndef AUDI_VECTORIZED_HPP
#define AUDI_VECTORIZED_HPP

#include <algorithm>
#include <cstddef>
#include <exception>
#include <initializer_list>
#include <type_traits>
#include <vector>

#include <boost/serialization/vector.hpp>

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
inline void negate(vectorized<T> &c)
{
    for (auto &x : c) {
        obake::negate(x);
    }
}

template <typename T>
inline std::size_t byte_size(const vectorized<T> &c)
{
    return sizeof(T) * c.size() + sizeof(c);
}

template <typename T>
inline void fma3(vectorized<T> &ret, const vectorized<T> &x, const vectorized<T> &y)
{
    if (x.size() == y.size()) {
        const auto size = x.size();
        ret.resize(size);
        std::transform(x.begin(), x.end(), y.begin(), ret.begin(), std::multiplies<T>());
    } else if (x.size() == 1u) {
        const auto size = y.size();
        ret.resize(size);
        std::transform(y.begin(), y.end(), ret.begin(), [scalar = x[0]](const auto &v) { return v * scalar; });
    } else if (y.size() == 1u) {
        const auto size = x.size();
        ret.resize(size);
        std::transform(x.begin(), x.end(), ret.begin(), [scalar = y[0]](const auto &v) { return v * scalar; });
    } else {
        throw std::invalid_argument("Coefficients of different sizes in fma3");
    }
}

} // namespace audi

#undef MAX_STREAMED_COMPONENTS
#endif
