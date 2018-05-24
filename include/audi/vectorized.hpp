#ifndef AUDI_VECTORIZED_HPP
#define AUDI_VECTORIZED_HPP

#include <algorithm>
#include <boost/serialization/vector.hpp>
#include <exception>
#include <vector>

#include <piranha/math.hpp>
#include <piranha/pow.hpp>
#include <piranha/s11n.hpp>

#include <audi/back_compatibility.hpp>

#include <audi/type_traits.hpp>

// The streaming operator will only output the first MAX_STREAMED_COMPONENTS elements of the vector
#define MAX_STREAMED_COMPONENTS 5u

namespace audi
{

// This class implements a vectorized coefficient to be used as the coefficient type in a piranha polynomial
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
    vectorized() : m_c({T(0.)}){};

    // Constructor from other, compatible, types
    template <typename T1>
    explicit vectorized(T1 a) : m_c({static_cast<T>(a)})
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
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(), [&d1](const T &x) { return x + d1.m_c[0]; });
            return *this;
        } else if (this->size() == 1u) {
            T scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(), [scalar](const T &x) { return x + scalar; });
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
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(), [&d1](const T &x) { return x - d1.m_c[0]; });
            return *this;
        } else if (this->size() == 1u) {
            T scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(), [scalar](const T &x) { return scalar - x; });
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
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(), [&d1](const T &x) { return x * d1.m_c[0]; });
            return *this;
        } else if (this->size() == 1u) {
            T scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(), [scalar](const T &x) { return scalar * x; });
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
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(), [&d1](const T &x) { return x / d1.m_c[0]; });
            return *this;
        } else if (this->size() == 1u) {
            T scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(), [scalar](const T &x) { return scalar / x; });
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
} // namespace audi

namespace piranha
{
namespace math
{

// ---------------------   impl functions needed for vectorized_double to pass piranha::is_cf type trait
template <typename T>
struct is_zero_impl<audi::vectorized<T>> {
    bool operator()(const audi::vectorized<T> &v) const
    {
        return std::all_of(v.begin(), v.end(), [](const T &x) { return x == 0.; });
    }
};

template <typename T>
struct mul3_impl<audi::vectorized<T>> {
    /// Call operator.
    /**
     * @param out the output value.
     * @param a the first operand.
     * @param b the second operand.
     *
     * @return the output of piranha::mp_integer::mul().
     */
    void operator()(audi::vectorized<T> &out, const audi::vectorized<T> &a, const audi::vectorized<T> &b) const
    {
        if (a.size() == b.size()) {
            if (out.size() != a.size()) {
                out.resize(a.size());
            }
            std::transform(a.begin(), a.end(), b.begin(), out.begin(), std::multiplies<T>());
            return;
        } else if (a.size() == 1u) {
            if (out.size() != b.size()) {
                out.resize(b.size());
            }
            std::transform(b.begin(), b.end(), out.begin(), [a](const T &item) { return item * a[0]; });
            return;
        } else if (b.size() == 1u) {
            if (out.size() != a.size()) {
                out.resize(a.size());
            }
            std::transform(a.begin(), a.end(), out.begin(), [b](const T &item) { return item * b[0]; });
            return;
        }
        throw std::invalid_argument("Coefficients of different sizes in mul3");
    }
};

// ------------------ impl functions needed to have the methods partial, integrate and subs
template <typename T>
struct partial_impl<audi::vectorized<T>> {
    /// Call operator.
    /**
     * @return an instance of piranha::mp_integer constructed from zero.
     */
    audi::vectorized<T> operator()(const audi::vectorized<T> &in, const std::string &) const
    {
        return audi::vectorized<T>(std::vector<T>(in.size(), 0.));
    }
};
template <typename T, typename U>
struct pow_impl<audi::vectorized<T>, U> {
    /// Call operator.
    /**
     * @param c the input
     * @param exp the exponent
     * @return the exp operator applied to all elements of the input
     */
    audi::vectorized<T> operator()(const audi::vectorized<T> &c, const U &exp) const
    {
        auto retval(c);
        std::transform(retval.begin(), retval.end(), retval.begin(), [exp](const T &x) { return piranha::math::pow(x, exp); });
        return retval;
    };
};

} // namespace math

template <typename Archive, typename T>
struct boost_save_impl<Archive, audi::vectorized<T>> : boost_save_via_boost_api<Archive, audi::vectorized<T>> {
};

template <typename Archive, typename T>
struct boost_load_impl<Archive, audi::vectorized<T>> : boost_load_via_boost_api<Archive, audi::vectorized<T>> {
};

template <typename T>
struct zero_is_absorbing<audi::vectorized<T>> : std::false_type {
};

} // namespace piranha

#undef MAX_STREAMED_COMPONENTS
#endif
