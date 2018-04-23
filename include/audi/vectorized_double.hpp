#ifndef AUDI_VECTORIZED_DOUBLE_HPP
#define AUDI_VECTORIZED_DOUBLE_HPP

#include <algorithm>
#include <boost/serialization/vector.hpp>
#include <exception>
#include <vector>

#include <piranha/math.hpp>
#include <piranha/pow.hpp>
#include <piranha/s11n.hpp>

// The streaming operator will only output the first MAX_STREAMED_COMPONENTS elements of the vector
#define MAX_STREAMED_COMPONENTS 5u

namespace audi
{

// This class implements a vectorized coefficient to be used as the coefficient type in a piranha polynomial
// The coefficient is, essentially a vector of doubles [a0,a1,a2,...an] on which all arithmetic operations and
// function calls operate element-wise.

class vectorized_double
{
public:
    // Default constructor. Constructs [0.]
    vectorized_double() : m_c({0.}){};
    // Constructor from int. Its mandatory for piranha::polynomial coefficient
    explicit vectorized_double(int a) : m_c({static_cast<double>(a)}){};
    // Constructor from double value a. Construct [a] TODO: should we make this explicit?
    vectorized_double(double a) : m_c({a}){};
    // Constructor from an std::vector
    explicit vectorized_double(const std::vector<double> &c) : m_c(c)
    {
        if (m_c.size() == 0) {
            throw std::invalid_argument("Cannot build an empty coefficient_v (lvalue)");
        }
    };
    // Constructor from an std::vector r value
    explicit vectorized_double(std::vector<double> &&c) : m_c(c)
    {
        if (m_c.size() == 0) {
            throw std::invalid_argument("Cannot build an empty coefficient_v (rvalue)");
        }
    };
    explicit vectorized_double(std::initializer_list<double> c) : m_c(c)
    {
        if (m_c.size() == 0) {
            throw std::invalid_argument("Cannot build an empty coefficient_v (initializer)");
        }
    };
    // ------------------- Binary arithmetic operators implemented using +=,-=, etc.
    friend vectorized_double operator+(const vectorized_double &d1, const vectorized_double &d2)
    {
        vectorized_double retval(d1);
        retval += d2;
        return retval;
    };
    friend vectorized_double operator-(const vectorized_double &d1, const vectorized_double &d2)
    {
        vectorized_double retval(d1);
        retval -= d2;
        return retval;
    };
    friend vectorized_double operator*(const vectorized_double &d1, const vectorized_double &d2)
    {
        vectorized_double retval(d1);
        retval *= d2;
        return retval;
    };
    friend vectorized_double operator*(int d1, const vectorized_double &d2)
    {
        vectorized_double retval(d1);
        retval *= d2;
        return retval;
    };
    friend vectorized_double operator/(const vectorized_double &d1, const vectorized_double &d2)
    {
        vectorized_double retval(d1);
        retval /= d2;
        return retval;
    };

    // ----------------- Juice implementation of the operators. It also deals with the case [b1] op [a1,a2,..an] to
    // take care of scalar multiplication/division etc.
    vectorized_double &operator+=(const vectorized_double &d1)
    {
        if (d1.size() == this->size()) {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(), std::plus<double>());
            return *this;
        } else if (d1.size() == 1u) {
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(),
                           [&d1](double x) { return x + d1.m_c[0]; });
            return *this;
        } else if (this->size() == 1u) {
            double scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(), [scalar](double x) { return x + scalar; });
            return *this;
        }
        throw std::invalid_argument("Coefficients of different sizes in +");
    }
    vectorized_double &operator-=(const vectorized_double &d1)
    {
        if (d1.size() == this->size()) {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(), std::minus<double>());
            return *this;
        } else if (d1.size() == 1u) {
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(),
                           [&d1](double x) { return x - d1.m_c[0]; });
            return *this;
        } else if (this->size() == 1u) {
            double scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(), [scalar](double x) { return scalar - x; });
            return *this;
        }
        throw std::invalid_argument("Coefficients of different sizes in -");
    }
    vectorized_double &operator*=(const vectorized_double &d1)
    {
        if (d1.size() == this->size()) {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(),
                           std::multiplies<double>());
            return *this;
        } else if (d1.size() == 1u) {
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(),
                           [&d1](double x) { return x * d1.m_c[0]; });
            return *this;
        } else if (this->size() == 1u) {
            double scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(), [scalar](double x) { return scalar * x; });
            return *this;
        }
        throw std::invalid_argument("Coefficients of different sizes in *");
    }
    vectorized_double &operator/=(const vectorized_double &d1)
    {
        if (d1.size() == this->size()) {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(),
                           std::divides<double>());
            return *this;
        } else if (d1.size() == 1u) {
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(),
                           [&d1](double x) { return x / d1.m_c[0]; });
            return *this;
        } else if (this->size() == 1u) {
            double scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(), [scalar](double x) { return scalar / x; });
            return *this;
        }
        throw std::invalid_argument("Coefficients of different sizes in /");
    }
    vectorized_double operator-() const
    {
        vectorized_double retval(m_c);
        transform(retval.m_c.begin(), retval.m_c.end(), retval.m_c.begin(), std::negate<double>());
        return retval;
    }
    friend bool operator==(const vectorized_double &d1, const vectorized_double &d2)
    {
        if (d1.size() == d2.size()) {
            return d1.m_c == d2.m_c;
        } else if (d1.size() == 1u) {
            return std::all_of(d2.begin(), d2.end(), [d1](double x) { return x == d1[0]; });
        } else if (d2.size() == 1u) {
            return std::all_of(d1.begin(), d1.end(), [d2](double x) { return x == d2[0]; });
        }
        return false;
    }
    friend bool operator!=(const vectorized_double &d1, const vectorized_double &d2)
    {
        return !(d1 == d2);
    }
    friend bool operator>(const vectorized_double &d1, const vectorized_double &d2)
    {
        if (d1.size() == d2.size()) {
            return d1.m_c > d2.m_c;
        } else if (d1.size() == 1u) {
            return std::all_of(d2.begin(), d2.end(), [d1](double x) { return d1[0] > x; });
        } else if (d2.size() == 1u) {
            return std::all_of(d1.begin(), d1.end(), [d2](double x) { return x > d2[0]; });
        }
        return false;
    }
    friend bool operator<(const vectorized_double &d1, const vectorized_double &d2)
    {
        if (d1.size() == d2.size()) {
            return d1.m_c < d2.m_c;
        } else if (d1.size() == 1u) {
            return std::all_of(d2.begin(), d2.end(), [d1](double x) { return d1[0] < x; });
        } else if (d2.size() == 1u) {
            return std::all_of(d1.begin(), d1.end(), [d2](double x) { return x < d2[0]; });
        }
        return false;
    }
    friend vectorized_double abs(vectorized_double in)
    {
        std::transform(in.m_c.begin(), in.m_c.end(), in.m_c.begin(), [](double x) { return std::abs(x); });
        return in;
    }
    friend std::ostream &operator<<(std::ostream &os, const vectorized_double &d)
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
    std::vector<double>::size_type size() const
    {
        return m_c.size();
    }
    std::vector<double>::const_iterator begin() const
    {
        return m_c.begin();
    }
    std::vector<double>::const_iterator end() const
    {
        return m_c.end();
    }
    std::vector<double>::iterator begin()
    {
        return m_c.begin();
    }
    std::vector<double>::iterator end()
    {
        return m_c.end();
    }
    void resize(std::vector<double>::size_type new_size)
    {
        m_c.resize(new_size);
    }
    void resize(std::vector<double>::size_type new_size, double val)
    {
        m_c.resize(new_size, val);
    }
    double operator[](const std::vector<double>::size_type idx) const
    {
        if (m_c.size() == 1u) {
            return m_c[0];
        }
        return m_c[idx];
    }
    void set_value(const std::vector<double>::size_type idx, double val)
    {
        m_c[idx] = val;
    }
    const std::vector<double> &get_v() const
    {
        return m_c;
    }

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar &m_c;
    }

    std::vector<double> m_c;
};

} // end of audi namespace

namespace piranha
{
namespace math
{

// ---------------------   impl functions needed for vectorized_double to pass piranha::is_cf type trait
template <typename T>
struct is_zero_impl<T, typename std::enable_if<std::is_same<T, audi::vectorized_double>::value>::type> {
    bool operator()(const T &v) const
    {
        return std::all_of(v.begin(), v.end(), [](double x) { return x == 0.; });
    }
};

template <typename T>
struct mul3_impl<T, typename std::enable_if<std::is_same<T, audi::vectorized_double>::value>::type> {
    /// Call operator.
    /**
     * @param out the output value.
     * @param a the first operand.
     * @param b the second operand.
     *
     * @return the output of piranha::mp_integer::mul().
     */
    void operator()(T &out, const T &a, const T &b) const
    {
        if (a.size() == b.size()) {
            if (out.size() != a.size()) {
                out.resize(a.size());
            }
            std::transform(a.begin(), a.end(), b.begin(), out.begin(), std::multiplies<double>());
            return;
        } else if (a.size() == 1u) {
            if (out.size() != b.size()) {
                out.resize(b.size());
            }
            std::transform(b.begin(), b.end(), out.begin(), [a](double item) { return item * a[0]; });
            return;
        } else if (b.size() == 1u) {
            if (out.size() != a.size()) {
                out.resize(a.size());
            }
            std::transform(a.begin(), a.end(), out.begin(), [b](double item) { return item * b[0]; });
            return;
        }
        throw std::invalid_argument("Coefficients of different sizes in mul3");
    }
};

// ------------------ impl functions needed to have the methods partial, integrate and subs
template <typename T>
struct partial_impl<T, typename std::enable_if<std::is_same<T, audi::vectorized_double>::value>::type> {
    /// Call operator.
    /**
     * @return an instance of piranha::mp_integer constructed from zero.
     */
    T operator()(const T &in, const std::string &) const
    {
        return T(std::vector<double>(in.size(), 0.));
    }
};
template <typename T, typename U>
struct pow_impl<T, U, typename std::enable_if<std::is_same<T, audi::vectorized_double>::value>::type> {
    /// Call operator.
    /**
     * @param c the input
     * @param exp the exponent
     * @return the exp operator applied to all elements of the input
     */
    T operator()(const T &c, const U &exp) const
    {
        auto retval(c);
        std::transform(retval.begin(), retval.end(), retval.begin(),
                       [exp](double x) { return piranha::math::pow(x, exp); });
        return retval;
    };
};

} // end of math namespace

template <typename Archive>
struct boost_save_impl<Archive, audi::vectorized_double> : boost_save_via_boost_api<Archive, audi::vectorized_double> {
};

template <typename Archive>
struct boost_load_impl<Archive, audi::vectorized_double> : boost_load_via_boost_api<Archive, audi::vectorized_double> {
};

template <>
struct zero_is_absorbing<audi::vectorized_double> : std::false_type {
};

} // end of piranha namespace

#undef MAX_STREAMED_COMPONENTS
#endif
