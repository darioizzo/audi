#include <piranha/polynomial.hpp>
#include <piranha/type_traits.hpp>
#include <piranha/pow.hpp>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <boost/timer/timer.hpp>

// The streaming operator will only output the first MAX_STREAMED_COMPONENTS elements of the vector
#define MAX_STREAMED_COMPONENTS 5u

namespace audi { namespace detail {

// This class implements a vectorized coefficient to be used as the coefficient type in a piranha polynomial
// The coefficient is, essentially a vector of doubles [a0,a1,a2,...an] on which all arithmetic operations and
// function calls operate element-wise.

class coefficient_v
{
public:
    // Default constructor. Constructs [0.]
    coefficient_v() : m_c({0.}) {};
    // Constructor from int. Its mandatory for piranha::polynomial coefficient
    coefficient_v(int a) : m_c({static_cast<double>(a)}) {};
    // Constructor from double value a. Construct [a]
    coefficient_v(double a) : m_c({a}) {};
    // Constructor from an std::vector
    coefficient_v(const std::vector<double> &c) : m_c(c) {
        if (m_c.size() == 0) {
            throw std::invalid_argument("Cannot build an empty coeffciient_v");
        }
    };
    // Constructor from an std::vector r value
    coefficient_v(std::vector<double> && c) : m_c(c) {
        if (m_c.size() == 0) {
            throw std::invalid_argument("Cannot build an empty coeffciient_v");
        }
    };
    coefficient_v(std::initializer_list<double> c) : m_c(c) {
        if (m_c.size() == 0) {
            throw std::invalid_argument("Cannot build an empty coeffciient_v");
        }
    };
    // ------------------- Binary arithmetic operators implemented using +=,-=, etc.
    friend coefficient_v operator+(const coefficient_v &d1, const coefficient_v &d2)
    {
        coefficient_v retval(d1);
        retval += d2;
        return retval;
    };
    friend coefficient_v operator-(const coefficient_v &d1, const coefficient_v &d2)
    {
        coefficient_v retval(d1);
        retval -= d2;
        return retval;
    };
    friend coefficient_v operator*(const coefficient_v &d1, const coefficient_v &d2)
    {
        coefficient_v retval(d1);
        retval *= d2;
        return retval;
    };
    friend coefficient_v operator/(const coefficient_v &d1, const coefficient_v &d2)
    {
        coefficient_v retval(d1);
        retval /= d2;
        return retval;
    };

    // ----------------- Juice implementation of the operators. It also deals with the case [b1] op [a1,a2,..an] to
    // take care of scalar multiplication/division etc.
    coefficient_v& operator+=(const coefficient_v &d1)
    {
        if (d1.size() == this->size())
        {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(), std::plus<double>());
            return *this;
        }
        else if (d1.size() == 1u) {
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(), [&d1](double x){return x + d1.m_c[0];});
            return *this;
        }
        else if (this->size() == 1u) {
            double scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(), [scalar](double x){return x + scalar;});
            return *this;
        }
        throw std::invalid_argument("Coefficients of different sizes in +");
    }
    coefficient_v& operator-=(const coefficient_v &d1)
    {
        if (d1.size() == this->size())
        {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(), std::minus<double>());
            return *this;
        }
        else if (d1.size() == 1u) {
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(), [&d1](double x){return x - d1.m_c[0];});
            return *this;
        }
        else if (this->size() == 1u) {
            double scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(), [scalar](double x){return scalar - x;});
            return *this;
        }
        throw std::invalid_argument("Coefficients of different sizes in -");
    }
    coefficient_v& operator*=(const coefficient_v &d1)
    {
        if (d1.size() == this->size())
        {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(), std::multiplies<double>());
            return *this;
        }
        else if (d1.size() == 1u) {
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(), [&d1](double x){return x * d1.m_c[0];});
            return *this;
        }
        else if (this->size() == 1u) {
            double scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(), [scalar](double x){return scalar * x;});
            return *this;
        }
        throw std::invalid_argument("Coefficients of different sizes in *");
    }
    coefficient_v& operator/=(const coefficient_v &d1)
    {
        if (d1.size() == this->size())
        {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(), std::divides<double>());
            return *this;
        }
        else if (d1.size() == 1u) {
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(), [&d1](double x){return x / d1.m_c[0];});
            return *this;
        }
        else if (this->size() == 1u) {
            double scalar = m_c[0];
            this->resize(d1.size());
            std::transform(d1.m_c.begin(), d1.m_c.end(), this->m_c.begin(), [scalar](double x){return scalar / x;});
            return *this;
        }
        throw std::invalid_argument("Coefficients of different sizes in /");
    }
    coefficient_v operator-() const
    {
        coefficient_v retval(m_c);
        transform (retval.m_c.begin(), retval.m_c.end(), retval.m_c.begin(), std::negate<double>());
        return retval;
    }
    friend bool operator==(const coefficient_v &d1, const coefficient_v &d2)
    {
        if (d1.size() == d2.size())
        {
            return d1.m_c == d2.m_c;
        }
        else if (d1.size() == 1u) {
            return std::all_of(d2.begin(),d2.end(),[d1](double x) {return x == d1[0];});
        }
        else if (d2.size() == 1u) {
            return std::all_of(d1.begin(),d1.end(),[d2](double x) {return x == d2[0];});
        }
        return false;
    }
    friend bool operator!=(const coefficient_v &d1, const coefficient_v &d2)
    {
        return !(d1 == d2);
    }
    friend std::ostream &operator<<(std::ostream &os, const coefficient_v &d) {
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
            os  << "... ]";
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
    void resize(std::vector<double>::size_type new_size) {
        m_c.resize(new_size);
    }
    double operator[] (const std::vector<double>::size_type idx) const {
        if (m_c.size() == 1u) {
            return m_c[0];
        }
        return m_c[idx];
    }
private:
    std::vector<double> m_c;
};

}} // end of audi::detail namespace

namespace piranha { namespace math {

// ---------------------   impl functions needed for coefficient_v to pass piranha::is_cf type trait
template <typename T>
struct is_zero_impl<T,typename std::enable_if<std::is_same<T,audi::detail::coefficient_v>::value>::type>
{
  bool operator()(const T &v) const
  {
      return std::all_of(v.begin(),v.end(),[](double x) {return x == 0.;});
  }
};

template <typename T>
struct mul3_impl<T, typename std::enable_if<std::is_same<T,audi::detail::coefficient_v>::value>::type> {
    /// Call operator.
    /**
     * @param[out] out the output value.
     * @param[in] a the first operand.
     * @param[in] b the second operand.
     *
     * @return the output of piranha::mp_integer::mul().
     */
    void operator()(T &out, const T &a, const T &b) const
    {
        if (a.size() == b.size())
        {
            if (out.size() != a.size()) {
                out.resize(a.size());
            }
            std::transform(a.begin(), a.end(), b.begin(), out.begin(), std::multiplies<double>());
            return;
        }
        else if (a.size() == 1u) {
            if (out.size() != b.size()) {
                out.resize(b.size());
            }
            std::transform(b.begin(), b.end(), out.begin(), [a](double item){return item * a[0];});
            return;
        }
        else if (b.size() == 1u) {
            if (out.size() != a.size()) {
                out.resize(a.size());
            }
            std::transform(a.begin(), a.end(), out.begin(), [b](double item){return item * b[0];});
            return;
        }
        throw std::invalid_argument("Coefficients of different sizes in mul3");
    }
};

// ------------------ impl functions needed to have the methods partial, integrate and subs
template <typename T>
struct partial_impl<T, typename std::enable_if<std::is_same<T,audi::detail::coefficient_v>::value>::type> {
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
struct pow_impl<T,U,typename std::enable_if<std::is_same<T,audi::detail::coefficient_v>::value>::type>
{
  audi::detail::coefficient_v operator()(const audi::detail::coefficient_v &c, const U &exp) const
  {
    auto retval(c);
    std::transform(retval.begin(),retval.end(),retval.begin(),[exp](double x) {return piranha::math::pow(x,exp);});
    return retval;
  };
};

}} // end of piranha::math namespace

#undef MAX_STREAMED_COMPONENTS
