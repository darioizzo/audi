#include <piranha/polynomial.hpp>
#include <piranha/type_traits.hpp>
#include <piranha/pow.hpp>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <boost/timer/timer.hpp>


class coefficient_v
{
public:
    coefficient_v() : m_c({0.}) {};
    coefficient_v(int a) : m_c({static_cast<double>(a)}) {};
    coefficient_v(double a) : m_c({a}) {};
    coefficient_v(const std::vector<double> &c) : m_c(c) {
        if (c.size() == 0) {
            throw std::invalid_argument("A coefficient of size zero is detected: NAAAA");
        }
    };

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
    coefficient_v& operator+=(const coefficient_v &d1)
    {
        if (d1.size() == this->size())
        {
            std::transform(this->m_c.begin(), this->m_c.end(), d1.m_c.begin(), this->m_c.begin(), std::plus<double>());
            return *this;
        }
        else if (d1.size() == 1u) {
            return *this;
            std::transform(this->m_c.begin(), this->m_c.end(), this->m_c.begin(), [&d1](double x){return x + d1.m_c[0];});
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
    friend bool operator==(const coefficient_v &d1, const coefficient_v &d2)
    {
        return d1.m_c == d2.m_c;
    }
    friend bool operator!=(const coefficient_v &d1, const coefficient_v &d2)
    {
        return d1.m_c != d2.m_c;
    }
    friend std::ostream &operator<<(std::ostream &os, const coefficient_v &d) {
        os << "[";
        for (auto i = 0u; i < d.m_c.size() - 1u; ++i) {
            os << d.m_c[i] << ", ";
        }
        os << d.m_c[d.m_c.size() - 1u] << "]";
        return os;
    }
    coefficient_v operator-() const
    {
        coefficient_v retval(m_c);
        transform (retval.m_c.begin(), retval.m_c.end(), retval.m_c.begin(), std::negate<double>());
        return retval;
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
        return m_c[idx];
    }
private:
    std::vector<double> m_c;
};

namespace piranha { namespace math {

template <typename T>
struct is_zero_impl<T,typename std::enable_if<std::is_same<T,coefficient_v>::value>::type>
{
  bool operator()(const T &v) const
  {
      return std::all_of(v.begin(),v.end(),[](double x) {return x == 0.;});
  }
};

template <typename T>
struct mul3_impl<T, typename std::enable_if<std::is_same<T,coefficient_v>::value>::type> {
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

template <typename T>
struct partial_impl<T, typename std::enable_if<std::is_same<T,coefficient_v>::value>::type> {
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
struct pow_impl<T,U,typename std::enable_if<std::is_same<T,coefficient_v>::value && is_exponentiable<double,U>::value>::type>
{
  coefficient_v operator()(const my_coefficient &c, const U &exp) const
  {
    auto retval(c);
    std::transform(retval.begin(),retval.end(),retval.begin(),[exp](double x) {return piranha::math::pow(x,exp);});
    return retval;
  };
};


}}
