#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <stdexcept>
#include <vector>

#include <audi/audi.hpp>

using namespace audi;

template<typename T>
void test_construction() {
    // Default Constructor
    {
        vectorized<T> x;
        BOOST_CHECK(x.size() == 1);
        BOOST_CHECK(*x.begin() == 0);
    }
    // Constructor from int.
    {
        vectorized<T> x{123};
        BOOST_CHECK(x.size() == 1);
        BOOST_CHECK(*x.begin() == 123);
    }
    // Constructor from T.
    {
        vectorized<T> x{T(123.)};
        BOOST_CHECK(x.size() == 1);
        BOOST_CHECK(*x.begin() == 123.);
    }
    // Constructor from an std::vector
    {
        std::vector<T> in{T(1.), T(2.), T(3.), T(4.), T(2.2)};
        vectorized<T> x(in);
        BOOST_CHECK(x.size() == 5);
        BOOST_CHECK(std::equal(in.begin(), in.end(), x.begin()));
        std::vector<T> empty;
        // rvalue
        BOOST_CHECK_THROW(vectorized<T>(std::vector<T>{}), std::invalid_argument);
        // lvalue
        BOOST_CHECK_THROW(vectorized<T>{empty}, std::invalid_argument);
        // initializer list
        BOOST_CHECK_THROW(vectorized<T>({}), std::invalid_argument);
    }
}

template<typename T>
void test_math() {
    {
        vectorized<T> x1{T(1.), T(1.), T(2.), T(-3.), T(4.)};
        vectorized<T> x2{T(1.), T(1.), T(2.), T(3.), T(4.)};
        vectorized<T> x3{T(-100.), T(-100.), T(-100.), T(-100.), T(-100.)};
        vectorized<T> x4{T(100.), T(100.), T(100.), T(100.), T(100.)};
        BOOST_CHECK(abs(x1) == x2);
        BOOST_CHECK(x3 < x1);
        BOOST_CHECK(x1 > x3);
        BOOST_CHECK(x2 < x4);
        BOOST_CHECK(x4 > x2);
        BOOST_CHECK(x1 < 100.);
        BOOST_CHECK(-100 < x1);
        BOOST_CHECK(100. > x1);
        BOOST_CHECK(x1 > -100);
        BOOST_CHECK(x1 != x2);
        BOOST_CHECK(x1 == x1);
        BOOST_CHECK(x3 == -100);
        BOOST_CHECK(x4 == 100);
        BOOST_CHECK(x3 != 32);
        BOOST_CHECK(x4 != 32);
        BOOST_CHECK(abs(x1) == x2);
        BOOST_CHECK(abs(x2) != x1);
        BOOST_CHECK(x1 - 1. == (vectorized<T>{T(0.), T(0.), T(1.), T(-4.), T(3.)}));
    }
    {
        vectorized<T> x1{T(1.), T(2.)};
        vectorized<T> x2{T(3.), T(3.)};
        // Interoperability against arith types
        BOOST_CHECK(x1 * 2. == (vectorized<T>{T(2.), T(4.)}));
        BOOST_CHECK(x1 / 2. == (vectorized<T>{T(0.5), T(1.)}));
        BOOST_CHECK(x1 + 2. == (vectorized<T>{T(3.), T(4.)}));
        BOOST_CHECK(x1 - 2. == (vectorized<T>{T(-1.), T(0.)}));
        BOOST_CHECK(x1 * 2 == (vectorized<T>{T(2.), T(4.)}));
        BOOST_CHECK(x1 / 2 == (vectorized<T>{T(0.5), T(1.)}));
        BOOST_CHECK(x1 + 2 == (vectorized<T>{T(3.), T(4.)}));
        BOOST_CHECK(x1 - 2 == (vectorized<T>{T(-1.), T(0.)}));
        BOOST_CHECK(2 * x1 == (vectorized<T>{T(2.), T(4.)}));
        BOOST_CHECK(2 / x1 == (vectorized<T>{T(2.), T(1.)}));
        BOOST_CHECK(2 + x1 == (vectorized<T>{T(3.), T(4.)}));
        BOOST_CHECK(2 - x1 == (vectorized<T>{T(1.), T(0.)}));
        BOOST_CHECK(2 * x1 == (vectorized<T>{T(2.), T(4.)}));
        BOOST_CHECK(2 / x1 == (vectorized<T>{T(2.), T(1.)}));
        BOOST_CHECK(2 + x1 == (vectorized<T>{T(3.), T(4.)}));
        BOOST_CHECK(2 - x1 == (vectorized<T>{T(1.), T(0.)}));
        BOOST_CHECK(x1 != 1.);
        BOOST_CHECK(x2 == 3.);
        BOOST_CHECK(x1 > -1.);
        BOOST_CHECK(x1 < 3.);
        // operators with vectorized-vectorized
        BOOST_CHECK(x1 + x2 == (vectorized<T>{T(4.), T(5.)}));
        BOOST_CHECK(x1 - x2 == (vectorized<T>{T(-2.), T(-1.)}));
        BOOST_CHECK(x1 / x2 == (vectorized<T>{T(1.)/T(3.), T(2.)/T(3.)}));
        BOOST_CHECK(x1 * x2 == (vectorized<T>{T(3.), T(6.)}));
    }
}

BOOST_AUTO_TEST_CASE(construction)
{
    test_construction<double>();
    test_construction<mppp::real128>();
}
BOOST_AUTO_TEST_CASE(math)
{
    test_math<double>();
    test_math<mppp::real128>();
}
