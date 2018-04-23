#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "helpers.hpp"
#include <audi/vectorized_double.hpp>

using namespace audi;

BOOST_AUTO_TEST_CASE(construction)
{
    // Default Constructor
    {
        vectorized_double x;
        BOOST_CHECK(x.size() == 1);
        BOOST_CHECK(*x.begin() == 0);
    }
    // Constructor from int.
    {
        vectorized_double x{123};
        BOOST_CHECK(x.size() == 1);
        BOOST_CHECK(*x.begin() == 123);
    }
    // Constructor from double.
    {
        vectorized_double x{123.};
        BOOST_CHECK(x.size() == 1);
        BOOST_CHECK(*x.begin() == 123.);
    }
    // Constructor from an std::vector
    {
        std::vector<double> in{1, 2, 3, 4, 2.2};
        vectorized_double x(in);
        BOOST_CHECK(x.size() == 5);
        BOOST_CHECK(std::equal(in.begin(), in.end(), x.begin()));
        std::vector<double> empty;
        // rvalue
        BOOST_CHECK_THROW(vectorized_double(std::vector<double>{}), std::invalid_argument);
        // lvalue
        BOOST_CHECK_THROW(vectorized_double{empty}, std::invalid_argument);
        // initializer list
        BOOST_CHECK_THROW(vectorized_double({}), std::invalid_argument);
    }
}
BOOST_AUTO_TEST_CASE(math)
{
    vectorized_double x1{-1., 1., 2., -3., 4.};
    vectorized_double x2{1., 1., 2., 3., 4.};
    vectorized_double x3{-100., -100., -100., -100., -100.};
    vectorized_double x4{100., 100., 100., 100., 100.};
    BOOST_CHECK(abs(x1)==x2);
    BOOST_CHECK(x3 < x1);
    BOOST_CHECK(x1 > x3);
    BOOST_CHECK(x2 < x4);
    BOOST_CHECK(x4 > x2);
    BOOST_CHECK(x1 < 100.);
    BOOST_CHECK(-100 < x1);
    BOOST_CHECK(100. > x1);
    BOOST_CHECK(x1 > -100);
}