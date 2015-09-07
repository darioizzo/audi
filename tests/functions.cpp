#include "helpers.hpp"
#include "../src/gdual.hpp"
#include "../src/functions.hpp"

#define BOOST_TEST_MODULE audi_functions_test
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

using namespace audi;

BOOST_AUTO_TEST_CASE(exponentiation)
{
    {
        gdual x("dx",3);
        gdual y("dy",3);

        auto p1 = x*x*y + x*y*x*x*x - 3*y*y*y*y*x*y*x + 3.2;
        BOOST_CHECK_EQUAL(pow(p1, 3), p1*p1*p1);									// calls pow(gdual, int)
        BOOST_CHECK(EPSILON_COMPARE(pow(p1, 3.), p1*p1*p1) == true); 				// calls pow(gdual, double)
        BOOST_CHECK(EPSILON_COMPARE(pow(p1, gdual(3,3)), p1*p1*p1) == true); 		// calls pow(gdual, gdual) with an integer exponent
        BOOST_CHECK(EPSILON_COMPARE(pow(p1, gdual(3.1,3)), pow(p1,3.1)) == true);	// calls pow(gdual, gdual) with an real exponent
    }

    gdual x("dx",3);
    gdual y("dy",3);
    gdual p1 = x+y-3*x*y+y*y;
    gdual p2 = p1 - 3.5;

    BOOST_CHECK(EPSILON_COMPARE(pow(3, gdual(3.2, 3)), std::pow(3, 3.2) * gdual(1, 3)) == true);
    BOOST_CHECK(EPSILON_COMPARE(pow(p2, 3), pow(p2, 3.)) == true);

    BOOST_CHECK_THROW(pow(p1, 3.1), std::domain_error); 			// zero p0
    BOOST_CHECK_THROW(pow(p2, 3.1), std::domain_error); 			// negative p0
    BOOST_CHECK_THROW(pow(p1, gdual(3.1,3)), std::domain_error); 	// zero p0
    BOOST_CHECK_THROW(pow(p2, gdual(3.1,3)), std::domain_error); 	// negative p0
}

BOOST_AUTO_TEST_CASE(square_root)
{
    gdual x("dx",3);
    gdual y("dy",3);
    x += 2.3;
    y += 1.5;
    auto p1 = x*x*y - x*y*x*x*x + 3*y*y*y*y*x*y*x;	// positive p0
    auto p2 = x*x*y - x*y*x*x*x - 3*y*y*y*y*x*y*x;	// negative coefficient
    BOOST_CHECK(EPSILON_COMPARE(sqrt(p1) * sqrt(p1), p1) == true);
    BOOST_CHECK_THROW(sqrt(p2), std::domain_error); // negative p0
}

BOOST_AUTO_TEST_CASE(cubic_root)
{
    gdual x("dx",3);
    gdual y("dy",3);
    x += 2.3;
    y += 1.5;
    auto p1 = x*x*y - x*y*x*x*x + 3*y*y*y*y*x*y*x;	// positive p0
    auto p2 = x*x*y - x*y*x*x*x - 3*y*y*y*y*x*y*x;	// negative coefficient
    BOOST_CHECK(EPSILON_COMPARE(cbrt(p1) * cbrt(p1) * cbrt(p1), p1) == true);
    BOOST_CHECK(EPSILON_COMPARE(cbrt(p2) * cbrt(p2) * cbrt(p2), p2) == true);
}

BOOST_AUTO_TEST_CASE(exponential)
{
    // we test the result of f = exp((1+dx)*dy)
    gdual x("dx",2);
    gdual y("dy",2);
    x += 1;
    y += 0;
    auto p1 = exp(x*y);
    BOOST_CHECK_EQUAL(p1.find_cf({0,0}), 1);
    BOOST_CHECK_EQUAL(p1.find_cf({1,0}), 0);
    BOOST_CHECK_EQUAL(p1.find_cf({0,1}), 1);
    BOOST_CHECK_EQUAL(p1.find_cf({1,1}), 1);
    BOOST_CHECK_EQUAL(p1.find_cf({2,0}), 0);
    BOOST_CHECK_EQUAL(p1.find_cf({0,2}), 0.5);
}

BOOST_AUTO_TEST_CASE(logarithm)
{
    {
        gdual x("dx",3);
        gdual y("dy",3);
        x += 2.3;
        y += 1.5;
        auto p1 = x*x*y - x*y*x*x*x + 3*y*y*y*y*x*y*x;
        BOOST_CHECK(EPSILON_COMPARE(exp(log(p1)),p1) == true);
    }

    gdual x("dx",3);
    gdual y("dy",3);
    gdual p1 = x+y-3*x*y+y*y;
    gdual p2 = p1 - 3.5;

    BOOST_CHECK_THROW(log(p1), std::domain_error); // zero p0
    BOOST_CHECK_THROW(log(p2), std::domain_error); // negative p0
}

BOOST_AUTO_TEST_CASE(sin_and_cos)
{
    gdual x("dx",3);
    gdual y("dy",3);
    x += 2.3;
    y += 1.5;
    auto p1 = x*x*y - x*y*x*x*x + 3*y*y*y*y*x*y*x;	// positive p0
	BOOST_CHECK(EPSILON_COMPARE(sin(2 * p1),2 * sin(p1) * cos(p1)) == true);
}