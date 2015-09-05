#include "helpers.hpp"
#include "../src/gdual.hpp"
#include "../src/functions.hpp"

#define BOOST_TEST_MODULE audi_functions_test
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

using namespace audi;

BOOST_AUTO_TEST_CASE(power)
{
	gdual x("dx",3);
	gdual y("dy",3);
	x += 2;
	y += 3;
	auto p1 = x*x*y + x*y*x*x*x - 3*y*y*y*y*x*y*x;
	BOOST_CHECK_EQUAL(pow(p1, 3), p1*p1*p1);						// calls pow(gdual, int)
	BOOST_CHECK(EPSILON_COMPARE(pow(p1, 3.), p1*p1*p1) == true); 	// calls pow(gdual, double)
}

BOOST_AUTO_TEST_CASE(square_root)
{
	gdual x("dx",3);
	gdual y("dy",3);
	x += 2.3;
	y += 1.5;
	auto p1 = x*x*y + x*y*x*x*x - 3*y*y*y*y*x*y*x;
	BOOST_CHECK(EPSILON_COMPARE(sqrt(p1) * sqrt(p1), p1) == true); 	
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
