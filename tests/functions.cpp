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
	BOOST_CHECK_EQUAL(pow(p1, 3), p1*p1*p1);
	BOOST_CHECK(EPSILON_COMPARE(pow(p1, 3.), p1*p1*p1) == true);
}

