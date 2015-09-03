#include "../src/gdual.hpp"

#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

using namespace audi;

BOOST_AUTO_TEST_CASE(addition)
{
	gdual p1(1, 4);
	gdual p2("x", 4);
	BOOST_CHECK(p1 + p2 == p2 + p1);
	BOOST_CHECK(1 + p2 == p2 + 1);
	BOOST_CHECK(1. + p2 == p2 + 1.);
}

BOOST_AUTO_TEST_CASE(subtraction)
{
	gdual p1(1, 4);
	gdual p2("x", 4);
	BOOST_CHECK(p1 - p1 == gdual(0,4));
	BOOST_CHECK(p2 - p2 == gdual(0,4));
	BOOST_CHECK(1 + p1 - p1 == gdual(1,4));
	BOOST_CHECK(1 + p2 - p2 == gdual(1,4));
	BOOST_CHECK(1. + p1 - p1 == gdual(1,4));
	BOOST_CHECK(1. + p2 - p2 == gdual(1,4));
	BOOST_CHECK((1 - p1) + (p1 + p2) == 1 + p2);
}

BOOST_AUTO_TEST_CASE(multiplication)
{
	// we test normalized derivatives computations of f = x + 3*x*y + y*y 
	gdual x("dx",2);
	gdual y("dy",2);
	x = x + 3;
	y = y + 7;
	auto result = x + 3*x*y + y*y;
	BOOST_CHECK_EQUAL(result.find_cf({0,0}), 115);
	BOOST_CHECK_EQUAL(result.find_cf({1,0}), 22);
	BOOST_CHECK_EQUAL(result.find_cf({0,1}), 23);
	BOOST_CHECK_EQUAL(result.find_cf({1,1}), 3);
	BOOST_CHECK_EQUAL(result.find_cf({2,0}), 0);
	BOOST_CHECK_EQUAL(result.find_cf({0,2}), 1);
}

BOOST_AUTO_TEST_CASE(division)
{
	// we test normalized derivatives computations of f = 1 / (x + 2*x*y + y*y)
	gdual x("dx",2);
	gdual y("dy",2);
	x = x + 0;
	y = y + 1;
	auto result = 1 / (x + 2*x*y + y*y);
	BOOST_CHECK_EQUAL(result.find_cf({0,0}), 1);
	BOOST_CHECK_EQUAL(result.find_cf({1,0}), -3);
	BOOST_CHECK_EQUAL(result.find_cf({0,1}), -2);
	BOOST_CHECK_EQUAL(result.find_cf({1,1}), 10);
	BOOST_CHECK_EQUAL(result.find_cf({2,0}), 9);
	BOOST_CHECK_EQUAL(result.find_cf({0,2}), 3);
}

BOOST_AUTO_TEST_CASE(find_cf)
{
	BOOST_CHECK_EQUAL(gdual("x",4).find_cf({5}),0.);
	BOOST_CHECK_THROW(gdual(1,4).find_cf({5}),std::invalid_argument);
	BOOST_CHECK_EQUAL(gdual("x",4).find_cf(std::vector<int>{5}),0.);
	BOOST_CHECK_THROW(gdual(1,4).find_cf(std::vector<int>{5}),std::invalid_argument);

	BOOST_CHECK_EQUAL((gdual(3,4) + gdual("x",4)).find_cf({0}),3);
	BOOST_CHECK_EQUAL((gdual(3,4) + gdual("x",4)).find_cf({1}),1);

	BOOST_CHECK_EQUAL((gdual("x",4) + gdual("y",4)).find_cf({0,1}),1);
	BOOST_CHECK_EQUAL((gdual("x",4) + gdual("y",4)).find_cf({1,0}),1);
}
