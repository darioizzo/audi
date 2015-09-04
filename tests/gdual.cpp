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
	BOOST_CHECK_EQUAL(p1 + p2, p2 + p1);
	BOOST_CHECK_EQUAL(1 + p2, p2 + 1);
	BOOST_CHECK_EQUAL(1. + p2, p2 + 1.);
}

BOOST_AUTO_TEST_CASE(subtraction)
{
	gdual p1(1, 4);
	gdual p2("x", 4);
	BOOST_CHECK_EQUAL(p1 - p1, gdual(0,4));
	BOOST_CHECK_EQUAL(p2 - p2, gdual(0,4));
	BOOST_CHECK_EQUAL(p1 - p2, - (p2 - p1));
	BOOST_CHECK_EQUAL(1 + p1 - p1, gdual(1,4));
	BOOST_CHECK_EQUAL(1 + p2 - p2, gdual(1,4));
	BOOST_CHECK_EQUAL(1. + p1 - p1, gdual(1,4));
	BOOST_CHECK_EQUAL(1. + p2 - p2, gdual(1,4));
	BOOST_CHECK_EQUAL((1 - p1) + (p1 + p2), 1 + p2);
}

BOOST_AUTO_TEST_CASE(multiplication)
{
	// we test normalized derivatives computations of f = x + 3*x*y + y*y 
	gdual x("dx",2);
	gdual y("dy",2);
	x += 3;
	y += 7;
	auto p1 = x + 3*x*y + y*y;
	BOOST_CHECK_EQUAL(p1.find_cf({0,0}), 115);
	BOOST_CHECK_EQUAL(p1.find_cf({1,0}), 22);
	BOOST_CHECK_EQUAL(p1.find_cf({0,1}), 23);
	BOOST_CHECK_EQUAL(p1.find_cf({1,1}), 3);
	BOOST_CHECK_EQUAL(p1.find_cf({2,0}), 0);
	BOOST_CHECK_EQUAL(p1.find_cf({0,2}), 1);

	// we test the truncation order (5) on a 3 variables polynomial p2 = (1 + x0 + x1 + x2)^10
	std::vector<gdual> vars;
	for (auto i = 0u; i<3; ++i)
	{
		vars.emplace_back("x" + std::to_string(i), 5);
	}
	gdual p2(1, 5);
	for (auto var : vars)
	{
		p2+=var;
	}
	p2 = p2*p2*p2*p2*p2*p2;
	BOOST_CHECK_EQUAL(p2.get_order(), 5);
}

BOOST_AUTO_TEST_CASE(division)
{
	// we test normalized derivatives computations of f = 1 / (x + 2*x*y + y*y)
	{
	gdual x("dx",2);
	gdual y("dy",2);
	x += 0;
	y += 1;
	auto p1 = 1 / (x + 2*x*y + y*y);
	BOOST_CHECK_EQUAL(p1.find_cf({0,0}), 1);
	BOOST_CHECK_EQUAL(p1.find_cf({1,0}), -3);
	BOOST_CHECK_EQUAL(p1.find_cf({0,1}), -2);
	BOOST_CHECK_EQUAL(p1.find_cf({1,1}), 10);
	BOOST_CHECK_EQUAL(p1.find_cf({2,0}), 9);
	BOOST_CHECK_EQUAL(p1.find_cf({0,2}), 3);
	}

	// we test that the division between polynomials
	gdual x("dx",4);
	gdual y("dy",4);
	x += 2;
	y += 3;
	auto p1 = x*x*y + x*y*x*x*x - 3*y*y*y*y*x*y*x;
	auto p2 = x*x*y+y*y*y+x-y;
	BOOST_CHECK_NO_THROW(p1/p2);
}

BOOST_AUTO_TEST_CASE(identities)
{
	// we test some trivial identities
	gdual x("dx",3);
	gdual y("dy",3);
	x += 2;
	y += 3;
	auto p1 = x*x+y-x*x*x*x*y-y*x*x;
	auto p2 = y*y-x+y*y*y*y*x-2*x;
	BOOST_CHECK_EQUAL((x + y)*(x + y), x*x + y*y + 2*x*y);
	BOOST_CHECK_EQUAL((p1 + p2)*(p1 + p2), p1*p1 + p2*p2 + 2*p1*p2);
	BOOST_CHECK_EQUAL(x*x*x*x-y*y*y*y, (x-y)*(x+y)*(x*x+y*y));
	BOOST_CHECK_EQUAL(p1*p1*p1*p1-p2*p2*p2*p2, (p1-p2)*(p1+p2)*(p1*p1+p2*p2));

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
