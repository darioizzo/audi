#define BOOST_TEST_MODULE gdual_test
#include <boost/test/unit_test.hpp>

#include "../src/audi.hpp"

using namespace audi;

BOOST_AUTO_TEST_CASE(add)
{
	piranha::environment env;
	gdual p1(1, 4);
	gdual p2("x", 4);
	BOOST_CHECK(p1 + p2 == p2 + p1);
	BOOST_CHECK(1 + p2 == p2 + 1);
	BOOST_CHECK(1. + p2 == p2 + 1.);
}

BOOST_AUTO_TEST_CASE(sub)
{
	piranha::environment env;
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

BOOST_AUTO_TEST_CASE(find_cf)
{
	piranha::environment env;
	BOOST_CHECK_THROW(gdual("x",4).find_cf({5}),std::invalid_argument);
	BOOST_CHECK_THROW(gdual(1,4).find_cf({5}),std::invalid_argument);
	BOOST_CHECK_THROW(gdual("x",4).find_cf(std::vector<int>{5}),std::invalid_argument);
	BOOST_CHECK_THROW(gdual(1,4).find_cf(std::vector<int>{5}),std::invalid_argument);

	BOOST_CHECK_EQUAL((gdual(3,4) + gdual("x",4)).find_cf({0}),3);
	BOOST_CHECK_EQUAL((gdual(3,4) + gdual("x",4)).find_cf({1}),1);

	BOOST_CHECK_EQUAL((gdual("x",4) + gdual("y",4)).find_cf({0,1}),1);
	BOOST_CHECK_EQUAL((gdual("x",4) + gdual("y",4)).find_cf({1,0}),1);
}