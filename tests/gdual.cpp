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