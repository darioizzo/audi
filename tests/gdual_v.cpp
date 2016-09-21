#include "../src/gdual_v.hpp"
#include "helpers.hpp"

#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/test/unit_test.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <stdexcept>
#include <vector>

using namespace audi;

BOOST_AUTO_TEST_CASE(construction)
{
    // Constructing a constant with a symbol
    gdual_v x({0.});
    std::cout << x << "\n";
    gdual_v y({1.2, 2.2});
    std::cout << y << "\n";
    gdual_v z({});
    std::cout << z << "\n";
    BOOST_CHECK_EQUAL(x.get_order(), 0);

    // TODO: other constructors
}
