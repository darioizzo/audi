#include "helpers.hpp"
#include "../src/gdual.hpp"
#include "../src/functions.hpp"
#include "../src/functions_from_d.hpp"

#define BOOST_TEST_MODULE audi_functions_floating_precision_test
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

#include <stdexcept>
#include <vector>

using namespace audi;

BOOST_AUTO_TEST_CASE(division_precision)
{

    {
    gdual x(0.1, "x",4);
    gdual y(0.13, "y",4);

    auto p1 = (x*x*y - x*y*x*x*x + 3*y*y*y*y*x*y*x);
    auto p2 = gdual(1.) / p1;
    BOOST_CHECK(EPSILON_COMPARE(p1 * p2, gdual(1), 1e-10) == true);
std::cout << p1 * p2 << std::endl;
    }

}
