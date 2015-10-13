#include "helpers.hpp"
#include "../src/gdual.hpp"
#include "../src/functions.hpp"
#include "../src/functions_from_d.hpp"

#define BOOST_TEST_MODULE audi_functions_from_d_test
#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>

using namespace audi;

BOOST_AUTO_TEST_CASE(functions_from_derivative)
{
    {
    unsigned int order = 5;
    gdual x(0.1, "x",order);
    gdual y(-0.2, "y",order);

    auto p1 = x + y;  
    BOOST_CHECK(EPSILON_COMPARE(atanh(p1), atanh_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 6;
    gdual x(0.1, "x",order);
    gdual y(-0.2, "y",order);

    auto p1 = x + y;  
    BOOST_CHECK(EPSILON_COMPARE(atanh(p1), atanh_d(p1), 1e-12) == true);
    }

    {
    unsigned int order = 5;
    gdual x(0.1, "x",order);
    gdual y(-0.2, "y",order);
    gdual z(-0.22, "z",order);

    auto p1 = (x + y) / (10. - y - z);  
    BOOST_CHECK(EPSILON_COMPARE(atanh(p1), atanh_d(p1), 1e-12) == true);
    }
    
    {
    unsigned int order = 6;
    gdual x(0.1, "x",order);
    gdual y(-0.2, "y",order);
    gdual z(-0.22, "z",order);

    auto p1 = (x + y) / (10. - y - z);  
    BOOST_CHECK(EPSILON_COMPARE(atanh(p1), atanh_d(p1), 1e-12) == true);
    }
}