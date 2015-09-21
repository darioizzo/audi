#include "../src/gdual.hpp"
#include "../src/functions.hpp"

#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <piranha/settings.hpp>

#include <vector>

using namespace audi;

void scalable_test_sin_and_cos(int m, int n)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual> variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back("x"+std::to_string(i), m);
    } 
    gdual p1(1, m);
    gdual sine(p1);
    gdual cosine(p1);
    for (int i = 0u; i < n; ++i) {p1 += variables[i];} // 1 + x1 + x2 + ...

    std::cout << "sin, cos: " << std::endl;   
    {
        boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
        sine=sin(p1);
        cosine=cos(p1);
    }

    std::cout << "sin_and_cos: " << std::endl;   
    {
        boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
        sin_and_cos(p1, sine, cosine);
    }
}

void scalable_test_tan(int m, int n)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual> variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back("x"+std::to_string(i), m);
    } 
    gdual p1(1, m);
    gdual tangent(p1);
    gdual sine(p1);
    gdual cosine(p1);
    for (int i = 0u; i < n; ++i) {p1 += variables[i];} // 1 + x1 + x2 + ...

    std::cout << "tan: " << std::endl;   
    {
        boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
        tangent=tan(p1);
    }

    std::cout << "sin / cos: " << std::endl;   
    {
        boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
        sin_and_cos(p1, sine, cosine);
std::cout << "HERE!" << std::endl;
        tangent = sine / cosine;
    }
}

BOOST_AUTO_TEST_CASE(trigonometry_perf)
{
    if (boost::unit_test::framework::master_test_suite().argc > 1) {
        piranha::settings::set_n_threads(boost::lexical_cast<unsigned>(boost::unit_test::framework::master_test_suite().argv[1u]));
    }

    // sin and cos
    std::cout << "Testing performance of sin_and_cos vs sin and cos on (1 + x1 + x2 + ...)" << std::endl;
    for (auto m = 10; m < 11; ++m) {
        for (auto n = 10; n < 11; ++n) {
            scalable_test_sin_and_cos(m,n);
        }
    }

    // tan
    std::cout << "Testing performance of tan vs sin/cos on (1 + x1 + x2 + ...)" << std::endl;
    for (auto m = 10; m < 11; ++m) {
        for (auto n = 10; n < 11; ++n) {
            scalable_test_tan(m,n);
        }
    }
}