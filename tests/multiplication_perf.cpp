#include "../src/gdual.hpp"

#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <piranha/settings.hpp>

#include <vector>

using namespace audi;

void scalable_mul(int m, int n)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual> variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back("x"+std::to_string(i), m);
    } 
    gdual p1(1, m);
    gdual p2(1, m);
    for (int i = 0u; i < n; ++i) {p1 += variables[i];} // 1 + x1 + x2 + ...
    for (int i = 0u; i < n; ++i) {p2 -= variables[i];} // 1 - x1 - x2 + ...

    auto result = p1 * p2;
    for (auto i = 1; i < m-1; ++i) {
        result*=result;
    }
    {
        boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
        result*=result;
    }
}

BOOST_AUTO_TEST_CASE(multiplication_performance)
{
    if (boost::unit_test::framework::master_test_suite().argc > 1) {
        piranha::settings::set_n_threads(boost::lexical_cast<unsigned>(boost::unit_test::framework::master_test_suite().argv[1u]));
    }
    std::cout << "Testing multiplication of (1 + x1 + .. + xn)^m * (1 - x1 - .. - xn)^m: " << std::endl;
    for (auto m = 5; m < 11; ++m) {
        for (auto n = 5; n < 11; ++n) {
            scalable_mul(m,n);
        }
    }
}