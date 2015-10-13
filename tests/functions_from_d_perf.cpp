#include "../src/gdual.hpp"
#include "../src/functions.hpp"
#include "../src/functions_from_d.hpp"

#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <piranha/settings.hpp>

#include <vector>

using namespace audi;

void scalable_test_atanh(int m, int n)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual> variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back("x"+std::to_string(i), m);
    } 
    gdual p1(0.1, m);
    gdual value(p1);
    for (int i = 0u; i < n; ++i) {p1 += variables[i];} // 1 + x1 + x2 + ...

    boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
    value=atanh(p1);
}

void scalable_test_atanh_d(int m, int n)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual> variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back("x"+std::to_string(i), m);
    } 
    gdual p1(0.1, m);
    gdual value(p1);
    for (int i = 0u; i < n; ++i) {p1 += variables[i];} // 1 + x1 + x2 + ...

    boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
    value=atanh_d(p1);
}

BOOST_AUTO_TEST_CASE(functions_from_derivative_vs_nilpotency)
{
    if (boost::unit_test::framework::master_test_suite().argc > 1) {
        piranha::settings::set_n_threads(boost::lexical_cast<unsigned>(boost::unit_test::framework::master_test_suite().argv[1u]));
    }

    unsigned int low=1, high=5;

    // we test the performance of atanh as computed by series expansion 
    std::cout << "Computing atanh(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test_atanh(m,n);
        }
    }

    // we test the performance of atanh as computed by exploting its derivative
    std::cout << "Computing atanh_d(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test_atanh_d(m,n);
        }
    }

    
}