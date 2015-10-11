#include "../src/gdual.hpp"
#include "../src/functions.hpp"

#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <piranha/settings.hpp>

#include <vector>

using namespace audi;

void scalable_test_sin_cos(int m, int n)
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

    boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
    sine=sin(p1);
    cosine=cos(p1);
}

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

    boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
    sin_and_cos(p1, sine, cosine);
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

    boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
    tangent=tan(p1);
}

void scalable_test_sin_over_cos(int m, int n)
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

    sin_and_cos(p1, sine, cosine);
    boost::timer::auto_cpu_timer t; // We only time the time cost of the following operation
    tangent = sine / cosine;
}

BOOST_AUTO_TEST_CASE(trigonometry_perf)
{
    if (boost::unit_test::framework::master_test_suite().argc > 1) {
        piranha::settings::set_n_threads(boost::lexical_cast<unsigned>(boost::unit_test::framework::master_test_suite().argv[1u]));
    }

    unsigned int low=10, high=12;

    // sin and cos
    std::cout << "Computing sin(1 + x1 + x2 + ...) and cos(1 + x1 + x2 + ...) separately: " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test_sin_cos(m,n);
        }
    }

    // sin and cos together
    std::cout << "\nComputing sin(1 + x1 + x2 + ...) and cos(1 + x1 + x2 + ...) at once: " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test_sin_and_cos(m,n);
        }
    }

    // tan
    std::cout << "\nTesting performance of tan(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test_tan(m,n);
        }
    }

    // sin / cos
    std::cout << "\nTesting division of sin(1 + x1 + x2 + ...) / cos(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test_sin_over_cos(m,n);
        }
    }
}