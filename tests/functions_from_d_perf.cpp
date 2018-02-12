#include <audi/gdual.hpp>
#include <audi/functions.hpp>
#include <audi/functions_from_d.hpp>

#define BOOST_TEST_MODULE audi_gdual_test
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <piranha/settings.hpp>

#include <vector>

using namespace audi;

void scalable_test(int m, int n, gdual<double>(*func)(const gdual<double>& d), double param)
{
    std::cout << "Testing for order, n_vars: " << m << ",\t" << n << std::endl;
    std::vector<gdual<double> > variables;
    for (auto i = 0; i < n; ++i) {
        variables.emplace_back(0., "x"+std::to_string(i), m);
    }
    gdual<double> p1(param);
    gdual<double> value(p1);
    for (int i = 0u; i < n; ++i) {p1 += variables[i];} // 1 + x1 + x2 + ...

    boost::timer::auto_cpu_timer t; // We only time the following operation
    value=(*func)(p1);
}

BOOST_AUTO_TEST_CASE(functions_from_derivative_vs_nilpotency)
{
    if (boost::unit_test::framework::master_test_suite().argc > 1) {
        piranha::settings::set_n_threads(boost::lexical_cast<unsigned>(boost::unit_test::framework::master_test_suite().argv[1u]));
    }

    unsigned int low=9, high=10;

    // we test the performance of atanh as computed by exploting the nilpotency
    std::cout << "Computing atanh(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test(m,n,&audi::atanh, 0.1);
        }
    }

    // we test the performance of atanh as computed by exploting its derivative
    std::cout << "Computing atanh_d(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test(m,n,&audi::atanh_d, 0.1);
        }
    }

    std::cout << std::endl;

    // we test the performance of atan as computed by exploting the nilpotency
    std::cout << "Computing atan(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test(m,n,&audi::atan, 0.1);
        }
    }

    // we test the performance of atanh as computed by exploting its derivative
    std::cout << "Computing atan_d(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test(m,n,&audi::atan_d, 0.1);
        }
    }

    std::cout << std::endl;

    // we test the performance of asinh as computed by exploting the identity
    std::cout << "Computing asinh(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test(m,n,&audi::asinh, 0.1);
        }
    }

    // we test the performance of asinh as computed by exploting its derivative
    std::cout << "Computing asinh_d(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test(m,n,&audi::asinh_d, 0.1);
        }
    }

    std::cout << std::endl;

    // we test the performance of asin as computed by exploting the identity
    std::cout << "Computing asin(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test(m,n,&audi::asin, 0.1);
        }
    }

    // we test the performance of asin as computed by exploting its derivative
    std::cout << "Computing asin_d(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test(m,n,&audi::asin_d, 0.1);
        }
    }

    std::cout << std::endl;

    // we test the performance of acosh as computed by exploting the identity
    std::cout << "Computing acosh(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test(m,n,&audi::acosh, 1.1);
        }
    }

    // we test the performance of acosh as computed by exploting its derivative
    std::cout << "Computing acosh_d(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test(m,n,&audi::acosh_d, 1.1);
        }
    }

    std::cout << std::endl;

    // we test the performance of acos as computed by exploting the identity
    std::cout << "Computing acos(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test(m,n,&audi::acos, 0.1);
        }
    }

    // we test the performance of acosh as computed by exploting its derivative
    std::cout << "Computing acos_d(1 + x1 + x2 + ...): " << std::endl;
    for (auto m = low; m < high; ++m) {
        for (auto n = low; n < high; ++n) {
            scalable_test(m,n,&audi::acos_d, 0.1);
        }
    }

}
