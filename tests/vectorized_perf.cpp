#include <piranha/polynomial.hpp>
#include <piranha/type_traits.hpp>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <boost/timer/timer.hpp>

#include "../src/audi.hpp"


std::vector<double> random_vector_double(unsigned int N, double lb, double ub)
{
    std::vector<double> retval(N);
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(lb, ub);
    for (auto &x : retval) {
        x = dist(mt);
    }
    return retval;
}

using namespace piranha;

int main() {
    unsigned int N = 10000u;
    auto lb = -0.5;
    auto ub = 0.5;
    // We generate random vector coefficients
    auto xcoeff = random_vector_double(N, lb, ub);
    auto ycoeff = random_vector_double(N, lb, ub);
    auto zcoeff = random_vector_double(N, lb, ub);
    auto coeff = random_vector_double(N, lb, ub);
    auto coeff_v = audi::gdual_v(coeff);
    const auto order = 2u;

    {
        boost::timer::auto_cpu_timer t;
        audi::gdual_v x{xcoeff, "x",order}, y{ycoeff, "y",order}, z{zcoeff, "z",order}, w{zcoeff, "w",order}, q{zcoeff, "q",order}, r{zcoeff, "r",order}, s{zcoeff, "s",order} ;
        auto foo = coeff_v + x + y + z - w + q + r + s;
        for (int i = 1; i < 25; ++i) {
            foo *= coeff_v + x + y + z - w + q + r + s;
        }
        foo *= foo;
        sin(foo);
        //std::cout << foo.get_derivative({0,1,0,0,0,0,0}) << std::endl;
    }
    {
        boost::timer::auto_cpu_timer t;
        for (auto i = 0u; i < N; ++i) {
            audi::gdual_d x{xcoeff[i], "x",order}, y{ycoeff[i], "y",order}, z{zcoeff[i], "z",order}, w{zcoeff[i], "w",order}, q{zcoeff[i], "q",order}, r{zcoeff[i], "r",order}, s{zcoeff[i], "s",order} ;
            auto foo = coeff[i] + x + y + z - w + q + r + s;
            for (int j = 1; j < 25; ++j) {
                foo *= coeff[i] + x + y + z - w + q + r + s;
            }
            foo *= foo;
            sin(foo);
            //if (i < 5) {
            //    std::cout << foo.get_derivative({0,1,0,0,0,0,0}) << ", ";
            //}
        }
    }

}
