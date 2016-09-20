#include <piranha/polynomial.hpp>
#include <piranha/type_traits.hpp>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <boost/timer/timer.hpp>
#include "src/gdual_v.hpp"



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
    using p_type_vector = piranha::polynomial<coefficient_v,piranha::monomial<char>>;
    using p_type = piranha::polynomial<double,piranha::monomial<char>>;

    auto c1  = 0.5;
    auto c2  = 0.23;
    auto c3  = -0.44;
    auto c4  = 0.1;


    unsigned int N = 100u;
    auto lb = -0.5;
    auto ub = 0.5;
    // We generate random vector coefficients
    auto xcoeff = random_vector_double(N, lb, ub);
    auto ycoeff = random_vector_double(N, lb, ub);
    auto zcoeff = random_vector_double(N, lb, ub);
    auto coeff = random_vector_double(N, lb, ub);

    {
        p_type_vector x{"x"}, y{"y"}, z{"z"};
        x *= xcoeff;
        y *= ycoeff;
        z *= zcoeff;

        auto foo = coeff + x + y + z;
        for (int i = 1; i < 15; ++i) {
            foo *= coeff + x + y + z;
        }

        boost::timer::auto_cpu_timer t;
        foo *= foo;
        // std::cout << foo << '\n';
    }
    {
        boost::timer::auto_cpu_timer t;
        for (auto i = 0u; i < N; ++i) {
            p_type x{"x"}, y{"y"}, z{"z"};
            x *= xcoeff[i];
            y *= ycoeff[i];
            z *= zcoeff[i];

            auto foo = coeff[i] + x + y + z;
            for (int i = 1; i < 15; ++i) {
                foo *= coeff[i] + x + y + z;
            }
            foo *= foo;
        }
    }

    audi::gdual_v dario({1.,2.}, "x", 3);
    std::cout << dario << std::endl;

}
