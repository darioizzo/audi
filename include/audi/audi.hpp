#ifndef AUDI_HPP
#define AUDI_HPP

#include <complex>

#include <mp++/real128.hpp>

#include <audi/vectorized_double.hpp>
#include <audi/gdual.hpp>
#include <audi/functions.hpp>
#include <audi/functions_from_d.hpp>
#include <audi/invert_map.hpp>
#include <audi/real128.hpp>

namespace audi
{
using gdual_d = audi::gdual<double>;
using gdual_v = audi::gdual<audi::vectorized_double>;
using gdual_c = audi::gdual<std::complex<double>>;
using gdual_mp = audi::gdual<mppp::real128>;
}

#endif
