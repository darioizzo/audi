#ifndef AUDI_HPP
#define AUDI_HPP

#include <complex>

#include <audi/vectorized_double.hpp>
#include <audi/gdual.hpp>

namespace audi
{
using gdual_d = audi::gdual<double>;
using gdual_v = audi::gdual<audi::vectorized_double>;
using gdual_c = audi::gdual<std::complex<double>>;
}

#endif
