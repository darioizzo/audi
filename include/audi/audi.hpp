#ifndef AUDI_HPP
#define AUDI_HPP

#include "detail/complex.hpp"
#include "functions.hpp"
#include "functions_from_d.hpp"
#include "gdual.hpp"

namespace audi
{
using gdual_d = audi::gdual<double>;
using gdual_v = audi::gdual<audi::vectorized_double>;
using gdual_c = audi::gdual<std::complex<double>>;
}

#endif
