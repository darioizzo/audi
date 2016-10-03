#include <cmath>
#include <vector>

#include "expose_gdual.hpp"
#include "../src/audi.hpp"

#if defined(__clang__) || defined(__GNUC__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wpedantic"
    #pragma GCC diagnostic ignored "-Wshadow"
    #pragma GCC diagnostic ignored "-Wsign-conversion"
    #pragma GCC diagnostic ignored "-Wdeprecated"
#endif

#include "pybind11/include/pybind11/pybind11.h"

#if defined(__clang__) || defined(__GNUC__)
    #pragma GCC diagnostic pop
#endif

using namespace audi;
using namespace pyaudi;
namespace py = pybind11; 

PYBIND11_PLUGIN(_core) {
    py::module m("_core", "pyaudi's core module");

    // We expose the gdual<double> using the expose_gdual defined in exposed_gdual.hpp as pyaudi::expose_gdual
    expose_gdual<double>(m, "double");

    // Similarly, we expose the gdual<vectorized_double> and we add two custom constructors to allow constructing it from lists
    auto a = expose_gdual<vectorized_double>(m, "vdouble");
    a.def(py::init<std::vector<double>>())
    .def(py::init<std::vector<double>, const std::string &, unsigned int>());

    // We expose the functions
    m.def("exp",[](const gdual_d &d) {return exp(d);},"Exponential (gdual_d).");
    m.def("exp",[](double x) {return std::exp(x);},"Exponential (double).");
    m.def("exp",[](const gdual_v &d) {return exp(d);},"Exponential (gdual_v).");

    m.def("log",[](const gdual_d &d) {return log(d);},"Natural logarithm (gdual_d).");
    m.def("log",[](double x) {return std::log(x);},"Natural logarithm (double).");
    m.def("log",[](const gdual_v &d) {return log(d);},"Natural logarithm (gdual_v).");

    m.def("sqrt",[](const gdual_d &d) {return sqrt(d);},"Square root (gdual_d).");
    m.def("sqrt",[](double x) {return std::sqrt(x);},"Square root (double).");
    m.def("sqrt",[](const gdual_v &d) {return sqrt(d);},"Square root (gdual_v).");

    m.def("cbrt",[](const gdual_d &d) {return cbrt(d);},"Cubic root (gdual_d).");
    m.def("cbrt",[](double x) {return std::cbrt(x);},"Cubic root (double).");
    m.def("cbrt",[](const gdual_d &d) {return cbrt(d);},"Cubic root (gdual_v).");

    m.def("sin",[](const gdual_d &d) {return sin(d);},"Sine (gdual_d).");
    m.def("sin",[](double x) {return std::sin(x);},"Sine (double).");
    m.def("sin",[](const gdual_v &d) {return sin(d);},"Sine (gdual_v).");
    // m.def("sin",py::vectorize([](double x) {return std::sin(x);}),"Sine (vectorized double).");

    m.def("asin",[](const gdual_d &d) {return asin(d);},"Arc sine (gdual_d).");
    m.def("asin",[](double x) {return std::asin(x);},"Arc sine (double).");
    m.def("asin",[](const gdual_v &d) {return asin(d);},"Arc sine (gdual_v).");

    m.def("cos",[](const gdual_d &d) {return cos(d);},"Cosine (gdual_d).");
    m.def("cos",[](double x) {return std::cos(x);},"Cosine (double).");
    m.def("cos",[](const gdual_v &d) {return cos(d);},"Cosine (gdual_v).");

    m.def("acos",[](const gdual_d &d) {return acos(d);},"Arc cosine (gdual_d).");
    m.def("acos",[](double x) {return std::acos(x);},"Arc cosine (double).");
    m.def("acos",[](const gdual_v &d) {return acos(d);},"Arc cosine (gdual_v).");

    m.def("sin_and_cos",[](const gdual_d &d) {return sin_and_cos(d);}, "Sine and Cosine at once (gdual_d).");
    m.def("sin_and_cos",[](const gdual_v &d) {return sin_and_cos(d);}, "Sine and Cosine at once (gdual_v).");

    m.def("tan",[](const gdual_d &d) {return tan(d);},"Tangent (gdual_d).");
    m.def("tan",[](double x) {return std::tan(x);},"Tangent (double).");
    m.def("tan",[](const gdual_v &d) {return tan(d);},"Tangent (gdual_v).");

    m.def("atan",[](const gdual_d &d) {return atan(d);},"Arc tangent (gdual_d).");
    m.def("atan",[](double x) {return std::atan(x);},"Arc tangent (double).");
    m.def("atan",[](const gdual_v &d) {return atan(d);},"Arc tangent (gdual_v).");

    m.def("sinh",[](const gdual_d &d) {return sinh(d);},"Hyperbolic sine (gdual_d).");
    m.def("sinh",[](double x) {return std::sinh(x);},"Hyperbolic sine (double).");
    m.def("sinh",[](const gdual_v &d) {return sinh(d);},"Hyperbolic sine (gdual_v).");

    m.def("asinh",[](const gdual_d &d) {return asinh(d);},"Inverse hyperbolic sine (gdual_d).");
    m.def("asinh",[](double x) {return std::asinh(x);},"Inverse hyperbolic sine (double).");
    m.def("asinh",[](const gdual_v &d) {return asinh(d);},"Inverse hyperbolic sine (gdual_v).");

    m.def("cosh",[](const gdual_d &d) {return cosh(d);},"Hyperbolic cosine (gdual_d).");
    m.def("cosh",[](double x) {return std::cosh(x);},"Hyperbolic cosine (double).");
    m.def("cosh",[](const gdual_v &d) {return cosh(d);},"Hyperbolic cosine (gdual_v).");

    m.def("acosh",[](const gdual_d &d) {return acosh(d);},"Inverse hyperbolic cosine (gdual_d).");
    m.def("acosh",[](double x) {return std::acosh(x);},"Inverse hyperbolic cosine (double).");
    m.def("acosh",[](const gdual_v &d) {return acosh(d);},"Inverse hyperbolic cosine (gdual_v).");

    m.def("sinh_and_cosh",[](const gdual_d &d) {return sinh_and_cosh(d);} ,"Hyperbolic sine and hyperbolic cosine at once (gdual_d).");
    m.def("sinh_and_cosh",[](const gdual_v &d) {return sinh_and_cosh(d);} ,"Hyperbolic sine and hyperbolic cosine at once (gdual_v).");

    m.def("tanh",[](const gdual_d &d) {return tanh(d);},"Hyperbolic tangent (gdual_d).");
    m.def("tanh",[](double x) {return std::tanh(x);},"Hyperbolic tangent (double).");
    m.def("tanh",[](const gdual_v &d) {return tanh(d);},"Hyperbolic tangent (gdual_v).");

    m.def("atanh",[](const gdual_d &d) {return atanh(d);},"Inverse hyperbolic arc tangent (gdual_d).");
    m.def("atanh",[](double x) {return std::atanh(x);},"Inverse hyperbolic arc tangent (double).");
    m.def("atanh",[](const gdual_v &d) {return atanh(d);},"Inverse hyperbolic arc tangent (gdual_v).");

    m.def("abs",[](const gdual_d &d) {return abs(d);},"Absolute value (gdual_d).");
    m.def("abs",[](double x) {return std::abs(x);},"Absolute value (double).");
    m.def("abs",[](const gdual_v &d) {return abs(d);},"Absolute value (gdual_v).");

    m.def("erf",[](const gdual_d &d) {return erf(d);},"Error function (gdual_d).");
    m.def("erf",[](double x) {return std::erf(x);},"Error function (double).");
    m.def("erf",[](const gdual_v &d) {return erf(d);},"Error function (gdual_v).");

    return m.ptr();
}
