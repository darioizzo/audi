#include <cmath>
#include <pybind11/pybind11.h>

#include "expose_gdual.hpp"

using namespace audi;
namespace py = pybind11;

PYBIND11_MODULE(core, m)
{
    // We expose simple mathematical functions to work on python floats
    m.def("exp", +[](double x) { return std::exp(x); }, "Exponential (double).");
    m.def("log", +[](double x) { return std::log(x); }, "Natural logarithm (double).");
    m.def("sqrt", +[](double x) { return std::sqrt(x); }, "Square root (double).");
    m.def("cbrt", +[](double x) { return std::cbrt(x); }, "Cubic root (double).");
    m.def("sin", +[](double x) { return std::sin(x); }, "Sine (double).");
    m.def("asin", +[](double x) { return std::asin(x); }, "Arc sine (double).");
    m.def("cos", +[](double x) { return std::cos(x); }, "Cosine (double).");
    m.def("acos", +[](double x) { return std::acos(x); }, "Arc cosine (double).");
    m.def("tan", +[](double x) { return std::tan(x); }, "Tangent (double).");
    m.def("atan", +[](double x) { return std::atan(x); }, "Arc tangent (double).");
    m.def("sinh", +[](double x) { return std::sinh(x); }, "Hyperbolic sine (double).");
    m.def("asinh", +[](double x) { return std::asinh(x); }, "Inverse hyperbolic sine (double).");
    m.def("cosh", +[](double x) { return std::cosh(x); }, "Hyperbolic cosine (double).");
    m.def("acosh", +[](double x) { return std::acosh(x); }, "Inverse hyperbolic cosine (double).");
    m.def("tanh", +[](double x) { return std::tanh(x); }, "Hyperbolic tangent (double).");
    m.def("atanh", +[](double x) { return std::atanh(x); }, "Inverse hyperbolic arc tangent (double).");
    m.def("abs", +[](double x) { return std::abs(x); }, "Absolute value (double).");
    m.def("erf", +[](double x) { return std::erf(x); }, "Error function (double).");

    // We expose the class gdual<double> or gdual_d
    pyaudi::expose_gdual<double>(m, "double");

}