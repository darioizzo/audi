#include <audi/audi.hpp>
#include <cmath>
#include <pybind11/pybind11.h>

#include "expose_gdual.hpp"

using namespace audi;
namespace py = pybind11;

PYBIND11_MODULE(core, m)
{
    // We expose simple mathematical functions to work on python floats
    m.def("exp", [](double x) { return std::exp(x); }, "Exponential (double).");
    m.def("log", [](double x) { return std::log(x); }, "Natural logarithm (double).");
    m.def("sqrt", [](double x) { return std::sqrt(x); }, "Square root (double).");
    m.def("cbrt", [](double x) { return std::cbrt(x); }, "Cubic root (double).");
    m.def("sin", [](double x) { return std::sin(x); }, "Sine (double).");
    m.def("asin", [](double x) { return std::asin(x); }, "Arc sine (double).");
    m.def("cos", [](double x) { return std::cos(x); }, "Cosine (double).");
    m.def("acos", [](double x) { return std::acos(x); }, "Arc cosine (double).");
    m.def("tan", [](double x) { return std::tan(x); }, "Tangent (double).");
    m.def("atan", [](double x) { return std::atan(x); }, "Arc tangent (double).");
    m.def("sinh", [](double x) { return std::sinh(x); }, "Hyperbolic sine (double).");
    m.def("asinh", [](double x) { return std::asinh(x); }, "Inverse hyperbolic sine (double).");
    m.def("cosh", [](double x) { return std::cosh(x); }, "Hyperbolic cosine (double).");
    m.def("acosh", [](double x) { return std::acosh(x); }, "Inverse hyperbolic cosine (double).");
    m.def("tanh", [](double x) { return std::tanh(x); }, "Hyperbolic tangent (double).");
    m.def("atanh", [](double x) { return std::atanh(x); }, "Inverse hyperbolic arc tangent (double).");
    m.def("abs", [](double x) { return std::abs(x); }, "Absolute value (double).");
    m.def("erf", [](double x) { return std::erf(x); }, "Error function (double).");

    // We expose the class gdual<double> or gdual_d and its arithmetic
    pyaudi::expose_gdual<double>(m, "double");

    m.def("exp", [](const gdual_d &d) { return exp(d); }, "Exponential (gdual_double).");
    m.def("log", [](const gdual_d &d) { return log(d); }, "Natural logarithm (gdual_d).");
    m.def("sqrt", [](const gdual_d &d) { return sqrt(d); }, "Square root (gdual_d).");
    m.def("cbrt", [](const gdual_d &d) { return cbrt(d); }, "Cubic root (gdual_d).");
    m.def("sin", [](const gdual_d &d) { return sin(d); }, "Sine (gdual_d).");
    m.def("asin", [](const gdual_d &d) { return asin(d); }, "Arc sine (gdual_d).");
    m.def("cos", [](const gdual_d &d) { return cos(d); }, "Cosine (gdual_d).");
    m.def("acos", [](const gdual_d &d) { return acos(d); }, "Arc cosine (gdual_d).");
    m.def("sin_and_cos", [](const gdual_d &d) { sin_and_cos(d); }, "Sine and Cosine at once (gdual_d).");
    m.def("tan", [](const gdual_d &d) { return tan(d); }, "Tangent (gdual_d).");
    m.def("atan", [](const gdual_d &d) { return atan(d); }, "Arc tangent (gdual_d).");
    m.def("sinh", [](const gdual_d &d) { return sinh(d); }, "Hyperbolic sine (gdual_d).");
    m.def("asinh", [](const gdual_d &d) { return asinh(d); }, "Inverse hyperbolic sine (gdual_d).");
    m.def("cosh", [](const gdual_d &d) { return cosh(d); }, "Hyperbolic cosine (gdual_d).");
    m.def("acosh", [](const gdual_d &d) { return acosh(d); }, "Inverse hyperbolic cosine (gdual_d).");
    m.def("sinh_and_cosh", [](const gdual_d &d) { return sinh_and_cosh(d); },
          "Hyperbolic sine and hyperbolic cosine at once (gdual_d).");
    m.def("tanh", [](const gdual_d &d) { return tanh(d); }, "Hyperbolic tangent (gdual_d).");
    m.def("atanh", [](const gdual_d &d) { return atanh(d); }, "Inverse hyperbolic arc tangent (gdual_d).");
    m.def("abs", [](const gdual_d &d) { return abs(d); }, "Absolute value (gdual_d).");
    m.def("erf", [](const gdual_d &d) { return erf(d); }, "Error function (gdual_d).");

    // We expose the class gdual<vectorized<double>> or gdual_v and its arithmetic. Two constructors are added to allow
    // construction from vector::double
    pyaudi::expose_gdual<audi::vectorized<double>>(m, "vdouble")
        .def(py::init<std::vector<double>>())
        .def(py::init<std::vector<double>, const std::string &, unsigned int>());

    m.def("exp", +[](const gdual_v &d) { return exp(d); }, "Exponential (gdual_vdouble).");
    m.def("log", +[](const gdual_v &d) { return log(d); }, "Natural logarithm (gdual_v).");
    m.def("sqrt", +[](const gdual_v &d) { return sqrt(d); }, "Square root (gdual_v).");
    m.def("cbrt", +[](const gdual_v &d) { return cbrt(d); }, "Cubic root (gdual_v).");
    m.def("sin", +[](const gdual_v &d) { return sin(d); }, "Sine (gdual_v).");
    m.def("asin", +[](const gdual_v &d) { return asin(d); }, "Arc sine (gdual_v).");
    m.def("cos", +[](const gdual_v &d) { return cos(d); }, "Cosine (gdual_v).");
    m.def("acos", +[](const gdual_v &d) { return acos(d); }, "Arc cosine (gdual_v).");
    m.def("sin_and_cos", +[](const gdual_v &d) { return sin_and_cos(d); }, "Sine and Cosine at once (gdual_v).");
    m.def("tan", +[](const gdual_v &d) { return tan(d); }, "Tangent (gdual_v).");
    m.def("atan", +[](const gdual_v &d) { return atan(d); }, "Arc tangent (gdual_v).");
    m.def("sinh", +[](const gdual_v &d) { return sinh(d); }, "Hyperbolic sine (gdual_v).");
    m.def("asinh", +[](const gdual_v &d) { return asinh(d); }, "Inverse hyperbolic sine (gdual_v).");
    m.def("sinh_and_cosh", +[](const gdual_v &d) { return sinh_and_cosh(d); },
          "Hyperbolic sine and hyperbolic cosine at once (gdual_v).");
    m.def("cosh", +[](const gdual_v &d) { return cosh(d); }, "Hyperbolic cosine (gdual_v).");
    m.def("acosh", +[](const gdual_v &d) { return acosh(d); }, "Inverse hyperbolic cosine (gdual_v).");
    m.def("tanh", +[](const gdual_v &d) { return tanh(d); }, "Hyperbolic tangent (gdual_v).");
    m.def("atanh", +[](const gdual_v &d) { return atanh(d); }, "Inverse hyperbolic arc tangent (gdual_v).");
    m.def("abs", +[](const gdual_v &d) { return abs(d); }, "Absolute value (gdual_v).");
    m.def("erf", +[](const gdual_v &d) { return erf(d); }, "Error function (gdual_v).");
}