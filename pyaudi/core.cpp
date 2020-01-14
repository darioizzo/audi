#include <audi/audi.hpp>
#include <cmath>
#include <pybind11/pybind11.h>
#include <sstream>

#if defined(AUDI_WITH_MPPP)
#include <audi/real128.hpp>
#include <mp++/extra/pybind11.hpp>
#include <mp++/mp++.hpp>
#endif

#include "common_utils.hpp"
//#include "expose_gdual.hpp"

using namespace audi;
namespace py = pybind11;

PYBIND11_MODULE(core, m)
{
    m.def(
        "exp", [](const gdual_d &d) { return exp(d); }, "Exponential (gdual_d).");
    m.def(
        "log", [](const gdual_d &d) { return log(d); }, "Natural logarithm (gdual_d).");
    m.def(
        "sqrt", [](const gdual_d &d) { return sqrt(d); }, "Square root (gdual_d).");
    m.def(
        "cbrt", [](const gdual_d &d) { return cbrt(d); }, "Cubic root (gdual_d).");
    m.def(
        "sin", [](const gdual_d &d) { return sin(d); }, "Sine (gdual_d).");
    m.def(
        "asin", [](const gdual_d &d) { return asin(d); }, "Arc sine (gdual_d).");
    m.def(
        "cos", [](const gdual_d &d) { return cos(d); }, "Cosine (gdual_d).");
    m.def(
        "acos", [](const gdual_d &d) { return acos(d); }, "Arc cosine (gdual_d).");
    m.def(
        "sin_and_cos", [](const gdual_d &d) { return pyaudi::v_to_l(sin_and_cos(d)); },
        "Sine and Cosine at once (gdual_d).");
    m.def(
        "tan", [](const gdual_d &d) { return tan(d); }, "Tangent (gdual_d).");
    m.def(
        "atan", [](const gdual_d &d) { return atan(d); }, "Arc tangent (gdual_d).");
    m.def(
        "sinh", [](const gdual_d &d) { return sinh(d); }, "Hyperbolic sine (gdual_d).");
    m.def(
        "asinh", [](const gdual_d &d) { return asinh(d); }, "Inverse hyperbolic sine (gdual_d).");
    m.def(
        "cosh", [](const gdual_d &d) { return cosh(d); }, "Hyperbolic cosine (gdual_d).");
    m.def(
        "acosh", [](const gdual_d &d) { return acosh(d); }, "Inverse hyperbolic cosine (gdual_d).");
    m.def(
        "sinh_and_cosh", [](const gdual_d &d) { return pyaudi::v_to_l(sinh_and_cosh(d)); },
        "Hyperbolic sine and hyperbolic cosine at once (gdual_d).");
    m.def(
        "tanh", [](const gdual_d &d) { return tanh(d); }, "Hyperbolic tangent (gdual_d).");
    m.def(
        "atanh", [](const gdual_d &d) { return atanh(d); }, "Inverse hyperbolic arc tangent (gdual_d).");
    m.def(
        "abs", [](const gdual_d &d) { return abs(d); }, "Absolute value (gdual_d).");
    m.def(
        "erf", [](const gdual_d &d) { return erf(d); }, "Error function (gdual_d).");
}