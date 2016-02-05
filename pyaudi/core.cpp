#include <cmath>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>

#include "../src/audi.hpp"
#include "pybind11/include/pybind11/operators.h"
#include "pybind11/include/pybind11/pybind11.h"
#include "pybind11/include/pybind11/stl.h"

using namespace audi;
namespace py = pybind11;

PYBIND11_PLUGIN(_core) {
    py::module m("_core", "pyaudi's core module");

    py::class_<gdual>(m,"gdual")
        .def(py::init<>())
        .def(py::init<const gdual &>())
        .def(py::init<gdual &&>())
        .def(py::init<const std::string &, unsigned int>())
        .def(py::init<double, unsigned int>())
        .def(py::init<double>())
        .def(py::init<double, const std::string &, unsigned int>())
        .def("__repr__",[](const gdual &g) -> std::string {
            std::ostringstream oss;
            oss << g;
            return oss.str();
        })
        .def_property_readonly("symbol_set",&gdual::get_symbol_set)
        .def_property_readonly("symbol_set_size",&gdual::get_symbol_set_size)
        .def("extend_symbol_set", &gdual::extend_symbol_set, "Extends the symbol set")
        .def("integrate", &gdual::integrate, "Integrate with respect to argument")
        .def("partial", &gdual::partial, "Partial derivative with respect to argument")
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * py::self)
        .def(py::self / py::self)
        .def(py::self + double())
        .def(py::self - double())
        .def(py::self * double())
        .def(py::self / double())
        .def(double() + py::self)
        .def(double() - py::self)
        .def(double() * py::self)
        .def(double() / py::self)
        .def(py::self + int())
        .def(py::self - int())
        .def(py::self * int())
        .def(py::self / int())
        .def(int() + py::self)
        .def(int() - py::self)
        .def(int() * py::self)
        .def(int() / py::self)
        .def("__pow__",[](const gdual &g, double x) {return pow(g,x);} ,"Exponentiation (gdual, double).")
        .def("__pow__",[](const gdual &g, int x) {return pow(g,x);} ,"Exponentiation (gdual, int).")
        .def("__pow__",[](const gdual & base, const gdual &g) {return pow(base,g);} ,"Exponentiation (gdual, gdual).")
        .def("__rpow__",[](const gdual &g, double x) {return pow(x,g);} ,"Exponentiation (double, gdual).")
        .def("__rpow__",[](const gdual &g, int x) {return pow(x,g);} ,"Exponentiation (int, gdual).")
    ;

    m.def("exp",[](const gdual &d) {return exp(d);},"Exponential (gdual).");
    m.def("exp",[](double x) {return std::exp(x);},"Exponential (double).");

    m.def("log",[](const gdual &d) {return log(d);},"Natural logarithm (gdual).");
    m.def("log",[](double x) {return std::log(x);},"Natural logarithm (double).");

    m.def("sqrt",[](const gdual &d) {return sqrt(d);},"Square root (gdual).");
    m.def("sqrt",[](double x) {return std::sqrt(x);},"Square root (double).");

    m.def("cbrt",[](const gdual &d) {return cbrt(d);},"Cubic root (gdual).");
    m.def("cbrt",[](double x) {return std::cbrt(x);},"Cubic root (double).");

    m.def("sin",[](const gdual &d) {return sin(d);},"Sine (gdual).");
    m.def("sin",[](double x) {return std::sin(x);},"Sine (double).");
    // m.def("sin",py::vectorize([](double x) {return std::sin(x);}),"Sine (vectorized double).");

    m.def("asin",[](const gdual &d) {return asin(d);},"Arc sine (gdual).");
    m.def("asin",[](double x) {return std::asin(x);},"Arc sine (double).");

    m.def("cos",[](const gdual &d) {return cos(d);},"Cosine (gdual).");
    m.def("cos",[](double x) {return std::cos(x);},"Cosine (double).");

    m.def("acos",[](const gdual &d) {return acos(d);},"Arc cosine (gdual).");
    m.def("acos",[](double x) {return std::acos(x);},"Arc cosine (double).");

    m.def("sin_and_cos",[](const gdual &d) -> std::vector<gdual> {
        auto arr = sin_and_cos(d);
        std::vector<gdual> retval;
        std::move(arr.begin(),arr.end(),std::back_inserter(retval));
        return retval;
    },"Sine and Cosine at once (gdual).");

    m.def("tan",[](const gdual &d) {return tan(d);},"Tangent (gdual).");
    m.def("tan",[](double x) {return std::tan(x);},"Tangent (double).");

    m.def("atan",[](const gdual &d) {return atan(d);},"Arc tangent (gdual).");
    m.def("atan",[](double x) {return std::atan(x);},"Arc tangent (double).");

    m.def("sinh",[](const gdual &d) {return sinh(d);},"Hyperbolic sine (gdual).");
    m.def("sinh",[](double x) {return std::sinh(x);},"Hyperbolic sine (double).");

    m.def("asinh",[](const gdual &d) {return asinh(d);},"Inverse hyperbolic sine (gdual).");
    m.def("asinh",[](double x) {return std::asinh(x);},"Inverse hyperbolic sine (double).");

    m.def("cosh",[](const gdual &d) {return cosh(d);},"Hyperbolic cosine (gdual).");
    m.def("cosh",[](double x) {return std::cosh(x);},"Hyperbolic cosine (double).");

    m.def("acosh",[](const gdual &d) {return acosh(d);},"Inverse hyperbolic cosine (gdual).");
    m.def("acosh",[](double x) {return std::acosh(x);},"Inverse hyperbolic cosine (double).");

    m.def("sinh_and_cosh",[](const gdual &d) -> std::vector<gdual> {
        auto arr = sinh_and_cosh(d);
        std::vector<gdual> retval;
        std::move(arr.begin(),arr.end(),std::back_inserter(retval));
        return retval;
    },"Hyperbolic sine and hyperbolic cosine at once (gdual).");

    m.def("tanh",[](const gdual &d) {return tanh(d);},"Hyperbolic tangent (gdual).");
    m.def("tanh",[](double x) {return std::tanh(x);},"Hyperbolic tangent (double).");

    m.def("atanh",[](const gdual &d) {return atanh(d);},"Inverse hyperbolic arc tangent (gdual).");
    m.def("atanh",[](double x) {return std::atanh(x);},"Inverse hyperbolic arc tangent (double).");

    m.def("abs",[](const gdual &d) {return abs(d);},"Absolute value (gdual).");
    m.def("abs",[](double x) {return std::abs(x);},"Absolute value (double).");

    m.def("erf",[](const gdual &d) {return erf(d);},"Error function (gdual).");
    m.def("erf",[](double x) {return std::erf(x);},"Error function (double).");

    return m.ptr();
}