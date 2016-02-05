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
        .def("__rpow__",[](double x, const gdual &g) {return pow(x,g);} ,"Exponentiation (double, gdual).")
        .def("__rpow__",[](int x, const gdual &g) {return pow(x,g);} ,"Exponentiation (int, gdual).")
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

    m.def("cos",[](const gdual &d) {return cos(d);},"Cosine (gdual).");
    m.def("cos",[](double x) {return std::cos(x);},"Cosine (double).");

    m.def("sin_and_cos",[](const gdual &d) -> std::vector<gdual> {
        auto arr = sin_and_cos(d);
        std::vector<gdual> retval;
        std::move(arr.begin(),arr.end(),std::back_inserter(retval));
        return retval;
    },"Sine and Cosine at once (gdual).");


    return m.ptr();
}