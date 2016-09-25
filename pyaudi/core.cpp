#include <boost/lexical_cast.hpp>
#include <cmath>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>

#include "../src/audi.hpp"

#if defined(__clang__) || defined(__GNUC__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wpedantic"
    #pragma GCC diagnostic ignored "-Wshadow"
    #pragma GCC diagnostic ignored "-Wsign-conversion"
    #pragma GCC diagnostic ignored "-Wdeprecated"
#endif

#include "pybind11/include/pybind11/operators.h"
#include "pybind11/include/pybind11/pybind11.h"
#include "pybind11/include/pybind11/stl.h"
#include "pybind11/include/pybind11/cast.h"

#if defined(__clang__) || defined(__GNUC__)
    #pragma GCC diagnostic pop
#endif

using namespace audi;
using gdual_d = gdual<double>;
namespace py = pybind11;

PYBIND11_PLUGIN(_core) {
    py::module m("_core", "pyaudi's core module");

    py::class_<gdual_d>(m,"gdual")
        .def(py::init<>())
        .def(py::init<const gdual_d &>())
        .def(py::init<double>())
        .def(py::init<double, const std::string &, unsigned int>())
        .def("__repr__",[](const gdual_d &g) -> std::string {
            std::ostringstream oss;
            oss << g;
            return oss.str();
        })
        .def("_repr_latex_",[](const gdual_d &g) -> std::string {
            std::ostringstream oss;
            g._poly().print_tex(oss);
            auto retval = oss.str();
            retval += std::string("+\\mathcal{O}\\left(")
                + boost::lexical_cast<std::string>(g.get_order() + 1) +  "\\right) \\]";
            return std::string("\\[ ") + retval;
        })
        .def("__getstate__", [](const gdual_d &p) {
            // Returns a tuple that contains the string
            // representation of a gdual_d as obtained
            // from the boost serialization library
            std::stringstream ss;
            boost::archive::text_oarchive oa(ss);
            oa << p;
            return py::make_tuple(ss.str());
        })
        .def("__setstate__", [](gdual_d &p, py::tuple t) {
            if (t.size() != 1)
                throw std::runtime_error("Invalid state!");
            // Invoke the default constructor.
            new (&p) gdual_d;
            // Reconstruct the gdual_d
            std::stringstream ss(t[0].cast<std::string>());
            boost::archive::text_iarchive ia(ss);
            ia >> p;
        })
        .def_property_readonly("symbol_set",&gdual_d::get_symbol_set)
        .def_property_readonly("symbol_set_size",&gdual_d::get_symbol_set_size)
        .def_property_readonly("degree",&gdual_d::degree)
        .def_property_readonly("order",&gdual_d::get_order)
        .def_property_readonly("constant_cf",&gdual_d::constant_cf)
        .def("extend_symbol_set", &gdual_d::extend_symbol_set, "Extends the symbol set")
        .def("integrate", &gdual_d::integrate<>, "Integrate with respect to argument")
        .def("partial", &gdual_d::partial<>, "Partial derivative with respect to argument")
        .def("evaluate",[](const gdual_d &g, const std::map< std::string, double> &dict) {return g.evaluate(std::unordered_map< std::string, double>(dict.begin(), dict.end()));} , "Evaluates the Taylor polynomial")
        .def("find_cf", [](const gdual_d &g, const std::vector<int> &v) {
            return g.find_cf(v);
        },"Find the coefficient of the Taylor expansion")
        .def("get_derivative", [](const gdual_d &g, const std::vector<int> &v) {
            return g.get_derivative(v);
        },"Finds the derivative (i.e. the coefficient of the Taylor expansion discounted of a factorial factor")
        .def("get_derivative", [](const gdual_d &g, const std::unordered_map<std::string, unsigned int> &dict) {
            return g.get_derivative(dict);
        },"Finds the derivative (i.e. the coefficient of the Taylor expansion discounted of a factorial factor")
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * py::self)
        .def(py::self / py::self)
        .def(py::self + double())
        .def(py::self - double())
        .def(py::self * double())
        .def(py::self / double())
        .def(-py::self)
        .def(+py::self)
        .def(double() + py::self)
        .def(double() - py::self)
        .def(double() * py::self)
        .def(double() / py::self)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__pow__",[](const gdual_d &g, double x) {return pow(g,x);} ,"Exponentiation (gdual_d, double).")
        .def("__pow__",[](const gdual_d &base, const gdual_d &g) {return pow(base,g);} ,"Exponentiation (gdual_d, gdual_d).")
        .def("__rpow__",[](const gdual_d &g, double x) {return pow(x,g);} ,"Exponentiation (double, gdual_d).")
    ;

    m.def("exp",[](const gdual_d &d) {return exp(d);},"Exponential (gdual_d).");
    m.def("exp",[](double x) {return std::exp(x);},"Exponential (double).");

    m.def("log",[](const gdual_d &d) {return log(d);},"Natural logarithm (gdual_d).");
    m.def("log",[](double x) {return std::log(x);},"Natural logarithm (double).");

    m.def("sqrt",[](const gdual_d &d) {return sqrt(d);},"Square root (gdual_d).");
    m.def("sqrt",[](double x) {return std::sqrt(x);},"Square root (double).");

    m.def("cbrt",[](const gdual_d &d) {return cbrt(d);},"Cubic root (gdual_d).");
    m.def("cbrt",[](double x) {return std::cbrt(x);},"Cubic root (double).");

    m.def("sin",[](const gdual_d &d) {return sin(d);},"Sine (gdual_d).");
    m.def("sin",[](double x) {return std::sin(x);},"Sine (double).");
    // m.def("sin",py::vectorize([](double x) {return std::sin(x);}),"Sine (vectorized double).");

    m.def("asin",[](const gdual_d &d) {return asin(d);},"Arc sine (gdual_d).");
    m.def("asin",[](double x) {return std::asin(x);},"Arc sine (double).");

    m.def("cos",[](const gdual_d &d) {return cos(d);},"Cosine (gdual_d).");
    m.def("cos",[](double x) {return std::cos(x);},"Cosine (double).");

    m.def("acos",[](const gdual_d &d) {return acos(d);},"Arc cosine (gdual_d).");
    m.def("acos",[](double x) {return std::acos(x);},"Arc cosine (double).");

    m.def("sin_and_cos",[](const gdual_d &d) {return sin_and_cos(d);}, "Sine and Cosine at once (gdual_d).");

    m.def("tan",[](const gdual_d &d) {return tan(d);},"Tangent (gdual_d).");
    m.def("tan",[](double x) {return std::tan(x);},"Tangent (double).");

    m.def("atan",[](const gdual_d &d) {return atan(d);},"Arc tangent (gdual_d).");
    m.def("atan",[](double x) {return std::atan(x);},"Arc tangent (double).");

    m.def("sinh",[](const gdual_d &d) {return sinh(d);},"Hyperbolic sine (gdual_d).");
    m.def("sinh",[](double x) {return std::sinh(x);},"Hyperbolic sine (double).");

    m.def("asinh",[](const gdual_d &d) {return asinh(d);},"Inverse hyperbolic sine (gdual_d).");
    m.def("asinh",[](double x) {return std::asinh(x);},"Inverse hyperbolic sine (double).");

    m.def("cosh",[](const gdual_d &d) {return cosh(d);},"Hyperbolic cosine (gdual_d).");
    m.def("cosh",[](double x) {return std::cosh(x);},"Hyperbolic cosine (double).");

    m.def("acosh",[](const gdual_d &d) {return acosh(d);},"Inverse hyperbolic cosine (gdual_d).");
    m.def("acosh",[](double x) {return std::acosh(x);},"Inverse hyperbolic cosine (double).");

    m.def("sinh_and_cosh",[](const gdual_d &d) {return sinh_and_cosh(d);} ,"Hyperbolic sine and hyperbolic cosine at once (gdual_d).");

    m.def("tanh",[](const gdual_d &d) {return tanh(d);},"Hyperbolic tangent (gdual_d).");
    m.def("tanh",[](double x) {return std::tanh(x);},"Hyperbolic tangent (double).");

    m.def("atanh",[](const gdual_d &d) {return atanh(d);},"Inverse hyperbolic arc tangent (gdual_d).");
    m.def("atanh",[](double x) {return std::atanh(x);},"Inverse hyperbolic arc tangent (double).");

    m.def("abs",[](const gdual_d &d) {return abs(d);},"Absolute value (gdual_d).");
    m.def("abs",[](double x) {return std::abs(x);},"Absolute value (double).");

    m.def("erf",[](const gdual_d &d) {return erf(d);},"Error function (gdual_d).");
    m.def("erf",[](double x) {return std::erf(x);},"Error function (double).");

    return m.ptr();
}
