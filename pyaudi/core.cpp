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
namespace py = pybind11;

template<typename T>
void expose_gdual(const py::module &m, std::string type)
{
    py::class_<gdual<T>>(m,("gdual_"+type).c_str())
    .def(py::init<>())
    .def(py::init<const gdual<T> &>())
    .def(py::init<double>())
    .def(py::init<double, const std::string &, unsigned int>())
    .def("__repr__",[](const gdual<T> &g) -> std::string {
        std::ostringstream oss;
        oss << g;
        return oss.str();
    })
    .def("_repr_latex_",[](const gdual<T> &g) -> std::string {
        std::ostringstream oss;
        g._poly().print_tex(oss);
        auto retval = oss.str();
        retval += std::string("+\\mathcal{O}\\left(")
            + boost::lexical_cast<std::string>(g.get_order() + 1) +  "\\right) \\]";
        return std::string("\\[ ") + retval;
    })
    .def("__getstate__", [](const gdual<T> &p) {
        // Returns a tuple that contains the string
        // representation of a gdual<T> as obtained
        // from the boost serialization library
        std::stringstream ss;
        boost::archive::text_oarchive oa(ss);
        oa << p;
        return py::make_tuple(ss.str());
    })
    .def("__setstate__", [](gdual<T> &p, py::tuple t) {
        if (t.size() != 1)
            throw std::runtime_error("Invalid state!");
        // Invoke the default constructor.
        new (&p) gdual<T>;
        // Reconstruct the gdual<T>
        std::stringstream ss(t[0].cast<std::string>());
        boost::archive::text_iarchive ia(ss);
        ia >> p;
    })
    .def_property_readonly("symbol_set",&gdual<T>::get_symbol_set)
    .def_property_readonly("symbol_set_size",&gdual<T>::get_symbol_set_size)
    .def_property_readonly("degree",&gdual<T>::degree)
    .def_property_readonly("order",&gdual<T>::get_order)
    .def_property_readonly("constant_cf",&gdual<T>::constant_cf)
    .def("extend_symbol_set", &gdual<T>::extend_symbol_set, "Extends the symbol set")
    .def("integrate", &gdual<T>::template integrate<>, "Integrate with respect to argument")
    .def("partial", &gdual<T>::template partial<>, "Partial derivative with respect to argument")
    .def("evaluate",[](const gdual<T> &g, const std::map< std::string, double> &dict) {return g.evaluate(std::unordered_map< std::string, double>(dict.begin(), dict.end()));} , "Evaluates the Taylor polynomial")
    .def("find_cf", [](const gdual<T> &g, const std::vector<int> &v) {
        return g.find_cf(v);
    },"Find the coefficient of the Taylor expansion")
    .def("get_derivative", [](const gdual<T> &g, const std::vector<int> &v) {
        return g.get_derivative(v);
    },"Finds the derivative (i.e. the coefficient of the Taylor expansion discounted of a factorial factor")
    .def("get_derivative", [](const gdual<T> &g, const std::unordered_map<std::string, unsigned int> &dict) {
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
    .def("__pow__",[](const gdual<double> &g, double x) {return pow(g,x);} ,("Exponentiation (gdual_"+type+", double).").c_str())
    .def("__pow__",[](const gdual<double> &base, const gdual<double> &g) {return pow(base,g);} ,("Exponentiation (gdual_"+type+", gdual_"+type+").").c_str())
    .def("__rpow__",[](const gdual<double> &g, double x) {return pow(x,g);} ,("Exponentiation (double, gdual_"+type+").").c_str())
    ;
}

PYBIND11_PLUGIN(_core) {
    py::module m("_core", "pyaudi's core module");

    expose_gdual<double>(m, "double");
    expose_gdual<vectorized_double>(m, "vdouble");

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

    m.def("exp",[](const gdual_v &d) {return exp(d);},"Exponential (gdual_v).");

    m.def("log",[](const gdual_v &d) {return log(d);},"Natural logarithm (gdual_v).");

    m.def("sqrt",[](const gdual_v &d) {return sqrt(d);},"Square root (gdual_v).");

    m.def("cbrt",[](const gdual_d &d) {return cbrt(d);},"Cubic root (gdual_v).");

    m.def("sin",[](const gdual_v &d) {return sin(d);},"Sine (gdual_v).");
    // m.def("sin",py::vectorize([](double x) {return std::sin(x);}),"Sine (vectorized double).");

    m.def("asin",[](const gdual_v &d) {return asin(d);},"Arc sine (gdual_v).");

    m.def("cos",[](const gdual_v &d) {return cos(d);},"Cosine (gdual_v).");

    m.def("acos",[](const gdual_v &d) {return acos(d);},"Arc cosine (gdual_v).");

    m.def("sin_and_cos",[](const gdual_v &d) {return sin_and_cos(d);}, "Sine and Cosine at once (gdual_v).");

    m.def("tan",[](const gdual_v &d) {return tan(d);},"Tangent (gdual_v).");

    m.def("atan",[](const gdual_v &d) {return atan(d);},"Arc tangent (gdual_v).");

    m.def("sinh",[](const gdual_v &d) {return sinh(d);},"Hyperbolic sine (gdual_v).");

    m.def("asinh",[](const gdual_v &d) {return asinh(d);},"Inverse hyperbolic sine (gdual_v).");

    m.def("cosh",[](const gdual_v &d) {return cosh(d);},"Hyperbolic cosine (gdual_v).");

    m.def("acosh",[](const gdual_v &d) {return acosh(d);},"Inverse hyperbolic cosine (gdual_v).");

    m.def("sinh_and_cosh",[](const gdual_v &d) {return sinh_and_cosh(d);} ,"Hyperbolic sine and hyperbolic cosine at once (gdual_v).");

    m.def("tanh",[](const gdual_v &d) {return tanh(d);},"Hyperbolic tangent (gdual_v).");

    m.def("atanh",[](const gdual_v &d) {return atanh(d);},"Inverse hyperbolic arc tangent (gdual_v).");

    m.def("abs",[](const gdual_v &d) {return abs(d);},"Absolute value (gdual_v).");

    m.def("erf",[](const gdual_v &d) {return erf(d);},"Error function (gdual_v).");

    return m.ptr();
}
