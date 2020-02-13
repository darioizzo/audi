#ifndef AUDI_EXPOSE_GDUAL_H
#define AUDI_EXPOSE_GDUAL_H

#if defined(AUDI_WITH_QUADMATH)
#include <audi/real128.hpp>
#endif

#include <obake/symbols.hpp>
#include <obake/tex_stream_insert.hpp>

#include <audi/audi.hpp>
#include <audi/back_compatibility.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <sstream>   //ostringstream, stringstream
#include <stdexcept> // stringstream
#include <string>
#include <type_traits>
#include <vector>

#include "common_utils.hpp"
#include "docstrings.hpp"

using namespace audi;
namespace py = pybind11;

namespace pyaudi
{

template <typename T>
py::class_<gdual<T>> expose_gdual(const py::module &m, std::string type)
{
    auto th
        = py::class_<gdual<T>>(m, ("gdual_" + type).c_str(), gdual_docstring().c_str())
              .def(py::init<>())
              .def(py::init<const gdual<T> &>())
              .def(py::init<T>())
              .def(py::init<T, const std::string &, unsigned int>())
              .def("__repr__",
                   [](const gdual<T> &g) -> std::string {
                       std::ostringstream oss;
                       oss << g;
                       return oss.str();
                   })
              .def("_repr_latex_",
                   [](const gdual<T> &g) -> std::string {
                       std::ostringstream oss;
                       obake::tex_stream_insert(oss, g._poly());
                       auto retval = oss.str();
                       retval += std::string("+\\mathcal{O}\\left(")
                                 + boost::lexical_cast<std::string>(g.get_order() + 1) + "\\right) \\]";
                       return std::string("\\[ ") + retval;
                   })
              .def(py::pickle(
                  [](const gdual<T> &p) { // __getstate__
                                          /* Return a tuple that fully encodes the state of the object */
                      std::stringstream ss;
                      boost::archive::text_oarchive oa(ss);
                      oa << p;
                      return py::make_tuple(ss.str());
                  },
                  [](py::tuple t) { // __setstate__
                      if (t.size() != 1) throw std::runtime_error("Invalid state!");
                      // Create a C++ instance
                      gdual<T> p;
                      // Reconstruct the gdual<T>
                      std::stringstream ss(t[0].cast<std::string>());
                      boost::archive::text_iarchive ia(ss);
                      ia >> p;
                      return p;
                  }))
              .def_property_readonly("symbol_set", &gdual<T>::get_symbol_set)
              .def_property_readonly("symbol_set_size", &gdual<T>::get_symbol_set_size)
              .def_property_readonly("degree", &gdual<T>::degree, gdual_degree_docstring().c_str())
              .def_property_readonly("order", &gdual<T>::get_order, "truncation order (>= degree)")
              .def("info", &gdual<T>::info, "Full information on the gdual")
              .def("extend_symbol_set", &gdual<T>::extend_symbol_set, "Extends the symbol set")
              .def("integrate", &gdual<T>::template integrate<>, "Integrate with respect to argument")
              .def("partial", &gdual<T>::template partial<>, "Partial derivative with respect to argument")
              .def("is_zero", &gdual<T>::is_zero, "checks if all coefficients of the gdual are zero within a tolerance")
              .def("trim", &gdual<T>::trim,
                   "returns a new gdual removing  all coefficients that are smaller than a tolerance")
              .def("extract_terms", &gdual<T>::extract_terms,
                   "returns a new gdual containing only terms of a given order")
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
              .def(py::self < py::self)
              .def(py::self > py::self)
              .def(
                  "__pow__", [](const gdual<T> &gd, double x) { return pow(gd, x); },
                  ("Exponentiation (gdual_" + type + ", double).").c_str())
              .def(
                  "__pow__", [](const gdual<T> &base, const gdual<T> &gd) { return pow(base, gd); },
                  ("Exponentiation (gdual_" + type + ", gdual_" + type + ").").c_str())
              .def(
                  "__rpow__", [](const gdual<T> &gd, double x) { return pow(x, gd); },
                  ("Exponentiation (double, gdual_" + type + ").").c_str())
              .def_property_readonly("constant_cf", &gdual<T>::constant_cf, "Constant term of the polynomial")
              .def(
                  "evaluate",
                  [](const gdual<T> &g, const py::dict &in) {
                      return g.evaluate(pyaudi::py_dict_to_obake_sm<double>(in));
                  },
                  "Evaluates the Taylor polynomial")
              .def(
                  "find_cf", [](const gdual<T> &g, const std::vector<int> &v) { return g.find_cf(v); },
                  "Find the coefficient of the Taylor expansion")
              .def(
                  "get_derivative",
                  [](const gdual<T> &g, const std::vector<unsigned> &v) { return g.get_derivative(v); },
                  "Finds the derivative (i.e. the coefficient of the Taylor expansion discounted of a factorial factor")
              .def(
                  "get_derivative",
                  [](const gdual<T> &g, const std::unordered_map<std::string, unsigned> &dict) {
                      return g.get_derivative(dict);
                  },
                  "Finds the derivative (i.e. the coefficient of the Taylor expansion discounted of a factorial "
                  "factor")
              .def(
                  "subs", [](gdual<T> &gd, const std::string &sym, const T &in) { return gd.subs(sym, in); },
                  "Substitutes a symbol with a value (does not remove symbol from symbol set)")
              .def(
                  "subs", [](gdual<T> &gd, const std::string &sym, const gdual<T> &in) { return gd.subs(sym, in); },
                  "Substitutes a symbol with a gdual");

    // Functions exposed as members of gduals so that numpy arithmetics (for example np.exp(x)) would also work
    // with gduals. We provide both arc- and a- inverse nomenclature as to be consistent to both numpy and our previous naming
    th.def(
          "exp", [](const gdual<T> &d) { return exp(d); }, "Exponential.")
        .def(
            "log", [](const gdual<T> &d) { return log(d); }, "Natural logarithm.")
        .def(
            "sqrt", [](const gdual<T> &d) { return sqrt(d); }, "Square root.")
        .def(
            "cbrt", [](const gdual<T> &d) { return cbrt(d); }, "Cubic root.")
        .def(
            "sin", [](const gdual<T> &d) { return sin(d); }, "Sine.")
        .def(
            "arcsin", [](const gdual<T> &d) { return asin(d); }, "Arc sine.")
        .def(
            "asin", [](const gdual<T> &d) { return asin(d); }, "Arc sine.")
        .def(
            "cos", [](const gdual<T> &d) { return cos(d); }, "Cosine.")
        .def(
            "arccos", [](const gdual<T> &d) { return acos(d); }, "Arc cosine.")
        .def(
            "acos", [](const gdual<T> &d) { return acos(d); }, "Arc cosine.")
        .def(
            "tan", [](const gdual<T> &d) { return tan(d); }, "Tangent.")
        .def(
            "arctan", [](const gdual<T> &d) { return atan(d); }, "Arc tangent.")
        .def(
            "atan", [](const gdual<T> &d) { return atan(d); }, "Arc tangent.")
        .def(
            "sinh", [](const gdual<T> &d) { return sinh(d); }, "Hyperbolic sine.")
        .def(
            "arcsinh", [](const gdual<T> &d) { return asinh(d); }, "Inverse hyperbolic sine.")
        .def(
            "asinh", [](const gdual<T> &d) { return asinh(d); }, "Inverse hyperbolic sine.")
        .def(
            "cosh", [](const gdual<T> &d) { return cosh(d); }, "Hyperbolic cosine.")
        .def(
            "arccosh", [](const gdual<T> &d) { return acosh(d); }, "Inverse hyperbolic cosine.")
        .def(
            "acosh", [](const gdual<T> &d) { return acosh(d); }, "Inverse hyperbolic cosine.")
        .def(
            "tanh", [](const gdual<T> &d) { return tanh(d); }, "Hyperbolic tangent.")
        .def(
            "arctanh", [](const gdual<T> &d) { return atanh(d); }, "Inverse hyperbolic arc tangent.")
        .def(
            "atanh", [](const gdual<T> &d) { return atanh(d); }, "Inverse hyperbolic arc tangent.")
        .def(
            "abs", [](const gdual<T> &d) { return abs(d); }, "Absolute value.")
        .def(
            "erf", [](const gdual<T> &d) { return erf(d); }, "Error function.")
        .def(
            "log10", [](const gdual<T> &d) { return log(d) / std::log(10.); }, "Base 10 Logarithm.")
        .def(
            "log2", [](const gdual<T> &d) { return log(d) / std::log(2); }, "Base 2 Logarithm.")
        .def(
            "log1p", [](const gdual<T> &d) { return log(d + 1.0); }, "log(x)+1, inverse of expm1")
        .def(
            "expm1", [](const gdual<T> &d) { return exp(d) - 1.0; }, "exp(x) - 1, inverse of log1p")
        .def(
            "erfc", [](const gdual<T> &d) { return 1.0 - erf(d); }, "Error Function Complement");
    return th;
}
} // namespace pyaudi
#endif
