#ifndef AUDI_EXPOSE_TAYLOR_MODEL_H
#define AUDI_EXPOSE_TAYLOR_MODEL_H

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
#include <boost/numeric/interval.hpp>
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

py::class_<int_d> expose_int_d(const py::module &m)
{
    auto th = py::class_<int_d>(m, "int_d")
                  .def(py::init<double, double>(), py::arg("lower"), py::arg("upper"))
                  .def_property_readonly("lower", &int_d::lower)
                  .def_property_readonly("upper", &int_d::upper);
    return th;
}

py::class_<taylor_model> expose_taylor_model(const py::module &m)
{

    auto th
        = py::class_<taylor_model>(m, "taylor_model")
              // default constructor
              .def(py::init<>())

              // copy constructor
              .def(py::init<const taylor_model &>())

              // constant constructor
              .def(py::init<double>(), py::arg("constant"))

              // full constructor
              .def(py::init<const gdual<double> &, const int_d &, const var_map_d &, const var_map_i &>(),
                   py::arg("tpol"), py::arg("rem_bound"), py::arg("exp_point"), py::arg("domain"))

              // static identity() functions
              .def_static("identity", py::overload_cast<>(&taylor_model::identity))
              .def_static("identity", py::overload_cast<int_d, var_map_d, var_map_i>(&taylor_model::identity),
                          py::arg("rem"), py::arg("exp"), py::arg("domain"))

              // Getters (and setters)
              .def_property("tpol", &taylor_model::get_tpol, &taylor_model::set_tpol, "Taylor polynomial")
              .def_property("rem_bound", &taylor_model::get_rem, &taylor_model::set_rem, "Remainder bound")
              .def_property("exp_point", &taylor_model::get_exp, &taylor_model::set_exp, "Expansion point(s)")
              .def_property("domain", &taylor_model::get_dom, &taylor_model::set_dom, "Domain")
              .def_property_readonly("order", &taylor_model::get_order)
              .def_property_readonly("ndim", &taylor_model::get_ndim)

              // bounds methods
              .def(
                  "get_bounds",
                  [](const taylor_model &tm) {
                      return tm.get_bounds(); // calls the instance method
                  },
                  "Compute the bounds of this Taylor model")

              .def_static(
                  "get_custom_bounds",
                  [](const gdual<double> &tpol, const var_map_d &exp_points, const var_map_i &domain) {
                      return taylor_model::get_bounds(tpol, exp_points, domain);
                  },
                  py::arg("tpol"), py::arg("exp_points"), py::arg("domain"),
                  "Compute the bounds of a gdual polynomial over a given domain")

              .def(py::self + py::self)
              .def(py::self + double())
              .def(double() + py::self)
              .def(int_d() + py::self)
              .def(py::self + int_d())
              .def(py::self - py::self)
              .def(py::self - double())
              .def(double() - py::self)
              .def(int_d() - py::self)
              .def(py::self - int_d())
              .def(-py::self)
              .def(py::self * py::self)
              .def(double() * py::self)
              .def(py::self * double())
              .def(int_d() * py::self)
              .def(py::self * int_d())
              .def(py::self / py::self)
              .def(py::self / double())
              .def(double() / py::self)
              .def(int_d() / py::self)
              .def(py::self / int_d())
              .def(py::self == py::self)
              .def(py::self != py::self)
              .def(
                  "__pow__", [](const taylor_model &tm, int x) { return pow(tm, x); },
                  "Exponentiation (taylor_model, int).")
              .def(
                  "__pow__", [](const taylor_model &tm, double x) { return pow(tm, x); },
                  "Exponentiation (taylor_model, double).")

              .def_static(
                  "interval_equal",
                  [](const int_d &a, const int_d &b, std::optional<double> tol = std::nullopt) {
                      return taylor_model::interval_equal(a, b, tol);
                  },
                  py::arg("a"), py::arg("b"), py::arg("tol") = std::nullopt,
                  "Compare two intervals for equality with optional tolerance")

              .def_static(
                  "map_interval_equal",
                  [](const std::unordered_map<std::string, int_d> &a, const std::unordered_map<std::string, int_d> &b,
                     std::optional<double> tol = std::nullopt) { return taylor_model::map_interval_equal(a, b, tol); },
                  py::arg("a"), py::arg("b"), py::arg("tol") = std::nullopt,
                  "Compare two maps of intervals for equality with optional tolerance")

              .def_static(
                  "map_equal",
                  [](const std::unordered_map<std::string, double> &a, const std::unordered_map<std::string, double> &b,
                     std::optional<double> tol = std::nullopt) { return taylor_model::map_equal(a, b, tol); },
                  py::arg("a"), py::arg("b"), py::arg("tol") = std::nullopt,
                  "Compare two maps of scalar values for equality with optional tolerance")

              .def("__repr__",
                   [](const taylor_model &tm) {
                       std::ostringstream oss;
                       oss << tm;
                       return oss.str();
                   })
              .def("__str__", [](const taylor_model &tm) {
                  std::ostringstream oss;
                  oss << tm;
                  return oss.str();
              });

    // Functions exposed as members of taylor_models so that numpy arithmetics (for example np.exp(x)) would also
    // work with taylor_models. We provide both arc- and a- inverse nomenclature as to be consistent to both numpy
    // and our previous naming
    th.def(
          "exp", [](const taylor_model &d) { return exp(d); }, "Exponential.")
        .def(
            "log", [](const taylor_model &d) { return log(d); }, "Natural logarithm.")
        .def(
            "sqrt", [](const taylor_model &d) { return sqrt(d); }, "Square root.")
        .def(
            "sin", [](const taylor_model &d) { return sin(d); }, "Sine.")
        .def(
            "arcsin", [](const taylor_model &d) { return asin(d); }, "Arc sine.")
        .def(
            "asin", [](const taylor_model &d) { return asin(d); }, "Arc sine.")
        .def(
            "cos", [](const taylor_model &d) { return cos(d); }, "Cosine.")
        .def(
            "arccos", [](const taylor_model &d) { return acos(d); }, "Arc cosine.")
        .def(
            "acos", [](const taylor_model &d) { return acos(d); }, "Arc cosine.")
        .def(
            "tan", [](const taylor_model &d) { return tan(d); }, "Tangent.")
        .def(
            "arctan", [](const taylor_model &d) { return atan(d); }, "Arc tangent.")
        .def(
            "atan", [](const taylor_model &d) { return atan(d); }, "Arc tangent.")
        .def(
            "sinh", [](const taylor_model &d) { return sinh(d); }, "Hyperbolic sine.")
        .def(
            "arcsinh", [](const taylor_model &d) { return asinh(d); }, "Inverse hyperbolic sine.")
        .def(
            "asinh", [](const taylor_model &d) { return asinh(d); }, "Inverse hyperbolic sine.")
        .def(
            "cosh", [](const taylor_model &d) { return cosh(d); }, "Hyperbolic cosine.")
        .def(
            "arccosh", [](const taylor_model &d) { return acosh(d); }, "Inverse hyperbolic cosine.")
        .def(
            "acosh", [](const taylor_model &d) { return acosh(d); }, "Inverse hyperbolic cosine.")
        .def(
            "tanh", [](const taylor_model &d) { return tanh(d); }, "Hyperbolic tangent.")
        .def(
            "arctanh", [](const taylor_model &d) { return atanh(d); }, "Inverse hyperbolic arc tangent.")
        .def(
            "atanh", [](const taylor_model &d) { return atanh(d); }, "Inverse hyperbolic arc tangent.")
        .def("abs", [](const taylor_model &d) { return abs(d); }, "Absolute value.");
    return th;
}
} // namespace pyaudi
#endif
