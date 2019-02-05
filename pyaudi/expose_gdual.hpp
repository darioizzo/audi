#ifndef AUDI_EXPOSE_GDUAL_H
#define AUDI_EXPOSE_GDUAL_H

#include <audi/functions.hpp>
#include <audi/gdual.hpp>

#include <boost/lexical_cast.hpp>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <sstream>
#include <string>

#include "docstrings.hpp"

using namespace audi;
namespace py = pybind11;

namespace pyaudi
{
// This is the interface common across types
template <typename T>
py::class_<gdual<T>> expose_gdual(const py::module &m, std::string type)
{
    auto th
        = py::class_<gdual<double>>(m, ("gdual_" + type).c_str(), gdual_docstring().c_str())
              .def(py::init<>())
              .def(py::init<const gdual<double> &>())
              .def(py::init<double>())
              .def(py::init<double, const std::string &, unsigned int>())
              .def("__repr__",
                   +[](const gdual<T> &g) -> std::string {
                       std::ostringstream oss;
                       oss << g;
                       return oss.str();
                   })
              .def("_repr_latex_",
                   +[](const gdual<T> &g) -> std::string {
                       std::ostringstream oss;
                       g._poly().print_tex(oss);
                       auto retval = oss.str();
                       retval += std::string("+\\mathcal{O}\\left(")
                                 + boost::lexical_cast<std::string>(g.get_order() + 1) + "\\right) \\]";
                       return std::string("\\[ ") + retval;
                   })
              .def_property_readonly("symbol_set", &gdual<T>::get_symbol_set)
              .def_property_readonly("symbol_set_size", &gdual<T>::get_symbol_set_size)
              .def_property_readonly("degree", &gdual<T>::degree, gdual_degree_docstring().c_str())
              .def_property_readonly("order", &gdual<T>::get_order, "truncation order (>= degree)")
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
              .def("__pow__", +[](const gdual<T> &gd, double x) { return pow(gd, x); },
                   ("Exponentiation (gdual_" + type + ", double).").c_str())
              .def("__pow__", +[](const gdual<T> &base, const gdual<T> &gd) { return pow(base, gd); },
                   ("Exponentiation (gdual_" + type + ", gdual_" + type + ").").c_str())
              .def("__rpow__", +[](const gdual<T> &gd, double x) { return pow(x, gd); },
                   ("Exponentiation (double, gdual_" + type + ").").c_str())
              .def_property_readonly("constant_cf", &gdual<T>::constant_cf, "Constant term of the polynomial")
              .def("evaluate", &gdual<T>::evaluate, "Evaluates the Taylor polynomial")
              .def("find_cf", &gdual<T>::find_cf, "Finds the coefficient of the Taylor expansion")
              .def("get_derivative", py::overload_cast<std::vector<unsigned>>(&gdual<T>::get_derivative),
                   "Finds the derivative (i.e. the coefficient of the Taylor expansion discounted by the factorial "
                   "factor");
              //.def("get_derivative",
               //    py::overload_cast<std::unordered_map<std::string, unsigned>>(&gdual<T>::get_derivative),
                //   "Finds the derivative (i.e. the coefficient of the Taylor expansion discounted by the factorial "
                  // "factor");
    return th;
}
} // namespace pyaudi
#endif
