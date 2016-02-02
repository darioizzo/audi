#include <sstream>
#include <string>

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
	;

	return m.ptr();
}
