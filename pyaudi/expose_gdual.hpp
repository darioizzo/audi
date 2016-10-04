#ifndef AUDI_EXPOSE_GDUAL_H
#define AUDI_EXPOSE_GDUAL_H

#include <boost/lexical_cast.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/python.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <string>
#include <sstream> //ostringstream, stringstream
#include <stdexcept> // stringstream
#include <vector>

#include "common_utils.hpp"
#include "python_includes.hpp"
#include "../src/audi.hpp"

namespace bp = boost::python;
using namespace audi;

namespace pyaudi {

template<typename T>
auto expose_gdual(std::string type)
{
     return bp::class_<gdual<T>>(("gdual_"+type).c_str())
    .def(bp::init<>())
    .def(bp::init<const gdual<T> &>())
    .def(bp::init<T>())
    .def(bp::init<T, const std::string &, unsigned int>())
    .def("__repr__", +[](const gdual<T> &g) -> std::string {
        std::ostringstream oss;
        oss << g;
        return oss.str();
    })
    .def("_repr_latex_",+[](const gdual<T> &g) -> std::string {
        std::ostringstream oss;
        g._poly().print_tex(oss);
        auto retval = oss.str();
        retval += std::string("+\\mathcal{O}\\left(")
            + boost::lexical_cast<std::string>(g.get_order() + 1) +  "\\right) \\]";
        return std::string("\\[ ") + retval;
    })
    .def("__getstate__", +[](const gdual<T> &p) {
        // Returns a tuple that contains the string
        // representation of a gdual<T> as obtained
        // from the boost serialization library
        std::stringstream ss;
        boost::archive::text_oarchive oa(ss);
        oa << p;
        auto s = ss.str();
        return bp::make_tuple(pyaudi::make_bytes(s.data(),boost::numeric_cast<Py_ssize_t>(s.size())));
    })
    .def("__setstate__", +[](gdual<T> &p, bp::tuple state) {
        if (len(state) != 1) {
            pyaudi_throw(PyExc_ValueError,"the state tuple must have a single element");
        }
        auto ptr = PyBytes_AsString(bp::object(state[0]).ptr());
        if (!ptr) {
            pyaudi_throw(PyExc_TypeError,"a bytes object is needed to deserialize a population");
        }
        const auto size = len(state[0]);
        std::string s(ptr,ptr + size);
        std::istringstream iss;
        iss.str(s);
        boost::archive::text_iarchive ia(iss);
        ia >> p;
    })
    .add_property("symbol_set", +[](const gdual<T> &gd){return pyaudi::v_to_l(gd.get_symbol_set());})
    .add_property("symbol_set_size",&gdual<T>::get_symbol_set_size)
    .add_property("degree",&gdual<T>::degree, "polynomial degree (<= order)")
    .add_property("order",&gdual<T>::get_order, "truncation order (>= degree)")
    .add_property("constant_cf",&gdual<T>::constant_cf, "Constant term of the polynomial")
    .def("extend_symbol_set", +[](gdual<T> &gd, const bp::object& in){return gd.extend_symbol_set(pyaudi::l_to_v<std::string>(in));}, "Extends the symbol set")
    .def("integrate", &gdual<T>::template integrate<>, "Integrate with respect to argument")
    .def("partial", &gdual<T>::template partial<>, "Partial derivative with respect to argument")
    .def("evaluate", +[](const gdual<T> &gd, const bp::dict &pydict) {return gd.evaluate(pydict_to_umap<std::string, double>(pydict));} , "Evaluates the Taylor polynomial")
    .def("find_cf", +[](const gdual<T> &gd, const bp::object& in) {
        return gd.find_cf(pyaudi::l_to_v<int>(in));
    }, "Find the coefficient of the Taylor expansion")
    .def("get_derivative", +[](const gdual<T> &gd, const bp::object& in) {
        return gd.get_derivative(pyaudi::l_to_v<int>(in));
    },"Finds the derivative (i.e. the coefficient of the Taylor expansion discounted of a factorial factor")
    .def("get_derivative", +[](const gdual<T> &g, const bp::dict &pydict) {
        return g.get_derivative(pydict_to_umap<std::string, unsigned int>(pydict));
    },"Finds the derivative (i.e. the coefficient of the Taylor expansion discounted of a factorial factor")
    .def(bp::self + bp::self)
    .def(bp::self - bp::self)
    .def(bp::self * bp::self)
    .def(bp::self / bp::self)
    .def(bp::self + double())
    .def(bp::self - double())
    .def(bp::self * double())
    .def(bp::self / double())
    .def(-bp::self)
    .def(+bp::self)
    .def(double() + bp::self)
    .def(double() - bp::self)
    .def(double() * bp::self)
    .def(double() / bp::self)
    .def(bp::self == bp::self)
    .def(bp::self != bp::self)
    .def("__pow__",+[](const gdual<T> &gd, double x) {return pow(gd,x);} ,("Exponentiation (gdual_"+type+", double).").c_str())
    .def("__pow__",+[](const gdual<T> &base, const gdual<T> &gd) {return pow(base,gd);} ,("Exponentiation (gdual_"+type+", gdual_"+type+").").c_str())
    .def("__rpow__",+[](const gdual<T> &gd, double x) {return pow(x,gd);} ,("Exponentiation (double, gdual_"+type+").").c_str())
    ;
}
}
#endif
