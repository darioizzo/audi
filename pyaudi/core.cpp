#include <boost/python.hpp>
#include <cmath>
#include <piranha/thread_pool.hpp>
#include <vector>

#include <audi/audi.hpp>
#include <audi/functions.hpp>
#include <audi/functions_from_d.hpp>
#include <audi/invert_map.hpp>

#include "expose_gdual.hpp"

using namespace audi;
namespace bp = boost::python;

BOOST_PYTHON_MODULE(core)
{
    // We register a converter between vectorized_double and a python list (handy to avoid a long separate interface in
    // expose gdual)
    bp::to_python_converter<audi::vectorized_double, pyaudi::vectorized_double_to_python_list>();

    // We expose the gdual<double> using the expose_gdual defined in exposed_gdual.hpp
    pyaudi::expose_gdual<double>("double");

    // Similarly, we expose the gdual<vectorized_double> and we add two custom constructors to allow constructing it
    // from lists
    pyaudi::expose_gdual<vectorized_double>("vdouble")
        .def("__init__", bp::make_constructor(
                             +[](const bp::object &value) { return ::new gdual_v(pyaudi::l_to_v<double>(value)); }))
        .def("__init__",
             bp::make_constructor(+[](const bp::object &value, const std::string &symbol, unsigned int order) {
                 return ::new gdual_v(pyaudi::l_to_v<double>(value), symbol, order);
             }));

    // We expose the functions
    bp::def("exp", +[](const gdual_d &d) { return exp(d); }, "Exponential (gdual_d).");
    bp::def("exp", +[](double x) { return std::exp(x); }, "Exponential (double).");
    bp::def("exp", +[](const gdual_v &d) { return exp(d); }, "Exponential (gdual_v).");

    bp::def("log", +[](const gdual_d &d) { return log(d); }, "Natural logarithm (gdual_d).");
    bp::def("log", +[](double x) { return std::log(x); }, "Natural logarithm (double).");
    bp::def("log", +[](const gdual_v &d) { return log(d); }, "Natural logarithm (gdual_v).");

    bp::def("sqrt", +[](const gdual_d &d) { return sqrt(d); }, "Square root (gdual_d).");
    bp::def("sqrt", +[](double x) { return std::sqrt(x); }, "Square root (double).");
    bp::def("sqrt", +[](const gdual_v &d) { return sqrt(d); }, "Square root (gdual_v).");

    bp::def("cbrt", +[](const gdual_d &d) { return cbrt(d); }, "Cubic root (gdual_d).");
    bp::def("cbrt", +[](double x) { return std::cbrt(x); }, "Cubic root (double).");
    bp::def("cbrt", +[](const gdual_v &d) { return cbrt(d); }, "Cubic root (gdual_v).");

    bp::def("sin", +[](const gdual_d &d) { return sin(d); }, "Sine (gdual_d).");
    bp::def("sin", +[](double x) { return std::sin(x); }, "Sine (double).");
    bp::def("sin", +[](const gdual_v &d) { return sin(d); }, "Sine (gdual_v).");
    // bp::def("sin",py::vectorize(+[](double x) {return std::sin(x);}),"Sine (vectorized double).");

    bp::def("asin", +[](const gdual_d &d) { return asin(d); }, "Arc sine (gdual_d).");
    bp::def("asin", +[](double x) { return std::asin(x); }, "Arc sine (double).");
    bp::def("asin", +[](const gdual_v &d) { return asin(d); }, "Arc sine (gdual_v).");

    bp::def("cos", +[](const gdual_d &d) { return cos(d); }, "Cosine (gdual_d).");
    bp::def("cos", +[](double x) { return std::cos(x); }, "Cosine (double).");
    bp::def("cos", +[](const gdual_v &d) { return cos(d); }, "Cosine (gdual_v).");

    bp::def("acos", +[](const gdual_d &d) { return acos(d); }, "Arc cosine (gdual_d).");
    bp::def("acos", +[](double x) { return std::acos(x); }, "Arc cosine (double).");
    bp::def("acos", +[](const gdual_v &d) { return acos(d); }, "Arc cosine (gdual_v).");

    bp::def("sin_and_cos", +[](const gdual_d &d) { return pyaudi::v_to_l(sin_and_cos(d)); },
            "Sine and Cosine at once (gdual_d).");
    bp::def("sin_and_cos", +[](const gdual_v &d) { return pyaudi::v_to_l(sin_and_cos(d)); },
            "Sine and Cosine at once (gdual_v).");

    bp::def("tan", +[](const gdual_d &d) { return tan(d); }, "Tangent (gdual_d).");
    bp::def("tan", +[](double x) { return std::tan(x); }, "Tangent (double).");
    bp::def("tan", +[](const gdual_v &d) { return tan(d); }, "Tangent (gdual_v).");

    bp::def("atan", +[](const gdual_d &d) { return atan(d); }, "Arc tangent (gdual_d).");
    bp::def("atan", +[](double x) { return std::atan(x); }, "Arc tangent (double).");
    bp::def("atan", +[](const gdual_v &d) { return atan(d); }, "Arc tangent (gdual_v).");

    bp::def("sinh", +[](const gdual_d &d) { return sinh(d); }, "Hyperbolic sine (gdual_d).");
    bp::def("sinh", +[](double x) { return std::sinh(x); }, "Hyperbolic sine (double).");
    bp::def("sinh", +[](const gdual_v &d) { return sinh(d); }, "Hyperbolic sine (gdual_v).");

    bp::def("asinh", +[](const gdual_d &d) { return asinh(d); }, "Inverse hyperbolic sine (gdual_d).");
    bp::def("asinh", +[](double x) { return std::asinh(x); }, "Inverse hyperbolic sine (double).");
    bp::def("asinh", +[](const gdual_v &d) { return asinh(d); }, "Inverse hyperbolic sine (gdual_v).");

    bp::def("cosh", +[](const gdual_d &d) { return cosh(d); }, "Hyperbolic cosine (gdual_d).");
    bp::def("cosh", +[](double x) { return std::cosh(x); }, "Hyperbolic cosine (double).");
    bp::def("cosh", +[](const gdual_v &d) { return cosh(d); }, "Hyperbolic cosine (gdual_v).");

    bp::def("acosh", +[](const gdual_d &d) { return acosh(d); }, "Inverse hyperbolic cosine (gdual_d).");
    bp::def("acosh", +[](double x) { return std::acosh(x); }, "Inverse hyperbolic cosine (double).");
    bp::def("acosh", +[](const gdual_v &d) { return acosh(d); }, "Inverse hyperbolic cosine (gdual_v).");

    bp::def("sinh_and_cosh", +[](const gdual_d &d) { return pyaudi::v_to_l(sinh_and_cosh(d)); },
            "Hyperbolic sine and hyperbolic cosine at once (gdual_d).");
    bp::def("sinh_and_cosh", +[](const gdual_v &d) { return pyaudi::v_to_l(sinh_and_cosh(d)); },
            "Hyperbolic sine and hyperbolic cosine at once (gdual_v).");

    bp::def("tanh", +[](const gdual_d &d) { return tanh(d); }, "Hyperbolic tangent (gdual_d).");
    bp::def("tanh", +[](double x) { return std::tanh(x); }, "Hyperbolic tangent (double).");
    bp::def("tanh", +[](const gdual_v &d) { return tanh(d); }, "Hyperbolic tangent (gdual_v).");

    bp::def("atanh", +[](const gdual_d &d) { return atanh(d); }, "Inverse hyperbolic arc tangent (gdual_d).");
    bp::def("atanh", +[](double x) { return std::atanh(x); }, "Inverse hyperbolic arc tangent (double).");
    bp::def("atanh", +[](const gdual_v &d) { return atanh(d); }, "Inverse hyperbolic arc tangent (gdual_v).");

    bp::def("abs", +[](const gdual_d &d) { return abs(d); }, "Absolute value (gdual_d).");
    bp::def("abs", +[](double x) { return std::abs(x); }, "Absolute value (double).");
    bp::def("abs", +[](const gdual_v &d) { return abs(d); }, "Absolute value (gdual_v).");

    bp::def("erf", +[](const gdual_d &d) { return erf(d); }, "Error function (gdual_d).");
    bp::def("erf", +[](double x) { return std::erf(x); }, "Error function (double).");
    bp::def("erf", +[](const gdual_v &d) { return erf(d); }, "Error function (gdual_v).");

    bp::def("invert_map",
            +[](const bp::object &map_in, bool verbose) {
                return pyaudi::v_to_l(invert_map(pyaudi::l_to_v<gdual_d>(map_in), verbose));
            },
            "invert a Taylor map (gdual_d)");
            
    // Define a cleanup functor to be run when the module is unloaded.
    struct audi_cleanup_functor {
        void operator()() const
        {
            std::cout << "Shutting down the thread pool.\n";
            piranha::thread_pool_shutdown<void>();
        }
    };
    // Expose it.
    bp::class_<audi_cleanup_functor> cl_c("_audi_cleanup_functor", bp::init<>());
    cl_c.def("__call__", &audi_cleanup_functor::operator());
    // Register it.
    bp::object atexit_mod = bp::import("atexit");
    atexit_mod.attr("register")(audi_cleanup_functor{});
}
