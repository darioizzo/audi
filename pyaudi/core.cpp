#include <boost/python.hpp>
#include <cmath>
#include <piranha/thread_pool.hpp>
#include <vector>

#include <audi/config.hpp>

#if defined(AUDI_WITH_MPPP)
#include <audi/real128.hpp>
#endif

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
    bp::to_python_converter<audi::vectorized<double>, pyaudi::vectorized_to_python_list<double>>();

    // We expose the gdual<double> arithmetic
    pyaudi::expose_gdual<double>("double");
    bp::def("exp", +[](const gdual_d &d) { return exp(d); }, "Exponential (gdual_double).");
    bp::def("log", +[](const gdual_d &d) { return log(d); }, "Natural logarithm (gdual_d).");
    bp::def("sqrt", +[](const gdual_d &d) { return sqrt(d); }, "Square root (gdual_d).");
    bp::def("cbrt", +[](const gdual_d &d) { return cbrt(d); }, "Cubic root (gdual_d).");
    bp::def("sin", +[](const gdual_d &d) { return sin(d); }, "Sine (gdual_d).");
    bp::def("asin", +[](const gdual_d &d) { return asin(d); }, "Arc sine (gdual_d).");
    bp::def("cos", +[](const gdual_d &d) { return cos(d); }, "Cosine (gdual_d).");
    bp::def("acos", +[](const gdual_d &d) { return acos(d); }, "Arc cosine (gdual_d).");
    bp::def("sin_and_cos", +[](const gdual_d &d) { return pyaudi::v_to_l(sin_and_cos(d)); },
            "Sine and Cosine at once (gdual_d).");
    bp::def("tan", +[](const gdual_d &d) { return tan(d); }, "Tangent (gdual_d).");
    bp::def("atan", +[](const gdual_d &d) { return atan(d); }, "Arc tangent (gdual_d).");
    bp::def("sinh", +[](const gdual_d &d) { return sinh(d); }, "Hyperbolic sine (gdual_d).");
    bp::def("asinh", +[](const gdual_d &d) { return asinh(d); }, "Inverse hyperbolic sine (gdual_d).");
    bp::def("cosh", +[](const gdual_d &d) { return cosh(d); }, "Hyperbolic cosine (gdual_d).");
    bp::def("acosh", +[](const gdual_d &d) { return acosh(d); }, "Inverse hyperbolic cosine (gdual_d).");
    bp::def("sinh_and_cosh", +[](const gdual_d &d) { return pyaudi::v_to_l(sinh_and_cosh(d)); },
            "Hyperbolic sine and hyperbolic cosine at once (gdual_d).");
    bp::def("tanh", +[](const gdual_d &d) { return tanh(d); }, "Hyperbolic tangent (gdual_d).");
    bp::def("atanh", +[](const gdual_d &d) { return atanh(d); }, "Inverse hyperbolic arc tangent (gdual_d).");
    bp::def("abs", +[](const gdual_d &d) { return abs(d); }, "Absolute value (gdual_d).");
    bp::def("erf", +[](const gdual_d &d) { return erf(d); }, "Error function (gdual_d).");

    // We expose the gdual<vectorized<double>> arithmetic
    pyaudi::expose_gdual<vectorized<double>>("vdouble")
        .def("__init__", bp::make_constructor(
                             +[](const bp::object &value) { return ::new gdual_v(pyaudi::l_to_v<double>(value)); }))
        .def("__init__",
             bp::make_constructor(+[](const bp::object &value, const std::string &symbol, unsigned int order) {
                 return ::new gdual_v(pyaudi::l_to_v<double>(value), symbol, order);
             }));
    bp::def("exp", +[](const gdual_v &d) { return exp(d); }, "Exponential (gdual_vdouble).");
    bp::def("log", +[](const gdual_v &d) { return log(d); }, "Natural logarithm (gdual_v).");
    bp::def("sqrt", +[](const gdual_v &d) { return sqrt(d); }, "Square root (gdual_v).");
    bp::def("cbrt", +[](const gdual_v &d) { return cbrt(d); }, "Cubic root (gdual_v).");
    bp::def("sin", +[](const gdual_v &d) { return sin(d); }, "Sine (gdual_v).");
    bp::def("asin", +[](const gdual_v &d) { return asin(d); }, "Arc sine (gdual_v).");
    bp::def("cos", +[](const gdual_v &d) { return cos(d); }, "Cosine (gdual_v).");
    bp::def("acos", +[](const gdual_v &d) { return acos(d); }, "Arc cosine (gdual_v).");
    bp::def("sin_and_cos", +[](const gdual_v &d) { return pyaudi::v_to_l(sin_and_cos(d)); },
            "Sine and Cosine at once (gdual_v).");
    bp::def("tan", +[](const gdual_v &d) { return tan(d); }, "Tangent (gdual_v).");
    bp::def("atan", +[](const gdual_v &d) { return atan(d); }, "Arc tangent (gdual_v).");
    bp::def("sinh", +[](const gdual_v &d) { return sinh(d); }, "Hyperbolic sine (gdual_v).");
    bp::def("asinh", +[](const gdual_v &d) { return asinh(d); }, "Inverse hyperbolic sine (gdual_v).");
    bp::def("cosh", +[](const gdual_v &d) { return cosh(d); }, "Hyperbolic cosine (gdual_v).");
    bp::def("acosh", +[](const gdual_v &d) { return acosh(d); }, "Inverse hyperbolic cosine (gdual_v).");
    bp::def("tanh", +[](const gdual_v &d) { return tanh(d); }, "Hyperbolic tangent (gdual_v).");
    bp::def("atanh", +[](const gdual_v &d) { return atanh(d); }, "Inverse hyperbolic arc tangent (gdual_v).");
    bp::def("abs", +[](const gdual_v &d) { return abs(d); }, "Absolute value (gdual_v).");
    bp::def("erf", +[](const gdual_v &d) { return erf(d); }, "Error function (gdual_v).");

    // We expose simple mathematical functions to work on python floats
    bp::def("exp", +[](double x) { return std::exp(x); }, "Exponential (double).");
    bp::def("log", +[](double x) { return std::log(x); }, "Natural logarithm (double).");
    bp::def("sqrt", +[](double x) { return std::sqrt(x); }, "Square root (double).");
    bp::def("cbrt", +[](double x) { return std::cbrt(x); }, "Cubic root (double).");
    bp::def("sin", +[](double x) { return std::sin(x); }, "Sine (double).");
    bp::def("asin", +[](double x) { return std::asin(x); }, "Arc sine (double).");
    bp::def("cos", +[](double x) { return std::cos(x); }, "Cosine (double).");
    bp::def("acos", +[](double x) { return std::acos(x); }, "Arc cosine (double).");
    bp::def("tan", +[](double x) { return std::tan(x); }, "Tangent (double).");
    bp::def("atan", +[](double x) { return std::atan(x); }, "Arc tangent (double).");
    bp::def("sinh", +[](double x) { return std::sinh(x); }, "Hyperbolic sine (double).");
    bp::def("asinh", +[](double x) { return std::asinh(x); }, "Inverse hyperbolic sine (double).");
    bp::def("cosh", +[](double x) { return std::cosh(x); }, "Hyperbolic cosine (double).");
    bp::def("acosh", +[](double x) { return std::acosh(x); }, "Inverse hyperbolic cosine (double).");
    bp::def("sinh_and_cosh", +[](const gdual_v &d) { return pyaudi::v_to_l(sinh_and_cosh(d)); },
            "Hyperbolic sine and hyperbolic cosine at once (gdual_v).");
    bp::def("tanh", +[](double x) { return std::tanh(x); }, "Hyperbolic tangent (double).");
    bp::def("atanh", +[](double x) { return std::atanh(x); }, "Inverse hyperbolic arc tangent (double).");
    bp::def("abs", +[](double x) { return std::abs(x); }, "Absolute value (double).");
    bp::def("erf", +[](double x) { return std::erf(x); }, "Error function (double).");

#if defined(AUDI_WITH_MPPP)
    // We expose the mppp::real128 class. NOTE: using boost::shared_ptr as holder is necessary for proper alignement of __float128 by boost.
    // Otherwise: segfault.
    // See http://webcache.googleusercontent.com/search?q=cache:xhGhABO-W58J:fhtagn.net/prog/2015/04/16/quaternion_boost_python.html&num=1&client=firefox-b&hl=en&gl=de&strip=1&vwsrc=0
    bp::class_<mppp::real128, boost::shared_ptr<mppp::real128>>("real128", "Multiple precision float (a wrapper around __float128)")
              .def(bp::init<const std::string &>())
              .def("__repr__",
                   +[](const mppp::real128 &g) -> std::string {
                       std::ostringstream oss;
                       oss << g;
                       return oss.str();
                   })
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
              .def("__pow__", +[](const mppp::real128 &gd, double x) { return audi::pow(gd, x); },
                   "Exponentiation (mppp::real128, double).")
              .def("__pow__", +[](const mppp::real128 &base, const mppp::real128 &gd) { return audi::pow(base, gd); },
                   "Exponentiation (mppp::real128, mppp::real128).")
              .def("__rpow__", +[](const mppp::real128 &gd, double x) { return audi::pow(x, gd); },
                   "Exponentiation (double, mppp::real128).");

    // We expose the gdual<mppp::real128> arithmetic
    pyaudi::expose_gdual<mppp::real128>("real128")
        .def(bp::init<double>())
        .def(bp::init<double, const std::string &, unsigned int>())
        .def(bp::init<const std::string &>())
        .def(bp::init<const std::string &, const std::string &, unsigned int>());
    bp::def("exp", +[](const gdual_mp &d) { return exp(d); }, "Exponential (gdual_real128).");
    bp::def("log", +[](const gdual_mp &d) { return log(d); }, "Natural logarithm (gdual_real128).");
    bp::def("sqrt", +[](const gdual_mp &d) { return sqrt(d); }, "Square root (gdual_real128).");
    bp::def("cbrt", +[](const gdual_mp &d) { return cbrt(d); }, "Cubic root (gdual_real128).");
    bp::def("sin", +[](const gdual_mp &d) { return sin(d); }, "Sine (gdual_real128).");
    bp::def("asin", +[](const gdual_mp &d) { return asin(d); }, "Arc sine (gdual_real128).");
    bp::def("cos", +[](const gdual_mp &d) { return cos(d); }, "Cosine (gdual_real128).");
    bp::def("acos", +[](const gdual_mp &d) { return acos(d); }, "Arc cosine (gdual_real128).");
    bp::def("sin_and_cos", +[](const gdual_mp &d) { return pyaudi::v_to_l(sin_and_cos(d)); },
            "Sine and Cosine at once (gdual_real128).");
    bp::def("tan", +[](const gdual_mp &d) { return tan(d); }, "Tangent (gdual_real128).");
    bp::def("atan", +[](const gdual_mp &d) { return atan(d); }, "Arc tangent (gdual_real128).");
    bp::def("sinh", +[](const gdual_mp &d) { return sinh(d); }, "Hyperbolic sine (gdual_real128).");
    bp::def("asinh", +[](const gdual_mp &d) { return asinh(d); }, "Inverse hyperbolic sine (gdual_real128).");
    bp::def("cosh", +[](const gdual_mp &d) { return cosh(d); }, "Hyperbolic cosine (gdual_real128).");
    bp::def("acosh", +[](const gdual_mp &d) { return acosh(d); }, "Inverse hyperbolic cosine (gdual_real128).");
    bp::def("sinh_and_cosh", +[](const gdual_mp &d) { return pyaudi::v_to_l(sinh_and_cosh(d)); },
            "Hyperbolic sine and hyperbolic cosine at once (gdual_real128).");
    bp::def("tanh", +[](const gdual_mp &d) { return tanh(d); }, "Hyperbolic tangent (gdual_real128).");
    bp::def("atanh", +[](const gdual_mp &d) { return atanh(d); }, "Inverse hyperbolic arc tangent (gdual_real128).");
    bp::def("abs", +[](const gdual_mp &d) { return abs(d); }, "Absolute value (gdual_real128).");
    bp::def("erf", +[](const gdual_mp &d) { return erf(d); }, "Error function (gdual_real128).");

    // We expose simple mathematical functions to work on mppp::real128
    bp::def("exp", +[](mppp::real128 x) { return audi::exp(x); }, "Exponential (mppp::real128).");
    bp::def("log", +[](mppp::real128 x) { return audi::log(x); }, "Natural logarithm (mppp::real128).");
    bp::def("sqrt", +[](mppp::real128 x) { return audi::sqrt(x); }, "Square root (mppp::real128).");
    bp::def("cbrt", +[](mppp::real128 x) { return audi::cbrt(x); }, "Cubic root (mppp::real128).");
    bp::def("sin", +[](mppp::real128 x) { return audi::sin(x); }, "Sine (mppp::real128).");
    bp::def("asin", +[](mppp::real128 x) { return audi::asin(x); }, "Arc sine (mppp::real128).");
    bp::def("cos", +[](mppp::real128 x) { return audi::cos(x); }, "Cosine (mppp::real128).");
    bp::def("acos", +[](mppp::real128 x) { return audi::acos(x); }, "Arc cosine (mppp::real128).");
    bp::def("tan", +[](mppp::real128 x) { return audi::tan(x); }, "Tangent (mppp::real128).");
    bp::def("atan", +[](mppp::real128 x) { return audi::atan(x); }, "Arc tangent (mppp::real128).");
    bp::def("sinh", +[](mppp::real128 x) { return audi::sinh(x); }, "Hyperbolic sine (mppp::real128).");
    bp::def("asinh", +[](mppp::real128 x) { return audi::asinh(x); }, "Inverse hyperbolic sine (mppp::real128).");
    bp::def("cosh", +[](mppp::real128 x) { return audi::cosh(x); }, "Hyperbolic cosine (mppp::real128).");
    bp::def("acosh", +[](mppp::real128 x) { return audi::acosh(x); }, "Inverse hyperbolic cosine (mppp::real128).");
    bp::def("sinh_and_cosh", +[](const gdual_v &d) { return pyaudi::v_to_l(sinh_and_cosh(d)); },
            "Hyperbolic sine and hyperbolic cosine at once (gdual_v).");
    bp::def("tanh", +[](mppp::real128 x) { return audi::tanh(x); }, "Hyperbolic tangent (mppp::real128).");
    bp::def("atanh", +[](mppp::real128 x) { return audi::atanh(x); }, "Inverse hyperbolic arc tangent (mppp::real128).");
    bp::def("abs", +[](mppp::real128 x) { return audi::abs(x); }, "Absolute value (mppp::real128).");
    bp::def("erf", +[](mppp::real128 x) { return audi::erf(x); }, "Error function (mppp::real128).");

#endif

    // Miscellanea functions
    bp::def("invert_map",
            +[](const bp::object &map_in, bool verbose) {
                return pyaudi::v_to_l(invert_map(pyaudi::l_to_v<gdual_d>(map_in), verbose));
            },
            "Inverts a Taylor map (gdual_d)", (bp::arg("map"), bp::arg("verbose") = false));

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
