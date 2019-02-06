#include <audi/audi.hpp>
#include <cmath>
#include <pybind11/pybind11.h>
#include <sstream>

#if defined(AUDI_WITH_MPPP)
#include <audi/real128.hpp>
#include <mp++/extra/pybind11.hpp>
#include <mp++/mp++.hpp>
#endif

#include "common_utils.hpp"
#include "expose_gdual.hpp"

using namespace audi;
namespace py = pybind11;

PYBIND11_MODULE(core, m)
{

    // We expose simple mathematical functions to work on python floats
    m.def("exp", [](double x) { return std::exp(x); }, "Exponential (double).");
    m.def("log", [](double x) { return std::log(x); }, "Natural logarithm (double).");
    m.def("sqrt", [](double x) { return std::sqrt(x); }, "Square root (double).");
    m.def("cbrt", [](double x) { return std::cbrt(x); }, "Cubic root (double).");
    m.def("sin", [](double x) { return std::sin(x); }, "Sine (double).");
    m.def("asin", [](double x) { return std::asin(x); }, "Arc sine (double).");
    m.def("cos", [](double x) { return std::cos(x); }, "Cosine (double).");
    m.def("acos", [](double x) { return std::acos(x); }, "Arc cosine (double).");
    m.def("tan", [](double x) { return std::tan(x); }, "Tangent (double).");
    m.def("atan", [](double x) { return std::atan(x); }, "Arc tangent (double).");
    m.def("sinh", [](double x) { return std::sinh(x); }, "Hyperbolic sine (double).");
    m.def("asinh", [](double x) { return std::asinh(x); }, "Inverse hyperbolic sine (double).");
    m.def("cosh", [](double x) { return std::cosh(x); }, "Hyperbolic cosine (double).");
    m.def("acosh", [](double x) { return std::acosh(x); }, "Inverse hyperbolic cosine (double).");
    m.def("tanh", [](double x) { return std::tanh(x); }, "Hyperbolic tangent (double).");
    m.def("atanh", [](double x) { return std::atanh(x); }, "Inverse hyperbolic arc tangent (double).");
    m.def("abs", [](double x) { return std::abs(x); }, "Absolute value (double).");
    m.def("erf", [](double x) { return std::erf(x); }, "Error function (double).");

    // We expose the class gdual<double> or gdual_d and its arithmetic
    pyaudi::expose_gdual<double>(m, "double");

    m.def("exp", [](const gdual_d &d) { return exp(d); }, "Exponential (gdual_double).");
    m.def("log", [](const gdual_d &d) { return log(d); }, "Natural logarithm (gdual_d).");
    m.def("sqrt", [](const gdual_d &d) { return sqrt(d); }, "Square root (gdual_d).");
    m.def("cbrt", [](const gdual_d &d) { return cbrt(d); }, "Cubic root (gdual_d).");
    m.def("sin", [](const gdual_d &d) { return sin(d); }, "Sine (gdual_d).");
    m.def("asin", [](const gdual_d &d) { return asin(d); }, "Arc sine (gdual_d).");
    m.def("cos", [](const gdual_d &d) { return cos(d); }, "Cosine (gdual_d).");
    m.def("acos", [](const gdual_d &d) { return acos(d); }, "Arc cosine (gdual_d).");
    m.def("sin_and_cos", [](const gdual_d &d) { return sin_and_cos(d); }, "Sine and Cosine at once (gdual_d).");
    m.def("tan", [](const gdual_d &d) { return tan(d); }, "Tangent (gdual_d).");
    m.def("atan", [](const gdual_d &d) { return atan(d); }, "Arc tangent (gdual_d).");
    m.def("sinh", [](const gdual_d &d) { return sinh(d); }, "Hyperbolic sine (gdual_d).");
    m.def("asinh", [](const gdual_d &d) { return asinh(d); }, "Inverse hyperbolic sine (gdual_d).");
    m.def("cosh", [](const gdual_d &d) { return cosh(d); }, "Hyperbolic cosine (gdual_d).");
    m.def("acosh", [](const gdual_d &d) { return acosh(d); }, "Inverse hyperbolic cosine (gdual_d).");
    m.def("sinh_and_cosh", [](const gdual_d &d) { return sinh_and_cosh(d); },
          "Hyperbolic sine and hyperbolic cosine at once (gdual_d).");
    m.def("tanh", [](const gdual_d &d) { return tanh(d); }, "Hyperbolic tangent (gdual_d).");
    m.def("atanh", [](const gdual_d &d) { return atanh(d); }, "Inverse hyperbolic arc tangent (gdual_d).");
    m.def("abs", [](const gdual_d &d) { return abs(d); }, "Absolute value (gdual_d).");
    m.def("erf", [](const gdual_d &d) { return erf(d); }, "Error function (gdual_d).");

    // We expose the class gdual<vectorized<double>> or gdual_v and its arithmetic. Two constructors are added to allow
    // construction from vector::double
    pyaudi::expose_gdual<audi::vectorized<double>>(m, "vdouble")
        .def(py::init<std::vector<double>>())
        .def(py::init<std::vector<double>, const std::string &, unsigned int>());

    m.def("exp", [](const gdual_v &d) { return exp(d); }, "Exponential (gdual_vdouble).");
    m.def("log", [](const gdual_v &d) { return log(d); }, "Natural logarithm (gdual_v).");
    m.def("sqrt", [](const gdual_v &d) { return sqrt(d); }, "Square root (gdual_v).");
    m.def("cbrt", [](const gdual_v &d) { return cbrt(d); }, "Cubic root (gdual_v).");
    m.def("sin", [](const gdual_v &d) { return sin(d); }, "Sine (gdual_v).");
    m.def("asin", [](const gdual_v &d) { return asin(d); }, "Arc sine (gdual_v).");
    m.def("cos", [](const gdual_v &d) { return cos(d); }, "Cosine (gdual_v).");
    m.def("acos", [](const gdual_v &d) { return acos(d); }, "Arc cosine (gdual_v).");
    m.def("sin_and_cos", [](const gdual_v &d) { return sin_and_cos(d); }, "Sine and Cosine at once (gdual_v).");
    m.def("tan", [](const gdual_v &d) { return tan(d); }, "Tangent (gdual_v).");
    m.def("atan", [](const gdual_v &d) { return atan(d); }, "Arc tangent (gdual_v).");
    m.def("sinh", [](const gdual_v &d) { return sinh(d); }, "Hyperbolic sine (gdual_v).");
    m.def("asinh", [](const gdual_v &d) { return asinh(d); }, "Inverse hyperbolic sine (gdual_v).");
    m.def("sinh_and_cosh", [](const gdual_v &d) { return sinh_and_cosh(d); },
          "Hyperbolic sine and hyperbolic cosine at once (gdual_v).");
    m.def("cosh", [](const gdual_v &d) { return cosh(d); }, "Hyperbolic cosine (gdual_v).");
    m.def("acosh", [](const gdual_v &d) { return acosh(d); }, "Inverse hyperbolic cosine (gdual_v).");
    m.def("tanh", [](const gdual_v &d) { return tanh(d); }, "Hyperbolic tangent (gdual_v).");
    m.def("atanh", [](const gdual_v &d) { return atanh(d); }, "Inverse hyperbolic arc tangent (gdual_v).");
    m.def("abs", [](const gdual_v &d) { return abs(d); }, "Absolute value (gdual_v).");
    m.def("erf", [](const gdual_v &d) { return erf(d); }, "Error function (gdual_v).");

#if defined(AUDI_WITH_MPPP)
    // Init the pybind11 integration of mppp.
    mppp_pybind11::init();

    // We expose simple mathematical functions to work on mppp::real128
    m.def("exp", [](mppp::real128 x) { return audi::exp(x); }, "Exponential (mppp::real128).");
    m.def("log", [](mppp::real128 x) { return audi::log(x); }, "Natural logarithm (mppp::real128).");
    m.def("sqrt", [](mppp::real128 x) { return audi::sqrt(x); }, "Square root (mppp::real128).");
    m.def("cbrt", [](mppp::real128 x) { return audi::cbrt(x); }, "Cubic root (mppp::real128).");
    m.def("sin", [](mppp::real128 x) { return audi::sin(x); }, "Sine (mppp::real128).");
    m.def("asin", [](mppp::real128 x) { return audi::asin(x); }, "Arc sine (mppp::real128).");
    m.def("cos", [](mppp::real128 x) { return audi::cos(x); }, "Cosine (mppp::real128).");
    m.def("acos", [](mppp::real128 x) { return audi::acos(x); }, "Arc cosine (mppp::real128).");
    m.def("tan", [](mppp::real128 x) { return audi::tan(x); }, "Tangent (mppp::real128).");
    m.def("atan", [](mppp::real128 x) { return audi::atan(x); }, "Arc tangent (mppp::real128).");
    m.def("sinh", [](mppp::real128 x) { return audi::sinh(x); }, "Hyperbolic sine (mppp::real128).");
    m.def("asinh", [](mppp::real128 x) { return audi::asinh(x); }, "Inverse hyperbolic sine (mppp::real128).");
    m.def("cosh", [](mppp::real128 x) { return audi::cosh(x); }, "Hyperbolic cosine (mppp::real128).");
    m.def("acosh", [](mppp::real128 x) { return audi::acosh(x); }, "Inverse hyperbolic cosine (mppp::real128).");
    m.def("tanh", [](mppp::real128 x) { return audi::tanh(x); }, "Hyperbolic tangent (mppp::real128).");
    m.def("atanh", [](mppp::real128 x) { return audi::atanh(x); }, "Inverse hyperbolic arc tangent (mppp::real128).");
    m.def("abs", [](mppp::real128 x) { return audi::abs(x); }, "Absolute value (mppp::real128).");
    m.def("erf", [](mppp::real128 x) { return audi::erf(x); }, "Error function (mppp::real128).");

    // We expose the class gdual<mppp::real128> and its arithmetic. Some constructors are added to allow
    // construction from std::string and doubles (with precision loss)
    pyaudi::expose_gdual<mppp::real128>(m, "real128")
        .def(py::init<double>())
        .def(py::init<double, const std::string &, unsigned int>())
        .def(py::init<const std::string &>())
        .def(py::init<const std::string &, const std::string &, unsigned int>());
    m.def("exp", [](const gdual_mp &d) { return exp(d); }, "Exponential (gdual_real128).");
    m.def("log", [](const gdual_mp &d) { return log(d); }, "Natural logarithm (gdual_real128).");
    m.def("sqrt", [](const gdual_mp &d) { return sqrt(d); }, "Square root (gdual_real128).");
    m.def("cbrt", [](const gdual_mp &d) { return cbrt(d); }, "Cubic root (gdual_real128).");
    m.def("sin", [](const gdual_mp &d) { return sin(d); }, "Sine (gdual_real128).");
    m.def("asin", [](const gdual_mp &d) { return asin(d); }, "Arc sine (gdual_real128).");
    m.def("cos", [](const gdual_mp &d) { return cos(d); }, "Cosine (gdual_real128).");
    m.def("acos", [](const gdual_mp &d) { return acos(d); }, "Arc cosine (gdual_real128).");
    m.def("sin_and_cos", [](const gdual_mp &d) { return sin_and_cos(d); }, "Sine and Cosine at once (gdual_real128).");
    m.def("tan", [](const gdual_mp &d) { return tan(d); }, "Tangent (gdual_real128).");
    m.def("atan", [](const gdual_mp &d) { return atan(d); }, "Arc tangent (gdual_real128).");
    m.def("sinh", [](const gdual_mp &d) { return sinh(d); }, "Hyperbolic sine (gdual_real128).");
    m.def("asinh", [](const gdual_mp &d) { return asinh(d); }, "Inverse hyperbolic sine (gdual_real128).");
    m.def("cosh", [](const gdual_mp &d) { return cosh(d); }, "Hyperbolic cosine (gdual_real128).");
    m.def("acosh", [](const gdual_mp &d) { return acosh(d); }, "Inverse hyperbolic cosine (gdual_real128).");
    m.def("sinh_and_cosh", [](const gdual_mp &d) { return sinh_and_cosh(d); },
          "Hyperbolic sine and hyperbolic cosine at once (gdual_real128).");
    m.def("tanh", [](const gdual_mp &d) { return tanh(d); }, "Hyperbolic tangent (gdual_real128).");
    m.def("atanh", [](const gdual_mp &d) { return atanh(d); }, "Inverse hyperbolic arc tangent (gdual_real128).");
    m.def("abs", [](const gdual_mp &d) { return abs(d); }, "Absolute value (gdual_real128).");
    m.def("erf", [](const gdual_mp &d) { return erf(d); }, "Error function (gdual_real128).");
#endif

    // Miscellanea functions
    m.def("invert_map", [](const std::vector<gdual_d> &map_in, bool verbose) { return invert_map(map_in, verbose); },
          "Inverts a Taylor map (gdual_d)", py::arg("map"), py::arg("verbose") = false);

    // Define a cleanup functor to be run when the module is unloaded.
    struct audi_cleanup_functor {
        void operator()() const
        {
            std::cout << "Shutting down the thread pool.\n";
            piranha::thread_pool_shutdown<void>();
        }
    };
    // Expose it.
    py::class_<audi_cleanup_functor>(m, "_audi_cleanup_functor")
        .def(py::init<>())
        .def("__call__", &audi_cleanup_functor::operator());
    // Register it.
    py::module atexit_mod = py::module::import("atexit");
    atexit_mod.attr("register")(audi_cleanup_functor{});
}