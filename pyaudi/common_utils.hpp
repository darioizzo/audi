#ifndef PYAUDI_COMMON_UTILS_HPP
#define PYAUDI_COMMON_UTILS_HPP

#include <algorithm>
#include <audi/vectorized.hpp>
#include <obake/math/safe_cast.hpp>
#include <obake/symbols.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace pyaudi
{
// Convert a dictionary into a symbol map of type T.
template <typename T>
inline ::obake::symbol_map<T> py_dict_to_obake_sm(const py::dict &d)
{
    typename ::obake::symbol_map<T>::sequence_type seq;
    seq.reserve(::obake::safe_cast<decltype(seq.size())>(py::len_hint(d)));

    for (const auto [k, v] : d) {
        seq.emplace_back(k.template cast<::std::string>(), v.template cast<T>());
    }

    ::std::sort(seq.begin(), seq.end(), [](const auto &p1, const auto &p2) { return p1.first < p2.first; });

    ::obake::symbol_map<T> retval;
    retval.adopt_sequence(::boost::container::ordered_unique_range_t{}, ::std::move(seq));

    return retval;
}
} // namespace pyaudi

// We add a specialization for the type audi::vectorized<T> so that
// we can easily have functions with this type as input and output
namespace pybind11
{
namespace detail
{
template <>
struct type_caster<audi::vectorized<double>> {
public:
    /**
     * This macro establishes the name 'vectorized<double>' in
     * function signatures and declares a local variable
     * 'value' of type vectorized<double>
     */
    PYBIND11_TYPE_CASTER(audi::vectorized<double>, _("vectorized<double>"));

    /**
     * Conversion part 1 (Python->C++): convert a PyObject into a vectorized<T>
     * instance or return false upon failure. The second argument
     * indicates whether implicit conversions should be applied.
     */
    bool load(handle src, bool)
    {
        if (!isinstance<list>(src)) return false;
        auto s = reinterpret_borrow<list>(src);
        value.m_c.clear();
        for (auto it : s) {
            value.m_c.push_back(it.cast<double>());
        }
        return true;
    }

    /**
     * Conversion part 2 (C++ -> Python): convert a vectorized<double> instance into
     * a Python object. The second and third arguments are used to
     * indicate the return value policy and parent object (for
     * ``return_value_policy::reference_internal``) and are generally
     * ignored by implicit casters.
     */
    static handle cast(audi::vectorized<double> src, return_value_policy /* policy */, handle /* parent */)
    {
        list l(src.size());
        for (decltype(src.size()) i = 0; i < src.size(); ++i) {
            PyList_SET_ITEM(l.ptr(), i, pybind11::cast(src[i]).release().ptr()); // steals a reference
        }
        return l.release();
    }
};
} // namespace detail
} // namespace pybind11
#endif
