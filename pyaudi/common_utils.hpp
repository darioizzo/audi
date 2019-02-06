#ifndef PYAUDI_COMMON_UTILS_HPP
#define PYAUDI_COMMON_UTILS_HPP

#include <audi/vectorized.hpp>
#include <pybind11/pybind11.h>

namespace pybind11
{
namespace detail
{
template <typename T>
struct type_caster<audi::vectorized<T>> {
public:
    /**
     * This macro establishes the name 'inty' in
     * function signatures and declares a local variable
     * 'value' of type inty
     */
    PYBIND11_TYPE_CASTER(audi::vectorized<T>, _("vectorized<T>"));

    /**
     * Conversion part 1 (Python->C++): convert a PyObject into a inty
     * instance or return false upon failure. The second argument
     * indicates whether implicit conversions should be applied.
     */
    bool load(handle src, bool)
    {
        if (!isinstance<list>(src))
            return false;
        auto s = reinterpret_borrow<list>(src);
        value.clear();
        for (auto it : s) {
            value.push_back(it.cast<T>());
        }
        return true;
    }

    /**
     * Conversion part 2 (C++ -> Python): convert an inty instance into
     * a Python object. The second and third arguments are used to
     * indicate the return value policy and parent object (for
     * ``return_value_policy::reference_internal``) and are generally
     * ignored by implicit casters.
     */
    static handle cast(audi::vectorized<T> src, return_value_policy /* policy */, handle /* parent */)
    {
        list l(src.size());
        size_t index = 0;
        for (decltype(src.size()) i = 0; i < src.size(); ++i) {
            PyList_SET_ITEM(l.ptr(), i, pybind11::cast(src[i]).release().ptr()); // steals a reference
        }
        return l.release();
    }
};
} // namespace detail
} // namespace pybind11
#endif
