#ifndef PYAUDI_COMMON_UTILS_HPP
#define PYAUDI_COMMON_UTILS_HPP

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <unordered_map>

#include <audi/audi.hpp>

// A throwing macro similar to audi_throw, only for Python. This will set the global
// error string of Python to "msg", the exception type to "type", and then invoke the Boost
// Python function to raise the Python exception.
#define pyaudi_throw(type, msg)                                                                                        \
    PyErr_SetString(type, msg);                                                                                        \
    boost::python::throw_error_already_set();                                                                          \
    throw

namespace bp = boost::python;

namespace pyaudi
{
// Wrapper around the CPython function to create a bytes object from raw data.
bp::object make_bytes(const char *ptr, Py_ssize_t len)
{
    PyObject *retval;
    if (len) {
        retval = PyBytes_FromStringAndSize(ptr, len);
    } else {
        retval = PyBytes_FromStringAndSize(nullptr, 0);
    }
    if (!retval) {
        pyaudi_throw(PyExc_RuntimeError, "unable to create a bytes object: the 'PyBytes_FromStringAndSize()' "
                                         "function returned NULL");
    }
    return bp::object(bp::handle<>(retval));
}

// Converts a C++ vector to a python list
template <typename T>
inline bp::list v_to_l(std::vector<T> vector)
{
    bp::list list;
    for (auto iter = vector.begin(); iter != vector.end(); ++iter) {
        list.append(*iter);
    }
    return list;
}

// Converts a C++ array to a python list
template <typename T, std::size_t dim>
inline bp::list v_to_l(std::array<T, dim> vector)
{
    bp::list list;
    for (auto iter = vector.begin(); iter != vector.end(); ++iter) {
        list.append(*iter);
    }
    return list;
}

// Converts a C++ unordered map to a python dict (TO BE TESTED)
template <typename K, typename V>
inline bp::dict umap_to_pydict(std::map<K, V> map)
{
    boost::python::dict dictionary;
    for (auto iter = map.begin(); iter != map.end(); ++iter) {
        dictionary[iter->first] = iter->second;
    }
    return dictionary;
}

// Converts a python dict to an std::unordered_map if the Value type is double or unsigned int
template <typename Key_type, typename Value_type>
inline std::unordered_map<Key_type, Value_type> pydict_to_umap(const bp::dict &dict)
{
    std::unordered_map<Key_type, Value_type> retval;
    bp::stl_input_iterator<Key_type> it(dict), end;
    for (; it != end; ++it) {
        retval[*it] = bp::extract<Value_type>(dict[*it])();
    }
    return retval;
}

// Converts a python iterable to an std::vector
template <typename T>
inline std::vector<T> l_to_v(const bp::object &iterable)
{
    return std::vector<T>(bp::stl_input_iterator<T>(iterable), bp::stl_input_iterator<T>());
}

// Used to register the vectorized to list converter
template <typename T>
struct vectorized_to_python_list {
    static PyObject *convert(const audi::vectorized<T> &vd)
    {
        bp::list list;
        for (auto iter = vd.begin(); iter != vd.end(); ++iter) {
            list.append(*iter);
        }
        return boost::python::incref(list.ptr());
    }
};
}

#endif
