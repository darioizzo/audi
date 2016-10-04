#ifndef PYAUDI_COMMON_UTILS_HPP
#define PYAUDI_COMMON_UTILS_HPP

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <unordered_map>

// A throwing macro similar to pagmo_throw, only for Python. This will set the global
// error string of Python to "msg", the exception type to "type", and then invoke the Boost
// Python function to raise the Python exception.
#define pyaudi_throw(type,msg) \
PyErr_SetString(type,msg); \
boost::python::throw_error_already_set(); \
throw

namespace bp = boost::python;

namespace pyaudi{
// Wrapper around the CPython function to create a bytes object from raw data.
bp::object make_bytes(const char *ptr, Py_ssize_t len)
{
    PyObject *retval;
    if (len) {
        retval = PyBytes_FromStringAndSize(ptr,len);
    } else {
        retval = PyBytes_FromStringAndSize(nullptr,0);
    }
    if (!retval) {
        pyaudi_throw(PyExc_RuntimeError,"unable to create a bytes object: the 'PyBytes_FromStringAndSize()' "
            "function returned NULL");
    }
    return bp::object(bp::handle<>(retval));
}

// Converts a C++ vector to a python list
template <typename T>
inline bp::list v_to_l(std::vector<T> vector) {
    typename std::vector<T>::iterator iter;
    bp::list list;
    for (iter = vector.begin(); iter != vector.end(); ++iter) {
        list.append(*iter);
    }
    return list;
}

// Converts a C++ unordered map to a python dict (TO BE TESTED)
template <typename K, typename V>
inline bp::dict umap_to_pydict(std::map<K, V> map) {
    typename std::unordered_map<K, V>::iterator iter;
	boost::python::dict dictionary;
	for (iter = map.begin(); iter != map.end(); ++iter) {
		dictionary[iter->first] = iter->second;
	}
	return dictionary;
}

// Converts a python dict to an std::unordered_map
template<typename Key_type, typename Value_type>
inline std::unordered_map<Key_type, Value_type> pydict_to_umap(const bp::dict& pydict)
{
    std::unordered_map<Key_type,Value_type> retval;
    bp::list keys = pydict.keys();
    for (int i = 0; i < len(keys); ++i) {
       bp::extract<Key_type> extracted_key(keys[i]);
       if(!extracted_key.check()){
            std::cout<<"Key invalid, map might be incomplete"<<std::endl;
            continue;
       }
       std::string key = extracted_key;
       bp::extract<Value_type> extracted_val(pydict[key]);
       if(!extracted_val.check()){
       std::cout<<"Value invalid, map might be incomplete"<<std::endl;
            continue;
       }
       Value_type value = extracted_val;
       retval[key] = value;
  }
   return retval;
}


// Converts a python iterable to an std::vector
template<typename T>
inline std::vector<T> l_to_v(const bp::object& iterable)
{
    return std::vector<T>( bp::stl_input_iterator<T>(iterable), bp::stl_input_iterator<T>());
}

}

#endif
