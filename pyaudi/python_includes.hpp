#ifndef PYAUDI_PYTHON_INCLUDES_HPP
#define PYAUDI_PYTHON_INCLUDES_HPP

// NOTE: the order of inclusion in the first two items here is forced by these two issues:
// http://mail.python.org/pipermail/python-list/2004-March/907592.html
// http://mail.python.org/pipermail/new-bugs-announce/2011-March/010395.html
#if defined(_WIN32)
#include <Python.h>
#include <cmath>
#else
#include <Python.h>
#include <cmath>
#endif

#endif
