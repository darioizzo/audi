.. _installation_pyaudi:

Python bindings
^^^^^^^^^^^^^^^^

The main functionalities of AuDi are exposed into a python module called pyaudi. To compile the module you need to have
the following dependencies installed in your system

* The `boost <http://www.boost.org/>`_ C++ libraries: the following boost libraries are necessary: boost_python, boost_system, boost_unit_test_framework, boost_timer, boost_chrono. boost headers must be found in the system
* `obake <https://github.com/bluescarni/obake>`_
* `audi <https://github.com/darioizzo/audi>`_: audi headers must be found in the system
* `Eigen linear algebra library <https://eigen.tuxfamily.org/>`_: The Eigen headers must be found in the system

You can then clone the `audi github repository <https://github.com/darioizzo/audi>`_  and configure it using cmake.
Make sure to activate the AUDI_BUILD_PYAUDI option (you will also need to install the audi c++ headers first, so AUDI_BUILD_AUDI option
will likely have to be used first).

Check carefully what python version is detected and what libraries are linked to. 

The CMAKE_INSTALL_PREFIX will be used to construct the final location of headers and python module after install.

When done, type (in your build directory):

.. code-block:: bash

   make install

To check that all went well fire-up your python console and try the example in :ref:`quick-start example <getting_started>`.