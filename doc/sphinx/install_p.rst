.. _installation_pyaudi:

Installation (Python package)
===========================================

Pre-compiled pyaudi binaries are available both from the Pyhton Package Index (`PyPi <https://pypi.python.org/pypi/>`_)
and from `conda-forge <https://conda-forge.org/>`_. Not all architectures are supported, namely only win64 (PyPi), 
linux 64 (PyPi and conda) and osx (only conda). The best is to try the following:

.. code-block:: bash

   conda config --add channels conda-forge 
   conda install pyaudi

or

.. code-block:: bash

   pip install pyaudi --user

if none work, then, I am afraid you will need to compile pyaudi.

Compiling (pyaudi)
------------------

The main functionalities of AuDi are exposed into a python module called pyaudi. To compile the module you need to have
the following dependencies installed in your system

* The `boost <http://www.boost.org/>`_ C++ libraries: the following boost libraries are necessary: boost_python, boost_system, boost_unit_test_framework, boost_timer, boost_chrono. boost headers must be found in the system
* `obake <https://github.com/bluescarni/obake>`_
* `audi <https://github.com/darioizzo/audi>`_: audi headers must be found in the system
* `Eigen linear algebra library <https://eigen.tuxfamily.org/>`_: The Eigen headers must be found in the system

You can then clone the `audi github repository <https://github.com/darioizzo/audi>`_  and configure it using cmake.
Make sure to activate the AUDI_BUILD_PYAUDI option (you will also need to install the audi c++ headers first, so AUDI_BUILD_AUDI option
will likely have to be used first).

Check carefully what python version is detected and what libraries are linked to. In particular select the
correct boost_python according to the python version (2 or 3) you want to compile the module for.

The CMAKE_INSTALL_PREFIX will be used to construct the final location of headers and python module after install.

When done, type (in your build directory):

.. code-block:: bash

   make install

To check that all went well fire-up your python console and try the example in :ref:`quick-start example <getting_started>`.