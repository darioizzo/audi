.. _installation_pyaudi:


Installation (pyaudi)
===========================

The main functionalities of AuDi are exposed into a python module called pyaudi. To create the module you need to have
the boost python libraries installed and activate the BUILD_PYAUDI option from within cmake.

Check carefully what python version is detected and what libraries are linked to. In particular select the correct boost_python
according to the python version (2 or 3) you want to compile the module for.

The CMAKE_INSTALL_PREFIX will be used to construct the final location of headers and python module after install.

When done, type (in your build directory):

.. code-block:: bash

   make install

To check that all went well fire-up your python console and try the example in :ref:`quick-start example <getting_started>`.

Binaries
--------

For some architectures (i.e. windows with 64 bit python 27, 34 and 35), we provide pre-built binaries for pyaudi via PyPi. If you have
a system compatible with those requirements, you can just type:

.. code-block:: bash

   pip install pyaudi
