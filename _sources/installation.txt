.. _installationguide:


Installation guide
==================

.. contents::


C++
---

AuDi is a header only library which has the following third party dependencies

* The `boost <http://www.boost.org/>`_ C++ libraries: the following boost libraries are necessary: boost_system, boost_unit_test_framework, boost_timer, boost_chrono. boost headers must be found in the system
* `piranha <http://bluescarni.github.io/piranha/index.html>`_: piranha headers must be found in the system
* `GNU MPFR library <http://www.mpfr.org/>`_: needed by piranha
* `The GMP multiprecision library <https://gmplib.org/>`_: needed by piranha

After making sure the dependencies above are installed in your system (most linux / osx package managers include them), you may download the latest AuDi version via git:

.. code-block:: bash

   git clone https://github.com/darioizzo/audi.git

cand onfigure your build using CMake. When done, type (in your build directory):

.. code-block:: bash

   make install

The headers will be installed in the CMAKE_INSTALL_PREFIX/include directory. To check that all went well compile the :ref:`quick-start example <getting_started>`.

-----------------------------------------------------------------------

Python
------
