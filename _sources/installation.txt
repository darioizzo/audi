.. AuDi installation guide


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

After making sure the dependencies above are installed in your system (most common package managers can find them), you can download the latest AuDi version via git:

.. code-block:: bash

   git clone https://github.com/darioizzo/audi.git

and just copy the headers in the directory src into your favorite path location. You can also let CMake do the dirty work by typing (after configuring)

.. code-block:: bash

   make install

The headers should now be found in your system. To check that all went well compile the quick-start example.

-----------------------------------------------------------------------

Python
------
