.. _installation_audi:

Packages
========

Conda
^^^^^

audi is available in the `conda <https://conda.io/en/latest/>`__ package manager from the
`conda-forge <https://conda-forge.org/>`__ channel. Two
packages are available:

* `audi <https://anaconda.org/conda-forge/audi>`__, which contains the header only library.
* `pyaudi <https://anaconda.org/conda-forge/pyaudi>`__, which contains the python bindings.

In order to install obake via audi, you just need
to add ``conda-forge`` to the channels:

.. code-block:: console

   $ conda config --add channels conda-forge
   $ conda install audi pyaudi

Please refer to the `conda documentation <https://conda.io/en/latest/>`__ for instructions on how to setup and manage
your conda installation.

Pip
^^^^^
audi is also available in the Python Package Index  but only for certain architectures.
You can check the `audi PyPi web page <https://pypi.org/project/pyaudi>`__ to see if yours is included.

From source
=============================================

C++ header only library
^^^^^^^^^^^^^^^^^^^^^^^

Audi is a header only library which has the following third party dependencies

* The `boost <http://www.boost.org/>`_ C++ libraries: the following boost libraries are necessary: boost_system, boost_unit_test_framework, boost_timer, boost_chrono. boost headers must be found in the system
* `obake <https://github.com/bluescarni/obake>`_
* `Eigen linear algebra library <https://eigen.tuxfamily.org/>`_: The Eigen headers must be found in the system

After making sure the dependencies above are installed in your system (most linux / osx package managers include them), you may download the latest Audi version via git:

.. code-block:: bash

   git clone https://github.com/darioizzo/audi.git

and configure your build using CMake.

.. note::

   The option AUDI_BUILD_AUDI should be selected, while the option AUDI_BUILD_PYAUDI should be deactivated. 
   You may also build the test by activating the option AUTI_BUILD_TESTS.
   
When done, type (in your build directory):

.. code-block:: bash

   make install

The headers will be installed in the CMAKE_INSTALL_PREFIX/include directory. 
To check that all went well compile the :ref:`quick-start example <getting_started>`.
