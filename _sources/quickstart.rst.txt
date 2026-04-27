.. quickstart examples


Quick start
===========

Python
------

My first program
^^^^^^^^^^^^^^^^^^

If you have successfully compiled and installed pyaudi following the :ref:`installation_pyaudi` you will be able to test its use typing the following script.

.. literalinclude:: ../../doc/examples/getting_started.py
   :language: python
   :linenos:

Place it into a getting_started.py text file and run it with 

.. code-block:: bash

   python getting_started.py

We recommend the use of Jupyter or ipython do enjoy pyaudi the most. 

------------------------------------------------------------------------------------

C++
---

My first program
^^^^^^^^^^^^^^^^^^

After following the :ref:`installation_audi` you will be able to compile and run your first C++ AuDi program:

.. _getting_started:

.. literalinclude:: ../../doc/examples/getting_started.cpp
   :language: c++
   :linenos:

Place it into a getting_started.cpp text file. Compiling such an example file with C++ manually
means making sure that all dependencies (i.e. headers and libraries) are on the search path of the
compiler you are using (e.g. g++). On possible example of compiling the program could be:

.. code-block:: bash

   g++ -std=c++11 getting_started.cpp -lmpfr -lgmp -pthread


another example with manual specification (and clang++) could be:

.. code-block:: bash

   clang++ audi_test.cpp \
     -std=c++20 \
     -I$CONDA_PREFIX/include \
     -L$CONDA_PREFIX/lib \
     -lmpfr -lgmp \
     -lboost_system \
     -lobake \
     -ltbb \
     -lmp++ \
     -lfmt \
     -pthread \
     -o audi_test



