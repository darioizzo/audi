.. quickstart examples


Quick start examples
====================

.. contents::


C++
---

After following the :ref:`installationguide` you will be able to compile and run your first C++ AuDi program:

.. _getting_started:

.. literalinclude:: ../../doc/examples/getting_started.cpp
   :language: c++
   :linenos:

Place it into a getting_started.cpp text file and compile it with:

.. code-block:: bash

   g++ -std=c++11 getting_started.cpp -lmpfr -lgmp -pthread

-----------------------------------------------------------------------

Python
------

If you have succesfully compiled and installed pyaudi following the :ref:`installationguide` you will be able to test its use typing the following script.

.. literalinclude:: ../../doc/examples/getting_started.py
   :language: python
   :linenos:

Place it into a getting_started.py text file and run it with 

.. code-block:: bash

   python getting_started.py

We reccomend the use of Jupyter or ipython do enjoy pyaudi the most. 

Notebooks
^^^^^^^^^

Follow the links below to visualize juypiter notebooks on the use of pyaudi.

Basic
"""""

- `The very basics <http://localhost:8888/notebooks/audi/examples/example00.ipynb>`_: (by Francesco Biscani and Dario Izzo)

Advanced
""""""""

- `Training an artificial neural network <http://localhost:8888/notebooks/audi/examples/example10.ipynb>`_: (by Carlos Sanchez)

