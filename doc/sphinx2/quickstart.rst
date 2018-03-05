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

Notebooks
^^^^^^^^^

Follow the links below to visualize juypiter notebooks on the use of pyaudi.

- `The very basics <https://github.com/darioizzo/audi/blob/master/examples/example00.ipynb>`_

- `Understanding gduals and floats <https://github.com/darioizzo/audi/blob/master/examples/example01.ipynb>`_:

- `Differential Intelligence <https://github.com/darioizzo/audi/blob/master/examples/example10.ipynb>`_

- `Training an artificial neural network <https://github.com/darioizzo/audi/blob/master/examples/example11.ipynb>`_

- Tutorial on High Order Taylor Maps <https://github.com/darioizzo/audi/blob/master/examples/example12.ipynb>`_

-----------------------------------------------------------------------


C++
---

My first program
^^^^^^^^^^^^^^^^^^

After following the :ref:`installation_audi` you will be able to compile and run your first C++ AuDi program:

.. _getting_started:

.. literalinclude:: ../../doc/examples/getting_started.cpp
   :language: c++
   :linenos:

Place it into a getting_started.cpp text file and compile it with:

.. code-block:: bash

   g++ -std=c++11 getting_started.cpp -lmpfr -lgmp -pthread




