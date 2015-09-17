.. quickstart examples


Quick start examples
====================

.. contents::


C++
---

After following the installation instructions you will be able to compile and run your first C++ AuDi program:

.. _getting_started:

.. literalinclude:: ../../doc/examples/getting_started.cpp
   :language: c++
   :linenos:

Place it into a getting_started.cpp text file and compile it with:

.. code-block:: bash

   g++ -std=c++11 getting_started.cpp -lmpfr -lgmp

Running the executable you should see something like:

.. code-block:: bash

   Taylor polynomial: -0.049579*dy**7-38.9124*dx**2*dy+29.2889*dx**4*dy**2+461.017*dx+25.0185*dx**3*dy**2+122.112+29.0545*dx**5*dy**2-53.6544*dx**3*dy+18.0337*dx**2*dy**2-11.6457*dx**4*dy**3+9.42776*dx*dy**2+1602.74*dx**3-7.36459*dy+1006.88*dx**2+3.42634*dy**2-20.6462*dx*dy+1.5498*dx*dy**4+2142.92*dx**6-1.47471*dy**3-1.18951*dx**2*dy**5+...
   Derivative value: -349.372

-----------------------------------------------------------------------

Python
------
