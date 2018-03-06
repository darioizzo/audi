Utilities
========================

*#include <audi/invert_map.hpp>*

.. cpp:function:: std::vector<gdual_d> invert_map(const std::vector<gdual_d> &map_in, bool verbose = false)

   Inversion of Taylor series.

   Consider a system of non linear equations in :math:`n` variables:
   
   .. math:: 
   
      \left\{
      \begin{array}{l}
      f_1(x_1, x_2, .., x_n) = 0 \\
      ...                        \\
      f_n(x_1, x_2, .., x_n) = 0 \\
      \end{array}
      \right.
   
   and represent each equation via a Taylor polynomial truncated at order :math:`m` (i.e. :math:`\in \mathcal P_{n,m}`) around
   some point :math:`\overline {\mathbf x}`:
   
   .. math::
   
      \left\{
      \begin{array}{l}
      T_{f_1}(dx_1, dx_2, .., dx_n) = 0 \\
      ...                        \\
      T_{f_n}(dx_1, dx_2, .., dx_n) = 0 \\
      \end{array}
      \right.
   
   The above system of, now, polynomial equations can be written as :math:`C + M + N = C + \mathcal M = 0`, where we have
   explicitly indicated the constant term :math:`C`, the linear term :math:`M` and the rest :math:`N`. The symbol
   :math:`\mathcal M` indicates what is called Taylor map and its inverse is found by this function by performing
   some Picard iterations (fixed point) over the following simple formal equalities:
     
   .. math::
   
      \mathcal M \circ \mathcal M^{-1} = \mathcal I \rightarrow  (M + N) \circ \mathcal M^{-1}  = \mathcal I \rightarrow M
      \circ \mathcal M^{-1} + N \circ \mathcal M^{-1} = \mathcal I 
      
   and hence: 
   
   .. math::
   
      \mathcal M^{-1} = M^{-1} \circ \left(\mathcal I - N \circ \mathcal M^{-1}\right)
      
   which can be used to make a fixed point iterative scheme once
   the inverse of the linear part :math:`M^{-1}` is found. Note that the inverse is guaranteed to exist if the linear part is invertible.
     
   :param map_in: The input map represented as an *std::vector* of *gdual<double>*. They can have a constant coefficient which will be neglected
   :param verbose: when true some output is shown during the iterations
   
   :return: The inverse map :math:`\mathcal M^{-1}` represented as an *std::vector* of *gduals* with symbol set
      [p0, p1, .. pn]
   
   :exception: std::invalid_argument if *map_in* is empty, if its linear part is not invertable, if the symbol set of the
      map components are not all equal, if the map is not square (i.e. it is of size n with n symbols) or if the order of
      the gduals are not all the same.

   .. code-block:: c++

      gdual<double> x(0., "x1", 9);
      auto sinx = sin(x);
      auto asinx = asin(x);
      auto inv_sinx = invert_map({sinx}); 
      std::cout << "sin x:\t\t" << sinx << std::endl;
      std::cout << "asin x:\t\t" << asinx << std::endl;
      std::cout << "inv sinx:\t" << inv_sinx[0] << std::endl;

   The code above will produce the output

   .. code-block:: c++

      sin x:        dx1+2.75573e-06*dx1**9-0.166667*dx1**3+0.00833333*dx1**5-0.000198413*dx1**7
      asin x:       dx1+0.0303819*dx1**9+0.166667*dx1**3+0.075*dx1**5+0.0446429*dx1**7
      inv sinx:     dp0+0.0303819*dp0**9+0.166667*dp0**3+0.075*dp0**5+0.0446429*dp0**7



