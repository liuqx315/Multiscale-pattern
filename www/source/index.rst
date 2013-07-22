.. _index:


The ARKode Solver 
===================

ARKode is a solver for stiff, nonstiff and multi-rate systems of
ordinary differential equations (ODEs) given in explicit form,

.. math::
   M \dot{y} = f_E(t,y) + f_I(t,y).

Here,

* :math:`t` is the independent variable, e.g. time,
* :math:`y` is the set of dependent variables (in :math:`\Re^N`),  
* :math:`M` is a user-specified, non-singular operator from
  :math:`\Re^N \to \Re^N`  (the "mass matrix", possibly time
  dependent, but independent of :math:`y`),
* :math:`f_E(t,y)` is the portion of the right-hand side containing
  the "slow" time scale components in the system that should be
  integrated explicitly, and 
* :math:`f_I(t,y)` is the portion of the right-hand side containing
  the "fast" time scale components to be integrated implicitly.

Either of the :math:`f_E` or :math:`f_I` operators may be disabled,
allowing for fully explicit, fully implicit, or combination
implicit-explicit (IMEX) time integration.

ARKode is a component of the `SUNDIALS
<https://computation.llnl.gov/casc/sundials/main.html>`_ suite of 
nonlinear and differential/algebraic equation solvers. 

ARKode is written in C, with C++ and Fortran interfaces.


News:
-------

:7/23/2013:
   ARKode is currently in beta-testing, undergoing optimization
   enhancements, documentation upgrades, and development of new
   features. An official release as part of the SUNDIALS suite is
   planned for Fall 2013, but the public source code repository is now
   open (see the :ref:`Downloads` page).


Contact
----------

`Daniel Reynolds <http://faculty.smu.edu/reynolds/>`_, 
`Department of Mathematics <http://smu.edu/math/>`_, 
`Southern Methodist University <http://www.smu.edu/>`_


Support
----------

This work is supported by the `U.S. Department of Energy
<http://science.energy.gov/ascr>`_ through the `FASTMath SciDAC
Institute <http://fastmath-scidac.org>`_, under subcontract B598130
from `Lawrence Livermore National Laboratory <http://www.llnl.gov>`_.

.. toctree::
   :hidden:

   Description
   Download
   Documentation
   Regression


.. raw:: html

   <table style="width: 100%; text-align: left; margin-left: auto; margin-right: auto;"
    border="0" cellpadding="10" cellspacing="10">
     <tbody>
       <tr valign="top">
         <td style="vertical-align: top; text-align: left;">
           <a href="http://science.energy.gov/ascr">
	     <img src="http://faculty.smu.edu/reynolds/arkode/Pics/doe_logo.jpg" alt="DOE ASCR logo" border="0" height="75"></a>
         </td>
         <td style="vertical-align: top; text-align: center;">
           <a href="http://fastmath-scidac.org/">
	     <img src="http://faculty.smu.edu/reynolds/arkode/Pics/fastmath_logo.png" alt="FASTMath logo" border="0" height="75"></a>
         </td>
         <td style="vertical-align: top; text-align: right;">
           <a href="https://www.llnl.gov/">
	     <img src="http://faculty.smu.edu/reynolds/arkode/Pics/llnl_logo.jpg" alt="LLNL logo" border="0" height="75"></a>
         </td>
         <td style="vertical-align: top; text-align: right;">
           <a href="https://www.smu.edu/">
	     <img src="http://faculty.smu.edu/reynolds/arkode/Pics/smu_logo.png" alt="SMU logo" border="0" height="75"></a>
         </td>
       </tr>
     </tbody>
   </table>

