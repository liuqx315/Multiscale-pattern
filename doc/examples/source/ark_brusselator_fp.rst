..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3



.. _ark_brusselator_fp:

ark_brusselator_fp
===================================================

This test problem is a duplicate the ``ark_brusselator`` problem
above, but with a few key changes in the methods used for time
integration and nonlinear solver.  As with the previous test, this
problem has 3 dependent variables :math:`u`, :math:`v` and :math:`w`,
that depend on the independent variable :math:`t` via the IVP system

.. math::

   \frac{du}{dt} &= a - (w+1)u + v u^2, \\
   \frac{dv}{dt} &= w u - v u^2, \\
   \frac{dw}{dt} &= \frac{b-w}{\varepsilon} - w u.

We integrate over the interval :math:`0 \le t \le 10`, with the
initial conditions :math:`u(0) = u_0`, :math:`v(0) = v_0`, :math:`w(0) = w_0`.
After each unit time interval, the solution is output to the screen.

We have 3 different testing scenarios:

Test 1:  :math:`u_0=3.9`,  :math:`v_0=1.1`,  :math:`w_0=2.8`,
:math:`a=1.2`, :math:`b=2.5`, and :math:`\varepsilon=10^{-5}` 

Test 2:  :math:`u_0=1.2`, :math:`v_0=3.1`, :math:`w_0=3`, :math:`a=1`,
:math:`b=3.5`, and :math:`\varepsilon=5\cdot10^{-6}` 

Test 3:  :math:`u_0=3`, :math:`v_0=3`, :math:`w_0=3.5`, :math:`a=0.5`,
:math:`b=3`, and :math:`\varepsilon=5\cdot10^{-4}` 

These tests are selected within the input file (test = {1,2,3}), 
with the default set to test 2 in case the input is invalid.
Also in the input file, we allow specification of the desired 
relative and absolute tolerances.



Numerical method
----------------

This program solves the problem with the ARK method, in which we have
split the right-hand side into stiff (:math:`f_i(t,y)`) and non-stiff
(:math:`f_e(t,y)`) components,

.. math::

   f_i(t,y) = \left[\begin{array}{c} 
      0 \\ 0 \\ \frac{b-w}{\varepsilon} 
   \end{array}\right]
   \qquad
   f_e(t,y) = \left[\begin{array}{c} 
      a - (w+1)u + v u^2 \\ w u - v u^2 \\ - w u
   \end{array}\right].

Also unlike the previous test problem, we solve the resulting implicit
stages using the available accelerated fixed-point solver, enabled
through a call to ``ARKodeSetFixedPoint``, with an acceleration
subspace of dimension 3.

100 outputs are printed at equal intervals, and run statistics 
are printed at the end.


