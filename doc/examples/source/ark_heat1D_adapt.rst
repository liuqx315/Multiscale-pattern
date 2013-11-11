..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3



.. _ark_heat1D_adapt:

ark_heat1D_adapt
===================================================

This problem is a clone of the ``ark_heat1D`` test problem except that
unlike the previous uniform-grid problem, this test problem allows a
dynamically-evolving spatial mesh.  The PDE under consideration is a
simple one-dimensional heat equation, 

.. math::

   \frac{\partial u}{\partial t} = k \frac{\partial^2 u}{\partial x^2} + f,

for :math:`t \in [0, 10]`, and :math:`x \in [0, 1]`, with initial
condition :math:`u(0,x) = 0`, stationary boundary conditions,

.. math::

   \frac{\partial u}{\partial t}(t,0) = \frac{\partial u}{\partial t}(t,1) = 0,

and a point-source heating term, 

.. math::

   f(t,x) = \begin{cases} 1 & \text{if}\;\; x=1/2, \\
                          0 & \text{otherwise}. \end{cases}

 

Numerical method
----------------

We again employ a method-of-lines discretization approach.  The
spatial derivatives are computed using a three-point centered stencil,
that is accurate to :math:`O(\Delta x_i^2)` if the neighboring points are
equidistant from the central point, i.e. :math:`x_{i+1} - x_i = x_i -
x_{i-1}`, though if these are unequal the approximation reduces to
first-order accuracy.  The spatial mesh is initially distributed
uniformly over 21 points in :math:`[0,1]`, but as the simulation
proceeds the mesh is [crudely] adapted to add points to the center of
subintervals bordering any node where 
:math:`\left|\frac{\partial^2 u}{\partial x^2}\right| > 3\times10^{-3}`.  

This program solves the problem with a DIRK method, utilizing a Newton
iteration and the PCG iterative linear solver.  Additionally, the test
problem utilizes ARKode's spatial adaptivity support (via
``ARKodeResize``), allowing retention of the major ARKode data
structures across vector length changes.
