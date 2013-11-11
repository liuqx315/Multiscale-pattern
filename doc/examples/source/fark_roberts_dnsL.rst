..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3



.. _fark_roberts_dnsL:

fark_roberts_dnsL
===================================================

This problem is an ARKode clone of the CVODE problem,
``fcv_roberts_dnsL``.  This test problem models the kinetics of a
three-species autocatalytic reaction.  This is an ODE system with 3
components, :math:`Y = [y_1,\, y_2,\, y_3]^T`, satisfying the equations,

.. math::

   \frac{d y_1}{dt} &= -0.04 y_1 + 10^4 y_2 y_3, \\
   \frac{d y_2}{dt} &= 0.04 y_1 - 10^4 y_2 y_3 - 3\cdot10^7 y_2^2, \\
   \frac{d y_3}{dt} &= 3\cdot10^7 y_2^2.

We integrate over the interval :math:`0\le t\le 4\cdot10^{10}`, with initial
conditions  :math:`Y(0) = [1,\, 0,\, 0]^T`. 

Additionally, we supply the following two root-finding equations:

.. math::

   g_1(u) = u - 10^{-4}, \\
   g_2(w) = w - 10^{-2}.

While these are not inherently difficult nonlinear equations, they
easily serve the purpose of determining the times at which our
solutions attain desired target values.


Numerical method
----------------

This program solves the problem with a DIRK method, using a Newton
iteration with the dense LAPACK linear solver module.

As with the :ref:`ark_robertson_root` problem, we enable ARKode's
rootfinding module to find the times at which
either :math:`u=10^{-4}` or :math:`w=10^{-2}`.

Performance data and solution values are printed at
selected output times, along with additional output at rootfinding
events.  All performance counters are printed on completion.

