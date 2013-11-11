..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3



.. _fark_diag_non_p:

fark_diag_non_p
===================================================

This problem is an ARKode clone of the CVODE problem,
``fcv_diag_non_p``.  This test problem models a nonstiff, linear,
diagonal ODE system,

.. math::

   \frac{\partial y_i}{\partial t} &= -\alpha i y_i, \quad i=1,\ldots N.


Here :math:`\alpha=\frac{10}{N}` and :math:`N=10 N_P`, where :math:`N_P` is the
number of MPI tasks used for the problem.  The problem has initial
conditions :math:`y_i=1` and evolves for the time interval :math:`t\in [0,1]`.




Numerical method
----------------

This program solves the problem with an ERK method, and hence does not
require either a nonlinear or linear solver for integration.

Performance data is printed at selected output times, and maximum
errors and final performance counters are printed on completion.
