..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3



.. _fark_diag_kry_bbd_p:

fark_diag_kry_bbd_p
===================================================

This problem is an ARKode clone of the CVODE problem,
``fcv_diag_kry_bbd_p``.  This test problem models a stiff, linear,
diagonal ODE system,

.. math::

   \frac{\partial y_i}{\partial t} &= -\alpha i y_i, \quad i=1,\ldots N.


Here :math:`\alpha=10` and :math:`N=10 N_P`, where :math:`N_P` is the
number of MPI tasks used for the problem.  The problem has initial
conditions :math:`y_i=1` and evolves for the time interval :math:`t\in [0,1]`.




Numerical method
----------------

This program solves the problem with a DIRK method, using a Newton
iteration with the preconditioned ARKSPGMR iterative linear solver.

A diagonal preconditioner matrix is used, formed automatically through
difference quotients within the ARKBBDPRE module.  Since ARKBBDPRE is
developed for use of a block-banded preconditioner, in this solver
each block is set to have half-bandwidths ``mudq = mldq = 0`` to
retain only the diagonal portion.

Two runs are made for this problem, first with left and then with
right preconditioning (``IPRE`` is first set to 1 and then to 2).

Performance data is printed at selected output times, and maximum
errors and final performance counters are printed on completion.
