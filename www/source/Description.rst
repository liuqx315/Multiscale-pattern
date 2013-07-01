Description
============

ARKode is designed to have a similar user experience to the `CVODE
<https://computation.llnl.gov/casc/sundials/description/description.html#descr_cvode>`_
solver that is already included in `SUNDIALS
<https://computation.llnl.gov/casc/sundials/main.html>`_, but it is
designed for multi-rate systems of ordinary differential equations.
As such, it admits user modes to allow adaptive integration to specified
output times, return after each internal step and root-finding
capabilities, for calculations both in serial and parallel (via
MPI). The default integration and solver options should apply to most
users, though complete control over all internal parameters and time
adaptivity algorithms is enabled through optional interface routines.

The methods used in ARKode are adaptive-step additive Runge Kutta
methods. Such methods are defined using a combination of two
complementary Runge-Kutta methods: one explicit (ERK) and the other
diagonally implicit (DIRK). Through appropriately  partitioning the
ODE system, such methods enable accurate and efficient time
integration of multi-rate systems of ordinary differential equations,
wherein only the components in :math:`f_I(t,y)` must be solved
implicitly (but the components in :math:`f_E(t,y)` are treated
explicitly). A variety of built-in Butcher tables are included,
allowing for adaptive explicit methods of order 2-6, adaptive implicit
methods of orders 3-5, and adaptive IMEX methods of orders 3-5.

For implicit and IMEX methods, the resulting nonlinear system is
solved approximately at each integration step, using modified or
inexact Newton methods.  When used with the serial NVECTOR module in
SUNDIALS, ARKode provides both direct (dense and band) and
preconditioned Krylov iterative (GMRES, BiCGStab, TFQMR) linear
solvers.  When used with the parallel NVECTOR module or an appropriate 
user-provided vector data structure, only the Krylov solvers are
available.

