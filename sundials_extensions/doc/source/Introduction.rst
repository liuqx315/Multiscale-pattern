Introduction
============

This is the documentation for ARKode, an adaptive step time
integration package for stiff, nonstiff and multi-rate systems of
ordinary differential equations (ODEs) given in explicit form

.. math::
   M y' = f_E(t,y) + f_I(t,y)
   :label: ODE

where :math:`t` is the independent variable, :math:`y` is the set of
dependent variables (in :math:`\Re^n`), :math:`M` is a
user-specified, non-singular operator from :math:`\Re^n` to
:math:`\Re^n` (possibly time dependent, but independent of
:math:`y`), and the right-hand side function is partitioned into two
components: 

- :math:`f_E(t,y)` contains the "slow" time scale components to be
  integrated explicitly, and 
- :math:`f_I(t,y)`  contains the "fast" time scale components to be
  integrated implicitly. 

Either of these operators may be disabled, allowing for fully
explicit, fully implicit, or combination implicit-explicit (IMEX) time
integration. 

The methods used in ARKode are adaptive-step additive Runge Kutta
methods. Such methods are defined with two complementary Runge-Kutta
methods: one explicit (ERK) and the other diagonally implicit
(DIRK). Through appropriately partitioning the ODE system :eq:`ODE`, such
methods enable accurate and efficient time integration of multi-rate
systems of ordinary differential equations, wherein only the
components in :math:`f_I(t,y)` must be solved implicitly. A variety of
built-in Butcher tables are included, allowing for adaptive explicit
methods of order 2-6, adaptive implicit methods of orders 3-5, and
adaptive IMEX methods of orders 3-5. 

For implicit and IMEX methods, the resulting nonlinear system is
solved approximately at each integration step, using modified or
inexact Newton methods. When used with the serial NVECTOR module in
SUNDIALS, ARKode provides both direct (dense and band) and
preconditioned Krylov iterative (GMRES, BiCGStab, TFQMR) linear
solvers. When used with the parallel NVECTOR module or a user-provided
vector data structure, only the Krylov solvers are available. 
