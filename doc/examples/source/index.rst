..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

.. ARKode_example documentation master file, created by
   sphinx-quickstart on Sat Dec 22 20:38:03 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==============================================
ARKode Example documentation
==============================================

This is the documentation for the ARKode examples.  ARKode is an
adaptive step time integration package for stiff, nonstiff and
multi-rate systems of ordinary differential equations (ODEs).  
The ARKode solver is a component of the `SUNDIALS
<https://computation.llnl.gov/casc/sundials/main.html>`_ suite of
nonlinear and differential/algebraic equation solvers. It is designed
to have a similar user experience to the `CVODE
<https://computation.llnl.gov/casc/sundials/description/description.html#descr_cvode>`_
solver, with user modes to allow adaptive integration to specified
output times, return after each internal step and root-finding
capabilities, for calculations both in serial and parallel (via
MPI). The default integration and solver options should apply to most
users, though complete control over all internal parameters and time
adaptivity algorithms is enabled through optional interface routines.  

ARKode is developed by `Southern Methodist University
<http://www.smu.edu>`_, with support by the `US Department of Energy
<http://www.doe.gov>`_ through the `FASTMath
<http://www.fastmath-scidac.org/>`_ SciDAC Institute, under subcontract
B598130 from `Lawrence Livermore National Laboratory
<http://www.llnl.gov>`_. 

Along with the ARKode solver, we have created a suite of example
problems demonstrating its usage on applications written in C, C++ and
Fortran 77 and Fortran 90.  These examples demonstrate a large variety
of ARKode solver options, including explicit, implicit and ImEx
solvers, root-finding, Newton and fixed-point nonlinear solvers,
direct and iterative linear solvers, adaptive resize capabilities, and
the Fortran solver interface.  While these examples are not an
exhaustive set of all possible usage scenarios, they are designed to
show a variety of exemplars, and can be used as templates for new
problems using ARKode's solvers.

The following table summarizes the salient features of each of the
following example problems.  Each example is designed to be relatively
self-contained, so that you need only study and/or emulate the problem
that is most closely related to your own.

.. cssclass:: table-bordered

====================================================  ==========  ================  =============  =============  ========  ===============================================================
Problem                                               Integrator  Nonlinear Solver  Linear Solver  Size           Language  Extras
====================================================  ==========  ================  =============  =============  ========  ===============================================================
:ref:`ark_analytic <ark_analytic>`                    DIRK        Newton            Dense          1              C         Analytical solution, variable stiffness
:ref:`ark_analytic_nonlin <ark_analytic_nonlin>`      ERK         N.A.              N.A.           1              C         Nonlinear, analytical solution
:ref:`ark_analytic_sys <ark_analytic_sys>`            DIRK        Newton            Dense          3              C++       ODE system, analytical solution, variable stiffness
:ref:`ark_brusselator <ark_brusselator>`              DIRK        Newton            Dense          3              C         Stiff, nonlinear, ODE system, "standard" test problem
:ref:`ark_bruss <ark_bruss>`                          ARK         Newton            Dense          3              F90       Stiff, nonlinear, ODE system, "standard" test problem
:ref:`ark_robertson <ark_robertson>`                  DIRK        Newton            Dense          3              C         Stiff, nonlinear, ODE system, "standard" test problem
:ref:`ark_robertson_root <ark_robertson_root>`        DIRK        Newton            Dense          3              C         Utilizes rootfinding capabilities
:ref:`ark_brusselator1D <ark_brusselator1D>`          DIRK        Newton            Band           3N             C         Stiff, nonlinear, reaction-diffusion PDE system
:ref:`ark_heat1D <ark_heat1D>`                        DIRK        Newton            PCG            N              C         Stiff, linear, diffusion PDE, iterative linear solver
:ref:`ark_heat2D <ark_heat2D>`                        DIRK        Newton            PCG            :math:`nx*ny`  C++       Parallel, stiff, linear, diffusion PDE, iterative linear solver
:ref:`ark_KrylovDemo_prec <ark_KrylovDemo_prec>`      DIRK        Newton            SPGMR          216            C         Stiff, nonlinear, rx-diff PDE system, different preconditioners
:ref:`ark_brusselator_fp <ark_brusselator_fp>`        ARK         Fixed-point       N.A.           3              C         Stiff, nonlinear, ODE system
:ref:`ark_diurnal_kry_bbd_p <ark_diurnal_kry_bbd_p>`  DIRK        Newton            SPGMR          200            C         Stiff, nonlinear, PDE system, parallel, BBD preconditioner
:ref:`ark_diurnal_kry_p <ark_diurnal_kry_p>`          DIRK        Newton            SPGMR          200            C         Stiff, nonlinear, PDE system, parallel, block-diagonal precond.
:ref:`ark_heat1D_adapt <ark_heat1D_adapt>`            DIRK        Newton            PCG            (dynamic)      C         Stiff, linear, diffusion, PCG solver, adaptive vector resizing
:ref:`fark_diag_kry_bbd_p <fark_diag_kry_bbd_p>`      DIRK        Newton            SPGMR          10*NProcs      F77       Stiff, linear, diagonal ODE system, BBD preconditioner
:ref:`fark_diag_non_p <fark_diag_non_p>`              ERK         N.A.              N.A.           10*NProcs      F77       Nonstiff, linear, diagonal ODE system
:ref:`fark_diurnal_kry_bp <fark_diurnal_kry_bp>`      DIRK        Newton            SPGMR          10             F77       Stiff, nonlinear, PDE system, banded preconditioner
:ref:`fark_heat2D <fark_heat2D>`                      DIRK        Newton            PCG            :math:`nx*ny`  F90       Parallel, stiff, linear, diffusion PDE, iterative linear solver
:ref:`fark_roberts_dnsL <fark_roberts_dnsL>`          DIRK        Newton            Dense          3              F77       Stiff, nonlinear, ODE system, LAPACK dense solver, rootfinding
====================================================  ==========  ================  =============  =============  ========  ===============================================================


Further details on each of the above-listed examples, including both
source code and plots of the computed results, are provided in the
following sub-sections:

.. toctree::
   :maxdepth: 1

   Simple linear example (ark_analytic) <ark_analytic>

   Simple nonlinear example (ark_analytic_nonlin) <ark_analytic_nonlin>

   Simple linear system example (ark_analytic_sys) <ark_analytic_sys>

   Stiff nonlinear system example (ark_brusselator) <ark_brusselator>

   Stiff nonlinear system, Fortran 90 example (ark_bruss) <ark_bruss>

   Stiff nonlinear system example (ark_robertson) <ark_robertson>

   Stiff nonlinear system with rootfinding example (ark_robertson_root) <ark_robertson_root>

   Stiff PDE system example (ark_brusselator1D) <ark_brusselator1D>

   PDE example with iterative linear solver (ark_heat1D) <ark_heat1D>

   Parallel PDE example with iterative linear solver (ark_heat2D) <ark_heat2D>

   Stiff nonlinear system, uses GMRES with different preconditioners (ark_KrylovDemo_prec) <ark_KrylovDemo_prec>

   Stiff nonlinear system, ImEx example, using the accelerated fixed-point nonlinear solver (ark_brusselator) <ark_brusselator_fp>

   Stiff nonlinear parallel C example with BBD preconditioner (ark_diurnal_kry_bbd_p) <ark_diurnal_kry_bbd_p>

   Stiff nonlinear parallel C example with block-diagonal preconditioner (ark_diurnal_kry_p) <ark_diurnal_kry_p>

   PDE example with PCG linear solver, and dynamic vector "resize" capabilities (ark_heat1D_adapt) <ark_heat1D_adapt>

   Nonstiff parallel diagonal ODE, Fortran 77 example using the BBD preconditioner (fark_diag_kry_bbd_p) <fark_diag_kry_bbd_p>

   Stiff parallel diagonal ODE, Fortran 77  example (fark_diag_non_p) <fark_diag_non_p>

   Stiff kinetics-transport Fortran 77 example using the banded preconditioner (fark_diurnal_kry_bp) <fark_diurnal_kry_bp>

   Parallel Fortran 90 example that replicates the ark_heat2D example (fark_heat2D) <fark_heat2D>

   Stiff chemical kinetics Fortran 77 example using Lapack solver interface and rootfinding module (fark_roberts_dnsL) <fark_roberts_dnsL>

   
.. only:: html

   * :ref:`search`


