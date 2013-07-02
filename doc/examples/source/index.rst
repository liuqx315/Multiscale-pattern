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
Fortran.  These examples demonstrate a large variety of ARKode solver
options, including explicit, implicit and ImEx solvers,
root-finding, direct and iterative linear solvers, and the Fortran
solver interface, FARKODE.  While these examples are not an exhaustive
set of all possible usage scenarios, they are designed to show a
variety of usage scenarios, and can be used as templates for new
problems using ARKode's solvers.

The following table summarizes the salient features of each of the
following example problems.  Each example is designed to be relatively
self-contained, so that you need only study and/or emulate the problem
that is most closely related to your own.

.. cssclass:: table-bordered

================================================  ==========  =============  =============  ========  ===============================================================
Problem                                           Integrator  Linear Solver  Size           Language  Extras
================================================  ==========  =============  =============  ========  ===============================================================
:ref:`ark_analytic <ark_analytic>`                DIRK        Dense          1              C         Analytical solution, variable stiffness
:ref:`ark_analytic_nonlin <ark_analytic_nonlin>`  ERK         N.A.           1              C         Nonlinear, analytical solution
:ref:`ark_analytic_sys <ark_analytic_sys>`        DIRK        Dense          3              C++       ODE system, analytical solution, variable stiffness
:ref:`ark_brusselator <ark_brusselator>`          DIRK        Dense          3              C         Stiff, nonlinear, ODE system, "standard" test problem
:ref:`ark_bruss <ark_bruss>`                      ARK         Dense          3              F90       Stiff, nonlinear, ODE system, "standard" test problem
:ref:`ark_robertson <ark_robertson>`              DIRK        Dense          3              C         Stiff, nonlinear, ODE system, "standard" test problem
:ref:`ark_robertson_root <ark_robertson_root>`    DIRK        Dense          3              C         Utilizes root-finding capabilities
:ref:`ark_brusselator1D <ark_brusselator1D>`      DIRK        Band           3N             C         Stiff, nonlinear, reaction-diffusion PDE system
:ref:`ark_heat1D <ark_heat1D>`                    DIRK        PCG            N              C         Stiff, linear, diffusion PDE, iterative linear solver
:ref:`ark_heat2D <ark_heat2D>`                    DIRK        PCG            :math:`nx*ny`  C++       Parallel, stiff, linear, diffusion PDE, iterative linear solver
================================================  ==========  =============  =============  ========  ===============================================================


Further details on each of the above-listed examples, including both
source code and plots of the computed results, are provided in the
following sub-sections:

.. toctree::
   :maxdepth: 1

   Simple linear example (ark_analytic) <analytic>
   Simple nonlinear example (ark_analytic_nonlin) <analytic_nonlin>
   Simple linear system example (ark_analytic_sys) <analytic_sys>
   Stiff nonlinear system example (ark_brusselator) <brusselator>
   Stiff nonlinear system, Fortran example (ark_bruss) <bruss>
   Stiff nonlinear system example (ark_robertson) <robertson>
   Stiff nonlinear system with root-finding example (ark_robertson_root) <robertson_root>
   Stiff PDE system example (ark_brusselator1D) <brusselator1D>
   PDE example with iterative linear solver (ark_heat1D) <heat1D>
   Parallel PDE example with iterative linear solver (ark_heat2D) <heat2D>
   
.. only:: html

   * :ref:`search`


