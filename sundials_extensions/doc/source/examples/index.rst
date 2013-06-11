:tocdepth: 3

.. _Examples:

=================
 ARKode Examples
=================

ARKode comes packaged with a variety of example problems, that
exercise options including explicit, implicit and ImEx solvers,
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

============================================  ==========  =============  ===========  ========  ===============================================================
Problem                                       Integrator  Linear Solver  Size         Language  Extras
============================================  ==========  =============  ===========  ========  ===============================================================
:ref:`analytic <ark_analytic>`                DIRK        Dense          1            C         Analytical solution, variable stiffness
:ref:`analytic_nonlin <ark_analytic_nonlin>`  ERK         N.A.           1            C         Nonlinear, analytical solution
:ref:`analytic_sys <ark_analytic_sys>`        DIRK        Dense          3            C++       ODE system, analytical solution, variable stiffness
:ref:`brusselator <ark_brusselator>`          DIRK        Dense          3            C         Stiff, nonlinear, ODE system, "standard" test problem
:ref:`bruss <ark_bruss>`                      ARK         Dense          3            F90       Stiff, nonlinear, ODE system, "standard" test problem
:ref:`robertson <ark_robertson>`              DIRK        Dense          3            C         Stiff, nonlinear, ODE system, "standard" test problem
:ref:`robertson_root <ark_robertson_root>`    DIRK        Dense          3            C         Utilizes root-finding capabilities
:ref:`brusselator1D <ark_brusselator1D>`      DIRK        Band           3N           C         Stiff, nonlinear, reaction-diffusion PDE system
:ref:`heat1D <ark_heat1D>`                    DIRK        PCG            N            C         Stiff, linear, diffusion PDE, iterative linear solver
:ref:`heat2D <ark_heat2D>`                    DIRK        PCG            :math:`N^2`  C++       Parallel, stiff, linear, diffusion PDE, iterative linear solver
============================================  ==========  =============  ===========  ========  ===============================================================



.. toctree::
   :maxdepth: 1

   analytic
   analytic_nonlin
   analytic_sys
   brusselator
   bruss
   robertson
   robertson_root
   brusselator1D
   heat1D
   heat2D
   
