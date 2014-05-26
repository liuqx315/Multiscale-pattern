..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _Introduction:

Introduction
============

The ARKode solver library provides an adaptive-step time integration
package for stiff, nonstiff and multi-rate systems of ordinary
differential equations (ODEs) given in explicit form

.. math::
   M \dot{y} = f_E(t,y) + f_I(t,y),  \qquad y(t_0) = y_0,
   :label: ODE

where :math:`t` is the independent variable, :math:`y` is the set of
dependent variables (in :math:`\mathbb{R}^N`), :math:`M` is a
user-specified, nonsingular operator from :math:`\mathbb{R}^N` to
:math:`\mathbb{R}^N` (possibly time dependent, but independent of
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
methods. Such methods are defined through combining two complementary
Runge-Kutta methods: one explicit (ERK) and the other diagonally implicit
(DIRK). Through appropriately partitioning the ODE system into
explicit and implicit components :eq:`ODE`, such methods have the
potential to enable accurate and efficient time integration of
multi-rate systems of ordinary differential equations.  A key
feature allowing for high efficiency of these methods is that only the
components in :math:`f_I(t,y)` must be solved implicitly, allowing for
splittings tuned for use with optimal implicit solvers.  

This framework allows for significant freedom over the constitutive
methods used for each component, and ARKode is packaged with a wide
array of built-in methods for use.  These built-in Butcher tables
include adaptive explicit methods of orders 2-6, adaptive implicit
methods of orders 2-5, and adaptive IMEX methods of orders 3-5. 

For problems that include nonzero implicit term :math:`f_I(t,y)`, the
resulting implicit system (assumed nonlinear) is solved approximately
at each integration step, using a modified Newton method, an Inexact
Newton method, or an accelerated fixed-point solver.  For implicit
problems using a Newton-based solver and the serial NVECTOR module in
SUNDIALS, ARKode provides both direct (dense, band and sparse) and
preconditioned Krylov iterative (GMRES, BiCGStab, TFQMR, FGMRES, PCG)
linear solvers.  When used with the parallel NVECTOR module or a
user-provided vector data structure, only the Krylov solvers are
available, although a user may supply their own linear solver for any
data structures if desired.

The guide is separated into sections focused on the major aspects of
the ARKode library.  In the next section we provide a thorough
presentation of the underlying :ref:`mathematics <Mathematics>` that
relate these algorithms together.  We follow this with overview of how 
the source code for ARKode is :ref:`organized <Organization>`.  The
largest section follows, providing a full account of the ARKode user
interface, including a description of all user-accessible functions
and outlines for ARKode usage for serial and parallel applications.
Since ARKode is written in C, we first present :ref:`the C and C++
interface <CInterface>`, followed with a separate section on
:ref:`using ARKode within Fortran applications <FortranInterface>`.  The
following three sections discuss shared features between ARKode and
the rest of the SUNDIALS library: :ref:`vector data structures <NVectors>`,
:ref:`linear solvers <LinearSolvers>`, and the :ref:`installation
procedure <Installation>`.  The final sections catalog the full set of
:ref:`ARKode constants <Constants>`, that are used for both input
specifications and return codes, and the full set of 
:ref:`Butcher tables <Butcher>` that are packaged with
ARKode. 
