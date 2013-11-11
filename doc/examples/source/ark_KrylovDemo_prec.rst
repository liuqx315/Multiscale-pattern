..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3



.. _ark_KrylovDemo_prec:

ark_KrylovDemo_prec
============================================

This problem is an ARKode clone of the CVODE problem,
``cv_KrylovDemo_prec``.  This is a demonstration program using the
GMRES linear solver.  The program solves a stiff ODE system that arises
from a system of PDEs modeling a six-species food web population
model, with predator-prey interaction and diffusion on the unit square
in two dimensions. We have a system with 6 components, :math:`C =
[c^1,\, c^2,\,\ldots, c^6]^T` that satisfy the equations, 

.. math::

   \frac{\partial c^i}{\partial t} &= d_i \left(\frac{\partial^2 c^i}{\partial
      x^2} + \frac{\partial^2 c^i}{\partial y^2}\right) +
      f_i(x,y,c),\quad i=1,\ldots,6.

where

.. math::

   f_i(x,y,c) = c^i\left( b_i + \sum_{j=1}^{ns} a_{i,j} c^j\right).

Here, the first three species are prey and the last three are
predators.  The coefficients :math:`a_{i,j}, b_i, d_i` are:

.. math::

   a_{i,j} = \begin{cases}
               -1, \quad & i=j,\\
	       -0.5\times10^{-6}, \quad & i\le 3, j>3, \\
	        10^4, \quad & i>3, j\le3
             \end{cases}
   b_i = \begin{cases}
            (1+xy), \quad & i\le 3,\\
	   -(1+xy), \quad & i>3
         \end{cases}
   d_i = \begin{cases}
            1, \quad & i\le 3,\\
	    \frac12, \quad & i>3
         \end{cases}

The spatial domain is :math:`(x,y) \in [0, 1]^2`; the time domain is
:math:`t \in [0,10]`, with initial conditions 

.. math::

   c^i(x,y) &=  10 + i \sqrt{4x(1-x)}\sqrt{4y(1-y)}

and with homogeneous Neumann boundary conditions, 
:math:`\nabla c^i \cdot \vec{n} = 0`.




Numerical method
----------------

We employ a method of lines approach, wherein we first
semi-discretize in space to convert the system of 6 PDEs into a larger
system of ODEs.  To this end, the spatial derivatives are computed
using second-order centered differences, with the data distributed
over :math:`Mx*My` points on a uniform spatial grid.  Resultingly, ARKode
approaches the problem as one involving :math:`6*Mx*My` coupled ODEs.

This program solves the problem with a DIRK method, using a Newton
iteration with the preconditioned ARKSPGMR iterative linear solver.
The preconditioner matrix used is the product of two matrices: 

1. A matrix, only defined implicitly, based on a fixed number of
   Gauss-Seidel iterations using the diffusion terms only. 

2. A block-diagonal matrix based on the partial derivatives of the
   interaction terms :math:`f` only, using block-grouping (computing
   only a subset of the :math:`3\times3` blocks). 

Four different runs are made for this problem.  The product
preconditoner is applied on the left and on the right.  In each case,
both the modified and classical Gram-Schmidt orthogonalization options
are tested.  In the series of runs, ``ARKodeInit`` and ``ARKSpgmr``
are called only for the first run, whereas ``ARKodeReInit``,
``ARKSpilsSetPrecType`` and ``ARKSpilsSetGSType`` are called for each
of the remaining three runs. 

A problem description, performance statistics at selected output
times, and final statistics are written to standard output.  On the
first run, solution values are also printed at output times.  Error
and warning messages are written to standard error, but there should
be no such messages. 



