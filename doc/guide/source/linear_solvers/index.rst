..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _LinearSolvers:

========================
Linear Solvers in ARKode
========================

In this section, we describe the generic linear solver code modules
from SUNDIALS that are included in ARKode.  While these may be used in
conjunction with ARKode, they may also be used separately as generic
packages in themselves.  These generic linear solver modules are
organized in three families of solvers, the DLS family, which includes
direct linear solvers appropriate for sequential computations; the SLS
family, which includes direct linear solvers for sparse matrices in
serial or shared-memory computations; and the 
SPILS family, which includes scaled preconditioned iterative (Krylov)
linear solvers. The solvers in each family share common data
structures and functions. 


The :ref:`DLS <LinearSolvers.DLS>` family contains the following two
generic linear solvers: 

* The DENSE package, a linear solver for dense matrices either
  specified through a matrix type (defined below) or as simple
  arrays. 
* The BAND package, a linear solver for banded matrices either
  specified through a matrix type (defined below) or as simple
  arrays. 

We further note that this family also includes the BLAS/LAPACK linear
solvers (dense and band) available to the SUNDIALS solvers, but these
are not discussed here. 


The :ref:`SLS <LinearSolvers.SLS>` family contains the following two
sparse direct linear solvers: 

* The KLU package, a linear solver for compressed-sparse-column
  matrices, [KLU]_.
* The SUPERLUMT package, a threaded linear solver for
  compressed-sparse-column matrices, [SuperLUMT]_.


The :ref:`SPILS <LinearSolvers.SPILS>` family contains the following
generic linear solvers: 

* The SPGMR package, a solver for the scaled preconditioned GMRES
  method. 
* The SPBCG package, a solver for the scaled preconditioned Bi-CGStab
  method. 
* The SPTFQMR package, a solver for the scaled preconditioned TFQMR
  method. 
* The SPFGMR package, a solver for the scaled preconditioned Flexible
  GMRES method. 
* The PCG package, a solver for the preconditioned conjugate gradient
  method. 

For reasons related to installation, the names of the files involved
in these generic solvers begin with the prefix SUNDIALS. But despite
this, each of the DLS and SPILS solvers is in fact generic, in that
they are usable completely independently of SUNDIALS. 

For the sake of space, the functions for the DENSE, BAND modules
that work with a matrix type and the functions in the SPGMR, SPBCG,
SPTFQMR, SPFGMR and PCG modules are only summarized briefly, since
they are less likely to be of direct use in connection with a SUNDIALS
solver.  However, the functions for dense matrices treated as simple
arrays and sparse matrices are fully described, because we anticipate
that they will be useful in the implementation of preconditioners used
with the combination of one of the SUNDIALS solvers and one of the
SPILS linear solvers. 

Lastly, it is possible to supply customized linear solvers to ARKode,
in that the ARKode solvers only require the existence of a minimal set
of generic routines.  Through attaching user-supplied routines for
these function pointers, it is possible to use arbitrary approaches
for solution to the implicit linear systems arising during an ARKode
solve. 

Specifics of these built-in linear solver packages, as well as the
generic linear solver interface, are provided in the following
sub-sections:

.. toctree::
   :maxdepth: 1

   DLS
   SLS
   SPILS
   custom
