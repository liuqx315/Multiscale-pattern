..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _LinearSolvers.SPILS:

The SPILS modules: SPGMR, SPBCG, SPTFQMR, SPFGMR and PCG
==========================================================

Due to their reliance on general vector operations, the iterative
linear solvers in the SPILS family can only be used in conjunction
with a complete NVECTOR implementation (i.e. one that provides a
nearly-full set of the generic NVECTOR module functions).  While both
the :ref:`NVECTOR_SERIAL <NVectors.NVSerial>` and
:ref:`NVECTOR_PARALLEL <NVectors.NVParallel>` modules provided with
SUNDIALS meet these criteria, these criteria may also be easily met
through a user-supplied vector implementation.



The SPGMR module
-----------------------------------------

The SPGMR package, in the files ``sundials_spgmr.h`` and
``sundials_spgmr.c``, includes an implementation of the scaled
preconditioned GMRES method. A separate code module, implemented in 
``sundials_iterative.h`` and ``sundials_iterative.c``, contains
auxiliary functions that support SPGMR, as well as the other Krylov
solvers in SUNDIALS (SPBCG, SPTFQMR and PCG). For full details, including
usage instructions, see the header files ``sundials_spgmr.h`` and
``sundials_iterative.h``. 

The files comprising the SPGMR generic linear solver, and their
locations in the SUNDIALS ``srcdir``, are as follows:

* header files (located in ``srcdir/include/sundials``)

  ``sundials_spgmr.h``, ``sundials_iterative.h``,
  ``sundials_nvector.h``, ``sundials_types.h``, ``sundials_math.h``,
  ``sundials_config.h``

* source files (located in ``srcdir/src/sundials``)

  ``sundials_spgmr.c``, ``sundials_iterative.c``, ``sundials_nvector.c``


Only two of the preprocessing directives in the header file
``sundials_config.h`` are required to use the SPGMR package by itself
(see the section :ref:`Installation` for details): 

* (required) definition of the precision of the SUNDIALS type
  ``realtype``. One of the following lines must be present:

  .. code-block:: c

     #define SUNDIALS_DOUBLE_PRECISION 1
     #define SUNDIALS_SINGLE_PRECISION 1
     #define SUNDIALS_EXTENDED_PRECISION 1

* (optional) use of generic math functions:

  .. code-block:: c

     #define SUNDIALS USE GENERIC MATH 1


The ``sundials_types.h`` header file defines the SUNDIALS ``realtype``
and ``booleantype`` types and the macro ``RCONST``, while the
``sundials_math.h`` header file is needed for the ``MAX`` and ``ABS``
macros and ``RAbs`` and ``RSqrt`` functions.

The generic NVECTOR files, ``sundials_nvector.h`` and
``sundials_nvector.c`` are needed for the definition of the generic
``N_Vector`` type and functions. The NVECTOR functions used by the
SPGMR module are: :c:func:`N_VDotProd()`, :c:func:`N_VLinearSum()`,
:c:func:`N_VScale()`, :c:func:`N_VProd()`, :c:func:`N_VDiv()`,
:c:func:`N_VConst()`, :c:func:`N_VClone()`,
:c:func:`N_VCloneVectorArray()`, :c:func:`N_VDestroy()`, and
:c:func:`N_VDestroyVectorArray()`. 

The nine files listed above can be extracted from the SUNDIALS
``srcdir`` and compiled by themselves into an SPGMR library or into a
larger user code. 

The following functions are available in the SPGMR package:

* ``SpgmrMalloc``: allocation of memory for ``SpgmrSolve``;
* ``SpgmrSolve``: solution of :math:`Ax = b` by the SPGMR method;
* ``SpgmrFree``: free memory allocated by ``SpgmrMalloc``.

The following functions are available in the support package
``sundials_iterative.h`` and ``sundials_iterative.c``:

* ``ModifiedGS``: performs modified Gram-Schmidt procedure;
* ``ClassicalGS``: performs classical Gram-Schmidt procedure;
* ``QRfact``: performs QR factorization of Hessenberg matrix;
* ``QRsol``: solves a least squares problem with a Hessenberg matrix
  factored by ``QRfact``. 




The SPBCG module
-----------------------------------------

The SPBCG package, in the files ``sundials_spbcgs.h`` and
``sundials_spbcgs.c``, includes an implementation of the scaled
preconditioned Bi-CGStab method. For full details, including usage
instructions, see the file ``sundials_spbcgs.h``.

The files needed to use the SPBCG module by itself are the same as for
the SPGMR module, but with ``sundials_spbcgs.h`` and
``sundials_spbcgs.c`` in place of ``sundials_spgmr.h`` and
``sundials_spgmr.c``. 

The following functions are available in the SPBCG package:

* ``SpbcgMalloc``: allocation of memory for ``SpbcgSolve``;
* ``SpbcgSolve``: solution of :math:`Ax = b` by the SPBCG method;
* ``SpbcgFree``: free memory allocated by ``SpbcgMalloc``.



The SPTFQMR module
-----------------------------------------


The SPTFQMR package, in the files ``sundials_sptfqmr.h`` and
``sundials_sptfqmr.c``, includes an implementation of the scaled
preconditioned TFQMR method. For full details, including usage
instructions, see the file ``sundials_sptfqmr.h``.

The files needed to use the SPTFQMR module by itself are the same as
for the SPGMR module, but with ``sundials_sptfqmr.h`` and
``sundials_sptfqmr.c`` in place of ``sundials_spgmr.h`` and
``sundials_spgmr.c``. 

The following functions are available in the SPTFQMR package:

* ``SptfqmrMalloc``: allocation of memory for ``SptfqmrSolve``;
* ``SptfqmrSolve``: solution of :math:`Ax = b` by the SPTFQMR method;
* ``SptfqmrFree``: free memory allocated by ``SptfqmrMalloc``.



The SPFGMR module
-----------------------------------------

The SPFGMR package, in the files ``sundials_spfgmr.h`` and
``sundials_spfgmr.c``, includes an implementation of the scaled
preconditioned Flexible Generalized Minimum Residual method. For full
details, including usage instructions, see the file
``sundials_spfgmr.h``. 

The files needed to use the SPFGMR module by itself are the same as for
the SPGMR module, but with ``sundials_spfgmr.h`` and
``sundials_spfgmr.c`` in place of ``sundials_spgmr.h`` and
``sundials_spgmr.c``. 

The following functions are available in the SPFGMR package:

* ``SpfgmrMalloc``: allocation of memory for ``SpfgmrSolve``;
* ``SpfgmrSolve``: solution of :math:`Ax = b` by the SPFGMR method;
* ``SpfgmrFree``: free memory allocated by ``SpfgmrMalloc``.



The PCG module
-----------------------------------------

The PCG package, in the files ``sundials_pcg.h`` and
``sundials_pcg.c``, includes an implementation of the 
preconditioned conjugate gradient method.  We note that due to the
requirement of symmetric linear systems for the conjugate gradient
method, this solver should only be used for problems with symmetric
linear operators.  Furthermore, aside from allowing a weight vector
for computing weighted convergence norms, no variable or equation
scaling is allowed for systems using this solver.  For full details,
including usage instructions, see the file ``sundials_pcg.h``.

The files needed to use the PCG module by itself are the same as for
the SPGMR module, but with ``sundials_pcg.h`` and
``sundials_pcs.c`` in place of ``sundials_spgmr.h`` and
``sundials_spgmr.c``. 

The following functions are available in the PCG package:

* ``PcgMalloc``: allocation of memory for ``PcgSolve``;
* ``PcgSolve``: solution of :math:`Ax = b` by the PCG method;
* ``PcgFree``: free memory allocated by ``PcgMalloc``.

