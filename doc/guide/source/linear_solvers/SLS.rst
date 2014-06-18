..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2014, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _LinearSolvers.SLS:

The SLS modules: KLU and SUPERLUMT
========================================

The files comprising the SPARSE linear solver package, and their
locations in the SUNDIALS ``srcdir``, are as follows:

* header files (located in ``srcdir/include/sundials``):

  ``sundials_sparse.h``, ``sundials_klu_impl.h``,
  ``sundials_superlumt_impl.h``, ``sundials_types.h``,
  ``sundials_math.h``, ``sundials_config.h``

* source files (located in ``srcdir/src/sundials``):

  ``sundials_sparse.c``, ``sundials_math.c``

Only two of the preprocessing directives in the header file
``sundials_config.h`` are relevant to the SPARSE package by
itself (see the section :ref:`Installation` for details): 

* (required) definition of the precision of the SUNDIALS type
  ``realtype``. One of the following lines must be present:

  .. code-block:: c
 
     #define SUNDIALS_DOUBLE_PRECISION 1
     #define SUNDIALS_SINGLE_PRECISION 1
     #define SUNDIALS_EXTENDED_PRECISION 1

* (optional) use of generic math functions: 

  .. code-block:: c

     #define SUNDIALS_USE_GENERIC_MATH 1

The ``sundials_types.h`` header file defines the SUNDIALS ``realtype``
and ``booleantype`` types and the macro ``RCONST``, while the
``sundials_math.h`` header file is needed for the ``MIN``, ``MAX``,
and ``ABS`` macros and ``RAbs`` and ``RSqrt`` functions.

The files listed above for either module can be extracted from the
SUNDIALS ``srcdir`` and compiled by themselves into a separate library
or into a larger user code.



SlsMat
--------------------

The type :c:type:`SlsMat`, defined in ``sundials_sparse.h`` is a
pointer to a structure defining a generic matrix, and is used with all
linear solvers in the SLS family: 

.. c:type:: SlsMat

   .. code-block:: c

      typedef struct _SlsMat {
        int M;
        int N;
        int NNZ;
        realtype *data;
        int *rowvals;
        int *colptrs;
      } *SlsMat;

The fields of this structure are as follows (note that a dense matrix
of type :c:type:`SlsMat` need not be square): 

:M: -- number of rows
:N: --  number of columns
:NNZ: -- maximum number of nonzero entries in the matrix (allocated
   length of **data** and **rowvals** arrays)
:data: -- pointer to a contiguous block of ``realtype`` variables (of
   length **NNZ**), containing the values of the nonzero entries in the
   matrix.
:rowvals: -- pointer to a contiguous block of ``int`` variables (of
   length **NNZ**), containing the row indices of each nonzero
   entry held in **data**.
:colptrs: -- pointer to a contiguous block of ``int`` variables (of
  length **N+1**).  Each entry provides the index of the first column
  entry into the **data** and **rowvals** arrays, e.g. if
  **colptr[3]=7**, then the first nonzero entry in the fourth column
  of the matrix is located in **data[7]**, and is located in row
  **rowvals[7]** of the matrix.  The last entry points just past the
  end of the active data in **data** and **rowvals**.

..
   .. _SLS_figure:

   .. figure:: figs/sls_diagram.png

      SLS Diagram: caption




Functions in the SPARSE module
-------------------------------------------

The SPARSE module defines functions that act on sparse matrices of
type :c:type:`SlsMat`.  For full details, see the header file
``sundials_sparse.h``.


.. c:function:: SlsMat NewSparseMat(int M, int N, int NNZ)
   
   Allocates a :c:type:`SlsMat` sparse matrix having *M* rows, *N*
   columns, and storage for *NNZ* nonzero entries.

.. c:function:: SlsMat SlsConvertDls(DlsMat A)

   Converts a dense matrix of type :c:type:`DlsMat` into a sparse
   matrix of type :c:type:`SlsMat` by retaining only the nonzero
   values of the dense matrix.

.. c:function:: void DestroySparseMat(SlsMat A)

   Frees memory for a :c:type:`SlsMat` matrix.

.. c:function:: void SlsSetToZero(SlsMat A)

   Zeros out a :c:type:`SlsMat` matrix (but retains its storage).

.. c:function:: void CopySparseMat(SlsMat A, SlsMat B)

   Copies one sparse matrix to another.  If *B* has insufficient
   storage, its data arrays are reallocated to match those from *A*.

.. c:function:: void ScaleSparseMat(realtype c, SlsMat A)

   Scales a sparse matrix by a scalar.

.. c:function:: void AddIdentitySparseMat(SlsMat A)

   Increments a sparse matrix by the identity matrix.  If *A* is not
   square, only the existing diagonal values are incremented.  Resizes
   the data arrays of *A* upon completion to exactly match the
   nonzero storage for the result.

.. c:function:: int SlsAddMat(SlsMat A, SlsMat B)

   Adds two sparse matrices: :math:`A = A+B`.  Resizes the data arrays
   of *A* upon completion to exactly match the nonzero storage for
   the result.  Upon successful completion, the return value is zero;
   otherwise 1 is returned.

.. c:function:: void ReallocSparseMat(SlsMat A)

   This function eliminates unused storage in *A* by reallocating
   the internal ``data`` and ``rowvals`` arrays to contain
   ``colptrs[N]`` nonzeros.

.. c:function:: int SlsMatvec(SlsMat A, realtype *x, realtype *y)

   Computes the sparse matrix-vector product, :math:`y=Ax`.  If *A*
   is a sparse matrix of dimension :math:`M\times N`, then it is assumed that *x*
   is a ``realtype`` array of  length :math:`N`, and *y* is a
   ``realtype`` array of length :math:`M`. Upon successful completion, the
   return value is zero; otherwise 1 is returned.

.. c:function:: void PrintSparseMat(DlsMat A)

   Prints a :c:type:`SlsMat` matrix to standard output.





The KLU solver
-------------------------------------------

KLU is a sparse matrix factorization and solver library written by Tim
Davis [KLU]_.  In order to use KLU-enabled SUNDIALS solvers, it is
assumed that KLU has been installed on the system prior to
installation of SUNDIALS, and that SUNDIALS has been configured
appropriately to link with KLU (see :ref:`Installation` for details).

Designed for serial calculations only, KLU is supported for
calculations employing SUNDIALS' serial or shared-memory parallel
``N_Vector`` modules (see :ref:`NVectors.NVSerial`,
:ref:`NVectors.OpenMP` and :ref:`NVectors.Pthreads`).



The SuperLU_MT solver
-------------------------------------------

SuperLU_MT is a threaded sparse matrix factorization and solver
library written by X. Sherry Li [SuperLUMT]_.  In order to use 
SuperLU_MT enabled SUNDIALS solvers, it is assumed that SuperLU_MT has
been installed on the system prior to installation of SUNDIALS, and
that SUNDIALS has been configured appropriately to link with
SuperLU_MT (see :ref:`Installation` for details).

Designed for serial and threaded calculations only, SuperLU_MT is
supported for calculations employing SUNDIALS' serial or shared-memory
parallel ``N_Vector`` modules (see :ref:`NVectors.NVSerial`,
:ref:`NVectors.OpenMP` and :ref:`NVectors.Pthreads`).
