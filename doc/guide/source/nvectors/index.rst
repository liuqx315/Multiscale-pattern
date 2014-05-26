..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _NVectors:

=======================
Vector Data Structures
=======================

The SUNDIALS library comes packaged with two NVECTOR implementations,
one designed for serial simulations and the second for
distributed-memory parallel simulations.  Both implementations assume
that the process-local data is stored contiguously, and they in turn
provide a variety of standard vector algebra operations that may be
performed on the data.  

In addition, SUNDIALS provides a simple interface for generic vectors
(akin to a C++ *abstract base class*).  All of the major SUNDIALS
solvers (CVODE, IDA, KINSOL, ARKODE) in turn are constructed to only
depend on these generic vector operations, making them immediately
extensible to new user-defined vector objects.  The only exceptions to
this rule relate to the dense and banded linear system solvers, since
they rely on particular data storage and access patterns in the
NVECTORS used.

Details on the generic NVECTOR module are below.  However, to jump to
specific descriptions of the serial and parallel vector modules
provided by SUNDIALS, or ARKode's requirements for routines comprising
a user-supplied NVECTOR module, the following links are provided:

.. toctree::
   :maxdepth: 1

   NVector_Serial
   NVector_Parallel
   NVector_OpenMP
   NVector_Pthreads
   ARKode_requirements




.. _NVectors.Description:

Description of the NVECTOR Modules
======================================

The SUNDIALS solvers are written in a data-independent manner. They
all operate on generic vectors (of type ``N_Vector``) through a set of
operations defined by, and specific to, the particular NVECTOR
implementation. Users can provide a custom implementation of the
NVECTOR module or use one of two provided within SUNDIALS, one serial
and the other an MPI parallel implementation.

The generic ``N_Vector`` type is a pointer to a structure that has an
implementation-dependent *content* field containing the description
and actual data of the vector, and an *ops* field pointing to a
structure with generic vector operations. The type ``N_Vector`` is
defined as:

.. code-block:: c

   typedef struct _generic_N_Vector *N_Vector;
   
   struct _generic_N_Vector { 
      void *content;
      struct _generic_N_Vector_Ops *ops;
   };

Here, the ``_generic_N_Vector_Op`` structure is essentially a list of
function pointers to the various actual vector operations, and is
defined as  

.. code-block:: c

   struct _generic_N_Vector_Ops { 
      N_Vector    (*nvclone)(N_Vector); 
      N_Vector    (*nvcloneempty)(N_Vector); 
      void        (*nvdestroy)(N_Vector); 
      void        (*nvspace)(N_Vector, long int *, long int *); 
      realtype*   (*nvgetarraypointer)(N_Vector); 
      void        (*nvsetarraypointer)(realtype *, N_Vector); 
      void        (*nvlinearsum)(realtype, N_Vector, realtype, N_Vector, N_Vector);
      void        (*nvconst)(realtype, N_Vector);
      void        (*nvprod)(N_Vector, N_Vector, N_Vector); 
      void 	  (*nvdiv)(N_Vector, N_Vector, N_Vector);
      void	  (*nvscale)(realtype, N_Vector, N_Vector);
      void	  (*nvabs)(N_Vector, N_Vector); 
      void	  (*nvinv)(N_Vector, N_Vector);
      void	  (*nvaddconst)(N_Vector, realtype, N_Vector);
      realtype	  (*nvdotprod)(N_Vector, N_Vector); 
      realtype	  (*nvmaxnorm)(N_Vector);
      realtype	  (*nvwrmsnorm)(N_Vector, N_Vector);
      realtype	  (*nvwrmsnormmask)(N_Vector, N_Vector, N_Vector);
      realtype	  (*nvmin)(N_Vector);
      realtype	  (*nvwl2norm)(N_Vector, N_Vector); 
      realtype	  (*nvl1norm)(N_Vector);
      void	  (*nvcompare)(realtype, N_Vector, N_Vector); 
      booleantype (*nvinvtest)(N_Vector, N_Vector); 
      booleantype (*nvconstrmask)(N_Vector, N_Vector, N_Vector); 
      realtype	  (*nvminquotient)(N_Vector, N_Vector);
   };


The generic NVECTOR module defines and implements the vector
operations acting on a ``N_Vector``. These routines are nothing but
wrappers for the vector operations defined by a particular NVECTOR
implementation, which are accessed through the *ops* field of the
``N_Vector`` structure. To illustrate this point we show below the
implementation of a typical vector operation from the generic NVECTOR
module, namely ``N_VScale``, which performs the scaling of a vector
``x`` by a scalar ``c``:

.. code-block:: c

   void N_VScale(realtype c, N_Vector x, N_Vector z) {
      z->ops->nvscale(c, x, z);
   }

The subsection :ref:`NVectors.Ops` contains a complete list of all
vector operations defined by the generic NVECTOR module.  Finally, we
note that the generic NVECTOR module defines the functions
``N_VCloneVectorArray`` and ``N_VCloneEmptyVectorArray``. Both
functions create (by cloning) an array of *count* variables of type
``N_Vector``, each of the same type as an existing ``N_Vector``. Their
prototypes are: 

.. code-block:: c

   N_Vector *N_VCloneVectorArray(int count, N_Vector w);
   N_Vector *N_VCloneEmptyVectorArray(int count, N_Vector w);

and their definitions are based on the implementation-specific
``N_VClone`` and ``N_VCloneEmpty`` operations, respectively. 

Similarly, an array of variables of type ``N_Vector`` can be destroyed
by calling ``N_VDestroyVectorArray``, whose prototype is 

.. code-block:: c
   
   void N_VDestroyVectorArray(N_Vector *vs, int count); 

and whose definition is based on the implementation-specific
``N_VDestroy`` operation. 



In particular, any implementation of the NVECTOR module **must**:

* Specify the *content* field of the ``N_Vector``.

* Define and implement the necessary vector operations. Note that the
  names of these routines should be unique to that implementation in
  order to permit using more than one NVECTOR module (each with
  different ``N_Vector`` internal data representations) in the same
  code.  We further note that not all of the defined operations are
  required for each solver in SUNDIALS.  The list of required
  operations for use with ARKode is given in the section
  :reF:`NVectors.ARKode`. 

* Define and implement user-callable constructor and destructor
  routines to create and free a ``N_Vector`` with the new *content*
  field and with *ops* pointing to the new vector operations. 

* Optionally, define and implement additional user-callable routines
  acting on the newly defined ``N_Vector`` (e.g., a routine to print the
  *content* for debugging purposes). 

* Optionally, provide accessor macros as needed for that particular
  implementation to be used to access different parts in the content
  field of the newly defined ``N_Vector``. 



.. _NVectors.Ops:

Description of the NVECTOR operations
=========================================

For each of the ``N_vector`` operations, we give the name, usage
of the function, and a description of its mathematical operations
below.


.. c:function:: N_Vector N_VClone(N_Vector w)

   Creates a new ``N_Vector`` of the same type as an existing vector
   *w* and sets the *ops* field.  It does not copy the vector, but
   rather allocates storage for the new vector.
   
   Usage:

   .. code-block:: c

      v = N_VClone(w);


.. c:function:: N_Vector N_VCloneEmpty(N_Vector w)

   Creates a new ``N_Vector`` of the same type as an existing vector
   *w* and sets the *ops* field.  It does not allocate storage for the
   new vector's data array. 

   Usage:

   .. code-block:: c

      v = N VCloneEmpty(w);


.. c:function:: void N_VDestroy(N_Vector v)
 
   Destroys the ``N_Vector`` *v* and frees memory allocated for its
   internal data.  

   Usage:

   .. code-block:: c

      N_VDestroy(v);


.. c:function:: void N_VSpace(N_Vector v, long int* lrw, long int* liw)

   Returns storage requirements for the ``N_Vector`` *v*: *lrw* contains
   the number of ``realtype`` words and *liw* contains the number of
   integer words.  This function is advisory only, for use in
   determining a user's total space requirements; it could be a dummy
   function in a user-supplied NVECTOR module if that information is
   not of interest.  
 
   Usage:

   .. code-block:: c
 
      N_VSpace(nvSpec, &lrw, &liw);


.. c:function:: realtype* N_VGetArrayPointer(N_Vector v)

   Returns a pointer to a ``realtype`` array from the ``N_Vector``
   *v*.  Note that this assumes that the internal data in the
   ``N_Vector`` is a contiguous array of ``realtype``. This routine is
   only used in the solver-specific interfaces to the dense and banded
   (serial) linear solvers, and in the interfaces to the banded
   (serial) and band-block-diagonal (parallel) preconditioner modules
   provided with SUNDIALS.  

   Usage:

   .. code-block:: c

      vdata = NVGetArrayPointer(v);


.. c:function:: void N_VSetArrayPointer(realtype* vdata, N_Vector v)

   Replaces the data array pointer in an ``N_Vector`` with a given
   array of ``realtype``.  Note that this assumes that the internal
   data in the ``N_Vector`` is a contiguous array of
   ``realtype``. This routine is only used in the interfaces to the
   dense (serial) linear solver, hence need not exist in a
   user-supplied NVECTOR module. 

   Usage:

   .. code-block:: c

      NVSetArrayPointer(vdata,v);


.. c:function:: void N_VLinearSum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)

   Performs the operation *z = ax + by*, where *a* and *b* are 
   scalars and *x* and *y* are of type ``N_Vector``:

   .. math::
      z_i = a x_i + b y_i, \quad i=1,\ldots,n. 

   Usage:

   .. code-block:: c

      N_VLinearSum(a, x, b, y, z);


.. c:function:: void N_VConst(realtype c, N_Vector z)

   Sets all components of the ``N_Vector`` *z* to the scalar value
   *c*:

   .. math::
      z_i = c, \quad i=1,\ldots,n. 

   Usage:

   .. code-block:: c

      N_VConst(c, z);


.. c:function:: void N_VProd(N_Vector x, N_Vector y, N_Vector z)

   Sets the ``N_Vector`` *z* to be the component-wise product of the 
   ``N_Vector`` inputs *x* and *y*:

   .. math::
      z_i = x_i y_i, \quad i=1,\ldots,n.
 
   Usage:
 
   .. code-block:: c

      N_VProd(x, y, z);


.. c:function:: void N_VDiv(N_Vector x, N_Vector y, N_Vector z)

   Sets the ``N_Vector`` *z* to be the component-wise ratio of the
   ``N_Vector`` inputs *x* and *y*: 

   .. math::
      z_i = \frac{x_i}{y_i}, \quad i=1,\ldots,n.

   The :math:`y_i` may not be tested for 0 values. It should only be
   called with a *y* that is guaranteed to have all nonzero components.  

   Usage: 

   .. code-block:: c

      N_VDiv(x, y, z);


.. c:function:: void N_VScale(realtype c, N_Vector x, N_Vector z)

   Scales the ``N_Vector`` *x* by the scalar *c* and returns the
   result in *z*:

   .. math::
      z_i = c x_i, \quad i=1,\ldots,n.

   Usage:

   .. code-block:: c

      N_VScale(c, x, z);


.. c:function:: void N_VAbs(N_Vector x, N_Vector z)
 
   Sets the components of the ``N_Vector`` *z* to be the absolute
   values of the components of the ``N_Vector`` *x*: 

   .. math::
      y_i = |x_i|, \quad i=1,\ldots,n.

   Usage:

   .. code-block:: c

      N_VAbs(x, z);


.. c:function:: void N_VInv(N_Vector x, N_Vector z)

   Sets the components of the ``N_Vector`` *z* to be the inverses of
   the components of the ``N_Vector`` *x*: 

   .. math::
      z_i = 1.0/x_i, \quad i=1,\ldots,n.

   This routine may not check for division by 0.  It should be called
   only with an *x* which is guaranteed to have all nonzero components.

   Usage:

   .. code-block:: c

      N_VInv(x, z);


.. c:function:: void N_VAddConst(N_Vector x, realtype b, N_Vector z)

   Adds the scalar *b* to all components of *x* and returns the
   result in the ``N_Vector`` *z*:

   .. math::
      z_i = x_i+b, \quad i=1,\ldots,n.

   Usage:

   .. code-block:: c

      N_VAddConst(x, b, z);


.. c:function:: realtype N_VDotProd(N_Vector x, N_Vector z)

   Returns the value of the dot-product of the ``N_Vectors`` *x* and *y*:

   .. math::
      d = \sum_{i=1}^{n} x_i y_i.

   Usage:

   .. code-block:: c

      d = N_VDotProd(x, y);


.. c:function:: realtype N_VMaxNorm(N_Vector x)
 
   Returns the value of the :math:`l_{\infty}` norm of the
   ``N_Vector`` *x*:

   .. math::
      m = \max_{1\le i\le n} |x_i|.

   Usage: 
 
   .. code-block:: c

      m = N_VMaxNorm(x);


.. c:function:: realtype N_VWrmsNorm(N_Vector x, N_Vector w)

   Returns the weighted root-mean-square norm of the ``N_Vector`` *x*
   with (positive) weight vector *w*: 
 
   .. math::
      m = \left( \frac1n \sum_{i=1}^{n} \left(x_i w_i\right)^2\right)^{1/2}.  

   Usage:

   .. code-block:: c

      m = N_VWrmsNorm(x, w);


.. c:function:: realtype N_VWrmsNormMask(N_Vector x, N_Vector w, N_Vector id)

   Returns the weighted root mean square norm of the ``N_Vector`` *x*
   with (positive) weight vector *w* built using only the elements of *x*
   corresponding to nonzero elements of the ``N_Vector`` *id*:
  
   .. math::
      m = \left( \frac1n \sum_{i=1}^{n} \left(x_i w_i \text{sign}(id_i)\right)^2 \right)^{1/2}. 

   .. code-block:: c

      m = N_VWrmsNormMask(x, w, id);

.. c:function:: realtype N_VMin(N_Vector x)
 
   Returns the smallest element of the ``N_Vector`` *x*: 

   .. math::
      m = \min_{1\le i\le n} x_i.

   Usage:

   .. code-block:: c

      m = N_VMin(x);

.. c:function:: realtype N_VWl2Norm(N_Vector x, N_Vector w)

   Returns the weighted Euclidean :math:`l_2` norm of the ``N_Vector``
   *x* with weight vector *w*: 

   .. math::
      m = \left(\sum_{i=1}^{n}\left(x_i w_i\right)^2\right)^{1/2}.  

   Usage:

   .. code-block:: c

      m = N_VWL2Norm(x, w);

.. c:function:: realtype N_VL1Norm(N_Vector x)

   Returns the :math:`l_1` norm of the ``N_Vector`` *x*: 

   .. math::
      m = \sum_{i=1}^{n} |x_i|. 

   Usage:

   .. code-block:: c

      m = N_VL1Norm(x);


.. c:function:: void N_VCompare(realtype c, N_Vector x, N_Vector z)

   Compares the components of the ``N_Vector`` *x* to the scalar *c*
   and returns an ``N_Vector`` *z* such that for all :math:`1\le i\le n`,

   .. math::
      z_i = \begin{cases} 1.0 &\quad\text{if}\; |x_i| \ge c,\\
                          0.0 &\quad\text{otherwise}\end{cases}.

   Usage:

   .. code-block:: c

      N_VCompare(c, x, z);

.. c:function:: booleantype N_VInvTest(N_Vector x, N_Vector z)
 
   Sets the components of the ``N_Vector`` *z* to be the inverses of
   the components of the ``N_Vector`` *x*, with prior testing for
   zero values: 

   .. math::
      z_i = 1.0/x_i, \quad i=1,\ldots,n.  

   This routine returns ``TRUE`` if all components of *x* are nonzero
   (successful inversion) and returns ``FALSE`` otherwise.

   Usage:

   .. code-block:: c

      t = N_VInvTest(x, z);

.. c:function:: booleantype N_VConstrMask(N_Vector c, N_Vector x, N_Vector m)
 
   Performs the following constraint tests based on the values in
   :math:`c_i`: 

   .. math::
      x_i > 0 \;\text{if}\; c_i = 2, \\
      x_i \ge 0 \;\text{if}\; c_i = 1, \\
      x_i < 0 \;\text{if}\; c_i = -2, \\
      x_i \le 0 \;\text{if}\; c_i = -1.

   There is no constraint on :math:`x_i` if :math:`c_i = 0`. This
   routine returns ``FALSE`` if any element failed the constraint test,
   ``TRUE`` if all passed. It also sets a mask vector *m*, with
   elements equal to 1.0 where the constraint test failed, and 0.0
   where the test passed. This routine is used only for constraint
   checking.

   Usage:

   .. code-block:: c

      t = N_VConstrMask(c, x, m);

.. c:function:: realtype N_VMinQuotient(N_Vector num, N_Vector denom)

   This routine returns the minimum of the quotients obtained by
   termwise dividing the elements of *n* by the elements in *d*:

   .. math::
      \min_{i=1,\ldots,n} \frac{\text{num}_i}{\text{denom}_i}. 

   A zero element in *denom* will be skipped.  If no such quotients
   are found, then the large value ``BIG_REAL`` (defined in the header
   file ``sundials_types.h``) is returned. 

   Usage:

   .. code-block:: c

      minq = N_VMinQuotient(num, denom);

