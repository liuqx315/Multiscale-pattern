.. _NVectors:

Description of the NVECTOR Modules
======================================

The SUNDIALS solvers are written in a data-independent manner. They
all operate on generic vectors (of type ``N_Vector``) through a set of
operations defined by the particular NVECTOR implementation. Users can
provide their own specific implementation of the NVECTOR module or use
one of two provided within SUNDIALS, a serial and an MPI parallel
implementation.

The generic ``N_Vector`` type is a pointer to a structure that has an
implementation-dependent `content` field containing the description
and actual data of the vector, and an `ops` field pointing to a
structure with generic vector operations. The type ``N_Vector`` is
defined as:

.. code-block: c

   typedef struct _generic_N_Vector *N_Vector;
   
   struct _generic_N_Vector { 
      void *content;
      struct _generic_N_Vector_Ops *ops;
   };

The ``_generic_N_Vector_Op``` structure is essentially a list of
pointers to the various actual vector operations, and is defined as 

.. code-block: c

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
implementation, which are accessed through the `ops` field of the
``N_Vector`` structure. To illustrate this point we show below the
implementation of a typical vector operation from the generic NVECTOR
module, namely ``N_VScale``, which performs the scaling of a vector
``x`` by a scalar ``c``:

.. code-block:: c

   void N_VScale(realtype c, N_Vector x, N_Vector z) 
   {
      z->ops->nvscale(c, x, z);
   }

The subsection :ref:`NVectors.Ops` contains a complete list of all
vector operations defined by the generic NVECTOR module. Finally, note
that the generic NVECTOR module defines the functions
``N_VCloneVectorArray`` and ``N_VCloneEmptyVectorArray``. Both
functions create (by cloning) an array of ``count`` variables of type
``N_Vector``, each of the same type as an existing ``N_Vector``. Their
prototypes are: 

.. code-block:: c

   N_Vector *N_VCloneVectorArray(int count, N_Vector w);
   N_Vector *N_VCloneEmptyVectorArray(int count, N_Vector w);

and their definitions are based on the implementation-specific
``N_VClone`` and ``N_VCloneEmpty`` operations, respectively. 

An array of variables of type ``N_Vector`` can be destroyed by calling
``N_VDestroyVectorArray``, whose prototype is 

.. code-block:: c
   
   void N_VDestroyVectorArray(N_Vector *vs, int count); 

and whose definition is based on the implementation-specific
``N_VDestroy`` operation. 

A particular implementation of the NVECTOR module must:

* Specify the `content` field of the ``N_Vector``.

* Define and implement the vector operations. Note that the names of
  these routines should be unique to that implementation in order to
  permit using more than one NVECTOR module (each with different
  ``N_Vector`` internal data representations) in the same code. 

* Define and implement user-callable constructor and destructor
  routines to create and free a ``N_Vector`` with the new `content`
  field and with `ops` pointing to the new vector operations. 

* Optionally, define and implement additional user-callable routines
  acting on the newly defined ``N_Vector`` (e.g., a routine to print the
  `content` for debugging purposes). 

* Optionally, provide accessor macros as needed for that particular
  implementation to be used to access different parts in the content
  field of the newly defined ``N_Vector``. 



.. _NVectors.Ops:

Description of the NVECTOR operations
-----------------------------------------

For each of the ``N_vector`` operations, we give the name, usage
of the function, and a description of its mathematical operations
below.

* N_VClone

  .. code-block:: c

     v = N_VClone(w);

  Creates a new ``N_Vector`` of the same type as an existing vector
  ``w`` and sets the `ops` field. It does not copy the vector, but
  rather allocates storage for the new vector.

* N_VCloneEmpty

  .. code-block:: c

     v = N VCloneEmpty(w);

  Creates a new ``N_Vector`` of the same type as an existing vector
  ``w`` and sets the `ops` field. It does not allocate storage for the
  data array. 

* N_VDestroy

  .. code-block:: c

     N_VDestroy(v);

  Destroys the ``N_Vector v`` and frees memory allocated for its
  internal data.  

* N_VSpace

  .. code-block:: c

     N_VSpace(nvSpec, &lrw, &liw);

  Returns storage requirements for one ``N_Vector``. ``lrw`` contains
  the number of ``realtype`` words and ``liw`` contains the number of
  integer words. This function is advisory only, for use in
  determining a user's total space requirements; it could be a dummy
  function in a user-supplied NVECTOR module if that information is
  not of interest.  

* N_VGetArrayPointer

  .. code-block:: c

     vdata = NVGetArrayPointer(v);

  Returns a pointer to a ``realtype`` array from the ``N_Vector
  v``. Note that this assumes that the internal data in the
  ``N_Vector`` is a contiguous array of ``realtype``. This routine is
  only used in the solver-specific interfaces to the dense and banded
  (serial) linear solvers, and in the interfaces to the banded
  (serial) and band-block-diagonal (parallel) preconditioner modules
  provided with SUNDIALS.  

* N_VSetArrayPointer

  .. code-block:: c

     NVSetArrayPointer(vdata,v);

  Overwrites the data in an ``N_Vector`` with a given array of
  ``realtype``. Note that this assumes that the internal data in the
  ``N_Vector`` is a contiguous array of ``realtype``. This routine is
  only used in the interfaces to the dense (serial) linear solver,
  hence need not exist in a user-supplied NVECTOR module.

* N_VLinearSum

  .. code-block:: c

     N_VLinearSum(a, x, b, y, z);

  Performs the operation ``z = ax + by``, where ``a`` and ``b`` are
  scalars and ``x`` and ``y`` are of type ``N_Vector``: :math:`z_i = a
  x_i + b y_i, \; i=0,\ldots,n-1`. 

* N_VConst

  .. code-block:: c

     N_VConst(c, z);

  Sets all components of the ``N_Vector`` ``z`` to ``c``: :math:`z_i =
  c, \; i=0,\ldots,n-1`. 

* N_VProd

  .. code-block:: c

     N_VProd(x, y, z);

  Sets the ``N_Vector z`` to be the component-wise product of the 
  ``N_Vector`` inputs ``x`` and ``y``: :math:`z_i = x_i y_i, \;
  i=0,\ldots,n-1`.

* N_VDiv

  .. code-block:: c

     N_VDiv(x, y, z);

  Sets the ``N_Vector`` ``z`` to be the component-wise ratio of the
  ``N_Vector`` inputs ``x`` and ``y``: :math:`z_i = x_i/y_i, \;
  i=0,\ldots,n-1`.  The yi may not be tested for 0 values. It should
  only be called with a ``y`` that is guaranteed to have all nonzero
  components.  

* N_VScale

  .. code-block:: c

     N_VScale(c, x, z);

  Scales the ``N_Vector`` ``x`` by the scalar ``c`` and returns the
  result in ``z``: :math:`z_i = c x_i, \; i=0,\ldots,n-1`.

* N_VAbs

  .. code-block:: c

     N_VAbs(x, z);

  Sets the components of the ``N_Vector`` ``z`` to be the absolute
  values of the components of the ``N_Vector`` ``x``: :math:`y_i =
  |x_i|, \; i=0,\ldots,n-1`.

* N_VInv

  .. code-block:: c

     N_VInv(x, z);

  Sets the components of the ``N_Vector`` ``z`` to be the inverses of
  the components of the ``N_Vector`` ``x``: :math:`z_i = 1.0/x_i, \;
  i=0,\ldots,n-1`.  This routine may not check for division by 0. It
  should be called only with an x which is guaranteed to have all
  nonzero components.

* N_VAddConst

  .. code-block:: c

     N_VAddConst(x, b, z);

  Adds the scalar ``b`` to all components of ``x`` and returns the
  result in the ``N_Vector`` ``z``: :math:`z_i = x_i+b, \;
  i=0,\ldots,n-1`.  

* N_VDotProd

  .. code-block:: c

     d = N_VDotProd(x, y);

  Returns the value of the ordinary dot product of ``x`` and ``y``:
  :math:`d = \sum_{i=0}^{n-1} x_i y_i`.

* N_VMaxNorm

  .. code-block:: c

     m = N_VMaxNorm(x);

  Returns the maximum norm of the ``N_Vector x``: :math:`m = \max_i
  |x_i|`.

* N_VWrmsNorm

  .. code-block:: c

     m = N_VWrmsNorm(x, w);

  Returns the weighted root-mean-square norm of the ``N_Vector`` ``x``
  with weight vector ``w``: 
 
  .. math::
     m = \left( \frac1n \sum_{i=0}^{n-1} \left(x_i w_i\right)^2\right)^{1/2}.  

* N_VWrmsNormMask

  .. code-block:: c

     m = N_VWrmsNormMask(x, w, id);

  Returns the weighted root mean square norm of the ``N_Vector`` ``x``
  with weight vector ``w`` built using only the elements of ``x``
  corresponding to nonzero elements of the ``N_Vector`` ``id``:
  
  .. math::
     m = \left( \frac1n \sum_{i=0}^{n-1} \left(x_i w_i \text{sign}(id_i)\right)^2 \right)^{1/2}. 

* N_VMin

  .. code-block:: c

     m = N_VMin(x);

  Returns the smallest element of the ``N_Vector x``: :math:`m =
  \min_i x_i`.

* N_VWl2Norm

  .. code-block:: c

     m = N_VWL2Norm(x, w);

  Returns the weighted Euclidean :math:`l_2` norm of the ``N_Vector
  x`` with weight vector ``w``: 

  .. math::
     m = \left(\sum_{i=0}^{n-1}\left(x_i w_i\right)^2\right)^{1/2}.  

* N_VL1Norm

  .. code-block:: c

     m = N_VL1Norm(x);

  Returns the :math:`l_1` norm of the ``N_Vector x``: :math:`m = \sum_{i=0}^{n-1} |x_i|`. 

* N_VCompare

  .. code-block:: c

     N_VCompare(c, x, z);

  Compares the components of the ``N_Vector x`` to the scalar ``c``
  and returns an ``N_Vector z`` such that: 

  .. math::
     z_i = \begin{cases} 1.0 &\;\text{if}\; |x_i| \ge c,\\
                         0.0 &\;\text{otherwise}\end{cases}.

* N_VInvTest

  .. code-block:: c

     t = N_VInvTest(x, z);

  Sets the components of the ``N_Vector`` ``z`` to be the inverses of
  the components of the ``N_Vector`` ``x``, with prior testing for
  zero values: :math:z_i = 1.0/x_i, \; i=0,\ldots,n-1`.  This routine
  returns ``TRUE`` if all components of ``x`` are nonzero (successful
  inversion) and returns ``FALSE`` otherwise.

* N_VConstrMask

  .. code-block:: c

     t = N_VConstrMask(c, x, m);

  Performs the following constraint tests based on the values in
  :math:`c_i`: :math:`x_i > 0 \;\text{if}\; c_i = 2,\quad`
  :math:`x_i \ge 0 \;\text{if}\; c_i = 1,\quad`
  :math:`x_i < 0 \;\text{if}\; c_i = -2,\quad`
  :math:`x_i \le 0 \;\text{if}\; c_i = -1.\quad`
  There is no constraint on :math:`x_i` if :math:`c_i = 0`. This
  routine returns ``FALSE`` if any element failed the constraint test,
  ``TRUE`` if all passed. It also sets a mask vector ``m``, with
  elements equal to 1.0 where the constraint test failed, and 0.0
  where the test passed. This routine is used only for constraint
  checking. 

* N_VMinQuotient

  .. code-block:: c

     minq = N_VMinQuotient(n, d);

  This routine returns the in ``minq`` the minimum of the quotients
  obtained by termwise dividing :math:`n_i/d_i`. A zero element in
  ``d`` will be skipped. If no such quotients are found, then the
  large value ``BIG_REAL`` (defined in the header file
  ``sundials_types.h``) is returned. 





.. _NVectors.NVSerial:

The NVECTOR_SERIAL implementation
-----------------------------------------

The serial implementation of the NVECTOR module provided with
SUNDIALS, NVECTOR_SERIAL, defines the `content` field of a
``N_Vector`` to be a structure containing the length of the vector, a
pointer to the beginning of a contiguous data array, and a boolean
flag `own_data` which specifies the ownership of data. 

.. code-block:: c

   struct _N_VectorContent_Serial { 
      long int length; 
      booleantype own_data; 
      realtype *data;
   };

The following five macros are provided to access the content of an
NVECTOR_SERIAL vector. The suffix ``_S`` in the names denotes serial
version. 

* ``NV_CONTENT_S``

  This routine gives access to the contents of the serial vector
  ``N_Vector``. 

  The assignment ``v_cont = NV_CONTENT_S(v)`` sets ``v_cont`` to be a
  pointer to the serial ``N_Vector`` `content` structure. 

  Implementation:
  
  .. code-block:: c

     #define NV_CONTENT_S(v) ( (N_VectorContent_Serial)(v->content) ) 

* ``NV_OWN_DATA_S``, ``NV_DATA_S``, ``NV_LENGTH_S``

  These macros give individual access to the parts of the content of a
  serial ``N_Vector``. 
  
  The assignment ``v_data = NV_DATA_S(v)`` sets ``v_data`` to be a
  pointer to the first component of the `data` for the ``N_Vector
  v``. 

  The assignment ``NV_DATA_S(v) = v_data`` sets the component
  array of ``v`` to be ``v_data`` by storing the pointer ``v_data``.

  The assignment ``v_len = NV_LENGTH_S(v)`` sets ``v_len`` to be the
  `length` of ``v``. On the other hand, the call ``NV_LENGTH_S(v) =
  len_v`` sets the `length` of ``v`` to be ``len_v``. 

  Implementation:

  .. code-block:: c
 
     #define NV_OWN_DATA_S(v) ( NV_CONTENT_S(v)->own_data ) 
     #define NV_DATA_S(v) ( NV_CONTENT_S(v)->data ) 
     #define NV_LENGTH_S(v) ( NV_CONTENT_S(v)->length )

* ``NV_Ith_S``

  This macro gives access to the individual components of the `data`
  array of an ``N_Vector``. 

  The assignment ``r = NV_Ith_S(v,i)`` sets ``r`` to be the value of
  the ``i``-th component of ``v``. 

  The assignment ``NV_Ith_S(v,i) = r`` sets the value of the ``i``-th
  component of ``v`` to be ``r``. 

  Here ``i`` ranges from 0 to :math:`n-1` for a vector of length
  :math:`n`. 

  Implementation: 

  .. code-block:: c

    #define NV_Ith_S(v,i) ( NV_DATA_S(v)[i] )

The NVECTOR_SERIAL module defines serial implementations of all vector
operations listed in the section :ref:`NVectors.Ops`. Their names are
obtained from those in that section by appending the suffix
``_Serial``. The module NVECTOR_SERIAL provides the following
additional user-callable routines: 

* ``N_VNew_Serial`` 

  This function creates and allocates memory for a serial
  ``N_Vector``. Its only argument is the vector length.

  .. code-block:: c

     N_Vector N_VNew_Serial(long int vec_length);

* ``N_VNewEmpty_Serial``

  This function creates a new serial ``N_Vector`` with an empty
  (``NULL``) data array. 

  .. code-block:: c

     N_Vector N_VNewEmpty_Serial(long int vec_length);

* ``N_VMake_Serial``

  This function creates and allocates memory for a serial vector with
  user-provided data array. 

  .. code-block:: c

     N_Vector N_VMake_Serial(long int vec_length, realtype *v_data); 

* ``N_VCloneVectorArray_Serial``

  This function creates (by cloning) an array of ``count`` serial
  vectors. 

  .. code-block:: c

     N_Vector *N_VCloneVectorArray_Serial(int count, N_Vector w);

* ``N_VCloneEmptyVectorArray_Serial``

  This function creates (by cloning) an array of ``count`` serial
  vectors, each with an empty (```NULL``) data array.

  .. code-block:: c

     N_Vector *N_VCloneEmptyVectorArray_Serial(int count, N_Vector w);

* ``N_VDestroyVectorArray_Serial``
  
  This function frees memory allocated for the array of ``count``
  variables of type ``N_Vector`` created with
  ``N_VCloneVectorArray_Serial`` or with
  ``N_VCloneEmptyVectorArray_Serial``. 

  .. code-block:: c

     void N_VDestroyVectorArray_Serial(N_Vector *vs, int count);

* ``N_VPrint_Serial``

  This function prints the content of a serial vector to ``stdout``.

  .. code-block:: c

     void N_VPrint_Serial(N_Vector v);

**Notes**

* When looping over the components of an ``N_Vector v``, it is more
  efficient to first obtain the component array via ``v_data =
  NV_DATA_S(v)`` and then access ``v_data[i]`` within the loop than it
  is to use ``NV_Ith_S(v,i)`` within the loop. 
* ``N_VNewEmpty_Serial``, ``N_VMake_Serial``, and
  ``N_VCloneEmptyVectorArray_Serial`` set the field `own_data` to
  ``FALSE``.  ``N_VDestroy_Serial`` and
  ``N_VDestroyVectorArray_Serial`` will not attempt to free the
  pointer data for any ``N_Vector`` with `own_data` set to ``FALSE``.
  In such a case, it is the user's responsibility to deallocate the
  data pointer. 
* To maximize efficiency, vector operations in the NVECTOR_SERIAL
  implementation that have more than one ``N_Vector`` argument do not
  check for consistent internal representation of these vectors. It is
  the user's responsibility to ensure that such routines are called
  with ``N_Vector`` arguments that were all created with the same
  internal representations. 





.. _NVectors.NVParallel:

The NVECTOR_PARALLEL implementation
-----------------------------------------

The parallel implementation of the NVECTOR module provided with
SUNDIALS, NVECTOR_PARALLEL, defines the `content` field of a
``N_Vector`` to be a structure containing the global and local lengths
of the vector, a pointer to the beginning of a contiguous local data
array, an MPI communicator, an a boolean flag `own_data` indicating
ownership of the data array `data`. 

.. code-block:: c

   struct _N_VectorContent_Parallel { 
      long int local_length; 
      long int global_length; 
      booleantype own_data;
      realtype *data;
      MPI_Comm comm; 
   };

The following seven macros are provided to access the content of a
NVECTOR_PARALLEL vector. The suffix ``_P`` in the names denotes
parallel version. 

* ``NV_CONTENT_P``
 
  This macro gives access to the contents of the parallel vector
  ``N_Vector``. 

  The assignment ``v_cont = NV_CONTENT_P(v)`` sets ``v_cont`` to be a
  pointer to the ``N_Vector`` `content` structure of type ``struct
  N_VectorParallelContent``. 

  Implementation:

  .. code-block:: c

     #define NV_CONTENT_P(v) ( (N_VectorContent_Parallel)(v->content) )

* ``NV_OWN_DATA_P``, ``NV_DATA_P``, ``NV_LOCLENGTH_P``,
  ``NV_GLOBLENGTH_P``

  These macros give individual access to the parts of the content of a
  parallel ``N_Vector``.
 
  The assignment ``v_data = NV_DATA_P(v)`` sets ``v_data`` to be a
  pointer to the first component of the `local_data` for the
  ``N_Vector v``. 

  The assignment ``NV_DATA_P(v) = v_data`` sets the component array of
  ``v`` to be ``v_data`` by storing the pointer ``v_data`` into
  `data`.

  The assignment ``v_llen = NV_LOCLENGTH_P(v)`` sets ``v_llen`` to be
  the length of the local part of ``v``. 

  The call ``NV_LENGTH_P(v) = llen_v`` sets the `local_length` of
  ``v`` to be ``llen_v``. 

  The assignment ``v_glen = NV_GLOBLENGTH_P(v)`` sets ``v_glen`` to be
  the `global_length` of the vector ``v``. The call
  ``NV_GLOBLENGTH_P(v) = glen_v`` sets the `global_length` of ``v`` to
  be ``glen_v``. 

  Implementation:
 
  .. code-block:: c

     #define NV_OWN_DATA_P(v)   ( NV_CONTENT_P(v)->own_data ) 
     #define NV_DATA_P(v)       ( NV_CONTENT_P(v)->data ) 
     #define NV_LOCLENGTH_P(v)  ( NV_CONTENT_P(v)->local_length ) 
     #define NV_GLOBLENGTH_P(v) ( NV_CONTENT_P(v)->global_length )

* ``NV_COMM_P``
 
  This macro provides access to the MPI communicator used by the
  NVECTOR_PARALLEL vectors. 

  Implementation: 

  .. code-block:: c

     #define NV_COMM_P(v) ( NV_CONTENT_P(v)->comm )

* ``NV_Ith_P``

  This macro gives access to the individual components of the
  `local_data` array of an ``N_Vector``. 

  The assignment ``r = NV_Ith_P(v,i)`` sets ``r`` to be the value of
  the ``i``-th component of the local part of ``v``. 

  The assignment ``NV_Ith_P(v,i) = r`` sets the value of the ``i``-th
  component of the local part of ``v`` to be ``r``.

  Here ``i`` ranges from 0 to :math:`n-1`, where :math:`n` is the
  `local_length`. 

  Implementation: 

  .. code-block:: c
  
     #define NV_Ith_P(v,i) ( NV_DATA_P(v)[i] )

The NVECTOR_PARALLEL module defines parallel implementations of all
vector operations listed in the section :ref:`NVectors.Ops`.  Their
names are obtained from those that section by appending the suffix
``_Parallel``. The module NVECTOR_PARALLEL provides the following
additional user-callable routines: 

* ``N_VNew_Parallel``

  This function creates and allocates memory for a parallel vector.

  .. code-block:: c

     N_Vector N_VNew_Parallel(MPI_Comm comm, long int local_length, 
                              long int global_length);

* ``N_VNewEmpty_Parallel``

  This function creates a new parallel ``N_Vector`` with an empty
  (``NULL``) data array. 
 
  .. code-block:: c

     N_Vector N_VNewEmpty_Parallel(MPI_Comm comm, long int local_length, 
                                   long int global_length); 

* ``N_VMake_Parallel``

  This function creates and allocates memory for a parallel vector
  with user-provided data array. 

  .. code-block:: c

     N_Vector N_VMake_Parallel(MPI_Comm comm, long int local_length,
                               long int global_length, realtype *v_data); 

* ``N_VCloneVectorArray_Parallel``

  This function creates (by cloning) an array of ``count`` parallel vectors.

  .. code-block:: c

     N_Vector *N_VCloneVectorArray_Parallel(int count, N_Vector w);

* ``N_VCloneEmptyVectorArray_Parallel``

  This function creates (by cloning) an array of ``count`` parallel
  vectors, each with an empty (``NULL``) data array. 

  .. code-block:: c

     N_Vector *N_VCloneEmptyVectorArray_Parallel(int count, N_Vector w);

* ``N_VDestroyVectorArray_Parallel``

  This function frees memory allocated for the array of ``count``
  variables of type ``N_Vector`` created with
  ``N_VCloneVectorArray_Parallel`` or with
  ``N_VCloneEmptyVectorArray_Parallel``. 

  .. code-block:: c

     void N_VDestroyVectorArray_Parallel(N_Vector *vs, int count);

* ``N_VPrint_Parallel``

  This function prints the content of a parallel vector to
  ``stdout``. 

  .. code-block:: c

     void N_VPrint_Parallel(N_Vector v);


**Notes**

* When looping over the components of an ``N_Vector`` ``v``, it is
  more efficient to first obtain the local component array via ``v_data
  = NV_DATA_P(v)`` and then access ``v_data[i]`` within the loop than it
  is to use ``NV_Ith_P(v,i)`` within the loop. 
* ``N_VNewEmpty_Parallel``, ``N_VMake_Parallel``, and
  ``N_VCloneEmptyVectorArray_Parallel`` set the field `own_data` to
  ``FALSE``. ``N_VDestroy_Parallel`` and
  ``N_VDestroyVectorArray_Parallel`` will not attempt to free the
  pointer data for any ``N_Vector`` with `own_data` set to
  ``FALSE``. In such a case, it is the user's responsibility to
  deallocate the data pointer. 
* To maximize efficiency, vector operations in the NVECTOR_PARALLEL
  implementation that have more than one ``N_Vector`` argument do not
  check for consistent internal representation of these vectors. It is
  the user's responsibility to ensure that such routines are called
  with ``N_Vector`` arguments that were all created with the same
  internal representations.



.. _NVectors.ARKode:

NVECTOR functions used by ARKode
-----------------------------------------

(to be added)
