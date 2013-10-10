..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _NVectors.NVParallel:

The NVECTOR_PARALLEL Module
================================

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
