..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2014, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _NVectors.OpenMP:

The NVECTOR_OPENMP Module
======================================

The OpenMP implementation of the NVECTOR module provided with SUNDIALS,
NVECTOR_OPENMP, defines the *content* field of ``N_Vector`` to be a structure 
containing the length of the vector, a pointer to the beginning of a contiguous 
data array, and a boolean flag *own_data* which specifies the ownership 
of *data*.  Operations on the vector are threaded using OpenMP, 
the number of threads used is based on the supplied argument in 
the vector constructor.

.. code-block:: c

   struct _N_VectorContent_OpenMP {
     long int length;
     booleantype own_data;
     realtype *data;
     int num_threads;
   };

The following six macros are provided to access the content of an NVECTOR_OPENMP
vector. The suffix ``_OMP`` in the names denotes OpenMP version.


.. c:macro:: NV_CONTENT_OMP(v)

   This macro gives access to the contents of the OpenMP vector
   ``N_Vector`` *v*.
  
   The assignment ``v_cont = NV_CONTENT_OMP(v)`` sets ``v_cont`` to be
   a pointer to the OpenMP ``N_Vector`` content structure.
  
   Implementation:
  
   .. code-block:: c

      #define NV_CONTENT_OMP(v) ( (N_VectorContent_OpenMP)(v->content) ) 


.. c:macro:: NV_OWN_DATA_OMP(v)

   Access the *own_data* component of the OpenMP ``N_Vector`` *v*.

   Implementation:

   .. code-block:: c
 
      #define NV_OWN_DATA_OMP(v) ( NV_CONTENT_OMP(v)->own_data ) 


.. c:macro:: NV_DATA_OMP(v)

   The assignment ``v_data = NV_DATA_OMP(v)`` sets ``v_data`` to be a
   pointer to the first component of the *data* for the ``N_Vector``
   ``v``. 

   Similarly, the assignment ``NV_DATA_OMP(v) = v_data`` sets the component
   array of ``v`` to be ``v_data`` by storing the pointer ``v_data``.

   Implementation:

   .. code-block:: c
 
      #define NV_DATA_OMP(v) ( NV_CONTENT_OMP(v)->data ) 


.. c:macro:: NV_LENGTH_OMP(v)

   Access the *length* component of the OpenMP ``N_Vector`` *v*.

   The assignment ``v_len = NV_LENGTH_OMP(v)`` sets ``v_len`` to be the
   *length* of ``v``. On the other hand, the call ``NV_LENGTH_OMP(v) =
   len_v`` sets the *length* of ``v`` to be ``len_v``. 

   Implementation:

   .. code-block:: c
 
      #define NV_LENGTH_OMP(v) ( NV_CONTENT_OMP(v)->length )


.. c:macro:: NV_NUM_THREADS_OMP(v)

   Access the *num_threads* component of the OpenMP ``N_Vector`` *v*.

   The assignment ``v_threads = NV_NUM_THREADS_OMP(v)`` sets
   ``v_threads`` to be the *num_threads* of ``v``. On the other hand,
   the call ``NV_NUM_THREADS_OMP(v) = num_threads_v`` sets the
   *num_threads* of ``v`` to be ``num_threads_v``.

   Implementation:

   .. code-block:: c
 
      #define NV_NUM_THREADS_OMP(v) ( NV_CONTENT_OMP(v)->num_threads )


.. c:macro:: NV_Ith_OMP(v,i)

   This macro gives access to the individual components of the *data*
   array of an ``N_Vector``, using standard 0-based C indexing. 

   The assignment ``r = NV_Ith_OMP(v,i)`` sets ``r`` to be the value of
   the ``i``-th component of ``v``. 

   The assignment ``NV_Ith_OMP(v,i) = r`` sets the value of the ``i``-th
   component of ``v`` to be ``r``. 

   Here ``i`` ranges from 0 to :math:`n-1` for a vector of length
   :math:`n`. 

   Implementation: 

   .. code-block:: c

      #define NV_Ith_OMP(v,i) ( NV_DATA_OMP(v)[i] )




The NVECTOR_OPENMP module defines OpenMP implementations of all vector
operations listed in the section :ref:`NVectors.Ops`.  Their names are
obtained from those in that section by appending the suffix
``_OpenMP``.  

In addition, the module NVECTOR_OPENMP provides the following
additional user-callable routines:


.. c:function:: N_Vector N_VNew_OpenMP(long int vec_length, int num_threads)

   This function creates and allocates memory for a OpenMP
   ``N_Vector``. Arguments are the vector length and number of threads.


.. c:function:: N_Vector N_VNewEmpty_OpenMP(long int vec_length, int num_threads)

   This function creates a new OpenMP ``N_Vector`` with an empty
   (``NULL``) data array. 


.. c:function:: N_Vector N_VMake_OpenMP(long int vec_length, realtype* v_data, int num_threads)

   This function creates and allocates memory for a OpenMP vector with
   user-provided data array, *v_data*. 


.. c:function:: N_Vector* N_VCloneVectorArray_OpenMP(int count, N_Vector w)

   This function creates (by cloning) an array of *count* OpenMP
   vectors. 


.. c:function:: N_Vector* N_VCloneEmptyVectorArray_OpenMP(int count, N_Vector w)

   This function creates (by cloning) an array of *count* OpenMP
   vectors, each with an empty (```NULL``) data array.


.. c:function:: void N_VDestroyVectorArray_OpenMP(N_Vector* vs, int count)
  
   This function frees memory allocated for the array of *count*
   variables of type ``N_Vector`` created with
   :c:func:`N_VCloneVectorArray_OpenMP()` or with
   :c:func:`N_VCloneEmptyVectorArray_OpenMP()`. 


.. c:function:: void N_VPrint_OpenMP(N_Vector v)

   This function prints the content of a OpenMP vector to ``stdout``.

    

**Notes**

* When looping over the components of an ``N_Vector v``, it is more
  efficient to first obtain the component array via ``v_data =
  NV_DATA_OMP(v)`` and then access ``v_data[i]`` within the loop than it 
  is to use ``NV_Ith_S(v,i)`` within the loop. 

* :c:func:`N_VNewEmpty_OpenMP()`, :c:func:`N_VMake_OpenMP()`, and
  :c:func:`N_VCloneEmptyVectorArray_OpenMP()` set the field *own_data*
  to ``FALSE``.  The functions :c:func:`N_VDestroy_OpenMP()` and
  :c:func:`N_VDestroyVectorArray_OpenMP()` will not attempt to free the
  pointer data for any ``N_Vector`` with *own_data* set to ``FALSE``.
  In such a case, it is the user's responsibility to deallocate the
  data pointer. 

* To maximize efficiency, vector operations in the NVECTOR_OPENMP
  implementation that have more than one ``N_Vector`` argument do not
  check for consistent internal representation of these vectors. It is
  the user's responsibility to ensure that such routines are called
  with ``N_Vector`` arguments that were all created with the same
  internal representations.
