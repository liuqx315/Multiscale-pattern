..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _NVectors.NVSerial:

The NVECTOR_SERIAL Module
======================================

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
  array of an ``N_Vector``, using standard 0-based C indexing. 

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
