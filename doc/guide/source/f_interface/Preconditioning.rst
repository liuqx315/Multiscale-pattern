..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _FInterface.Preconditioning:

Usage of the FARKODE interface to built-in preconditioners
============================================================

The FARKODE interface enables usage of the two built-in
preconditioning modules ARKBANDPRE and ARKBBDPRE.  Details on how
these preconditioners work are provided in the section
:ref:`CInterface.PreconditionerModules`.  In this section, we focus
specifically on the Fortran interface to these modules.



.. _FInterface.BandPre:

Usage of the FARKBP interface to ARKBANDPRE
-----------------------------------------------

The FARKBP interface module is a package of C functions which,
as part of the FARKODE interface module, support the use of the
ARKode solver with the serial or threaded NVector modules
(:ref:`NVectors.NVSerial`, :ref:`NVectors.OpenMP` or
:ref:`NVectors.Pthreads`), and the combination of the ARKBANDPRE
preconditioner module (see the section :ref:`CInterface.BandPre`) with
any of the Krylov iterative linear solvers. 

The two user-callable functions in this package, with the
corresponding ARKode function around which they wrap, are: 

* :f:func:`FARKBPINIT()` interfaces to :c:func:`ARKBandPrecInit()`.

* :f:func:`FARKBPOPT()` interfaces to the ARKBANDPRE optional output
  functions, :c:func:`ARKBandPrecGetWorkSpace()` and
  :c:func:`ARKBandPrecGetNumRhsEvals()`.

As with the rest of the FARKODE routines, the names of the
user-supplied routines are mapped to actual values through a series of
definitions in the header file ``farkbp.h``. 

The following is a summary of the usage of this module.  Steps that
are unchanged from the main program described in the section
:ref:`FInterface.Usage` are *italicized*.

1. *Right-hand side specification*

2. *NVECTOR module initialization*

3. *Problem specification*

4. *Set optional inputs*

5. Linear solver specification 

   First, specify one of the ARKSPILS iterative linear solvers, by
   calling one of :f:func:`FARKSPGMR()`, :f:func:`FARKSPBCG()`, 
   :f:func:`FARKSPTFQMR()`, :f:func:`FARKSPFGMR()` or
   :f:func:`FARKPCG()`.

   Optionally, to specify that SPGMR, SPBCG, SPTFQMR, SPFGMR or PCG
   should use the supplied :f:func:`FARKJTIMES()` routine, the user
   should call :f:func:`FARKSPILSSETJAC()` with FLAG :math:`\ne 0`,
   as described in the section :ref:`FInterface.SpilsUserSupplied`.

   Then, to initialize the ARKBANDPRE preconditioner, call the
   routine :f:func:`FARKBPINIT()`, as follows:


   .. f:subroutine:: FARKBPINIT(NEQ, MU, ML, IER)
   
      Interfaces with the :c:func:`ARKBandPrecInit()`
      function to allocate memory and initialize data associated
      with the ARKBANDPRE preconditioner.
   
      **Arguments:** 
         * *NEQ* (``long int``, input) -- problem size. 
         * *MU* (``long int``, input) -- upper half-bandwidth of the
	   band matrix that is retained as an approximation of the
	   Jacobian. 
         * *ML*  (``long int``, input) -- lower half-bandwidth of the
	   band matrix approximation to the Jacobian.
         * *IER*  (``int``, output) -- return flag  (0 if success, -1
	   if a memory failure). 


6. *Problem solution*

7. ARKBANDPRE optional outputs 

   Optional outputs specific to the SPGMR, SPBCG, SPTFQMR, SPFGMR or
   PCG solver are listed in :ref:`FInterface.SpilsIOUTTable`.  To
   obtain the optional outputs associated with the ARKBANDPRE module,
   the user should call the :f:func:`FARKBPOPT()`, as specified below: 


   .. f:subroutine:: FARKBPOPT(LENRWBP, LENIWBP, NFEBP)
      
      Interfaces with the ARKBANDPRE optional output functions.
         
      **Arguments:** 
         * *LENRWBP* (``long int``, output) -- length of real
	   preconditioner work space (from
	   :c:func:`ARKBandPrecGetWorkSpace()`). 
         * *LENIWBP* (``long int``, output) -- length of integer
	   preconditioner work space, in integer words (from
	   :c:func:`ARKBandPrecGetWorkSpace()`). 
         * *NFEBP* (``long int``, output) -- number of
	   :math:`f_I(t,y)` evaluations (from 
	   :c:func:`ARKBandPrecGetNumRhsEvals()`)  

8. *Additional solution output*

9. *Problem reinitialization*

10. *Memory deallocation* 

    (The memory allocated for the FARKBP module is deallocated
    automatically by :f:func:`FARKFREE()`)




.. _FInterface.BBDPre:

Usage of the FARKBBD interface to ARKBBDPRE
-----------------------------------------------

The FARKBBD interface module is a package of C functions which, as
part of the FARKODE interface module, support the use of the ARKode
solver with the parallel vector module (:ref:`NVectors.NVParallel`),
and the combination of the ARKBBDPRE preconditioner module (see the
section :ref:`CInterface.BBDPre`) with any of the Krylov iterative
linear solvers. 

The user-callable functions in this package, with the corresponding
ARKode and ARKBBDPRE functions, are as follows:

* :f:func:`FARKBBDINIT()` interfaces to :c:func:`ARKBBDPrecInit()`.

* :f:func:`FARKBBDREINIT()` interfaces to :c:func:`ARKBBDPrecReInit()`.

* :f:func:`FARKBBDOPT()` interfaces to the ARKBBDPRE optional output
  functions.

In addition to the functions required for general FARKODE usage, the
user-supplied functions required by this package are listed in the
table below, each with the corresponding interface function which
calls it (and its type within ARKBBDPRE or ARKode).


*Table: FARKBBD function mapping*

.. cssclass:: table-bordered

+--------------------------+------------------------+-----------------------------------+
| FARKBBD routine          | ARKode routine         | ARKode interface                  |
| (FORTRAN, user-supplied) | (C, interface)         | function type                     |
+==========================+========================+===================================+
| :f:func:`FARKJTIMES()`   | FARKJtimes             | :c:func:`ARKSpilsJacTimesVecFn()` |
+--------------------------+------------------------+-----------------------------------+
| :f:func:`FARKGLOCFN()`   | FARKgloc               | :c:func:`ARKLocalFn()`            |
+--------------------------+------------------------+-----------------------------------+
| :f:func:`FARKCOMMFN()`   | FARKcfn                | :c:func:`ARKCommFn()`             |
+--------------------------+------------------------+-----------------------------------+

As with the rest of the FARKODE routines, the names of all
user-supplied routines here are fixed, in order to maximize
portability for the resulting mixed-language program.  Additionally,
based on flags discussed above in the section :ref:`FInterface.Routines`,
the names of the user-supplied routines are mapped to actual values
through a series of definitions in the header file ``farkbbd.h``. 

The following is a summary of the usage of this module. Steps that are
unchanged from the main program described in the section
:ref:`FInterface.Usage` are *italicized*. 

1. *Right-hand side specification*

2. *NVECTOR module initialization*

3. *Problem specification*

4. *Set optional inputs*

5. Linear solver specification 

   First, specify one of the ARKSPILS iterative linear solvers, by
   calling one of :f:func:`FARKSPGMR()`, :f:func:`FARKSPBCG()`, 
   :f:func:`FARKSPTFQMR()`, :f:func:`FARKSPFGMR()` or
   :f:func:`FARKPCG()`.

   Optionally, to specify that SPGMR, SPBCG, SPTFQMR, SPFGMR or PCG
   should use the supplied :f:func:`FARKJTIMES()` routine, the user
   should call :f:func:`FARKSPILSSETJAC()` with FLAG :math:`\ne 0`,
   as described in the section :ref:`FInterface.SpilsUserSupplied`.

   Then, to initialize the ARKBBDPRE preconditioner, call the function
   :f:func:`FARKBBDINIT()`, as described below:


   .. f:subroutine:: FARKBBDINIT(NLOCAL, MUDQ, MLDQ, MU, ML, DQRELY, IER)
      
      Interfaces with the :c:func:`ARKBBDPrecInit()`
      routine to initialize the ARKBBDPRE preconditioning module.
         
      **Arguments:** 
	 * *NLOCAL* (``long int``, input) -- local vector size on this
	   process. 
   	 * *MUDQ* (``long int``, input) -- upper half-bandwidth to be
	   used in the computation of the local Jacobian blocks by
	   difference quotients.  These may be smaller than the 
   	   true half-bandwidths of the Jacobian of the local block
   	   of :math:`g`, when smaller values may provide greater
	   efficiency.
	 * *MLDQ* (``long int``, input) -- lower half-bandwidth to be
	   used in the computation of the local Jacobian blocks by
	   difference quotients.
	 * *MU* (``long int``, input) -- upper half-bandwidth of the
	   band matrix that is retained as an approximation of the
	   local Jacobian block (may be smaller than *MUDQ*).
	 * *ML* (``long int``, input) -- lower half-bandwidth of the
	   band matrix that is retained as an approximation of the
	   local Jacobian block (may be smaller than *MLDQ*). 
	 * *DQRELY* (``realtype``, input) -- relative increment factor
	   in :math:`y` for difference quotients (0.0 indicates to use
	   the default).
         * *IER*  (``int``, output) -- return flag (0 if success, -1
	   if a memory failure).


6. *Problem solution*

7. ARKBBDPRE optional outputs

   Optional outputs specific to the SPGMR, SPBCG, SPTFQMR, SPFGMR or
   PCG solver are listed in :ref:`FInterface.SpilsIOUTTable`.  To
   obtain the optional outputs associated with the ARKBBDPRE module,
   the user should call the :f:func:`FARKBBDOPT()`, as specified below:


   .. f:subroutine:: FARKBBDOPT(LENRWBBD, LENIWBBD, NGEBBD)
      
      Interfaces with the ARKBBDPRE optional output functions.
         
      **Arguments:** 
	 * *LENRWBP* (``long int``, output) -- length of real
	   preconditioner work space on this process (from
	   :c:func:`ARKBBDPrecGetWorkSpace()`). 
         * *LENIWBP* (``long int``, output) -- length of integer
	   preconditioner work space on this process (from
	   :c:func:`ARKBBDPrecGetWorkSpace()`).
         * *NGEBBD* (``long int``, output) -- number of :math:`g(t,y)`
	   evaluations (from :c:func:`ARKBBDPrecGetNumGfnEvals()`) so
	   far.

8. *Additional solution output*

9. Problem reinitialization

   If a sequence of problems of the same size is being solved using
   the same linear solver (SPGMR, SPBCG, SPTFQMR, SPFGMR or PCG) in
   combination with the ARKBBDPRE preconditioner, then the ARKode
   package can be re-initialized for the second and subsequent
   problems by calling :f:func:`FARKREINIT()`, following which a call
   to :f:func:`FARKBBDREINIT()` may or may not be needed. If the input
   arguments are the same, no :f:func:`FARKBBDREINIT()` call is
   needed.

   If there is a change in input arguments other than *MU* or
   *ML*, then the user program should call :f:func:`FARKBBDREINIT()`
   as specified below: 


   .. f:subroutine:: FARKBBDREINIT(NLOCAL, MUDQ, MLDQ, DQRELY, IER)
      
      Interfaces with the
      :c:func:`ARKBBDPrecReInit()` function to reinitialize the
      ARKBBDPRE module.
         
      **Arguments:**  The arguments of the same names have the same
      meanings as in :f:func:`FARKBBDINIT()`.


   However, if the value of *MU* or *ML* is being changed, then a call
   to :f:func:`FARKBBDINIT()` must be made instead. 

   Finally, if there is a change in any of the linear solver inputs,
   then a call to :f:func:`FARKSPGMR()`, :f:func:`FARKSPBCG()`,
   :f:func:`FARKSPTFQMR()`, :f:func:`FARKSPFGMR()` or
   :f:func:`FARKPCG()` must also be made; in this case the linear
   solver memory is reallocated.  


10. Problem resizing

    If a sequence of problems of different sizes (but with similar
    dyanamical time scales) is being solved using the same linear
    solver (SPGMR, SPBCG, SPTFQMR, SPFGMR or PCG) in combination with
    the ARKBBDPRE preconditioner, then the ARKode package can be
    re-initialized for the second and subsequent problems by calling
    :f:func:`FARKRESIZE()`, following which a call to
    :f:func:`FARKBBDINIT()` is required to delete and re-allocate the 
    preconditioner memory of the correct size.


    .. f:subroutine:: FARKBBDREINIT(NLOCAL, MUDQ, MLDQ, DQRELY, IER)
      
       Interfaces with the
       :c:func:`ARKBBDPrecReInit()` function to reinitialize the
       ARKBBDPRE module.
         
       **Arguments:**  The arguments of the same names have the same
       meanings as in :f:func:`FARKBBDINIT()`.


    However, if the value of MU or ML is being changed, then a call to
    :f:func:`FARKBBDINIT()` must be made instead. 

    Finally, if there is a change in any of the linear solver inputs,
    then a call to :f:func:`FARKSPGMR()`, :f:func:`FARKSPBCG()`,
    :f:func:`FARKSPTFQMR()`, :f:func:`FARKSPFGMR()` or
    :f:func:`FARKPCG()` must also be made; in this case the linear
    solver memory is reallocated.  


11. `Memory deallocation` 

    (The memory allocated for the FARKBBD module is deallocated
    automatically by :f:func:`FARKFREE()`).

12. User-supplied routines 

    The following two routines must be supplied for use with the
    ARKBBDPRE module:


    .. f:subroutine:: FARKGLOCFN(NLOC, T, YLOC, GLOC, IPAR, RPAR, IER)
      
       User-supplied routine (of type :c:func:`ARKLocalFn()`) that
       computes a processor-local approximation :math:`g(t,y)` to
       the right-hand side function :math:`f_I(t,y)`.
         
       **Arguments:** 
          * *NLOC* (``long int``, input) -- local problem size. 
          * *T* (``realtype``, input) -- current value of the
	    independent variable. 
	  * *YLOC* (``realtype``, input) -- array containing local
	    dependent state variables. 
	  * *GLOC* (``realtype``, output) -- array containing local
	    dependent state derivatives. 
          * *IPAR* (``long int``, input/output) -- array containing
	    integer user data that was passed to
	    :f:func:`FARKMALLOC()`. 
          * *RPAR* (``realtype``, input/output) -- array containing
	    real user data that was passed to :f:func:`FARKMALLOC()`.
          * *IER* (``int``, output) -- return flag (0 if success, >0
	    if a recoverable error occurred, <0 if an unrecoverable
	    error occurred).


    .. f:subroutine:: FARKCOMMFN(NLOC, T, YLOC, IPAR, RPAR, IER)
      
       User-supplied routine (of type :c:func:`ARKCommFn()`) that
       performs all interprocess communication necessary for the
       executation of the :f:func:`FARKGLOCFN()` function above, using
       the input vector *YLOC*.
         
       **Arguments:** 
          * *NLOC* (``long int``, input) -- local problem size. 
	  * *T* (``realtype``, input) -- current value of the
	    independent variable. 
	  * *YLOC* (``realtype``, input) -- array containing local
	    dependent state variables. 
          * *IPAR* (``long int``, input/output) -- array containing
	    integer user data that was passed to
	    :f:func:`FARKMALLOC()`. 
          * *RPAR* (``realtype``, input/output) -- array containing
	    real user data that was passed to :f:func:`FARKMALLOC()`.
          * *IER* (``int``, output) -- return flag (0 if success, >0
	    if a recoverable error occurred, <0 if an unrecoverable
	    error occurred).

       **Notes:**
       This subroutine must be supplied even if it is not needed, and
       must return *IER = 0*.  



