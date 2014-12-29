..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _FInterface.Usage:

Usage of the FARKODE interface module
==========================================

The usage of FARKODE requires calls to five or more interface
functions, depending on the method options selected, and two or more
user-supplied routines which define the problem to be solved.  These 
function calls and user routines are summarized individually below.
Some details on specific argument options, and the user is referred to
the description of the corresponding C interface ARKode functions for
complete information on the arguments of any given user-callable
interface routine.  The usage of FARKODE for rootfinding and with
preconditioner modules is described in later subsections.

We note that at present, support for non-identity mass matrices,
:math:`M\ne I` is only provided in the C/C++ interface to ARKode,
although support within the Fortran interface is planned for the near
future.

In the instructions below, steps marked [**S**] apply to the serial 
NVECTOR implementation (NVECTOR_SERIAL) only, steps marked [**O**]
apply to the OpenMP NVECTOR implementation (NVECTOR_OPENMP) only,
steps marked [**T**] apply to the Pthreads NVECTOR implementation
(NVECTOR_PTHREADS) only, while those marked with a [**P**] apply to
NVECTOR_PARALLEL.  Some steps will be marked with a combination of the
above, e.g.  [**S**, **O**, **T**].  Steps not marked apply to all
supplied NVECTOR implementations.



.. _FInterface.RHS:

Right-hand side specification
--------------------------------------

The user must in all cases supply the following Fortran routines:

.. f:subroutine:: FARKIFUN(T, Y, YDOT, IPAR, RPAR, IER)
   
   Sets the *YDOT* array to :math:`f_I(t,y)`, the implicit portion of
   the right-hand side of the ODE system, as function of the
   independent variable *T* :math:`=t` and the array of dependent state
   variables *Y* :math:`=y`. 
      
   **Arguments:**
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *Y* (``realtype``, input) -- array containing state variables.
      * *YDOT* (``realtype``, output) -- array containing state derivatives.
      * *IPAR* (``long int``, input) -- array containing integer user
        data that was passed to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user
        data that was passed to :f:func:`FARKMALLOC()`.
      * *IER* (``int``, output) -- return flag (0 success, >0
        recoverable error, <0 unrecoverable error).
   

.. f:subroutine:: FARKEFUN(T, Y, YDOT, IPAR, RPAR, IER)
   
   Sets the *YDOT* array to :math:`f_E(t,y)`, the explicit portion of
   the right-hand side of the ODE system, as function of the
   independent variable *T* :math:`=t` and the array of dependent state
   variables *Y* :math:`=y`. 
      
   **Arguments:**
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *Y* (``realtype``, input) -- array containing state variables.
      * *YDOT* (``realtype``, output) -- array containing state derivatives.
      * *IPAR* (``long int``, input) -- array containing integer user
        data that was passed to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user
        data that was passed to :f:func:`FARKMALLOC()`.
      * *IER* (``int``, output) -- return flag (0 success, >0
        recoverable error, <0 unrecoverable error).

For purely explicit problems, although the routine
:f:func:`FARKIFUN()` must exist, it will never be called, and may
remain empty.  Similarly, for purely implicit problems, 
:f:func:`FARKEFUN()` will never be called and must exist and may
remain empty.



.. _FInterface.NVector:

NVECTOR module initialization
--------------------------------------

[**S**] To initialize the serial NVECTOR module, the user must
call the following function with the argument *KEY* = 4.

.. f:subroutine:: FNVINITS(KEY, NEQ, IER)
   
   Initializes the Fortran interface to the serial NVECTOR module.
      
   **Arguments:** 
      * *KEY* (``int``, input) -- integer flag denoting which solver
	is to be used (1 is CVODE, 2 is IDA, 3 is KINSOL and 4 is ARKode).
      * *NEQ* (``long int``, input) -- size of the ODE system.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


[**O**] To initialize the OpenMP NVECTOR module, the user must
call the following function with the argument *KEY* = 4.

.. f:subroutine:: FNVINITS_OPENMP(KEY, NEQ, NUM_THREADS, IER)
   
   Initializes the Fortran interface to the OpenMP NVECTOR module.
      
   **Arguments:** 
      * *KEY* (``int``, input) -- integer flag denoting which solver
	is to be used (1 is CVODE, 2 is IDA, 3 is KINSOL and 4 is ARKode).
      * *NEQ* (``long int``, input) -- size of the ODE system.
      * *NUM_THREADS* (``int``, input) -- number of threads to use in
	parallelized regions.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


[**T**] To initialize the Pthreads NVECTOR module, the user must
call the following function with the argument *KEY* = 4.

.. f:subroutine:: FNVINITS_PTHREADS(KEY, NEQ, NUM_THREADS, IER)
   
   Initializes the Fortran interface to the Pthreads NVECTOR module.
      
   **Arguments:** 
      * *KEY* (``int``, input) -- integer flag denoting which solver
	is to be used (1 is CVODE, 2 is IDA, 3 is KINSOL and 4 is ARKode).
      * *NEQ* (``long int``, input) -- size of the ODE system.
      * *NUM_THREADS* (``int``, input) -- number of threads to use in
	parallelized regions.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


[**P**] To initialize the parallel NVECTOR module, the user must
call the following function with the argument *KEY* = 4.

.. f:subroutine:: FNVINITP(COMM, KEY, NLOCAL, NGLOBAL, IER)
   
   Initializes the Fortran interface to the parallel NVECTOR module.
      
   **Arguments:** 
      * *COMM* (``int``, input) -- the MPI communicator.
      * *KEY* (``int``, input) -- integer flag denoting which solver is to be
        used (1 is CVODE, 2 is IDA, 3 is KINSOL and 4 is ARKode).
      * *NLOCAL* (``long int``, input) -- local vector size on this processor.
      * *NGLOBAL* (``long int``, input) -- the size of the ODE system,
	and the global size of vectors (the sum of all values of *NLOCAL*).
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).
      
   **Notes:** If the header file ``sundials_config.h`` defines
   ``SUNDIALS_MPI_COMM_F2C`` to be 1 (meaning the MPI implementation 
   used to build SUNDIALS includes the :c:func:`MPI_Comm_f2c()` function),
   then COMM can be any valid MPI communicator.  Otherwise,
   ``MPI_COMM_WORLD`` will be used, so the user can just pass an
   integer value as a placeholder.



.. _FInterface.Problem:

Problem specification
--------------------------------------

To set various problem and solution parameters and allocate internal
memory, the user must call :f:func:`FARKMALLOC()`.


.. f:subroutine:: FARKMALLOC(T0, Y0, IMEX, IATOL, RTOL, ATOL, IOUT, ROUT, IPAR, RPAR, IER)
   
   Initializes the Fortran interface to the ARKode solver, providing
   interfaces to the C routines :c:func:`ARKodeCreate()`,
   :c:func:`ARKodeSetUserData()`, and :c:func:`ARKodeInit()`, as well
   as one of :c:func:`ARKodeSStolerances()` or
   :c:func:`ARKodeSVtolerances()`.
      
   **Arguments:** 
      * *T0* (``realtype``, input) -- initial value of :math:`t`.
      * *Y0* (``realtype``, input) -- array of initial conditions. 
      * *IMEX* (``int``, input) -- flag denoting basic integration
	method: 0 = implicit, 1 = explicit, 2 = ImEx. 
      * *IATOL* (``int``, input) -- type for absolute tolerance input
	*ATOL*: 1 = scalar, 2 = array, 3 = user-supplied function; the
	user must subsequently call :f:func:`FARKEWTSET()` and supply
	a routine :f:func:`FARKEWT()` to compute the error weight vector.
      * *RTOL* (``realtype``, input) -- scalar relative tolerance.
      * *ATOL* (``realtype``, input) -- scalar or array absolute tolerance.
      * *IOUT* (``long int``, input/output) -- array of length 22 for integer optional outputs.
      * *ROUT* (``realtype``, input/output) -- array of length 6 for real optional outputs.
      * *IPAR* (``long int``, input/output) -- array of user integer data, which will be passed
        unmodified to all user-provided routines.
      * *RPAR* (``realtype``, input/output) -- array with user real data, which will be passed
        unmodified to all user-provided routines.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).
      
   **Notes:** Modifications to the user data arrays *IPAR* and *RPAR*
   inside a user-provided routine will be propagated to all
   subsequent calls to such routines. The optional outputs
   associated with the main ARKode integrator are listed in
   :ref:`FInterface.IOUTTable` and :ref:`FInterface.ROUTTable`, in
   the section :ref:`FInterface.OptionalOutputs`. 


As an alternative to providing tolerances in the call to
:f:func:`FARKMALLOC()`, the user may provide a routine to compute the
error weights used in the WRMS norm evaluations.  If supplied, it must
have the following form:

.. f:subroutine:: FARKEWT(Y, EWT, IPAR, RPAR, IER)
   
   It must set the positive components of the error weight
   vector *EWT* for the calculation of the WRMS norm of *Y*.
      
   **Arguments:** 
      * *Y* (``realtype``, input) -- array containing state variables.
      * *EWT* (``realtype``, output) -- array containing the error weight vector.
      * *IPAR* (``long int``, input) -- array containing the integer user data that was passed
        to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing the real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).

   
If the :f:func:`FARKEWT()` routine is provided, then, following the
call to :f:func:`FARKMALLOC()`, the user must call the function
:f:func:`FARKEWTSET()`. 

.. f:subroutine:: FARKEWTSET(FLAG, IER)
 
   Informs FARKODE to use the user-supplied :f:func:`FARKEWT()` function.
      
   **Arguments:** 
      * *FLAG* (``int``, input) -- flag, use "1" to denoting to use :f:func:`FARKEWT()`.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).




.. _FInterface.OptionalInputs:

Setting optional inputs
--------------------------------------

Unlike ARKode's C interface, that provides separate functions for
setting each optional input, FARKODE uses only two functions, that
accept keywords to specify which optional input should be set to the
provided value.  These routines are :f:func:`FARKSETIIN()` and
:f:func:`FARKSETRIN()`, and are further described below. 


.. f:subroutine:: FARKSETIIN(KEY, IVAL, IER)
   
   Specification routine to pass optional integer inputs
   to the :f:func:`FARKODE()` solver.
      
   **Arguments:** 
      * *KEY* (quoted string, input) -- which optional input
        is set (see :ref:`FInterface.IINOptionTable`).
      * *IVAL* (``long int``, input) -- the integer input value to be used.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).



.. _FInterface.IINOptionTable:

Table: Keys for setting FARKODE integer optional inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

=======================  =========================================
Key                      ARKode routine
=======================  =========================================
``ORDER``                :c:func:`ARKodeSetOrder()`
``DENSE_ORDER``          :c:func:`ARKodeSetDenseOrder()`
``LINEAR``               :c:func:`ARKodeSetLinear()`
``NONLINEAR``            :c:func:`ARKodeSetNonlinear()`
``FIXEDPOINT``           :c:func:`ARKodeSetFixedPoint()`
``NEWTON``               :c:func:`ARKodeSetNewton()`
``EXPLICIT``             :c:func:`ARKodeSetExplicit()`
``IMPLICIT``             :c:func:`ARKodeSetImplicit()`
``IMEX``                 :c:func:`ARKodeSetImEx()`
``IRK_TABLE_NUM``        :c:func:`ARKodeSetIRKTableNum()`
``ERK_TABLE_NUM``        :c:func:`ARKodeSetERKTableNum()`
``ARK_TABLE_NUM`` *(a)*  :c:func:`ARKodeSetARKTableNum()`      
``MAX_NSTEPS``           :c:func:`ARKodeSetMaxNumSteps()`
``HNIL_WARNS``           :c:func:`ARKodeSetMaxHnilWarns()`
``PREDICT_METHOD``       :c:func:`ARKodeSetPredictorMethod()`
``MAX_ERRFAIL``          :c:func:`ARKodeSetMaxErrTestFails()`
``MAX_CONVFAIL``         :c:func:`ARKodeSetMaxConvFails()`
``MAX_NITERS``           :c:func:`ARKodeSetMaxNonlinIters()`
``ADAPT_SMALL_NEF``      :c:func:`ARKodeSetSmallNumEFails()`
``LSETUP_MSBP``          :c:func:`ARKodeSetMaxStepsBetweenLSet()`
=======================  =========================================

*(a)* When setting ``ARK_TABLE_NUM``, pass in *IVAL* as an array of
length 2, specifying the IRK table number first, then the ERK table
number.
      

.. f:subroutine:: FARKSETRIN(KEY, RVAL, IER)
  
   Specification routine to pass optional real inputs
   to the :f:func:`FARKODE()` solver.
      
   **Arguments:** 
      * *KEY* (quoted string, input) -- which optional input
        is set (see :ref:`FInterface.RINOptionTable`).
      * *RVAL* (``realtype``, input) -- the real input value to be used.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).



.. _FInterface.RINOptionTable:

Table: Keys for setting FARKODE real optional inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

=================  =========================================
Key                ARKode routine
=================  =========================================
``INIT_STEP``      :c:func:`ARKodeSetInitStep()`
``MAX_STEP``       :c:func:`ARKodeSetMaxStep()`
``MIN_STEP``       :c:func:`ARKodeSetMinStep()`
``STOP_TIME``      :c:func:`ARKodeSetStopTime()`
``NLCONV_COEF``    :c:func:`ARKodeSetNonlinConvCoef()`
``ADAPT_CFL``      :c:func:`ARKodeSetCFLFraction()`
``ADAPT_SAFETY``   :c:func:`ARKodeSetSafetyFactor()`
``ADAPT_BIAS``     :c:func:`ARKodeSetErrorBias()`
``ADAPT_GROWTH``   :c:func:`ARKodeSetMaxGrowth()`
``ADAPT_ETAMX1``   :c:func:`ARKodeSetMaxFirstGrowth()`
``ADAPT_BOUNDS``   :c:func:`ARKodeSetFixedStepBounds()`
``ADAPT_ETAMXF``   :c:func:`ARKodeSetMaxEFailGrowth()`
``ADAPT_ETACF``    :c:func:`ARKodeSetMaxCFailGrowth()`
``NONLIN_CRDOWN``  :c:func:`ARKodeSetNonlinCRDown()`
``NONLIN_RDIV``    :c:func:`ARKodeSetNonlinRDiv()`
``LSETUP_DGMAX``   :c:func:`ARKodeSetDeltaGammaMax()`
``FIXED_STEP``     :c:func:`ARKodeSetFixedStep()`
=================  =========================================




If a user wishes to reset all of the options to their default values,
they may call the routine :f:func:`FARKSETDEFAULTS()`. 

.. f:subroutine:: FARKSETDEFAULTS(IER)
   
   Specification routine to reset all FARKODE optional
   inputs to their default values.
      
   **Arguments:** 
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).
   



Optional advanced FARKODE inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FARKODE supplies additional routines to specify optional advanced
inputs to the :c:func:`ARKode()` solver.  These are summarized below,
and the user is referred to their C routine counterparts for more
complete information. 



.. f:subroutine:: FARKSETERKTABLE(S, Q, P, C, A, B, BEMBED, IER)
   
   Interface to the routine :c:func:`ARKodeSetERKTable()`.
      
   **Arguments:** 
      * *S* (``int``, input) -- number of stages in the table.
      * *Q* (``int``, input) -- global order of accuracy of the method.
      * *P* (``int``, input) -- global order of accuracy of the embedding.
      * *C* (``realtype``, input) -- array of length *S* containing the stage times.
      * *A* (``realtype``, input) -- array of length *S*S* containing the ERK coefficients
        (stored in row-major, "C", order).
      * *B* (``realtype``, input) -- array of length *S* containing the solution coefficients.
      * *BEMBED* (``realtype``, input) -- array of length *S* containing the embedding
        coefficients.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


.. f:subroutine:: FARKSETIRKTABLE(S, Q, P, C, A, B, BEMBED, IER)
   
   Interface to the routine :c:func:`ARKodeSetIRKTable()`.
      
   **Arguments:** 
      * *S* (``int``, input) -- number of stages in the table.
      * *Q* (``int``, input) -- global order of accuracy of the method.
      * *P* (``int``, input) -- global order of accuracy of the embedding.
      * *C* (``realtype``, input) -- array of length *S* containing the stage times.
      * *A* (``realtype``, input) -- array of length *S*S* containing the IRK coefficients
        (stored in row-major, "C", order).
      * *B* (``realtype``, input) -- array of length *S* containing the solution coefficients.
      * *BEMBED* (``realtype``, input) -- array of length *S* containing the embedding
        coefficients.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


.. f:subroutine:: FARKSETARKTABLES(S, Q, P, C, AI, AE, B, BEMBED, IER)
   
   Interface to the routine :c:func:`ARKodeSetARKTables()`.
   
   **Arguments:** 
      * *S* (``int``, input) -- number of stages in the table.
      * *Q* (``int``, input) -- global order of accuracy of the method.
      * *P* (``int``, input) -- global order of accuracy of the embedding.
      * *C* (``realtype``, input) -- array of length *S* containing the stage times.
      * *AI* (``realtype``, input) -- array of length *S*S* containing the IRK coefficients
        (stored in row-major, "C", order) 
      * *AE* (``realtype``, input) -- array of length *S*S* containing the ERK coefficients
        (stored in row-major, "C", order) 
      * *B* (``realtype``, input) -- array of length *S* containing the solution coefficients 
      * *BEMBED* (``realtype``, input) -- array of length *S* containing the embedding
        coefficients 
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure) 
   


Additionally, a user may set the accuracy-based step size adaptivity
strategy (and it's associated parameters) through a call to
:f:func:`FARKSETADAPTIVITYMETHOD()`, as described below. 

.. f:subroutine:: FARKSETADAPTIVITYMETHOD(IMETHOD, IDEFAULT, IPQ, PARAMS, IER)
   
   Specification routine to set the step size adaptivity strategy and
   parameters within the :f:func:`FARKODE()` solver.  Interfaces with
   the C routine :c:func:`ARKodeSetAdaptivityMethod()`.
      
   **Arguments:** 
      * *IMETHOD* (``int``, input) -- choice of adaptivity method.
      * *IDEFAULT* (``int``, input) -- flag denoting whether to use
	default parameters (1) or that customized parameters will be
	supplied (1).
      * *IPQ* (``int``, input) -- flag denoting whether to use
	the embedding order of accuracy (0) or the method order of
	accuracy (1) within step adaptivity algorithm.
      * *PARAMS* (``realtype``, input) -- array of 3 parameters to be
	used within the adaptivity strategy.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


Lastly, the user may provide functions to aid/replace those within
ARKode for handling adaptive error control and explicit stability.
The former of these is designed for advanced users who wish to
investigate custom step adaptivity approaches as opposed to using any
of those built-in to ARKode.  In ARKode's C/C++ interface, this would be
provided by a function of type :c:func:`ARKAdaptFn()`; in the Fortran
interface this is provided through the user-supplied function:

.. f:subroutine:: FARKADAPT(Y, T, H1, H2, H3, E1, E2, E3, Q, P, HNEW, IPAR, RPAR, IER)
   
   It must set the new step size *HNEW* based on the three previous
   steps (*H1*, *H2*, *H3*) and the three previous error estimates
   (*E1*, *E2*, *E3*). 
      
   **Arguments:** 
      * *Y* (``realtype``, input) -- array containing state variables.
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *H1* (``realtype``, input) -- current step size.
      * *H2* (``realtype``, input) -- previous step size.
      * *H3* (``realtype``, input) -- previous-previous step size.
      * *E1* (``realtype``, input) -- estimated temporal error in current step.
      * *E2* (``realtype``, input) -- estimated temporal error in previous step.
      * *E3* (``realtype``, input) -- estimated temporal error in previous-previous step.
      * *Q* (``int``, input) -- global order of accuracy for RK method.
      * *P* (``int``, input) -- global order of accuracy for RK embedding.
      * *HNEW* (``realtype``, output) -- array containing the error weight vector.
      * *IPAR* (``long int``, input) -- array containing the integer
	user data that was passed to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing the real user
	data that was passed to :f:func:`FARKMALLOC()`.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).


This routine is enabled by a call to the activation routine:

.. f:subroutine:: FARKADAPTSET(FLAG, IER)
   
   Informs FARKODE to use the user-supplied :f:func:`FARKADAPT()` function.
      
   **Arguments:** 
      * *FLAG* (``int``, input) -- flag, use "1" to denoting to use
	:f:func:`FARKADAPT()`, or use "0" to denote a return to the
        default adaptivity strategy.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne
	0` failure).

   Note: The call to :f:func:`FARKADAPTSET()` must occur *after* the call
   to :f:func:`FARKMALLOC()`.


Similarly, if either an explicit or mixed implicit-explicit
integration method is to be employed, the user may specify a function
to provide the maximum explicitly-stable step for their problem.
Again, in the C/C++ interface this would be a function of type
:c:func:`ARKExpStabFn()`, while in ARKode's Fortran interface this
must be given through the user-supplied function:

.. f:subroutine:: FARKEXPSTAB(Y, T, HSTAB, IPAR, RPAR, IER)
   
   It must set the maximum explicitly-stable step size, *HSTAB*, based
   on the current solution, *Y*.
      
   **Arguments:** 
      * *Y* (``realtype``, input) -- array containing state variables.
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *HSTAB* (``realtype``, output) -- maximum explicitly-stable step size.
      * *IPAR* (``long int``, input) -- array containing the integer user data that was passed
        to :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing the real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).
 

This routine is enabled by a call to the activation routine:

.. f:subroutine:: FARKEXPSTABSET(FLAG, IER)
   
   Informs FARKODE to use the user-supplied :f:func:`FARKEXPSTAB()` function.
      
   **Arguments:** 
      * *FLAG* (``int``, input) -- flag, use "1" to denoting to use
	:f:func:`FARKEXPSTAB()`, or use "0" to denote a return to the 
        default error-based stability strategy.
      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne
	0` failure).

   Note: The call to :f:func:`FARKEXPSTABSET()` must occur *after* the call
   to :f:func:`FARKMALLOC()`.



   
.. _FInterface.LinearSolver:

Linear solver specification
---------------------------------

In the case of using either an implicit or ImEx method, the solution
of each Runge-Kutta stage may involve the solution of linear systems
related to the Jacobian :math:`J = \frac{\partial f_I}{\partial y}` of
the implicit portion of the ODE system. ARKode presently includes
seven choices for the treatment of these systems, and the user of
FARKODE must call a routine with a specific name to make the
desired choice. 


[**S**, **O**, **T**] Dense treatment of the linear system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use the direct dense linear solver based on the internal SUNDIALS
implementation, the user must call the :f:func:`FARKDENSE()` routine: 


.. f:subroutine:: FARKDENSE(NEQ, IER)
   
   Interfaces with the :c:func:`ARKDense()` function to
   specify use of the dense direct linear solver.
      
   **Arguments:** 
      * *NEQ* (``long int``, input) -- size of the ODE system.
      * *IER* (``int``, output) -- return flag (0 if success, -1 if a memory allocation
        error occurred, -2 for an illegal input).


Alteratively, to use the LAPACK-based direct dense linear solver, a
user must call the similar :f:func:`FARKLAPACKDENSE()` routine:

.. f:subroutine:: FARKLAPACKDENSE(NEQ, IER)
   
   Interfaces with the :c:func:`ARKLapackDense()` function
   to specify use of the LAPACK the dense direct linear solver.
      
   **Arguments:** 
      * *NEQ* (``int``, input) -- size of the ODE system.
      * *IER* (``int``, output) -- return flag (0 if success, -1 if a memory allocation
        error occurred, -2 for an illegal input).


As an option when using either of these dense linear solvers, the user
may supply a routine that computes a dense approximation of the system
Jacobian :math:`J = \frac{\partial f_I}{\partial y}`.  If supplied, it
must have the following form:


.. f:subroutine:: FARKDJAC(NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)
   
   Interface to provide a user-supplied dense Jacobian approximation
   function (of type :c:func:`ARKDlsDenseJacFn()`), to be used by the
   :f:func:`FARKDENSE()` solver. 
      
   **Arguments:** 
      * *NEQ* (``long int``, input) -- size of the ODE system.
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *Y* (``realtype``, input) -- array containing values of the dependent state variables.
      * *FY* (``realtype``, input) -- array containing values of the dependent state derivatives.
      * *DJAC* (``realtype`` of size (NEQ,NEQ), output) -- 2D array containing the Jacobian entries.
      * *H* (``realtype``, input) -- current step size.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *WK1*, *WK2*, *WK3*  (``realtype``, input) -- array containing temporary workspace
        of same size as *Y*.
      * *IER* (``int``, output) -- return flag (0 if success, >0 if a recoverable error
        occurred, <0 if an unrecoverable error occurred).
      
   **Notes:** Typically this routine will use only *NEQ*, *T*, *Y*, and
   *DJAC*. It must compute the Jacobian and store it column-wise in *DJAC*. 
  

   
If the above routine uses difference quotient approximations, it may
need to access the error weight array *EWT* in the calculation of
suitable increments. The array *EWT* can be obtained by calling
:f:func:`FARKGETERRWEIGHTS()` using one of the work arrays as
temporary storage for *EWT*. It may also need the unit roundoff, which
can be obtained as the optional output *ROUT(6)*, passed from the
calling program to this routine using either *RPAR* or a common block. 

If the :f:func:`FARKDJAC()` routine is provided, then, following the
call to :f:func:`FARKDENSE()` or :f:func:`FARKLAPACKDENSE()`, the user
must call the routine :f:func:`FARKDENSESETJAC()`:


.. f:subroutine:: FARKDENSESETJAC(FLAG, IER)
   
   Interface to the :c:func:`ARKDlsSetDenseJacFn()` function, specifying
   to use the user-supplied routine :f:func:`FARKDJAC()` for the
   Jacobian approximation. 
      
   **Arguments:** 
      * *FLAG* (``int``, input) -- any nonzero value specifies to use
	:f:func:`FARKDJAC()`. 
      * *IER* (``int``, output) -- return flag (0 if success,
	:math:`\ne 0` if an error occurred).
   

   


[**S**, **O**, **T**] Band treatment of the linear system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use the direct band linear solver that is based on the internal
SUNDIALS implementation, the user must call the :f:func:`FARKBAND()`
routine.


.. f:subroutine:: FARKBAND(NEQ, MU, ML, IER)
   
   Interfaces with the :c:func:`ARKBand()` function to
   specify use of the dense banded linear solver.
      
   **Arguments:** 
      * *NEQ* (``long int``, input) -- size of the ODE system.
      * *MU* (``long int``, input) -- upper half-bandwidth.
      * *ML* (``long int``, input) -- lower half-bandwidth.
      * *IER* (``int``, output) -- return flag (0 if success, -1 if a memory allocation
        error occurred, -2 for an illegal input).


Alteratively, to use the LAPACK-based direct banded linear solver, a
user must call the similar :f:func:`FARKLAPACKBAND()` routine:


.. f:subroutine:: FARKLAPACKBAND(NEQ, MU, ML, IER)
   
   Interfaces with the :c:func:`ARKLapackBand()` function
   to specify use of the dense banded linear solver.
      
   **Arguments:** 
      * *NEQ* (``int``, input) -- size of the ODE system.
      * *MU* (``int``, input) -- upper half-bandwidth.
      * *ML* (``int``, input) -- lower half-bandwidth.
      * *IER* (``int``, output) -- return flag (0 if success, -1 if a memory allocation
        error occurred, -2 for an illegal input).
   

   
As an option when using either of these banded linear solvers, the user
may supply a routine that computes a banded approximation of the
linear system Jacobian :math:`J = \frac{\partial f_I}{\partial y}`. If
supplied, it must have the following form:

.. f:subroutine:: FARKBJAC(NEQ, MU, ML, MDIM, T, Y, FY, BJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)
   
   Interface to provide a user-supplied band Jacobian approximation
   function (of type :c:func:`ARKDlsBandJacFn()`), to be used by the
   :f:func:`FARKBAND()` solver. 
     
   **Arguments:** 
      * *NEQ* (``long int``, input) -- size of the ODE system.
      * *MU*   (``long int``, input) -- upper half-bandwidth.
      * *ML*   (``long int``, input) -- lower half-bandwidth.
      * *MDIM* (``long int``, input) -- leading dimension of *BJAC* array.
      * *T*    (``realtype``, input) -- current value of the independent variable.
      * *Y*    (``realtype``, input) -- array containing dependent state variables.
      * *FY*   (``realtype``, input) -- array containing dependent state derivatives.
      * *BJAC* (``realtype`` of size *(MDIM,NEQ)*, output) -- 2D array
	containing the Jacobian entries. 
      * *H*    (``realtype``, input) -- current step size.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *WK1*, *WK2*, *WK3*  (``realtype``, input) -- array containing temporary workspace
        of same size as *Y*.
      * *IER* (``int``, output) -- return flag (0 if success, >0 if a recoverable error
        occurred, <0 if an unrecoverable error occurred).
      
   **Notes:**
   Typically this routine will use only *NEQ*, *MU*, *ML*, *T*, *Y*, and
   *BJAC*. It must load the *MDIM* by *N* array *BJAC* with the Jacobian
   matrix at the current :math:`(t,y)` in band form.  Store in
   *BJAC(k,j)* the Jacobian element :math:`J_{i,j}` with 
   *k = i - j + MU + 1* (or *k = 1, ..., ML+MU+1*) and *j = 1, ..., N*. 


If the above routine uses difference quotient approximations, it may
need to use the error weight array *EWT* in the calculation of
suitable increments. The array *EWT* can be obtained by calling
:f:func:`FARKGETERRWEIGHTS()` using one of the work 
arrays as temporary storage for *EWT*. It may also need the unit
roundoff, which can be obtained as the optional output *ROUT(6)*,
passed from the calling program to this routine using either *RPAR*
or a common block. 

If the :f:func:`FARKBJAC()` routine is provided, then, following the
call to either :f:func:`FARKBAND()` or :f:func:`FARKLAPACKBAND()`, the
user must call the routine :f:func:`FARKBANDSETJAC()`. 


.. f:subroutine:: FARKBANDSETJAC(FLAG, IER)
   
   Interface to the :c:func:`ARKDlsSetBandJacFn()` function, specifying
   to use the user-supplied routine :f:func:`FARKBJAC()` for the
   Jacobian approximation. 
      
   **Arguments:** 
      * *FLAG* (``int``, input) -- any nonzero value specifies to use
        :f:func:`FARKBJAC()`.
      * *IER* (``int``, output) -- return flag (0 if success, 
	:math:`\ne 0` if an error occurred).




[**S**, **O**, **T**] Sparse treatment of the linear system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use the sparse direct linear solver interface to the KLU library,
the user must call the :f:func:`FARKKLU()` routine:  


.. f:subroutine:: FARKKLU(NEQ, NNZ, IER)
   
   Interfaces with the :c:func:`ARKKLU()` function to
   specify use of the sparse direct linear solver.
      
   **Arguments:** 
      * *NEQ* (``int``, input) -- size of the ODE system.
      * *NNZ* (``int``, input) -- maximum number of nonzeros in
	the sparse Jacobian.
      * *IER* (``int``, output) -- return flag (0 if success, -1 if a
	memory allocation error occurred, -2 for an illegal input).


Alteratively, to use the SuperLU_MT-based threaded sparse direct
linear solver, a user must call the similar :f:func:`FARKSUPERLUMT()`
routine:

.. f:subroutine:: FARKSUPERLUMT(NTHREADS, NEQ, NNZ, IER)
   
   Interfaces with the :c:func:`ARKSuperLUMT()` function
   to specify use of the SuperLU_MT threaded sparse direct linear solver.
      
   **Arguments:** 
      * *NTHREADS* (``int``, input) -- number of threads to use in
	factorization and solution of the Jacobian systems.
      * *NEQ* (``int``, input) -- size of the ODE system.
      * *NNZ* (``int``, input) -- maximum number of nonzeros in
	the sparse Jacobian.
      * *IER* (``int``, output) -- return flag (0 if success, -1 if a
	memory allocation error occurred, -2 for an illegal input).


When using either of these sparse direct linear solvers, the user must
supply a routine that computes a compressed-sparse-column
approximation of the system Jacobian :math:`J = \frac{\partial
f_I}{\partial y}`, having the following form:


.. f:subroutine:: FARKSPJAC(T, Y, FY, N, NNZ, JDATA, JRVALS, JCPTRS, H, IPAR, RPAR, WK1, WK2, WK3, IER)
   
   Interface to provide a user-supplied sparse Jacobian approximation
   function (of type :c:func:`ARKSlsSparseJacFn()`), to be used by the
   :f:func:`FARKKLU()` or :f:func:`FARKSUPERLUMT()` solver. 
      
   **Arguments:** 
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *Y* (``realtype``, input) -- array containing values of the dependent state variables.
      * *FY* (``realtype``, input) -- array containing values of the dependent state derivatives.
      * *N* (``int``, input) -- number of matrix rows in Jacobian.
      * *NNZ* (``int``, input) -- allocated length of nonzero storage in Jacobian.
      * *JDATA* (``realtype`` of size NNZ, output) -- nonzero values in Jacobian.
      * *JRVALS* (``int`` of size NNZ, output) -- row indices for each
	nonzero Jacobian entry.
      * *JCPTRS* (``int`` of size N+1, output) -- indices of where
	each column's nonzeros begin in data array; last entry points
	just past end of data values.
      * *H* (``realtype``, input) -- current step size.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *WK1*, *WK2*, *WK3*  (``realtype``, input) -- array containing temporary workspace
        of same size as *Y*.
      * *IER* (``int``, output) -- return flag (0 if success, >0 if a recoverable error
        occurred, <0 if an unrecoverable error occurred).
      
   **Notes:** Due to the format of both the KLU and SuperLU_MT
   solvers, the number of matrix rows, number of matrix nonzeros, and
   row index array are all of type ``int`` and not ``long int``.
  

If the above routine uses difference quotient approximations to
compute the nonzero entries, it may need to access the error weight
array *EWT* in the calculation of suitable increments. The array *EWT*
can be obtained by calling :f:func:`FARKGETERRWEIGHTS()` using one of
the work arrays as temporary storage for *EWT*.  It may also need the
unit roundoff, which can be obtained as the optional output *ROUT(6)*,
passed from the calling program to this routine using either *RPAR* or
a common block.


   


SPGMR treatment of the linear systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the Scaled Preconditioned GMRES solution of the linear systems,
the user must call the :f:func:`FARKSPGMR()` routine:


.. f:subroutine:: FARKSPGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)
   
   Interfaces with the :c:func:`ARKSpgmr()` and ARKSpilsSet* routines
   to specify use of the SPGMR iterative linear solver.
      
   **Arguments:** 
      * *IPRETYPE* (``int``, input) -- preconditioner type: 0 = none,
	1 = left only, 2 = right only, 3 = both sides.
      * *IGSTYPE* (``int``, input) -- Gram-schmidt orthogonalization
	process: 1 = modified G-S, 2 = classical G-S.
      * *MAXL* (``int``; input) -- maximum Krylov subspace dimension
	(0 for default). 
      * *DELT* (``realtype``, input) -- linear convergence tolerance
	factor (0.0 for default). 
      * *IER* (``int``, output) -- return flag (0 if success, -1 if a
	memory allocation error occurred, -2 for an illegal input).


For descriptions of the optional user-supplied routines for use with
:f:func:`FARKSPGMR()` see the section :ref:`FInterface.SpilsUserSupplied`.



SPBCG treatment of the linear systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the Scaled Preconditioned Bi-CGStab solution of the linear systems,
the user must call the :f:func:`FARKSPBCG()` routine:


.. f:subroutine:: FARKSPBCG(IPRETYPE, MAXL, DELT, IER)
   
   Interfaces with the :c:func:`ARKSpbcg()` and
   ARKSpilsSet* routines to specify use of the SPBCG iterative
   linear solver.
      
   **Arguments:**  The arguments are the same as those with the
   same names for :f:func:`FARKSPGMR()`. 


For descriptions of the optional user-supplied routines for use with
:f:func:`FARKSPBCG()` see the section :ref:`FInterface.SpilsUserSupplied`.





SPTFQMR treatment of the linear systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the Scaled Preconditioned TFQMR solution of the linear systems,
the user must call the :f:func:`FARKSPTFQMR()` routine:


.. f:subroutine:: FARKSPTFQMR(IPRETYPE, MAXL, DELT, IER)
   
   Interfaces with the :c:func:`ARKSptfqmr()` and
   ARKSpilsSet* routines to specify use of the SPTFQMR iterative
   linear solver.
      
   **Arguments:**  The arguments are the same as those with the same names
   for :f:func:`FARKSPGMR()`.


For descriptions of the optional user-supplied routines for use with
:f:func:`FARKSPTFQMR()` see the next section.



SPFGMR treatment of the linear systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the Scaled Preconditioned Flexible Generalized Minimum Residual
solution of the linear systems, the user must call the
:f:func:`FARKSPFGMR()` routine:


.. f:subroutine:: FARKSPFGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)
   
   Interfaces with the :c:func:`ARKSpfgmr()` and
   ARKSpilsSet* routines to specify use of the SPFGMR iterative
   linear solver.
      
   **Arguments:**  The arguments are the same as those for
   :f:func:`FARKSPGMR()`.


For descriptions of the optional user-supplied routines for use with
:f:func:`FARKSPFGMR()` see the section :ref:`FInterface.SpilsUserSupplied`.





PCG treatment of the linear systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the Preconditioned Conjugate Gradient solution of symmetric linear
systems, the user must call the :f:func:`FARKPCG()` routine:


.. f:subroutine:: FARKPCG(IPRETYPE, MAXL, DELT, IER)
 
   Interfaces with the :c:func:`ARKPcg()` and
   ARKSpilsSet* routines to specify use of the PCG iterative
   linear solver.
      
   **Arguments:**  The arguments are the same as those with the
   same names for :f:func:`FARKSPGMR()`. 


For descriptions of the optional user-supplied routines for use with
:f:func:`FARKPCG()` see the section :ref:`FInterface.SpilsUserSupplied`.





.. _FInterface.SpilsUserSupplied:

User-supplied routines for SPGMR/SPBCG/SPTFQMR/SPFGMR/PCG
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With treatment of the linear systems by any of the Krylov iterative
solvers, there are three optional user-supplied routines --
:f:func:`FARKJTIMES()`, :f:func:`FARKPSET()` and :f:func:`FARKPSOL()`.
The specifications of these functions are given below.

The first of these optional routines when using a Krylov iterative
solver is a routine to compute the product of the system Jacobian
:math:`J = \frac{\partial f_I}{\partial y}` and a given vector
:math:`v`.  If supplied, it must have the following form: 


.. f:subroutine:: FARKJTIMES(V, FJV, T, Y, FY, H, IPAR, RPAR, WORK, IER)
   
   Interface to provide a user-supplied Jacobian-times-vector product
   approximation function (corresponding to a C interface routine of
   type :c:func:`ARKSpilsJacTimesVecFn()`), to be used by one of the
   Krylov iterative linear solvers.
      
   **Arguments:** 
      * *V*    (``realtype``, input) -- array containing the vector to multiply.
      * *FJV*  (``realtype``, output) -- array containing resulting product vector.
      * *T*    (``realtype``, input) -- current value of the independent variable.
      * *Y*    (``realtype``, input) -- array containing dependent state variables.
      * *FY*   (``realtype``, input) -- array containing dependent state derivatives.
      * *H*    (``realtype``, input) -- current step size.
      * *IPAR* (``long int``, input) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *WORK* (``realtype``, input) -- array containing temporary workspace of same size as
        *Y*.
      * *IER*  (``int``, output) -- return flag  (0 if success, :math:`\ne 0` if an error).
         
   **Notes:**
   Typically this routine will use only *NEQ*, *T*, *Y*, *V*, and
   *FJV*.  It must compute the product vector :math:`Jv`, where
   :math:`v` is given in *V*, and the product is stored in *FJV*. 
   

If this routine has been supplied by the user, then, following the
call to :f:func:`FARKSPGMR()`, :f:func:`FARKSPBCG()`,
:f:func:`FARKSPTFQMR()`, :f:func:`FARKSPFGMR()` or
:f:func:`FARKPCG()`, the user must call the routine
:f:func:`FARKSPILSSETJAC()` with *FLAG* :math:`\ne 0` to specify use
of the user-supplied Jacobian-times-vector function:


.. f:subroutine:: FARKSPILSSETJAC(FLAG, IER)
   
   Interface to the function :c:func:`ARKSpilsSetJacTimesVecFn()` to
   specify use of the user-supplied Jacobian-times-vector function
   :f:func:`FARKJTIMES()`. 
      
   **Arguments:** 
      * *FLAG* (``int``, input) -- flag denoting to use
	:f:func:`FARKJTIMES()` routine. 
      * *IER*  (``int``, output) -- return flag  (0 if success,
	:math:`\ne 0` if an error).


If preconditioning is to be performed during the Krylov solver 
(i.e. the solver was set up with *IPRETYPE* :math:`\ne 0`), then the 
user must also call the routine :f:func:`FARKSPILSSETPREC()` with
*FLAG* :math:`\ne 0`:


.. f:subroutine:: FARKSPILSSETPREC(FLAG, IER)
   
   Interface to the function :c:func:`ARKSpilsSetPreconditioner()` to
   specify use of the user-supplied preconditioner setup and solve
   functions, :f:func:`FARKPSET()` and :f:func:`FARKPSOL()`, respectively.
      
   **Arguments:** 
      * *FLAG* (``int``, input) -- flag denoting use of user-supplied
        preconditioning routines.
      * *IER*  (``int``, output) -- return flag  (0 if success,
	:math:`\ne 0` if an error).
         

In addition, the user must provide the following two routines to
implement the preconditioner setup and solve functions to be used
within the solve.


.. f:subroutine:: FARKPSET(T,Y,FY,JOK,JCUR,GAMMA,H,IPAR,RPAR,V1,V2,V3,IER)
   
   User-supplied preconditioner setup routine (of type
   :c:func:`ARKSpilsPrecSetupFn()`). 
      
   **Arguments:** 
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *Y* (``realtype``, input) -- current dependent state variable array.
      * *FY* (``realtype``, input) -- current dependent state variable derivative array.
      * *JOK* (``int``, input) -- flag indicating whether Jacobian-related data needs to be 
        recomputed: 0 = recompute, 1 = reuse with the current value of *GAMMA*.
      * *JCUR* (``realtype``, output) -- return flag to denote if
	Jacobian data was recomputed (1=yes, 0=no).
      * *GAMMA* (``realtype``, input) -- Jacobian scaling factor.
      * *H* (``realtype``, input) -- current step size.
      * *IPAR* (``long int``, input/output) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input/output) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *V1*, *V2*, *V3* (``realtype``, input) -- arrays containing temporary workspace of
        same size as *Y*. 
      * *IER*  (``int``, output) -- return flag  (0 if success, >0 if a recoverable
        failure, <0 if a non-recoverable failure).
      
   **Notes:**
   This routine must set up the preconditioner :math:`P` to be used in
   the subsequent call to :f:func:`FARKPSOL()`.  The preconditioner (or
   the product of the left and right preconditioners if using both)
   should be an approximation to the matrix  :math:`M - \gamma J`,
   where :math:`M` is the system mass matrix, :math:`\gamma` is the
   input *GAMMA*, and :math:`J = \frac{\partial f_I}{\partial y}`. 
   
   
.. f:subroutine:: FARKPSOL(T,Y,FY,R,Z,GAMMA,DELTA,LR,IPAR,RPAR,VT,IER)
   
   User-supplied preconditioner solve routine (of type
   :c:func:`ARKSpilsPrecSolveFn()`).
      
   **Arguments:** 
      * *T* (``realtype``, input) -- current value of the independent variable.
      * *Y* (``realtype``, input) -- current dependent state variable array.
      * *FY* (``realtype``, input) -- current dependent state variable derivative array.
      * *R* (``realtype``, input) -- right-hand side array.
      * *Z* (``realtype``, output) -- solution array.
      * *GAMMA* (``realtype``, input) -- Jacobian scaling factor.
      * *DELTA* (``realtype``, input) -- desired residual tolerance.
      * *LR* (``int``, input) -- flag denoting to solve the right or left preconditioner
        system: 1 = left preconditioner, 2 = right preconditioner.
      * *IPAR* (``long int``, input/output) -- array containing integer user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *RPAR* (``realtype``, input/output) -- array containing real user data that was passed to
        :f:func:`FARKMALLOC()`.
      * *VT* (``realtype``, input) -- array containing temporary workspace of same size as *Y*.
      * *IER*  (``int``, output) -- return flag  (0 if success, >0 if a recoverable
        failure, <0 if a non-recoverable failure).
      
   **Notes:**
   Typically this routine will use only *NEQ*, *T*, *Y*, *GAMMA*, *R*,
   *LR*, and *Z*.  It must solve the preconditioner linear system 
   :math:`Pz = r`.  The preconditioner (or the product of the left and
   right preconditioners if both are nontrivial) should be an
   approximation to the matrix  :math:`M - \gamma J`, where
   :math:`M` is the system mass matrix, :math:`\gamma` is the input
   GAMMA, and :math:`J = \frac{\partial f_I}{\partial y}`. 


Notes:

(a) If the user's :f:func:`FARKJTIMES()` or :f:func:`FARKPSET()` routine
    uses difference quotient approximations, it may need to use the
    error weight array *EWT* and/or the unit roundoff, in the
    calculation of suitable increments. Also, if :f:func:`FARKPSOL()`
    uses an iterative method in its solution, the residual vector
    :math:`\rho = r - Pz` of the system should be made less than
    :math:`\delta =` *DELTA* in the weighted l2 norm, i.e. 
    
    .. math::
       \left(\sum_i \left(\rho_i\, EWT_i\right)^2 \right)^{1/2} < \delta.

(b) If needed in :f:func:`FARKJTIMES()`, :f:func:`FARKPSOL()`, or
    :f:func:`FARKPSET()`, the error weight array *EWT* can be
    obtained by calling :f:func:`FARKGETERRWEIGHTS()` using one of the
    work arrays as temporary storage for *EWT*. 

(c) If needed in :f:func:`FARKJTIMES()`, :f:func:`FARKPSOL()`, or
    :f:func:`FARKPSET()`, the unit roundoff can be obtained as the
    optional output *ROUT(6)* (available after the call to
    :f:func:`FARKMALLOC()`) and can be passed using either the *RPAR*
    user data array or a common block. 




.. _FInterface.Solution:

Problem solution
-----------------------

Carrying out the integration is accomplished by making calls to
:f:func:`FARKODE()`.


.. f:subroutine:: FARKODE(TOUT, T, Y, ITASK, IER)
   
   Fortran interface to the C routine :c:func:`ARKode()`
   for performing the solve, along with many of the ARK*Get*
   routines for reporting on solver statistics.
   
   **Arguments:** 
      * *TOUT* (``realtype``, input) -- next value of :math:`t` at
	which a solution is desired.

      * *T* (``realtype``, output) -- value of independent variable
	that corresponds to the output *Y*

      * *Y* (``realtype``, output) -- array containing dependent state
	variables on output. 

      * *ITASK* (``int``, input) -- task indicator :

        * 1 = normal mode (overshoot *TOUT* and interpolate the solution)

        * 2 = one-step mode (return after each internal step taken)

        * 3 = normal 'tstop' mode (like 1, but integration never
          proceeds past *TSTOP*, which must be specified through a
          preceding call to :f:func:`FARKSETRIN()` using the key
          *STOP_TIME*)

        * 4 = one step 'tstop' mode (like 2, but integration never
	  goes past *TSTOP*).

      * *IER* (int, output) -- completion flag: 

	* 0 = success, 

	* 1 = tstop return, 

	* 2 = root return, 

	* values -1, ..., -10 are failure modes (see :c:func:`ARKode()` and
          :ref:`Constants`).
      
   **Notes:**
   The current values of the optional outputs are immediately
   available in *IOUT* and *ROUT* upon return from this function (see
   :ref:`FInterface.IOUTTable` and :ref:`FInterface.ROUTTable`). 

   A full description of error flags and output behavior of the solver
   (values filled in for *T* and *Y*) is provided in the description
   of :c:func:`ARKode()`.
   



.. _FInterface.AdditionalOutput:

Additional solution output
---------------------------------------

After a successful return from :f:func:`FARKODE()`, the routine
:f:func:`FARKDKY()` may be used to obtain a derivative of the solution,
of order up to 3, at any :math:`t` within the last step taken. 


.. f:subroutine:: FARKDKY(T, K, DKY, IER)
   
   Fortran interface to the C routine :f:func:`ARKDKY()` for
   interpolating output of the solution or its derivatives at any
   point within the last step taken.
      
   **Arguments:** 
      * *T* (``realtype``, input) -- time at which solution derivative
	is desired, within the interval :math:`[t_n-h,t_n]`.
      * *K* (``int``, input) -- derivative order :math:`(0 \le k \le 3)`.
      * *DKY* (``realtype``, output) -- array containing the computed
	*K*-th derivative of :math:`y`.
      * *IER* (``int``, output) -- return flag (0 if success, <0 if an
	illegal argument).



.. _FInterface.ReInit:

Problem reinitialization
---------------------------------------

To re-initialize the ARKode solver for the solution of a new
problem of the same size as one already solved, the user must call
:f:func:`FARKREINIT()`:


.. f:subroutine:: FARKREINIT(T0, Y0, IMEX, IATOL, RTOL, ATOL, IER)
   
   Re-initializes the Fortran interface to the ARKode solver.
      
   **Arguments:**  The arguments have the same names and meanings as those of
   :f:func:`FARKMALLOC()`.
      
   **Notes:**
   This routine performs no memory allocation, instead using the
   existing memory created by the previous :f:func:`FARKMALLOC()`
   call.  The call to specify the linear system solution method may
   or may not be needed. 


Following a call to :f:func:`FARKREINIT()`, a call to specify the
linear system solver must be made if the choice of linear solver is
being changed.  Otherwise, a call to reinitialize the linear solver
last used is only needed if linear solver input parameters need
modification. 

In the case of the BAND solver, for any change in the
half-bandwidth parameters, call :f:func:`FARKBAND()` (or
:f:func:`FARKLAPACKBAND()`) again, as described above.

In the case of SPGMR, for a change of inputs other than *MAXL*,
the user may call the routine :f:func:`FARKSPGMRREINIT()` to
reinitialize SPGMR without reallocating its memory, as follows: 



.. f:subroutine:: FARKSPGMRREINIT(IPRETYPE, IGSTYPE, DELT, IER)
   
   Re-initializes the Fortran interface to the SPGMR linear solver.
      
   **Arguments:**  The arguments have the same names and meanings as those of
   :f:func:`FARKSPGMR()`.
   


However, if *MAXL* is being changed, then the user should call
:f:func:`FARKSPGMR()` instead, since memory will need to be
deallocated/reallocated by the solver.

In the case of SPBCG, for a change in any inputs, the user can
reinitialize SPBCG without reallocating its memory by calling
:f:func:`FARKSPBCGREINIT()`, as follows:


.. f:subroutine:: FARKSPBCGREINIT(IPRETYPE, MAXL, DELT, IER)
   
   Re-initializes the Fortran interface to the SPBCG
   linear solver.
      
   **Arguments:**  The arguments have the same names and meanings as
   those of :f:func:`FARKSPBCG()`.



In the case of SPTFQMR, for a change in any inputs, the user can
reinitialize SPTFQMR without reallocating its memory by calling
:f:func:`FARKSPTFQMRREINIT()`, as follows:


.. f:subroutine:: FARKSPTFQMRREINIT(IPRETYPE, MAXL, DELT, IER)
   
   Re-initializes the Fortran interface to the SPBTFQMR
   linear solver.
      
   **Arguments:**  The arguments have the same names and meanings as
   those of :f:func:`FARKSPTFQMR()`.


In the case of SPFGMR, for a change of inputs other than *MAXL*,
the user may call the routine :f:func:`FARKSPFGMRREINIT()` to
reinitialize SPFGMR without reallocating its memory, as follows: 


.. f:subroutine:: FARKSPFGMRREINIT(IPRETYPE, IGSTYPE, DELT, IER)
   
   Re-initializes the Fortran interface to the SPFGMR
   linear solver.
      
   **Arguments:**  The arguments have the same names and meanings as
   those of :f:func:`FARKSPFGMR()`.
   
However, if *MAXL* is being changed, then the user should call
:f:func:`FARKSPFGMR()` instead, since memory will need to be
deallocated/reallocated by the solver.


In the case of PCG, for a change in any inputs, the user can
reinitialize PCG without reallocating its memory by calling
:f:func:`FARKPCGREINIT()`, as follows:


.. f:subroutine:: FARKPCGREINIT(IPRETYPE, MAXL, DELT, IER)
   
   Re-initializes the Fortran interface to the PCG
   linear solver.
      
   **Arguments:**  The arguments have the same names and meanings as
   those of :f:func:`FARKPCG()`.





.. _FInterface.Resize:

Resizing the ODE system
-----------------------------

For simulations involving changes to the number of equations and
unknowns in the ODE system (e.g. when solving a spatially-adaptive
PDE), the :f:func:`FARKODE()` integrator may be "resized" between
integration steps, through calls to the :f:func:`FARKRESIZE()`
function, that interfaces with the C routine :c:func:`ARKodeResize()`.
This function modifies ARKode's internal memory structures to use the
new problem size, without destruction of the temporal adaptivity
heuristics.  It is assumed that the dynamical time scales before and
after the vector resize will be comparable, so that all time-stepping
heuristics prior to calling :c:func:`FARKRESIZE` remain valid after
the call.  If instead the dynamics should be re-calibrated, the
FARKODE memory structure should be deleted with a call to
:f:func:`FARKFREE()`, and re-created with a call to
:f:func:`FARKMALLOC()`.


.. f:subroutine:: FARKRESIZE(T0, Y0, HSCALE, ITOL, RTOL, ATOL, IER)
   
   Re-initializes the Fortran interface to the ARKode solver for a
   differently-sized ODE system.
      
   **Arguments:** 
      * *T0* (``realtype``, input) -- initial value of the independent
	variable :math:`t`.

      * *Y0* (``realtype``, input) -- array of dependent-variable
	initial conditions.  

      * *HSCALE* (``realtype``, input) -- desired step size scale factor:

        * 1.0 is the default,

        * any value <= 0.0 results in the default.

      * *ITOL* (``int``, input) -- flag denoting that a new relative
	tolerance and vector of absolute tolerances are supplied in
	the *RTOL* and *ATOL* arguments: 

        * 0 = retain the current scalar-valued relative and absolute
	  tolerances, or the user-supplied error weight function,
	  :f:func:`FARKEWT()`. 

        * 1 = *RTOL* contains the new scalar-valued relative tolerance 
          and *ATOL* contains a new array of absolute tolerances.

      * *RTOL* (``realtype``, input) -- scalar relative tolerance.

      * *ATOL* (``realtype``, input) -- array of absolute tolerances.

      * *IER* (``int``, output) -- return flag (0 success, :math:`\ne 0` failure).
      
   **Notes:**
   This routine performs the opposite set of of operations as
   :f:func:`FARKREINIT()`: it does not reinitialize any of the
   time-step heuristics, but it does perform memory reallocation.  


Following a call to :f:func:`FARKRESIZE()`, a call to specify the
linear system solver must be made **after** the call to
:f:func:`FARKRESIZE()`, since the internal data structures for the
linear solver will also be the incorrect size.  

If any user-supplied linear solver helper routines were used (Jacobian 
evaluation, Jacobian-vector product, preconditioning, etc.), then the
relevant "set" routines to specify their usage must be called again
**following** the re-specification of the linear solver module.





.. _FInterface.Deallocation:

Memory deallocation
---------------------------------------

To free the internal memory created by :f:func:`FARKMALLOC()`, the user
may call :f:func:`FARKFREE()`, as follows:


.. f:subroutine:: FARKFREE()
   
   Frees the internal memory created by :f:func:`FARKMALLOC()`.
      
   **Arguments:** None.
