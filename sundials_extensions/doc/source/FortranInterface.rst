.. _FortranInterface:

FARKode, an Interface Module for FORTRAN Applications
=====================================================

The FARKODE interface module is a package of C functions which
support the use of the ARKODE solver, for the solution of ODE
systems 

.. math::
   M \dot{y} = f_E(t,y) + f_I(t,y),

in a mixed Fortran/C setting. While ARKODE is written in C, it is
assumed here that the user's calling program and user-supplied
problem-defining routines are written in Fortran. This package
provides the necessary interface to ARKODE for both the serial and
the parallel NVECTOR implementations.


.. _FInterface.Routines:

FARKODE routines
----------------

The user-callable functions, with the corresponding ARKODE
functions, are as follows:

- Interface to the NVECTOR modules

  - :c:func:`FNVINITS()` (defined by NVECTOR_SERIAL) interfaces to
    :c:func:`N_VNewEmpty_Serial()`.

  - :c:func:`FNVINITP()` (defined by NVECTOR_PARALLEL) interfaces to
    :c:func:`N_VNewEmpty_Parallel()`. 

- Interface to the main ARKODE module

  - :c:func:`FARKMALLOC()` interfaces to :c:func:`ARKodeCreate()`,
    :c:func:`ARKodeSetUserData()`, and :c:func:`ARKodeInit()`, as well
    as one of :c:func:`ARKodeSStolerances()` or :c:func:`ARKodeSVtolerances()`.

  - :c:func:`FARKREINIT()` interfaces to :c:func:`ARKodeReInit()`.

  - :c:func:`FARKSETIIN()` and :c:func:`FARKSETRIN()` interface to the
    ARKodeSet* functions (see :ref:`CInterface.OptionalInputs`).

  - :c:func:`FARKEWTSET()` interfaces to :c:func:`ARKodeWFtolerances()`.

  - :c:func:`FARKODE()` interfaces to :c:func:`ARKode()`, the
    ARKodeGet* functions (see :ref:`CInterface.OptionalOutputs`), 
    and to the optional output functions for the selected linear
    solver module (see :ref:`CInterface.OptionalOutputs`). 

  - :c:func:`FARKDKY()` interfaces to the interpolated output function
    :c:func:`ARKodeGetDky()`.

  - :c:func:`FARKGETERRWEIGHTS()` interfaces to
    :c:func:`ARKodeGetErrWeights()`.

  - :c:func:`FARKGETESTLOCALERR()` interfaces to
    :c:func:`ARKodeGetEstLocalErrors()`.

  - :c:func:`FARKFREE()` interfaces to :c:func:`ARKodeFree()`.

- Interface to the linear solver modules

  - :c:func:`FARKDENSE()` interfaces to :c:func:`ARKDense()`.

  - :c:func:`FARKDENSESETJAC()` interfaces to :c:func:`ARKDlsSetDenseJacFn()`.

  - :c:func:`FARKLAPACKDENSE()` interfaces to :c:func:`ARKLapackDense()`.

  - :c:func:`FARKLAPACKDENSESETJAC()` interfaces to :c:func:`ARKDlsSetDenseJacFn()`.

  - :c:func:`FARKBAND()` interfaces to :c:func:`ARKBand()`.

  - :c:func:`FARKBANDSETJAC()` interfaces to :c:func:`ARKDlsSetBandJacFn()`.

  - :c:func:`FARKLAPACKBAND()` interfaces to :c:func:`ARKLapackBand()`.

  - :c:func:`FARKLAPACKBANDSETJAC()` interfaces to :c:func:`ARKDlsSetBandJacFn()`.

  - :c:func:`FARKSPGMR()` interfaces to :c:func:`ARKSpgmr()` and the SPGMR optional input
    functions (see :ref:`CInterface.ARKSpilsInputTable`).

  - :c:func:`FARKSPGMRREINIT()` interfaces to the SPGMR optional input
    functions (see :ref:`CInterface.ARKSpilsInputTable`).

  - :c:func:`FARKSPBCG()` interfaces to :c:func:`ARKSpbcg()` and the SPBCG optional input
    functions (see :ref:`CInterface.ARKSpilsInputTable`).

  - :c:func:`FARKSPBCGREINIT()` interfaces to the SPBCG optional input
    functions.

  - :c:func:`FARKSPTFQMR()` interfaces to :c:func:`ARKSptfqmr()` and the SPTFQMR optional
    input functions.

  - :c:func:`FARKSPTFQMRREINIT()` interfaces to the SPTFQMR optional input
    functions.

  - :c:func:`FARKSPILSSETJAC()` interfaces to :c:func:`ARKSpilsSetJacTimesVecFn()`.

  - :c:func:`FARKSPILSSETPREC()` interfaces to :c:func:`ARKSpilsSetPreconditioner()`.


The user-supplied functions, each listed with the corresponding
internal interface function which calls it (and its type within
ARKode), are as follows:

.. cssclass:: table-bordered

+--------------------------+------------------------+-----------------------------------+
| FARKODE routine          | ARKode routine         | ARKode interface                  |
| (FORTRAN, user-supplied) | (C, interface)         | function type                     |
+==========================+========================+===================================+
| :c:func:`FARKIFUN()`     | FARKfi                 | :c:func:`ARKRhsFn()`              |
+--------------------------+------------------------+-----------------------------------+
| :c:func:`FARKEFUN()`     | FARKfe                 | :c:func:`ARKRhsFn()`              |
+--------------------------+------------------------+-----------------------------------+
| :c:func:`FARKDJAC()`     | FARKDenseJac           | :c:func:`ARKDlsDenseJacFn()`      |
+--------------------------+------------------------+-----------------------------------+
| :c:func:`FARKLDJAC()`    | FARKLapackDenseJac     | :c:func:`ARKDlsDenseJacFn()`      |
+--------------------------+------------------------+-----------------------------------+
| :c:func:`FARKBJAC()`     | FARKBandJac            | :c:func:`ARKDlsBandJacFn()`       |
+--------------------------+------------------------+-----------------------------------+
| :c:func:`FARKLBJAC()`    | FARKLapackBandJac      | :c:func:`ARKDlsBandJacFn()`       |
+--------------------------+------------------------+-----------------------------------+
| :c:func:`FARKPSET()`     | FARKPSet               | :c:func:`ARKSpilsPrecSetupFn()`   |
+--------------------------+------------------------+-----------------------------------+
| :c:func:`FARKPSOL()`     | FARKPSol               | :c:func:`ARKSpilsPrecSolveFn()`   |
+--------------------------+------------------------+-----------------------------------+
| :c:func:`FARKJTIMES()`   | FARKJtimes             | :c:func:`ARKSpilsJacTimesVecFn()` |
+--------------------------+------------------------+-----------------------------------+
| :c:func:`FARKEWT()`      | FARKEwtSet             | :c:func:`ARKEwtFn()`              |
+--------------------------+------------------------+-----------------------------------+

In contrast to the case of direct use of ARKode, and of most
Fortran ODE solvers, the names of all user-supplied routines here are
fixed, in order to maximize portability for the resulting
mixed-language program. 


.. _FInterface.Portability:

Important notes on portability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this package, the names of the interface functions, and the names
of the Fortran user routines called by them, appear as dummy names
which are mapped to actual values by a series of definitions in the
header files ``farkode.h``, ``farkroot.h``, ``farkbp.h``, and
``farkbbd.h``. By default, those mapping definitions depend in turn on
the C macro ``F77_FUNC`` defined in the header file
``sundials_config.h`` and decided upon at configuration time (see
:ref:`Installation`). 


.. _FInterface.DataTypes:

Data Types
"""""""""""

Throughout this documentation, we will refer to data types according
to their usage in SUNDIALS.  The equivalent types to these may
vary, depending on your computer architecture and on how SUNDIALS
was compiled (see :ref:`Installation`).  A Fortran user should take
care that all arguments passed through this Fortran/C interface  are
declared of the appropriate type. 

**Integers**: SUNDIALS uses both ``int`` and ``long int`` types:

   ``int`` -- equivalent to an ``INTEGER`` or ``INTEGER*4`` in Fortran

   ``long int`` -- this will depend on the computer architecture:
   
      32-bit -- equivalent to an ``INTEGER`` or ``INTEGER*4`` in Fortran

      64-bit -- equivalent to an ``INTEGER*8`` in Fortran
	      
**Real numbers**:  As discussed in :ref:`Installation`, at compilation
SUNDIALS allows the configuration option  ``--with-precision``,
that accepts values of ``single``, ``double`` or ``extended`` (the
default is ``double``).  This choice dictates the size of a
``realtype`` variable.  The corresponding Fortran types for these
``realtype`` sizes are: 

   ``single`` -- equivalent to a ``REAL`` or ``REAL*4`` in Fortran

   ``double`` -- equivalent to a ``DOUBLE PRECISION`` or ``REAL*8``
   in Fortran
 
   ``extended`` -- equivalent to a ``REAL*16`` in Fortran


.. _FInterface.Usage:

Usage of the FARKODE interface module
-------------------------------------

The usage of FARKODE requires calls to five or more interface
functions, depending on the method options selected, and one or more
user-supplied routines which define the problem to be solved.  These 
function calls and user routines are summarized separately below.
Some details are omitted, and the user is referred to the description
of the corresponding ARKode functions for complete information on
the arguments of any given user-callable interface routine.  The usage
of FARKODE for rootfinding and with preconditioner modules is
described in later subsections.

Steps marked [**S**] in the instructions below apply to the serial
NVECTOR implementation (NVECTOR_SERIAL) only, while those
marked with a [**P**] apply to NVECTOR_PARALLEL.


.. _FInterface.RHS:

Right-hand side specification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user must in all cases supply at least one of the following Fortran 
routines:



   .. c:function:: SUBROUTINE FARKIFUN(T, Y, YDOT, IPAR, RPAR, IER)
   
      Sets the YDOT array to :math:`f_I(t,y)`, the
      implicit portion of the right-hand side of the ODE system, as
      function of the independent variable T :math:`=t` and the array
      of dependent state variables Y :math:`=y`.
      
      **Arguments:**
         * T (``realtype``, input) -- current value of the independent variable
         * Y (``realtype``, input) -- array containing state variables  
         * YDOT (``realtype``, output) -- array containing state derivatives 
	 * IPAR (``long int``, input) -- array containing integer user
           data that was passed to :c:func:`FARKMALLOC()` 
         * RPAR (``realtype``, input) -- array containing real user
           data that was passed to :c:func:`FARKMALLOC()` 
         * IER (``int``, output) -- return flag (0 success, >0
           recoverable error, <0 unrecoverable error)  
   


   .. c:function:: SUBROUTINE FARKEFUN(T, Y, YDOT, IPAR, RPAR, IER)
   
      Sets the YDOT array to :math:`f_E(t,y)`, the
      explicit portion of the right-hand side of the ODE system, as
      function of the independent variable T :math:`=t` and the array
      of dependent state variables Y :math:`=y`.
      
      **Arguments:**
         * T (``realtype``, input) -- current value of the independent variable
         * Y (``realtype``, input) -- array containing state variables  
         * YDOT (``realtype``, output) -- array containing state derivatives 
	 * IPAR (``long int``, input) -- array containing integer user
           data that was passed to :c:func:`FARKMALLOC()` 
         * RPAR (``realtype``, input) -- array containing real user
           data that was passed to :c:func:`FARKMALLOC()` 
         * IER (``int``, output) -- return flag (0 success, >0
           recoverable error, <0 unrecoverable error)  


.. _FInterface.NVector:

NVECTOR module initialization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[**S**] To initialize the serial NVECTOR module, the user must
call the function FNVINITS with the argument KEY = 4.



   .. c:function:: SUBROUTINE FNVINITS(KEY, NEQ, IER)
   
      Initializes the Fortran interface to the serial
      NVECTOR module.
      
      **Arguments:** 
         * KEY (``int``, input) -- integer flag denoting which solver is to be
           used (1 is CVODE, 2 is IDA, 3 is KINSOL and 4 is
           ARKode) 
         * NEQ (``long int``, input) -- size of the ODE system 
         * IER (``int``, output) -- return flag (0 success, :math:`\ne 0` failure) 




[**P**] To initialize the parallel NVECTOR module, the user must
call the function FNVINITP with the argument KEY = 4.



   .. c:function:: SUBROUTINE FNVINITP(COMM, KEY, NLOCAL, NGLOBAL, IER)
   
      Initializes the Fortran interface to the parallel
      NVECTOR module.
      
      **Arguments:** 
         * COMM (``int``, input) -- the MPI communicator 
         * KEY (``int``, input) -- integer flag denoting which solver is to be
           used (1 is CVODE, 2 is IDA, 3 is KINSOL and 4 is
           ARKode) 
         * NLOCAL (``long int``, input) -- local size of vectors on this processor 
         * NGLOBAL (``long int``, input) -- the size of the ODE system, and the global size of
           vectors (the sum of all values of NLOCAL) 
         * IER (``int``, output) -- return flag (0 success, :math:`\ne 0` failure) 
      
      **Notes:** If the header file ``sundials_config.h`` defines
      ``SUNDIALS_MPI_COMM_F2C`` to be 1 (meaning the MPI implementation 
      used to build SUNDIALS includes the ``MPI_Comm_f2c`` function),
      then COMM can be any valid MPI communicator.  Otherwise,
      ``MPI_COMM_WORLD`` will be used, so the user can just pass an
      integer value as a placeholder.



.. _FInterface.Problem:

Problem specification
^^^^^^^^^^^^^^^^^^^^^^^

To set various problem and solution parameters and allocate internal
memory, the user must call FARKMALLOC.



   .. c:function:: SUBROUTINE FARKMALLOC(T0, Y0, IMEX, IATOL, RTOL, ATOL, IOUT, ROUT, IPAR, RPAR, IER)
   
      Initializes the Fortran interface to the ARKode
      solver, providing interfaces to the C routines :c:func:`ARKodeCreate()`,
      :c:func:`ARKodeSetUserData()`, and :c:func:`ARKodeInit()`, as well
      as one of :c:func:`ARKodeSStolerances()` or
      :c:func:`ARKodeSVtolerances()`.
      
      **Arguments:** 
         * T0 (``realtype``, input) -- initial value of :math:`t` 
         * Y0 (``realtype``, input) -- array of initial conditions 
         * IMEX (``int``, input) -- flag denoting basic integration method:
            * 0 = implicit, 
	    * 1 = explicit, 
	    * 2 = imex.
         * IATOL (``int``, input) -- type for absolute tolerance input ATOL:
            * 1 = scalar, 
            * 2 = array,
      	    * 3 = user-supplied function; the user must subsequently call
              :c:func:`FARKEWTSET()` and supply a routine :c:func:`FARKEWT()` to
              compute the error weight vector.
         * RTOL (``realtype``, input) -- scalar relative tolerance 
         * ATOL (``realtype``,
           input) -- scalar or array absolute tolerance 
         * IOUT (``long
           int``, input/output) -- array of length 22 for integer optional outputs 
         * ROUT (``realtype``, input/output) -- array of length 6 for real optional outputs
         * IPAR (``long int``, input/output) -- array of user integer data, which will be passed
           unmodified to all user-provided routines 
         * RPAR (``realtype``, input/output) -- array with user real data, which will be passed
           unmodified to all user-provided routines 
         * IER (``int``, output) -- return flag (0 success, :math:`\ne 0` failure) 
      
      **Notes:** Modifications to the user data arrays IPAR and RPAR
      inside a user-provided routine will be propagated to all
      subsequent calls to such routines. The optional outputs
      associated with the main ARKode integrator are listed in
      :ref:`FInterface.IOUTTable` and :ref:`FInterface.ROUTTable`, in
      the section :ref:`FInterface.OptionalOutputs`. 




As an alternative to providing tolerances in the call to
:c:func:`FARKMALLOC()`, the user may provide a routine to compute the
error weights used in the WRMS norm evaluations.  If supplied, it must
have the following form:



   .. c:function:: SUBROUTINE FARKEWT(Y, EWT, IPAR, RPAR, IER)
   
      It must set the positive components of the error weight
      vector EWT for the calculation of the WRMS norm of Y.
      
      **Arguments:** 
         * Y (``realtype``, input) -- array containing state variables  
         * EWT (``realtype``, output) -- array containing the error weight vector  
         * IPAR (``long int``, input) -- array containing the integer user data that was passed
           to :c:func:`FARKMALLOC()` 
         * RPAR (``realtype``, input) -- array containing the real user data that was passed to
           :c:func:`FARKMALLOC()` 
         * IER (``int``, output) -- return flag (0 success, :math:`\ne 0` failure) 



   
If the FARKEWT routine is provided, then, following the call to
:c:func:`FARKMALLOC()`, the user must call the function FARKEWTSET.



   .. c:function:: SUBROUTINE FARKEWTSET(FLAG, IER)
   
      Informs FARKODE to use the user-supplied
      :c:func:`FARKEWT()` function.
      
      **Arguments:** 
         * FLAG (``int``, input) -- flag, use "1" to denoting to use FARKEWT.
         * IER (``int``, output) -- return flag (0 success, :math:`\ne 0` failure) 



.. _FInterface.OptionalInputs:

Set optional inputs
^^^^^^^^^^^^^^^^^^^^^^^

To set desired optional inputs, the user can call the routines
:c:func:`FARKSETIIN()` and :c:func:`FARKSETRIN()`, as described below.



   .. c:function:: SUBROUTINE FARKSETIIN(KEY, IVAL, IER)
   
      Specification routine to pass optional integer inputs
      to the :c:func:`FARKODE()` solver.
      
      **Arguments:** 
         * KEY (quoted string, input) -- which optional input
           is set (see :ref:`FInterface.IINOptionTable`).
         * IVAL (``long int``, input) -- the integer input value to be used 
         * IER (``int``, output) -- return flag (0 success, :math:`\ne 0` failure) 


.. _FInterface.IINOptionTable:

Table: Keys for setting FARKODE integer optional inputs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

=================  =========================================
Key                ARKode routine
=================  =========================================
ORDER              :c:func:`ARKodeSetOrder()`
DENSE_ORDER        :c:func:`ARKodeSetDenseOrder()`
LINEAR             :c:func:`ARKodeSetLinear()`
NONLINEAR          :c:func:`ARKodeSetNonlinear()`
EXPLICIT           :c:func:`ARKodeSetExplicit()`
IMPLICIT           :c:func:`ARKodeSetImplicit()`
IMEX               :c:func:`ARKodeSetImEx()`
IRK_TABLE_NUM      :c:func:`ARKodeSetIRKTableNum()`
ERK_TABLE_NUM      :c:func:`ARKodeSetERKTableNum()`
ARK_TABLE_NUM `*`  :c:func:`ARKodeSetARKTableNum()`      
MAX_NSTEPS         :c:func:`ARKodeSetMaxNumSteps()`
HNIL_WARNS         :c:func:`ARKodeSetMaxHnilWarns()`
PREDICT_METHOD     :c:func:`ARKodeSetPredictorMethod()`
MAX_ERRFAIL        :c:func:`ARKodeSetMaxErrTestFails()`
MAX_NITERS         :c:func:`ARKodeSetMaxNonlinIters()`
MAX_CONVFAIL       :c:func:`ARKodeSetMaxConvFails()`
ADAPT_METHOD       :c:func:`ARKodeSetAdaptivityMethod()`
ADAPT_SMALL_NEF    :c:func:`ARKodeSetAdaptivityConstants()`
LSETUP_MSBP        :c:func:`ARKodeSetLSetupConstants()`
=================  =========================================

`*` When setting ARK_TABLE_NUM, pass in IVAL as an array of
length 2, specifying the IRK table number first, then the ERK table
number. 


      
   .. c:function:: SUBROUTINE FARKSETRIN(KEY, RVAL, IER)
   
      Specification routine to pass optional real inputs
      to the :c:func:`FARKODE()` solver.
      
      **Arguments:** 
         * KEY (quoted string, input) -- which optional input
           is set (see :ref:`FInterface.RINOptionTable`).
         * RVAL (``realtype``, input) -- the real input value to be used 
         * IER (``int``, output) -- return flag (0 success, :math:`\ne 0` failure) 


.. _FInterface.RINOptionTable:

Table: Keys for setting FARKODE real optional inputs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

============  =========================================
Key           ARKode routine
============  =========================================
INIT_STEP     :c:func:`ARKodeSetInitStep()`
MAX_STEP      :c:func:`ARKodeSetMaxStep()`
MIN_STEP      :c:func:`ARKodeSetMinStep()`
STOP_TIME     :c:func:`ARKodeSetStopTime()`
NLCONV_COEF   :c:func:`ARKodeSetNonlinConvCoef()`
ADAPT_CFL     :c:func:`ARKodeSetAdaptivityMethod()`
ADAPT_SAFETY  :c:func:`ARKodeSetAdaptivityMethod()`
ADAPT_BIAS    :c:func:`ARKodeSetAdaptivityMethod()`
ADAPT_GROWTH  :c:func:`ARKodeSetAdaptivityMethod()`
ADAPT_LB      :c:func:`ARKodeSetAdaptivityMethod()`
ADAPT_UB      :c:func:`ARKodeSetAdaptivityMethod()`
ADAPT_K1      :c:func:`ARKodeSetAdaptivityMethod()`
ADAPT_K2      :c:func:`ARKodeSetAdaptivityMethod()`
ADAPT_K3      :c:func:`ARKodeSetAdaptivityMethod()`
ADAPT_ETAMX1  :c:func:`ARKodeSetAdaptivityConstants()`
ADAPT_ETAMXF  :c:func:`ARKodeSetAdaptivityConstants()`
ADAPT_ETACF   :c:func:`ARKodeSetAdaptivityConstants()`
NEWT_CRDOWN   :c:func:`ARKodeSetNewtonConstants()`
NEWT_RDIV     :c:func:`ARKodeSetNewtonConstants()`
LSETUP_DGMAX  :c:func:`ARKodeSetLSetupConstants()`
============  =========================================


Alternatively, if a user wishes to reset all of the options to their
default values, they may call the routine FARKSETDEFAULTS.



   .. c:function:: SUBROUTINE FARKSETDEFAULTS(IER)
   
      Specification routine to reset all FARKODE optional
      inputs to their default values.
      
      **Arguments:** 
         * IER (``int``, output) -- return flag (0 success, :math:`\ne 0` failure) 
   


FARKODE supplies additional routines to specify optional advanced
inputs to the :c:func:`ARKode()` solver.  These are summarized below,
and the user is referred to their C routine counterparts for more
complete information. 



   .. c:function:: SUBROUTINE FARKSETERKTABLE(S, Q, P, C, A, B, BEMBED, IER)
   
      Interface to the routine :c:func:`ARKodeSetERKTable()`.
      
      **Arguments:** 
         * S (``int``, input) -- number of stages in the table 
         * Q (``int``, input) -- global order of accuracy of the method 
         * P (``int``, input) -- global order of accuracy of the embedding 
         * C (``realtype``, input) -- array of length S containing the stage times
         * A (``realtype``, input) -- array of length S*S containing the ERK coefficients
           (stored in row-major, "C", order) 
         * B (``realtype``, input) -- array of length S containing the solution coefficients 
         * BEMBED (``realtype``, input) -- array of length S containing the embedding
           coefficients 
         * IER (``int``, output) -- return flag (0 success, :math:`\ne 0` failure) 



   .. c:function:: SUBROUTINE FARKSETIRKTABLE(S, Q, P, C, A, B, BEMBED, IER)
   
      Interface to the routine :c:func:`ARKodeSetIRKTable()`.
      
      **Arguments:** 
         * S (``int``, input) -- number of stages in the table 
         * Q (``int``, input) -- global order of accuracy of the method 
         * P (``int``, input) -- global order of accuracy of the embedding 
         * C (``realtype``, input) -- array of length S containing the stage times
         * A (``realtype``, input) -- array of length S*S containing the IRK coefficients
           (stored in row-major, "C", order) 
         * B (``realtype``, input) -- array of length S containing the solution coefficients 
         * BEMBED (``realtype``, input) -- array of length S containing the embedding
           coefficients 
         * IER (``int``, output) -- return flag (0 success, :math:`\ne 0` failure) 
   

   
   .. c:function:: SUBROUTINE FARKSETARKTABLES(S, Q, P, C, AI, AE, B, BEMBED, IER)
   
      Interface to the routine :c:func:`ARKodeSetARKTables()`.
      
      **Arguments:** 
         * S (``int``, input) -- number of stages in the table 
         * Q (``int``, input) -- global order of accuracy of the method 
         * P (``int``, input) -- global order of accuracy of the embedding 
         * C (``realtype``, input) -- array of length S containing the stage times
         * AI (``realtype``, input) -- array of length S*S containing the IRK coefficients
           (stored in row-major, "C", order) 
         * AE (``realtype``, input) -- array of length S*S containing the ERK coefficients
           (stored in row-major, "C", order) 
         * B (``realtype``, input) -- array of length S containing the solution coefficients 
         * BEMBED (``realtype``, input) -- array of length S containing the embedding
           coefficients 
         * IER (``int``, output) -- return flag (0 success, :math:`\ne 0` failure) 
   

   
.. _FInterface.LinearSolver:

Linear solver specification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the case of using either an implicit or ImEx method, the solution
of each Runge-Kutta stage may involve the solution of linear systems
related to the Jacobian :math:`J = \frac{\partial f_I}{\partial y}` of
the implicit portion of the ODE system. ARKode presently includes
seven choices for the treatment of these systems, and the user of
FARKODE must call a routine with a specific name to make the
desired choice. 


[**S**] Dense treatment of the linear system
"""""""""""""""""""""""""""""""""""""""""""""""

To use the direct dense linear solver based on the internal ARKode
implementation, the user must call the FARKDENSE routine.



   .. c:function:: SUBROUTINE FARKDENSE(NEQ, IER)
   
      Interfaces with the :c:func:`ARKDense()` function to
      specify use of the dense direct linear solver.
      
      **Arguments:** 
         * NEQ (``long int``, input) -- size of the ODE system 
         * IER (``int``, output) -- return flag (0 if success, -1 if a memory allocation
           error occurred, -2 for an illegal input) 



Alteratively, to use the LAPACK-based direct dense linear solver, a
user must call the similar FARKLAPACKDENSE routine.



   .. c:function:: SUBROUTINE FARKLAPACKDENSE(NEQ, IER)
   
      Interfaces with the :c:func:`ARKLapackDense()` function
      to specify use of the LAPACK the dense direct linear solver.
      
      **Arguments:** 
         * NEQ (``int``, input) -- size of the ODE system 
         * IER (``int``, output) -- return flag (0 if success, -1 if a memory allocation
           error occurred, -2 for an illegal input) 



As an option when using either of these dense linear solvers, the user
may supply a routine that computes a dense approximation of the system
Jacobian :math:`J = \frac{\partial f_I}{\partial y}`. If supplied, it
must have one of the following forms:



   .. c:function:: SUBROUTINE FARKDJAC(NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)
   
      Interface to provide a user-supplied dense Jacobian
      approximation function (of type :c:func:`ARKDenseJacFn()`), to be
      used by the :c:func:`FARKDENSE()` solver.
      
      **Arguments:** 
         * NEQ (``long int``, input) -- size of the ODE system 
         * T (``realtype``, input) -- current value of the independent variable 
         * Y (``realtype``, input) -- array containing values of the dependent state variables 
         * FY (``realtype``, input) -- array containing values of the dependent state derivatives 
         * DJAC (``realtype`` of size (NEQ,NEQ), output) -- 2D array containing the Jacobian entries 
         * H (``realtype``, input) -- current step size 
         * IPAR (``long int``, input) -- array containing integer user data that was passed to
           :c:func:`FARKMALLOC()` 
         * RPAR (``realtype``, input) -- array containing real user data that was passed to
           :c:func:`FARKMALLOC()` 
         * WK1, WK2, WK3  (``realtype``, input) -- array containing temporary workspace
           of same size as Y 
         * IER (``int``, output) -- return flag (0 if success, >0 if a recoverable error
           occurred, <0 if an unrecoverable error occurred) 
      
      **Notes:** Typically this routine will use only NEQ, T, Y, and
      DJAC. It must compute the Jacobian and store it column-wise in DJAC. 
   

   
   .. c:function:: SUBROUTINE FARKLDJAC(NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)
   
      Interface to provide a user-supplied dense Jacobian
      approximation function (of type :c:func:`ARKLapackJacFn()`), to be
      used by the :c:func:`FARKLAPACKDENSE()` solver.
      
      **Arguments:** these all match those for :c:func:`FARKDJAC()`.
   

   
If either of the above routines (:c:func:`FARKDJAC()` or
:c:func:`FARKLDJAC()`) uses difference quotient approximations, it may
need to use the error weight array EWT and current stepsize H
in the calculation of suitable increments. The array EWT can be
obtained by calling :c:func:`FARKGETERRWEIGHTS()` using one of the work
arrays as temporary storage for EWT. It may also need the unit
roundoff, which can be obtained as the optional output ROUT(6),
passed from the calling program to this routine using either RPAR
or a common block. 

If the :c:func:`FARKDJAC()` routine is provided, then, following the
call to :c:func:`FARKDENSE()`, the user must call the routine
FARKDENSESETJAC. 



   .. c:function:: SUBROUTINE FARKDENSESETJAC(FLAG, IER)
   
      Interface to the :c:func:`ARKDenseSetJacFn()` function,
      specifying to use the user-supplied routine :c:func:`FARKDJAC()` for
      the Jacobian approximation.
      
      **Arguments:** 
         * FLAG (``int``, input) -- any nonzero value specifies to use :c:func:`FARKDJAC()` 
         * IER (``int``, output) -- return flag (0 if success, :math:`\ne 0` if an error
           occurred) 
   

   
Similarly, if the :c:func:`FARKLDJAC()` routine is provided, then,
following the call to :c:func:`FARKLAPACKDENSE()`, the user must call
the routine FARKLAPACKDENSESETJAC. 



   .. c:function:: SUBROUTINE FARKLAPACKDENSESETJAC(FLAG, IER)
   
      Interface to the :c:func:`ARKLapackSetJacFn()` function,
      specifying to use the user-supplied routine :c:func:`FARKLDJAC()` for
      the Jacobian approximation.
      
      **Arguments:** 
         * FLAG (``int``, input) -- any nonzero value specifies to use
           :c:func:`FARKLDJAC()` 
         * IER (``int``, output) -- return flag (0 if success, :math:`\ne 0` if an error
           occurred) 




[**S**] Band treatment of the linear system
"""""""""""""""""""""""""""""""""""""""""""""""

To use the direct band linear solver based on the internal ARKode
implementation, the user must call the FARKBAND routine.



   .. c:function:: SUBROUTINE FARKBAND(NEQ, MU, ML, IER)
   
      Interfaces with the :c:func:`ARKBand()` function to
      specify use of the dense banded linear solver.
      
      **Arguments:** 
         * NEQ (``long int``, input) -- size of the ODE system 
         * MU (``long int``, input) -- upper half-bandwidth 
         * ML (``long int``, input) -- lower half-bandwidth 
         * IER (``int``, output) -- return flag (0 if success, -1 if a memory allocation
           error occurred, -2 for an illegal input) 



Alteratively, to use the LAPACK-based direct banded linear solver, a
user must call the similar FARKLAPACKBAND routine.



   .. c:function:: SUBROUTINE FARKLAPACKBAND(NEQ, MU, ML, IER)
   
      Interfaces with the :c:func:`ARKLapackBand()` function
      to specify use of the dense banded linear solver.
      
      **Arguments:** 
         * NEQ (``int``, input) -- size of the ODE system 
         * MU (``int``, input) -- upper half-bandwidth 
         * ML (``int``, input) -- lower half-bandwidth 
         * IER (``int``, output) -- return flag (0 if success, -1 if a memory allocation
           error occurred, -2 for an illegal input) 
   

   
As an option when using either of these banded linear solvers, the user
may supply a routine that computes a banded approximation of the
linear system Jacobian :math:`J = \frac{\partial f_I}{\partial y}`. If
supplied, it must have one of the following forms:


   .. c:function:: SUBROUTINE FARKBJAC(NEQ, MU, ML, MDIM, T, Y, FY, BJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)
   
      Interface to provide a user-supplied band Jacobian
      approximation function (of type :c:func:`ARKBandJacFn()`), to be
      used by the :c:func:`FARKBAND()` solver.
      
      **Arguments:** 
         * NEQ (``long int``, input) -- size of the ODE system 
         * MU   (``long int``, input) -- upper half-bandwidth 
         * ML   (``long int``, input) -- lower half-bandwidth 
         * MDIM (``long int``, input) -- leading dimension of BJAC array 
         * T    (``realtype``, input) -- current value of the independent variable 
         * Y    (``realtype``, input) -- array containing dependent state variables 
         * FY   (``realtype``, input) -- array containing dependent state derivatives 
         * BJAC (``realtype`` of size
           (MDIM,NEQ), output) -- 2D array containing the Jacobian entries 
         * H    (``realtype``, input) -- current step size 
         * IPAR (``long int``, input) -- array containing integer user data that was passed to
           :c:func:`FARKMALLOC()` 
         * RPAR (``realtype``, input) -- array containing real user data that was passed to
           :c:func:`FARKMALLOC()` 
         * WK1, WK2, WK3  (``realtype``, input) -- array containing temporary workspace
           of same size as Y 
         * IER (``int``, output) -- return flag (0 if success, >0 if a recoverable error
           occurred, <0 if an unrecoverable error occurred) 
      
      **Notes:**
      Typically this routine will use only NEQ, MU, ML, T, Y, and
      BJAC. It must load the MDIM by N array BJAC with the Jacobian
      matrix at the current :math:`(t,y)` in band form.  Store in
      BJAC(k,j) the Jacobian element :math:`J_{i,j}` with :math:`k = i
      - j + MU + 1` (or :math:`k = 1, \ldots ML+MU+1`) and :math:`j =
      1, \ldots, N`. 



   .. c:function:: SUBROUTINE FARKLBJAC(NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)
   
      Interface to provide a user-supplied banded Jacobian
      approximation function (of type :c:func:`ARKLapackJacFn()`), to be
      used by the :c:func:`FARKLAPACKBAND()` solver.
      
      **Arguments:** these all match those for :c:func:`FARKBJAC()`.
   


If either of the above routines (:c:func:`FARKBJAC()` or
:c:func:`FARKLBJAC()`) uses difference quotient approximations, it may
need to use the error weight array EWT and current stepsize H
in the calculation of suitable increments. The array EWT can be
obtained by calling :c:func:`FARKGETERRWEIGHTS()` using one of the work
arrays as temporary storage for EWT. It may also need the unit
roundoff, which can be obtained as the optional output ROUT(6),
passed from the calling program to this routine using either RPAR
or a common block. 

If the :c:func:`FARKBJAC()` routine is provided, then, following the
call to :c:func:`FARKBAND()`, the user must call the routine
FARKBANDSETJAC. 



   .. c:function:: SUBROUTINE FARKBANDSETJAC(FLAG, IER)
   
      Interface to the :c:func:`ARKBandSetJacFn()` function,
      specifying to use the user-supplied routine :c:func:`FARKBJAC()` for
      the Jacobian approximation.
      
      **Arguments:** 
         * FLAG (``int``, input) -- any nonzero value specifies to use
           :c:func:`FARKBJAC()`  
         * IER (``int``, output) -- return flag (0 if success, :math:`\ne 0` if an error
           occurred) 



Similarly, if the :c:func:`FARKLBJAC()` routine is provided, then,
following the call to :c:func:`FARKLAPACKBAND()`, the user must call
the routine FARKLAPACKBANDSETJAC. 



   .. c:function:: SUBROUTINE FARKLAPACKBANDSETJAC(FLAG, IER)
   
      Interface to the :c:func:`ARKLapackSetJacFn()` function,
      specifying to use the user-supplied routine :c:func:`FARKLBJAC()` for
      the Jacobian approximation.
      
      **Arguments:** 
         * FLAG (``int``, input) -- any nonzero value specifies to use
           :c:func:`FARKLBJAC()` 
         * IER (``int``, output) -- return flag (0 if success, :math:`\ne 0` if an error
           occurred) 





[**S**][**P**] SPGMR treatment of the linear systems
"""""""""""""""""""""""""""""""""""""""""""""""""""""

For the Scaled Preconditioned GMRES solution of the linear systems,
the user must call the FARKSPGMR routine.



   .. c:function:: SUBROUTINE FARKSPGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)
   
      Interfaces with the :c:func:`ARKSpgmr()` and
      ARKSpilsSet* routines to specify use of the SPGMR iterative
      linear solver.
      
      **Arguments:** 
         * IPRETYPE (``int``, input) -- preconditioner type : 
            * 0 = none 
	    * 1 = left only
	    * 2 = right only
      	    * 3 = both sides
         * IGSTYPE (``int``, input) -- Gram-schmidt process type : 
            * 1 = modified G-S
     	    * 2 = classical G-S
         * MAXL (``int``; input) -- maximum Krylov subspace dimension (0 for default) .
         * DELT (``realtype``, input) -- linear convergence tolerance factor (0.0 for default) .
         * IER (``int``, output) -- return flag (0 if success, -1 if a memory allocation
           error occurred, -2 for an illegal input) 



For descriptions of the optional user-supplied routines for use with
:c:func:`FARKSPGMR()` see the section :ref:`FInterface.SpilsUserSupplied`.





[**S**][**P**] SPBCG treatment of the linear systems
"""""""""""""""""""""""""""""""""""""""""""""""""""""

For the Scaled Preconditioned Bi-CGStab solution of the linear systems,
the user must call the FARKSPBCG routine.



   .. c:function:: SUBROUTINE FARKSPBCG(IPRETYPE, MAXL, DELT, IER)
   
      Interfaces with the :c:func:`ARKSpbcg()` and
      ARKSpilsSet* routines to specify use of the SPBCG iterative
      linear solver.
      
      **Arguments:**  The arguments are the same as those with the
      same names for :c:func:`FARKSPGMR()`. 



For descriptions of the optional user-supplied routines for use with
:c:func:`FARKSPBCG()` see the section :ref:`FInterface.SpilsUserSupplied`.





[**S**][**P**] SPTFQMR treatment of the linear systems
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For the Scaled Preconditioned TFQMR solution of the linear systems,
the user must call the FARKSPTFQMR routine.



   .. c:function:: SUBROUTINE FARKSPTFQMR(IPRETYPE, MAXL, DELT, IER)
   
      Interfaces with the :c:func:`ARKSptfqmr()` and
      ARKSpilsSet* routines to specify use of the SPTFQMR iterative
      linear solver.
      
      **Arguments:**  The arguments are the same as those with the same names
      for :c:func:`FARKSPGMR()`.
   


For descriptions of the optional user-supplied routines for use with
:c:func:`FARKSPTFQMR()` see the next section.



.. _FInterface.SpilsUserSupplied:

[**S**][**P**] User-supplied routines for SPGMR/SPBCG/SPTFQMR
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

With treatment of the linear systems by any of the Krylov iterative
solvers, there are three optional user-supplied routines --
:c:func:`FARKJTIMES()`, :c:func:`FARKPSET()` and :c:func:`FARKPSOL()`.
The specifications of these functions are given below.

As an option when using the SPGMR, SPBCG or SPTFQMR linear
solvers, the user may supply a routine that computes the product of
the system Jacobian :math:`J = \frac{\partial f_I}{\partial y}` and a
given vector :math:`v`.  If supplied, it must have the following form:



   .. c:function:: SUBROUTINE FARKJTIMES(V, FJV, T, Y, FY, H, IPAR, RPAR, WORK, IER)
   
      Interface to provide a user-supplied
      Jacobian-times-vector product approximation function (of type
      :c:func:`ARKSpilsJacTimesVecFn()`), to be used by one of the Krylov
      iterative linear solvers.
      
      **Arguments:** 
         * V (``realtype``, input) -- array containing the vector to multiply
         * FJV  (``realtype``, output) -- array containing resulting product vector
         * T    (``realtype``, input) -- current value of the independent variable
         * Y    (``realtype``, input) -- array containing dependent state variables
         * FY   (``realtype``, input) -- array containing dependent state derivatives
         * H    (``realtype``, input) -- current step size 
         * IPAR (``long int``, input) -- array containing integer user data that was passed to
           :c:func:`FARKMALLOC()` 
         * RPAR (``realtype``, input) -- array containing real user data that was passed to
           :c:func:`FARKMALLOC()` 
         * WORK (``realtype``, input) -- array containing temporary workspace of same size as
           Y   
         * IER  (``int``, output) -- return flag  (0 if success, :math:`\ne 0` if an error)
         
      **Notes:**
      Typically this routine will use only NEQ, T, Y, V, and FJV.  It
      must compute the product vector :math:`Jv`, where :math:`v` is
      given in V, and the product is stored in FJV. 
   


If this routine has been supplied by the user, then, following the
call to :c:func:`FARKSPGMR()`, :c:func:`FARKSPBCG()` or
:c:func:`FARKSPTFQMR()`, the user must call the routine
FARKSPILSSETJAC with FLAG :math:`\ne 0` to specify use of the
user-supplied Jacobian-times-vector function.



   .. c:function:: SUBROUTINE FARKSPILSSETJAC(FLAG, IER)
   
      Interface to the function :c:func:`ARKSpilsSetJacTimesVecFn()` to specify use of the
      user-supplied Jacobian-times-vector function :c:func:`FARKJTIMES()`.
      
      **Arguments:** 
         * FLAG (``int``, input) -- flag denoting to use FARKJTIMES routine 
         * IER  (``int``, output) -- return flag  (0 if success, :math:`\ne 0` if an error)



If preconditioning is to be performed during the Krylov solver
(i.e. the solver was set up with IPRETYPE :math:`\ne 0`), then the
user must also call the routine FARKSPILSSETPREC with FLAG
:math:`\ne 0`. 



   .. c:function:: SUBROUTINE FARKSPILSSETJAC(FLAG, IER)
   
      Interface to the function :c:func:`ARKSpilsSetPreconditioner()` to specify use of the
      user-supplied preconditioner setup and solve functions,
      :c:func:`FARKPSET()` and :c:func:`FARKPSOL()`, respectively.
      
      **Arguments:** 
         * FLAG (``int``, input) -- flag denoting use of user-supplied
           preconditioning routines  
         * IER  (``int``, output) -- return flag  (0 if success, :math:`\ne 0` if an error)
         


In addition, the user must provide the following two routines to
implement the preconditioner setup and solve functions to be used
within the solve.



   .. c:function:: SUBROUTINE FARKPSET(T,Y,FY,JOK,JCUR,GAMMA,H,IPAR,RPAR,V1,V2,V3,IER)
   
      User-supplied preconditioner setup routine (of type
      :c:func:`ARKSpilsPrecSetupFn()`). 
      
      **Arguments:** 
         * T (``realtype``, input) -- current value of the independent variable
         * Y (``realtype``, input) -- current dependent state variable array 
         * FY (``realtype``, input) -- current dependent state variable derivative array 
         * JOK (``int``, input) -- flag indicating whether Jacobian-related data needs to be 
           recomputed:

  	    * 0 = recompute, 
	    * 1 = reuse with the current value of GAMMA.

         * JCUR (``realtype``, output) -- return flag to denote if Jacobian data was recomputed
           (1=yes, 0=no)  
         * GAMMA (``realtype``, input) -- Jacobian scaling factor 
         * H (``realtype``, input) -- current step size 
         * IPAR (``long int``, input/output) -- array containing integer user data that was passed to
           :c:func:`FARKMALLOC()` 
         * RPAR (``realtype``, input/output) -- array containing real user data that was passed to
           :c:func:`FARKMALLOC()` 
         * V1, V2, V3 (``realtype``, input) -- arrays containing temporary workspace of
           same size as Y 
         * IER  (``int``, output) -- return flag  (0 if success, >0 if a recoverable
           failure, <0 if a non-recoverable failure) 
      
      **Notes:**
      This routine must set up the preconditioner P to be used in the
      subsequent call to :c:func:`FARKPSOL()`.  The preconditioner (or
      the product of the left and right preconditioners if using both)
      should be an approximation to the matrix  :math:`M - \gamma J`,
      where :math:`M` is the system mass matrix, :math:`\gamma` is the
      input GAMMA, and :math:`J = \frac{\partial f_I}{\partial y}`. 
   

   
   .. c:function:: SUBROUTINE FARKPSOL(T,Y,FY,R,Z,GAMMA,DELTA,LR,IPAR,RPAR,VT,IER)
   
      User-supplied preconditioner solve routine (of type
      :c:func:`ARKSpilsPrecSolveFn()`). 
      
      **Arguments:** 
         * T (``realtype``, input) -- current value of the independent variable
         * Y (``realtype``, input) -- current dependent state variable array 
         * FY (``realtype``, input) -- current dependent state variable derivative array 
         * R (``realtype``, input) -- right-hand side array 
         * Z (``realtype``, output) -- solution array 
         * GAMMA (``realtype``, input) -- Jacobian scaling factor 
         * DELTA (``realtype``, input) -- desired residual tolerance 
         * LR (``int``, input) -- flag denoting to solve the right or left preconditioner
           system:

            * 1 = left preconditioner
	    * 2 = right preconditioner

         * IPAR (``long int``, input/output) -- array containing integer user data that was passed to
           :c:func:`FARKMALLOC()` 
         * RPAR (``realtype``, input/output) -- array containing real user data that was passed to
           :c:func:`FARKMALLOC()` 
         * VT (``realtype``, input) -- array containing temporary workspace of same size as Y  
         * IER  (``int``, output) -- return flag  (0 if success, >0 if a recoverable
           failure, <0 if a non-recoverable failure) 
      
      **Notes:**
      Typically this routine will use only NEQ, T, Y, GAMMA, R, LR,
      and Z.  It must solve the preconditioner linear system :math:`Pz
      = r`.  The preconditioner (or the product of the left and right
      preconditioners if both are nontrivial) should be an
      approximation to the matrix  :math:`M - \gamma J`, where
      :math:`M` is the system mass matrix, :math:`\gamma` is the input
      GAMMA, and :math:`J = \frac{\partial f_I}{\partial y}`. 



Notes:

(a) If the user's :c:func:`FARKJTIMES()` or :c:func:`FARKPSET()` routine
    uses difference quotient approximations, it may need to use the
    error weight array EWT, the current stepsize H, and/or the
    unit roundoff, in the calculation of suitable increments. Also, If
    :c:func:`FARKPSOL()` uses an iterative method in its solution, the
    residual vector :math:`\rho = r - Pz` of the system should be made
    less than :math:`\delta =` DELTA in the weighted l2 norm, i.e.
    
    .. math::
       \left(\sum_i \left(\rho_i * EWT_i\right)^2 \right)^{1/2} < \delta.

(b) If needed in :c:func:`FARKJTIMES()`, :c:func:`FARKPSOL()`, or
    :c:func:`FARKPSET()`, the error weight array EWT can be
    obtained by calling :c:func:`FARKGETERRWEIGHTS()` using one of the
    work arrays as temporary storage for EWT. 

(c) If needed in :c:func:`FARKJTIMES()`, :c:func:`FARKPSOL()`, or
    :c:func:`FARKPSET()`, the unit roundoff can be obtained as the
    optional output ROUT(6) (available after the call to
    :c:func:`FARKMALLOC()`) and can be passed using either the RPAR
    user data array or a common block. 




.. _FInterface.Solution:

Problem solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Carrying out the integration is accomplished by making calls to
:c:func:`FARKODE()`.



   .. c:function:: SUBROUTINE FARKODE(TOUT, T, Y, ITASK, IER)
   
      Fortran interface to the C routine :c:func:`ARKode()`
      for performing the solve, along with many of the ARK*Get*
      routines for reporting on solver statistics.
      
      **Arguments:** 
         * TOUT (``realtype``, input) -- next value of :math:`t` at which a solution is
           desired 
         * T (``realtype``, output) -- current value of independent variable reached by the solver
         * Y (``realtype``, output) -- array containing dependent state variables on output
         * ITASK (``int``, input) -- task indicator :
	    * 1 = normal mode (overshoot TOUT and interpolate)
	    * 2 = one-step mode (return after each internal step taken)
      	    * 3 = normal `tstop` mode (like 1, but integration never
              proceeds past TSTOP, which must be specified through a
              preceding call to :c:func:`FARKSETRIN()` using the key
              STOP_TIME)
      	    * 4 = one step `tstop` (like 2, but integration never goes past
              TSTOP) 
         * IER (int, output) -- completion flag : 
	    * 0 = success, 
	    * 1 = tstop return, 
	    * 2 = root return, 
	    * values -1 ... -10 are failure modes (see :c:func:`ARKode()` and
              :ref:`Constants`).
      
      **Notes:**
      The current values of the optional outputs are immediately
      available in IOUT and ROUT upon return from this function (see
      :ref:`FInterface.IOUTTable` and :ref:`FInterface.ROUTTable`). 
   


.. _FInterface.AdditionalOutput:

Additional solution output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After a successful return from :c:func:`FARKODE()`, the routine
:c:func:`FARKDKY()` may be used to obtain a derivative of the solution,
of order up to 3, at any :math:`t` within the last step taken. 



   .. c:function:: SUBROUTINE FARKDKY(T, K, DKY, IER)
   
      Fortran interface to the C routine :c:func:`ARKDKY()`
      for interpolating output of the solution or its derivatives at any
      point within the last step taken.
      
      **Arguments:** 
         * T (``realtype``, input) -- time at which solution derivative is desired,
           within the interval :math:`[t_n-h,t_n]`, .
         * K (``int``, input) -- derivative order :math:`(0 \le k \le 3)` 
         * DKY (``realtype``, output) -- array containing the computed K-th derivative of
           :math:`y`  
         * IER (``int``, output) -- return flag (0 if success, <0 if an illegal argument)



.. _FInterface.ReInit:

Problem reinitialization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To re-initialize the ARKode solver for the solution of a new
problem of the same size as one already solved, the user must call
:c:func:`FARKREINIT()`. 



   .. c:function:: SUBROUTINE FARKREINIT(T0, Y0, IMEX, IATOL, RTOL, ATOL, IER)
   
      Re-initializes the Fortran interface to the ARKode solver.
      
      **Arguments:**  The arguments have the same names and meanings as those of
      :c:func:`FARKMALLOC()`.
      
      **Notes:**
      This routine performs no memory allocation, instead using the
      existing memory created by the previous :c:func:`FARKMALLOC()`
      call.  The call to specify the linear system solution method may
      or may not be needed. 


Following a call to :c:func:`FARKREINIT()`, a call to specify the
linear system solver must be made if the choice of linear solver is
being changed. Otherwise, a call to reinitialize the linear solver
last used may or may not be needed, depending on changes in the inputs
to it. 

In the case of the BAND solver, for any change in the
half-bandwidth parameters, call :c:func:`FARKBAND()` (or
:c:func:`FARKLAPACKBAND()`) again described above.

In the case of SPGMR, for a change of inputs other than MAXL,
the user may call the routine :c:func:`FARKSPGMRREINIT()` to
reinitialize SPGMR without reallocating its memory, as follows: 



   .. c:function:: SUBROUTINE FARKSPGMRREINIT(IPRETYPE, IGSTYPE, DELT, IER)
   
      Re-initializes the Fortran interface to the SPGMR
      linear solver.
      
      **Arguments:**  The arguments have the same names and meanings as those of
      :c:func:`FARKSPGMR()`.
   


However, if MAXL is being changed, then the user should call
:c:func:`FARKSPGMR()` instead.

In the case of SPBCG, for a change in any inputs, the user can
reinitialize SPBCG without reallocating its memory by calling
:c:func:`FARKSPBCGREINIT()`, as follows:



   .. c:function:: SUBROUTINE FARKSPBCGREINIT(IPRETYPE, MAXL, DELT, IER)
   
      Re-initializes the Fortran interface to the SPBCG
      linear solver.
      
      **Arguments:**  The arguments have the same names and meanings as those of
      :c:func:`FARKSPBCG()`.



In the case of SPTFQMR, for a change in any inputs, the user can
reinitialize SPTFQMR without reallocating its memory by calling
:c:func:`FARKSPTFQMRREINIT()`, as follows:



   .. c:function:: SUBROUTINE FARKSPTFQMRREINIT(IPRETYPE, MAXL, DELT, IER)
   
      Re-initializes the Fortran interface to the SPBTFQMR
      linear solver.
      
      **Arguments:**  The arguments have the same names and meanings as those of
      :c:func:`FARKSPTFQMR()`.





.. _FInterface.Deallocation:

Memory deallocation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To free the internal memory created by :c:func:`FARKMALLOC()`, the user
may call :c:func:`FARKFREE()`, as follows:



   .. c:function:: SUBROUTINE FARKFREE()
   
      Frees the internal memory created by :c:func:`FARKMALLOC()`.
      
      **Arguments:** None.



.. _FInterface.OptionalOutputs:

FARKODE optional output
-----------------------------

The optional inputs to FARKODE have already been described in the
section :ref:`FInterface.OptionalInputs`.  


IOUT and ROUT arrays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The optional outputs from the :c:func:`ARKode()` solver are accessed
not through individual functions, but rather through a pair of arrays,
IOUT (``long int`` type) of dimension at least 22, and ROUT
(``realtype`` type) of dimension at least 6. These arrays are owned
(and allocated) by the user and are passed as arguments to
:c:func:`FARKMALLOC()`. 

:ref:`FInterface.IOUTTable` and
:ref:`FInterface.ROUTTable` list the entries in these
arrays associated with the main ARKode solver, along with the
relevant ARKode function that is actually called to extract the
optional output.  Similarly,
:ref:`FInterface.DlsIOUTTable` lists the IOUT
entries associated with the main ARKDENSE and ARKBAND direct
linear solvers, and :ref:`FInterface.SpilsIOUTTable`
lists the IOUT entries associated with the main ARKSPGMR,
ARKSPBCG and ARKSPTFQMR iterative linear solvers.

For more details on the optional inputs and outputs to ARKode, see
the sections :ref:`CInterface.OptionalInputs` and
:ref:`CInterface.OptionalOutputs`.



.. _FInterface.IOUTTable:

Table: Optional FARKODE integer outputs
""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

==============  ===============  =====================================================
IOUT Index      Optional output  ARKode function
==============  ===============  =====================================================
1               LENRW            :c:func:`ARKodeGetWorkSpace()`
2               LENIW            :c:func:`ARKodeGetWorkSpace()`
3               NST              :c:func:`ARKodeGetNumSteps()`
4               NST_STB          :c:func:`ARKodeGetNumExpSteps()`
5               NST_ACC          :c:func:`ARKodeGetNumAccSteps()`
6               NST_CNV          :c:func:`ARKodeGetNumConvSteps()`
7               NFE              :c:func:`ARKodeGetNumRhsEvals()` (:math:`f_E` calls)
8               NFI              :c:func:`ARKodeGetNumRhsEvals()` (:math:`f_I` calls)
9               NSETUPS          :c:func:`ARKodeGetNumLinSolvSetups()`
10              NETF             :c:func:`ARKodeGetNumErrTestFails()`
11              NNI              :c:func:`ARKodeGetNumNonlinSolvIters()`
12              NCFN             :c:func:`ARKodeGetNumNonlinSolvConvFails()`
13              NGE              :c:func:`ARKodeGetNumGEvals()`
==============  ===============  =====================================================



.. _FInterface.ROUTTable:

Table: Optional FARKODE real outputs 
""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

==============  ===============  ===================================================
ROUT Index      Optional output  ARKode function
==============  ===============  ===================================================
1               H0U              :c:func:`ARKodeGetActualInitStep()`
2               HU               :c:func:`ARKodeGetLastStep()`
3               HCUR             :c:func:`ARKodeGetCurrentStep()`
4               TCUR             :c:func:`ARKodeGetCurrentTime()`
5               TOLSF            :c:func:`ARKodeGetTolScaleFactor()`
6               UROUND           ``UNIT_ROUNDOFF`` (see :ref:`CInterface.DataTypes`)
==============  ===============  ===================================================



.. _FInterface.DlsIOUTTable:

Table: Optional ARKDENSE and ARKBAND outputs
""""""""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

==============  ===============  ===================================================
IOUT Index      Optional output  ARKode function
==============  ===============  ===================================================
14              LENRWLS          :c:func:`ARKDlsGetWorkSpace()`
15              LENIWLS          :c:func:`ARKDlsGetWorkSpace()`
16              LSTF             :c:func:`ARKDlsGetLastFlag()`
17              NFELS            :c:func:`ARKDlsGetNumRhsEvals()`
18              NJE              :c:func:`ARKDlsGetNumJacEvals()`
==============  ===============  ===================================================



.. _FInterface.SpilsIOUTTable:

Table: Optional ARKSPGMR, ARKSPBCG and ARKSPTFQMR outputs 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

==============  ===============  ===================================================
IOUT Index      Optional output  ARKode function
==============  ===============  ===================================================
14              LENRWLS          :c:func:`ARKSpilsGetWorkSpace()`
15              LENIWLS          :c:func:`ARKSpilsGetWorkSpace()`
16              LSTF             :c:func:`ARKSpilsGetLastFlag()`
17              NFELS            :c:func:`ARKSpilsGetNumRhsEvals()`
18              NJTV             :c:func:`ARKSpilsGetNumJtimesEvals()`
19              NPE              :c:func:`ARKSpilsGetNumPrecEvals()`
20              NPS              :c:func:`ARKSpilsGetNumPrecSolves()`
21              NLI              :c:func:`ARKSpilsGetNumLinIters()`
22              NCFL             :c:func:`ARKSpilsGetNumConvFails()`
==============  ===============  ===================================================



Additional optional output routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


In addition to the optional inputs communicated through FARKSET*
calls and the optional outputs extracted from IOUT and ROUT,
the following user-callable routines are available: 

To obtain the error weight array EWT, containing the
multiplicative error weights used the WRMS norms, the user may call
the routine :c:func:`FARKGETERRWEIGHTS()` as follows:



   .. c:function:: SUBROUTINE FARKGETERRWEIGHTS(EWT, IER)
   
      Retrieves the current error weight vector (interfaces
      with :c:func:`ARKodeGetErrWeights()`).
      
      **Arguments:** 
         * EWT (``realtype``, output) -- array containing the error weight vector
         * IER  (``int``, output) -- return flag  (0 if success, :math:`\ne 0` if an error)
      
      **Notes:**
      The array EWT, of length NEQ if using NVECTOR_SERIAL or NLOCAL
      if using NVECTOR_PARALLEL, must already have been declared by
      the user.



Similarly, to obtain the estimated local errors, following a
successful call to :c:func:`FARKODE()`, the user may call the routine
:c:func:`FARKGETESTLOCALERR()` as follows:



   .. c:function:: SUBROUTINE FARKGETESTLOCALERR(ELE, IER)
   
      Retrieves the current local truncation error estimate
      vector (interfaces with :c:func:`ARKodeGetEstLocalErrors()`).
      
      **Arguments:** 
         * ELE (``realtype``, output) -- array with the estimated local error vector
         * IER  (``int``, output) -- return flag  (0 if success, :math:`\ne 0` if an error)
      
      **Notes:**
      The array ELE, of length NEQ if using NVECTOR_SERIAL or NLOCAL
      if using NVECTOR_PARALLEL, must already have been declared by
      the user.  








.. _FInterface.RootFinding:

Usage of the FARKROOT interface to rootfinding
-----------------------------------------------

(to be added)








.. _FInterface.BandPre:

Usage of the FARKBP interface to ARKBANDPRE
-----------------------------------------------

The FARKBP interface sub-module is a package of C functions which,
as part of the FARKODE interface module, support the use of the
ARKode solver with the serial NVECTOR_SERIAL module, and the
combination of the ARKBANDPRE preconditioner module (see the
section :ref:`CInterface.BandPre`) with any of the Krylov iterative
linear solvers. 

The two user-callable functions in this package, with the
corresponding ARKode function around which they wrap, are: 

* :c:func:`FARKBPINIT()` interfaces to :c:func:`ARKBandPrecInit()`.

* :c:func:`FARKBPOPT()` interfaces to the ARKBANDPRE optional output
  functions, :c:func:`ARKBandPrecGetWorkSpace()` and
  :c:func:`ARKBandPrecGetNumRhsEvals()`. 

As with the rest of the FARKODE routines, the names of the
user-supplied routines are mapped to actual values through a series of
definitions in the header file ``farkbp.h``. 

The following is a summary of the usage of this module. Steps that are
unchanged from the main program described in the section
:ref:`FInterface.Usage` are `italicized`.


1. `Right-hand side specification`

2. `NVECTOR module initialization`

3. `Problem specification`

4. `Set optional inputs`

5. Linear solver specification 

   First, specify one of the ARKSPILS iterative linear solvers, by
   calling one of :c:func:`FARKSPGMR()`, :c:func:`FARKSPBCG()`, or
   :c:func:`FARKSPTFQMR()`. 

   Optionally, to specify that SPGMR, SPBCG, or SPTFQMR
   should use the supplied :c:func:`FARKJTIMES()` routine, the user
   should call :c:func:`FARKSPILSSETJAC()` with FLAG :math:`\ne 0`,
   as described in the section :ref:`FInterface.SpilsUserSupplied`.

   Then, to initialize the ARKBANDPRE preconditioner, call the
   routine :c:func:`FARKBPINIT()`, as follows:



      .. c:function:: SUBROUTINE FARKBPINIT(NEQ, MU, ML, IER)
   
         Interfaces with the :c:func:`ARKBandPrecInit()`
         function to allocates memory and initialize data associated
         with the ARKBANDPRE preconditioner.
   
         **Arguments:** 
	    * NEQ (``long int``, input) -- problem size 
            * MU (``long int``, input) -- upper half-bandwidth of the band matrix that is 
              retained as an approximation of the Jacobian 
            * ML  (``long int``, input) -- lower half-bandwidth of the band matrix approximant 
              to the Jacobian 	  
            * IER  (``int``, output) -- return flag  (0 if success, -1 if a memory failure)
            


6. `Problem solution`

7. ARKBANDPRE optional outputs 

   Optional outputs specific to the SPGMR, SPBCG, or
   SPTFQMR solver are listed in :ref:`FInterface.SpilsIOUTTable`. 
   To obtain the optional outputs associated with the ARKBANDPRE
   module, the user should call the :c:func:`FARKBPOPT()`, as specified
   below: 



      .. c:function:: SUBROUTINE FARKBPOPT(LENRWBP, LENIWBP, NFEBP)
      
         Interfaces with the ARKBANDPRE optional output
         functions.
         
         **Arguments:** 
	    * LENRWBP (``long int``, output) -- length of real preconditioner work
              space (from :c:func:`ARKBandPrecGetWorkSpace()`)  
            * LENIWBP (``long int``, output) -- length of integer preconditioner work space, in 
              integer words (from :c:func:`ARKBandPrecGetWorkSpace()`)  
            * NFEBP (``long int``, output) -- number of :math:`f_I(t,y)` evaluations (from
              :c:func:`ARKBandPrecGetNumRhsEvals()`)  



8. `Memory deallocation` 

   (The memory allocated for the FARKBP module is deallocated
   automatically by :c:func:`FARKFREE()`)






.. _FInterface.BBDPre:

Usage of the FARKBBD interface to ARKBBDPRE
-----------------------------------------------

The FARKBBD interface sub-module is a package of C functions which, as
part of the FARKODE interface module, support the use of the ARKode
solver with the parallel NVECTOR_PARALLEL module, and the combination
of the ARKBBDPRE preconditioner module (see the section
:ref:`CInterface.BBDPre`) with any of the Krylov iterative linear
solvers. 

The user-callable functions in this package, with the corresponding
ARKode and ARKBBDPRE functions, are as follows:

* :c:func:`FARKBBDINIT()` interfaces to :c:func:`ARKBBDPrecInit()`.

* :c:func:`FARKBBDREINIT()` interfaces to :c:func:`ARKBBDPrecReInit()`.

* :c:func:`FARKBBDOPT()` interfaces to the ARKBBDPRE optional output
  functions.

In addition to the Fortran right-hand side function
:c:func:`FARKFUN()`, the user-supplied functions used by this package
are listed in the table below, each with the
corresponding interface function which calls it (and its type within
ARKBBDPRE or ARKode).


*Table: FARKBBD function mapping*

.. cssclass:: table-bordered

+--------------------------+------------------------+-----------------------------------+
| FARKBBD routine          | ARKode routine         | ARKode interface                  |
| (FORTRAN, user-supplied) | (C, interface)         | function type                     |
+==========================+========================+===================================+
| :c:func:`FARKJTIMES()`   | FARKJtimes             | :c:func:`ARKSpilsJacTimesVecFn()` |
+--------------------------+------------------------+-----------------------------------+
| :c:func:`FARKLOCFN()`    | FARKgloc               | :c:func:`ARKLocalFn()`            |
+--------------------------+------------------------+-----------------------------------+
| :c:func:`FARKCOMMF()`    | FARKcfn                | :c:func:`ARKCommFn()`             |
+--------------------------+------------------------+-----------------------------------+

As with the rest of the FARKODE routines, the names of all
user-supplied routines here are fixed, in order to maximize
portability for the resulting mixed-language program. Additionally,
based on flags discussed above in the section :ref:`FInterface.Routines`,
the names of the user-supplied routines are mapped to actual values
through a series of definitions in the header file ``farkbbd.h``. 

The following is a summary of the usage of this module. Steps that are
unchanged from the main program described in the section
:ref:`FInterface.Usage` are `italicized`. 

1. `Right-hand side specification`

2. `NVECTOR module initialization`

3. `Problem specification`

4. `Set optional inputs`

5. Linear solver specification 

   First, specify one of the ARKSPILS iterative linear solvers, by
   calling one of :c:func:`FARKSPGMR()`, :c:func:`FARKSPBCG()`, or
   :c:func:`FARKSPTFQMR()`.  

   Optionally, to specify that SPGMR, SPBCG, or SPTFQMR
   should use the supplied :c:func:`FARKJTIMES()` routine, the user
   should call :c:func:`FARKSPILSSETJAC()` with FLAG :math:`\ne 0`,
   as described in the section :ref:`FInterface.SpilsUserSupplied`.

   Then, to initialize the ARKBBDPRE preconditioner, call the function
   :c:func:`FARKBBDINIT()`, as described below:



      .. c:function:: SUBROUTINE FARKBBDINIT(NLOCAL, MUDQ, MLDQ, MU, ML, DQRELY, IER)
      
         Interfaces with the :c:func:`ARKBBDPrecInit()`
         routine to initialize the ARKBBDPRE preconditioning module.
         
         **Arguments:** 
	    * NLOCAL (``long int``, input) -- local vector size on this process
   	    * MUDQ (``long int``, input) -- upper half-bandwidth to be
   	      used in the computation of the local Jacobian blocks by
   	      difference quotients.  These may be smaller than the
   	      true half-bandwidths of the Jacobian of the local block
   	      of :math:`g`, when smaller values may provide greater efficiency  
	    * MLDQ (``long int``, input) -- lower half-bandwidth to be used in the computation
              of the local Jacobian blocks by difference quotients
	    * MU (``long int``, input) -- upper half-bandwidth of the band matrix that is
              retained as an approximation of the local Jacobian block (may be smaller than MUDQ)  
	    * ML (``long int``, input) -- lower half-bandwidth of the band matrix that is
              retained as an approximation of the local Jacobian block (may be smaller than MLDQ)  
	    * DQRELY (``realtype``, input) -- relative increment factor in :math:`y` for
              difference quotients (0.0 indicates to use the default)
            * IER  (``int``, output) -- return flag  (0 if success, -1 if a memory
              failure) 



6. `Problem solution`

7. ARKBBDPRE optional outputs

   Optional outputs specific to the SPGMR, SPBCG, or SPTFQMR solver
   are listed in :ref:`FInterface.SpilsIOUTTable`.  To obtain the
   optional outputs associated with the ARKBBDPRE module, the user
   should call the :c:func:`FARKBBDOPT()`, as specified below:



      .. c:function:: SUBROUTINE FARKBBDOPT(LENRWBBD, LENIWBBD, NGEBBD)
      
         Interfaces with the ARKBBDPRE optional output
         functions.
         
         **Arguments:** 
	    * LENRWBP (``long int``, output) -- length of real preconditioner work
              space on this process (from :c:func:`ARKBBDPrecGetWorkSpace()`)  
            * LENIWBP (``long int``, output) -- length of integer preconditioner work space on
              this process (from :c:func:`ARKBBDPrecGetWorkSpace()`)
            * NGEBBD (``long int``, output) -- number of :math:`g(t,y)` evaluations (from
              :c:func:`ARKBBDPrecGetNumGfnEvals()`) so far  



8. Problem reinitialization

   If a sequence of problems of the same size is being solved using
   the same linear solver (SPGMR, SPBCG, or SPTFQMR) in combination
   with the ARKBBDPRE preconditioner, then the ARKode package can be
   re-initialized for the second and subsequent problems by calling
   :c:func:`FARKREINIT()`, following which a call to
   :c:func:`FARKBBDREINIT()` may or may not be needed. If the input
   arguments are the same, no :c:func:`FARKBBDREINIT()` call is
   needed.

   If there is a change in input arguments other than MU or
   ML, then the user program should call :c:func:`FARKBBDREINIT()` as
   specified beloe: 



      .. c:function:: SUBROUTINE FARKBBDREINIT(NLOCAL, MUDQ, MLDQ, DQRELY, IER)
      
         Interfaces with the
         :c:func:`ARKBBDPrecReInit()` function to reinitialize the
         ARKBBDPRE module.
         
         **Arguments:**  The arguments of the same names have the same
	 meanings as in :c:func:`FARKBBDINIT()`.



   However, if the value of MU or ML is being changed, then a call to
   :c:func:`FARKBBDINIT()` must be made instead. 

   Finally, if there is a change in any of the linear solver inputs,
   then a call to FARKSPGMR, FARKSPBCG, or FARKSPTFQMR must also be
   made; in this case the linear solver memory is reallocated. 

9. `Memory deallocation` 

   (The memory allocated for the FARKBBD module is deallocated
   automatically by :c:func:`FARKFREE()`) 

10. User-supplied routines 

    The following two routines must be supplied for use with the
    ARKBBDPRE module:



      .. c:function:: SUBROUTINE FARKGLOCFN(NLOC, T, YLOC, GLOC, IPAR, RPAR, IER)
      
         User-supplied routine (of type :c:func:`ARKLocalFn()`) that
	 computes a processor-local approximation :math:`g(t,y)` to
	 the right-hand side function :math:`f_I(t,y)`.
         
         **Arguments:** 
	    * NLOC (``long int``, input) -- local problem size 
            * T (``realtype``, input) -- current value of the independent variable
	    * YLOC (``realtype``, input) -- array containing local dependent state variables
	    * GLOC (``realtype``, output) -- array containing local dependent state derivatives
            * IPAR (``long int``, input/output) -- array containing integer user data that was passed to
              :c:func:`FARKMALLOC()` 
            * RPAR (``realtype``, input/output) -- array containing real user data that was passed to
              :c:func:`FARKMALLOC()` 
            * IER (``int``, output) -- return flag (0 if success, >0 if a recoverable error
              occurred, <0 if an unrecoverable error occurred) 



      .. c:function:: SUBROUTINE FARKCOMMFN(NLOC, T, YLOC, IPAR, RPAR, IER)
      
         User-supplied routine (of type
	 :c:func:`ARKCommFn()`) that performs all interprocess
         communication necessary for the executation of the
	 :c:func:`FARKGLOCFN()` function above, using the input vector
	 YLOC.
         
         **Arguments:** 
            * NLOC (``long int``, input) -- local problem size 
	    * T (``realtype``, input) -- current value of the independent variable
	    * YLOC (``realtype``, input) -- array containing local dependent state variables
            * IPAR (``long int``, input/output) -- array containing integer user data that was passed to
              :c:func:`FARKMALLOC()` 
            * RPAR (``realtype``, input/output) -- array containing real user data that was passed to
              :c:func:`FARKMALLOC()` 
            * IER (``int``, output) -- return flag (0 if success, >0 if a recoverable error
              occurred, <0 if an unrecoverable error occurred) 

         **Notes:**
	 The subroutine FARKCOMMFN must be supplied even if it is not
	 needed and must return IER=0.  



