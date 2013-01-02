.. _FortranInterface:

FARKode, an Interface Module for FORTRAN Applications
=====================================================

The ``FARKODE`` interface module is a package of C functions which
support the use of the ``ARKODE`` solver, for the solution of ODE
systems 

.. math::
   M \dot{y} = f_E(t,y) + f_I(t,y),

in a mixed Fortran/C setting. While ``ARKODE`` is written in C, it is
assumed here that the user's calling program and user-supplied
problem-defining routines are written in Fortran. This package
provides the necessary interface to ``ARKODE`` for both the serial and
the parallel ``NVECTOR`` implementations.


.. _FInterface.Routines:

FARKODE routines
----------------

The user-callable functions, with the corresponding ``ARKODE``
functions, are as follows:

- Interface to the ``NVECTOR`` modules

  - :ref:`Fcn.FNVINITS` (defined by ``NVECTOR_SERIAL``) interfaces to
    :ref:`Fcn.N_VNewEmpty_Serial`.

  - :ref:`Fcn.FNVINITP` (defined by ``NVECTOR_PARALLEL``) interfaces to
    :ref:`Fcn.N_VNewEmpty_Parallel`. 

- Interface to the main ``ARKODE`` module

  - :ref:`Fcn.FARKMALLOC` interfaces to :ref:`Fcn.ARKodeCreate`,
    :ref:`Fcn.ARKodeSetUserData`, and :ref:`Fcn.ARKodeInit`, as well
    as one of :ref:`Fcn.ARKodeSStolerances` or :ref:`Fcn.ARKodeSVtolerances`.

  - :ref:`Fcn.FARKREINIT` interfaces to :ref:`Fcn.ARKodeReInit`.

  - :ref:`Fcn.FARKSETIIN` and :ref:`Fcn.FARKSETRIN` interface to the
    ``ARKodeSet*`` functions (see :ref:`CInterface.OptionalInputs`).

  - :ref:`Fcn.FARKEWTSET` interfaces to :ref:`Fcn.ARKodeWFtolerances`.

  - :ref:`Fcn.FARKODE` interfaces to :ref:`Fcn.ARKode`, the
    ``ARKodeGet*`` functions (see :ref:`CInterface.OptionalOutputs`), 
    and to the optional output functions for the selected linear
    solver module (see :ref:`CInterface.OptionalOutputs`). 

  - :ref:`Fcn.FARKDKY` interfaces to the interpolated output function
    :ref:`Fcn.ARKodeGetDky`.

  - :ref:`Fcn.FARKGETERRWEIGHTS` interfaces to
    :ref:`Fcn.ARKodeGetErrWeights`.

  - :ref:`Fcn.FARKGETESTLOCALERR` interfaces to
    :ref:`Fcn.ARKodeGetEstLocalErrors`.

  - :ref:`Fcn.FARKFREE` interfaces to :ref:`Fcn.ARKodeFree`.

- Interface to the linear solver modules

  - :ref:`Fcn.FARKDENSE` interfaces to :ref:`Fcn.ARKDense`.

  - :ref:`Fcn.FARKDENSESETJAC` interfaces to :ref:`Fcn.ARKDlsSetDenseJacFn`.

  - :ref:`Fcn.FARKLAPACKDENSE` interfaces to :ref:`Fcn.ARKLapackDense`.

  - :ref:`Fcn.FARKLAPACKDENSESETJAC` interfaces to :ref:`Fcn.ARKDlsSetDenseJacFn`.

  - :ref:`Fcn.FARKBAND` interfaces to :ref:`Fcn.ARKBand`.

  - :ref:`Fcn.FARKBANDSETJAC` interfaces to :ref:`Fcn.ARKDlsSetBandJacFn`.

  - :ref:`Fcn.FARKLAPACKBAND` interfaces to :ref:`Fcn.ARKLapackBand`.

  - :ref:`Fcn.FARKLAPACKBANDSETJAC` interfaces to :ref:`Fcn.ARKDlsSetBandJacFn`.

  - :ref:`Fcn.FARKSPGMR` interfaces to :ref:`Fcn.ARKSpgmr` and the ``SPGMR`` optional input
    functions (see :ref:`CInterface.ARKSpilsInputTable`).

  - :ref:`Fcn.FARKSPGMRREINIT` interfaces to the ``SPGMR`` optional input
    functions (see :ref:`CInterface.ARKSpilsInputTable`).

  - :ref:`Fcn.FARKSPBCG` interfaces to :ref:`Fcn.ARKSpbcg` and the ``SPBCG`` optional input
    functions (see :ref:`CInterface.ARKSpilsInputTable`).

  - :ref:`Fcn.FARKSPBCGREINIT` interfaces to the ``SPBCG`` optional input
    functions.

  - :ref:`Fcn.FARKSPTFQMR` interfaces to :ref:`Fcn.ARKSptfqmr` and the ``SPTFQMR`` optional
    input functions.

  - :ref:`Fcn.FARKSPTFQMRREINIT` interfaces to the ``SPTFQMR`` optional input
    functions.

  - :ref:`Fcn.FARKSPILSSETJAC` interfaces to :ref:`Fcn.ARKSpilsSetJacTimesVecFn`.

  - :ref:`Fcn.FARKSPILSSETPREC` interfaces to :ref:`Fcn.ARKSpilsSetPreconditioner`.


The user-supplied functions, each listed with the corresponding
internal interface function which calls it (and its type within
``ARKode``), are as follows:

+--------------------------+------------------------+----------------------------------+
| ``FARKODE`` routine      | ``ARKode`` routine     | ``ARKode`` interface             |
| (FORTRAN, user-supplied) | (C, interface)         | function type                    |
+==========================+========================+==================================+
| :ref:`Fcn.FARKIFUN`      | ``FARKfi``             | :ref:`Fcn.ARKRhsFn`              |
+--------------------------+------------------------+----------------------------------+
| :ref:`Fcn.FARKEFUN`      | ``FARKfe``             | :ref:`Fcn.ARKRhsFn`              |
+--------------------------+------------------------+----------------------------------+
| :ref:`Fcn.FARKDJAC`      | ``FARKDenseJac``       | :ref:`Fcn.ARKDlsDenseJacFn`      |
+--------------------------+------------------------+----------------------------------+
| :ref:`Fcn.FARKLDJAC`     | ``FARKLapackDenseJac`` | :ref:`Fcn.ARKDlsDenseJacFn`      |
+--------------------------+------------------------+----------------------------------+
| :ref:`Fcn.FARKBJAC`      | ``FARKBandJac``        | :ref:`Fcn.ARKDlsBandJacFn`       |
+--------------------------+------------------------+----------------------------------+
| :ref:`Fcn.FARKLBJAC`     | ``FARKLapackBandJac``  | :ref:`Fcn.ARKDlsBandJacFn`       |
+--------------------------+------------------------+----------------------------------+
| :ref:`Fcn.FARKPSET`      | ``FARKPSet``           | :ref:`Fcn.ARKSpilsPrecSetupFn`   |
+--------------------------+------------------------+----------------------------------+
| :ref:`Fcn.FARKPSOL`      | ``FARKPSol``           | :ref:`Fcn.ARKSpilsPrecSolveFn`   |
+--------------------------+------------------------+----------------------------------+
| :ref:`Fcn.FARKJTIMES`    | ``FARKJtimes``         | :ref:`Fcn.ARKSpilsJacTimesVecFn` |
+--------------------------+------------------------+----------------------------------+
| :ref:`Fcn.FARKEWT`       | ``FARKEwtSet``         | :ref:`Fcn.ARKEwtFn`              |
+--------------------------+------------------------+----------------------------------+

In contrast to the case of direct use of ``ARKode``, and of most
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
to their usage in ``SUNDIALS``.  The equivalent types to these may
vary, depending on your computer architecture and on how ``SUNDIALS``
was compiled (see :ref:`Installation`).  A Fortran user should take
care that all arguments passed through this Fortran/C interface  are
declared of the appropriate type. 

**Integers**: ``SUNDIALS`` uses both ``int`` and ``long int`` types:

   ``int`` -- equivalent to an ``INTEGER`` or ``INTEGER*4`` in Fortran

   ``long int`` -- this will depend on the computer architecture:
   
      32-bit -- equivalent to an ``INTEGER`` or ``INTEGER*4`` in Fortran

      64-bit -- equivalent to an ``INTEGER*8`` in Fortran
	      
**Real numbers**:  As discussed in :ref:`Installation`, at compilation
``SUNDIALS`` allows the configuration option  ``--with-precision``,
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

The usage of ``FARKODE`` requires calls to five or more interface
functions, depending on the method options selected, and one or more
user-supplied routines which define the problem to be solved.  These 
function calls and user routines are summarized separately below.
Some details are omitted, and the user is referred to the description
of the corresponding ``ARKode`` functions for complete information on
the arguments of any given user-callable interface routine.  The usage
of ``FARKODE`` for rootfinding and with preconditioner modules is
described in later subsections.

Steps marked [**S**] in the instructions below apply to the serial
``NVECTOR`` implementation (``NVECTOR_SERIAL``) only, while those
marked with a [**P**] apply to ``NVECTOR_PARALLEL``.


.. _FInterface.RHS:

Right-hand side specification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user must in all cases supply at least one of the following Fortran 
routines:

.. _Fcn.FARKIFUN:

``FARKIFUN``
"""""""""""""

:Definition: ``SUBROUTINE FARKIFUN(T, Y, YDOT, IPAR, RPAR, IER)``

:Description:  Sets the ``YDOT`` array to :math:`f_I(t,y)``, the
   implicit portion of the right-hand side of the ODE system, as
   function of the independent variable ``T`` :math:`=t` and the array
   of dependent state variables ``Y`` :math:`=y`.

:Arguments: ``Y`` -- array containing state variables [``realtype``,
   input] 
       
   ``YDOT`` -- array containing state derivatives [``realtype``,
   output]

   ``IPAR`` -- array containing integer user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``long int``, input] 

   ``RPAR`` -- array containing real user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``realtype``, input] 

   ``IER`` -- return flag (0 success, >0 recoverable error, <0
   unrecoverable error) [``int``, output]


.. _Fcn.FARKEFUN:

``FARKEFUN``
"""""""""""""

:Definition: ``SUBROUTINE FARKEFUN(T, Y, YDOT, IPAR, RPAR, IER)``

:Description:  Sets the ``YDOT`` array to :math:`f_E(t,y)``, the
   explicit portion of the right-hand side of the ODE system, as
   function of the independent variable ``T`` :math:`=t` and the array
   of dependent state variables ``Y`` :math:`=y`.

:Arguments: ``Y`` -- array containing state variables [``realtype``,
   input] 
       
   ``YDOT`` -- array containing state derivatives [``realtype``,
   output]

   ``IPAR`` -- array containing integer user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``long int``, input] 

   ``RPAR`` -- array containing real user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``realtype``, input] 

   ``IER`` -- return flag (0 success, >0 recoverable error, <0
   unrecoverable error)  [``int``, output]


.. _FInterface.NVector:

``NVECTOR`` module initialization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[**S**] To initialize the serial ``NVECTOR`` module, the user must
call the function ``FNVINITS`` with the argument ``KEY = 4``.

.. _Fcn.FNVINITS:

``FNVINITS``
"""""""""""""

:Definition: ``SUBROUTINE FNVINITS(KEY, NEQ, IER)``

:Description:  Initializes the Fortran interface to the serial
   ``NVECTOR`` module.

:Arguments: ``KEY`` -- integer flag denoting which solver is to be
   used (1 is ``CVODE``, 2 is ``IDA``, 3 is ``KINSOL`` and 4 is
   ``ARKode``) [``int``, input]
       
   ``NEQ`` -- size of the ODE system [``long int``, input]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]



[**P**] To initialize the parallel ``NVECTOR`` module, the user must
call the function ``FNVINITP`` with the argument ``KEY = 4``.

.. _Fcn.FNVINITP:

``FNVINITP``
"""""""""""""

:Definition: ``SUBROUTINE FNVINITP(COMM, KEY, NLOCAL, NGLOBAL, IER)``

:Description:  Initializes the Fortran interface to the parallel
   ``NVECTOR`` module.

:Arguments: ``COMM`` -- the MPI communicator

   ``KEY`` -- integer flag denoting which solver is to be
   used (1 is ``CVODE``, 2 is ``IDA``, 3 is ``KINSOL`` and 4 is
   ``ARKode``) [``int``, input]
       
   ``NLOCAL`` -- local size of vectors on this processor [``long int``, input]

   ``NGLOBAL`` -- the size of the ODE system, and the global size of
   vectors (the sum of all values of ``NLOCAL``) [``long int``, input]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]

:Notes: If the header file ``sundials_config.h`` defines
   ``SUNDIALS_MPI_COMM_F2C`` to be 1 (meaning the MPI implementation
   used to build ``SUNDIALS`` includes the ``MPI_Comm_f2c`` function),
   then ``COMM`` can be any valid MPI communicator.  Otherwise,
   ``MPI_COMM_WORLD`` will be used, so the user can just pass an
   integer value as a placeholder.



.. _FInterface.Problem:

Problem specification
^^^^^^^^^^^^^^^^^^^^^^^

To set various problem and solution parameters and allocate internal
memory, the user must call ``FARKMALLOC``.

.. _Fcn.FARKMALLOC:

``FARKMALLOC``
""""""""""""""""

:Definition: ``SUBROUTINE FARKMALLOC(T0, Y0, IMEX, IATOL, RTOL, ATOL,
   IOUT, ROUT, IPAR, RPAR, IER)`` 

:Description:  Initializes the Fortran interface to the ``ARKode``
   solver, providing interfaces to the C routines :ref:`Fcn.ARKodeCreate`,
   :ref:`Fcn.ARKodeSetUserData`, and :ref:`Fcn.ARKodeInit`, as well
   as one of :ref:`Fcn.ARKodeSStolerances` or
   :ref:`Fcn.ARKodeSVtolerances`.

:Arguments: ``T0`` -- initial value of :math:`t` [``realtype``, input]

   ``Y0`` -- array of initial conditions [``realtype``, input]

   ``IMEX`` -- flag denoting basic integration method [``int``, input]:

      0 = implicit, 

      1 = explicit, 

      2 = imex.

   ``IATOL`` -- type for absolute tolerance input ``ATOL`` [``int``, input]:

      1 = scalar, 

      2 = array,

      3 = user-supplied function; the user must subsequently call
      :ref:`Fcn.FARKEWTSET` and supply a routine :ref:`Fcn.FARKEWT` to
      compute the error weight vector.

   ``RTOL`` -- scalar relative tolerance [``realtype``, input]

   ``ATOL`` -- scalar or array absolute tolerance [``realtype``,
   input]

   ``IOUT`` -- array of length 22 for integer optional outputs [``long
   int``, input/output]

   ``ROUT`` -- array of length 6 for real optional outputs
   [``realtype``, input/output] 

   ``IPAR`` -- array of user integer data, which will be passed
   unmodified to all user-provided routines [``long int``, input/output]

   ``RPAR`` -- array with user real data, which will be passed
   unmodified to all user-provided routines [``realtype``, input/output]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]

:Notes: Modifications to the user data arrays ``IPAR`` and ``RPAR``
   inside a user-provided routine will be propagated to all subsequent
   calls to such routines. The optional outputs associated with the
   main ``ARKode`` integrator are listed in
   :ref:`FInterface.OptionalIntegerOutputTable` and
   :ref:`FInterface.OptionalRealOutputTable`, in the section
   :ref:`FInterface.OptionalOutputs`.


As an alternative to providing tolerances in the call to
:ref:`Fcn.FARKMALLOC`, the user may provide a routine to compute the
error weights used in the WRMS norm evaluations.  If supplied, it must
have the following form:

.. _Fcn.FARKEWT:

``FARKEWT``
"""""""""""""

:Definition: ``SUBROUTINE FARKEWT(Y, EWT, IPAR, RPAR, IER)``

:Description:  It must set the positive components of the error weight
   vector ``EWT`` for the calculation of the WRMS norm of ``Y``.

:Arguments: ``Y`` -- array containing state variables [``realtype``,
   input] 
       
   ``EWT`` -- array containing the error weight vector [``realtype``,
   output] 

   ``IPAR`` -- array containing the integer user data that was passed
   to :ref:`Fcn.FARKMALLOC` [``long int``, input]

   ``RPAR`` -- array containing the real user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``realtype``, input]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]


If the ``FARKEWT`` routine is provided, then, following the call to
:ref:`Fcn.FARKMALLOC`, the user must call the function ``FARKEWTSET``.

.. _Fcn.FARKEWTSET:

``FARKEWTSET``
""""""""""""""""

:Definition: ``SUBROUTINE FARKEWTSET(FLAG, IER)``

:Description:  Informs ``FARKODE`` to use the user-supplied
   :ref:`Fcn.FARKEWT` function.

:Arguments: ``FLAG`` -- integer flag, use "1" to denoting to use ``FARKEWT``.

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]



.. _FInterface.OptionalInputs:

Set optional inputs
^^^^^^^^^^^^^^^^^^^^^^^

To set desired optional inputs, the user can call the routines
:ref:`Fcn.FARKSETIIN` and :ref:`Fcn.FARKSETRIN`, as described below.

.. _Fcn.FARKSETIIN:

``FARKSETIIN``
""""""""""""""""

:Definition: ``SUBROUTINE FARKSETIIN(KEY, IVAL, IER)``

:Description:  Specification routine to pass optional integer inputs
   to the :ref:`Fcn.FARKODE` solver.

:Arguments: ``KEY`` -- quoted string indicating which optional input
   is set (see :ref:`FInterface.IINOptionTable`).

   ``IVAL`` -- the integer input value to be used [``long int``, input]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]


.. _FInterface.IINOptionTable:

Table: Keys for setting ``FARKODE`` integer optional inputs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

=====================  ===================================
Key                    ``ARKode`` routine
=====================  ===================================
``ORDER``              :ref:`Fcn.ARKodeSetOrder`
``DENSE_ORDER``        :ref:`Fcn.ARKodeSetDenseOrder`
``LINEAR``             :ref:`Fcn.ARKodeSetLinear`
``NONLINEAR``          :ref:`Fcn.ARKodeSetNonlinear`
``EXPLICIT``           :ref:`Fcn.ARKodeSetExplicit`
``IMPLICIT``           :ref:`Fcn.ARKodeSetImplicit`
``IMEX``               :ref:`Fcn.ARKodeSetImEx`
``IRK_TABLE_NUM``      :ref:`Fcn.ARKodeSetIRKTableNum`
``ERK_TABLE_NUM``      :ref:`Fcn.ARKodeSetERKTableNum`
``ARK_TABLE_NUM`` `*`  :ref:`Fcn.ARKodeSetARKTableNum`      
``MAX_NSTEPS``         :ref:`Fcn.ARKodeSetMaxNumSteps`
``HNIL_WARNS``         :ref:`Fcn.ARKodeSetMaxHnilWarns`
``PREDICT_METHOD``     :ref:`Fcn.ARKodeSetPredictorMethod`
``MAX_ERRFAIL``        :ref:`Fcn.ARKodeSetMaxErrTestFails`
``MAX_NITERS``         :ref:`Fcn.ARKodeSetMaxNonlinIters`
``MAX_CONVFAIL``       :ref:`Fcn.ARKodeSetMaxConvFails`
=====================  ===================================

`*` When setting ``ARK_TABLE_NUM``, pass in ``IVAL`` as an array of
length 2, specifying the IRK table number first, then the ERK table
number. 


.. _Fcn.FARKSETRIN:

``FARKSETRIN``
""""""""""""""""

:Definition: ``SUBROUTINE FARKSETRIN(KEY, RVAL, IER)``

:Description:  Specification routine to pass optional real inputs
   to the :ref:`Fcn.FARKODE` solver.

:Arguments: ``KEY`` -- quoted string indicating which optional input
   is set (see :ref:`FInterface.RINOptionTable`).

   ``RVAL`` -- the real input value to be used [``realtype``, input]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]


.. _FInterface.RINOptionTable:

Table: Keys for setting ``FARKODE`` real optional inputs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

===============  ===================================
Key              ``ARKode`` routine
===============  ===================================
``INIT_STEP``    :ref:`Fcn.ARKodeSetInitStep`
``MAX_STEP``     :ref:`Fcn.ARKodeSetMaxStep`
``MIN_STEP``     :ref:`Fcn.ARKodeSetMinStep`
``STOP_TIME``    :ref:`Fcn.ARKodeSetStopTime`
``NLCONV_COEF``  :ref:`Fcn.ARKodeSetNonlinConvCoef`
===============  ===================================



Alternatively, if a user wishes to reset all of the options to their
default values, they may call the routine ``FARKSETDEFAULTS``.

.. _Fcn.FARKSETDEFAULTS:

``FARKSETDEFAULTS``
""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSETDEFAULTS(IER)``

:Description:  Specification routine to reset all ``FARKODE`` optional
   inputs to their default values.

:Arguments: ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]


``FARKODE`` supplies additional routines to specify optional advanced
inputs to the :ref:`Fcn.ARKode` solver.  These are summarized below,
and the user is referred to their C routine counterparts for more
complete information. 


.. _Fcn.FARKSETERKTABLE:

``FARKSETERKTABLE``
""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSETERKTABLE(S, Q, P, C, A, B, BEMBED, IER)``

:Description:  Interface to the routine :ref:`Fcn.ARKodeSetERKTable`.

:Arguments: ``S`` -- number of stages in the table [``int``, input]

   ``Q`` -- global order of accuracy of the method [``int``, input]

   ``P`` -- global order of accuracy of the embedding [``int``, input]

   ``C`` -- array of length ``S`` containing the stage times
   [``realtype``, input] 

   ``A`` -- array of length ``S*S`` containing the ERK coefficients
   (stored in row-major, "C", order) [``realtype``, input]

   ``B`` -- array of length ``S`` containing the solution coefficients 
   [``realtype``, input]

   ``BEMBED`` -- array of length ``S`` containing the embedding
   coefficients [``realtype``, input]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]



.. _Fcn.FARKSETIRKTABLE:

``FARKSETIRKTABLE``
""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSETIRKTABLE(S, Q, P, C, A, B, BEMBED, IER)``

:Description:  Interface to the routine :ref:`Fcn.ARKodeSetIRKTable`.

:Arguments: ``S`` -- number of stages in the table [``int``, input]

   ``Q`` -- global order of accuracy of the method [``int``, input]

   ``P`` -- global order of accuracy of the embedding [``int``, input]

   ``C`` -- array of length ``S`` containing the stage times
   [``realtype``, input] 

   ``A`` -- array of length ``S*S`` containing the IRK coefficients
   (stored in row-major, "C", order) [``realtype``, input]

   ``B`` -- array of length ``S`` containing the solution coefficients 
   [``realtype``, input]

   ``BEMBED`` -- array of length ``S`` containing the embedding
   coefficients [``realtype``, input]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]



.. _Fcn.FARKSETARKTABLES:

``FARKSETARKTABLES``
""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSETARKTABLES(S, Q, P, C, AI, AE, B, BEMBED, IER)``

:Description:  Interface to the routine :ref:`Fcn.ARKodeSetARKTables`.

:Arguments: ``S`` -- number of stages in the table [``int``, input]

   ``Q`` -- global order of accuracy of the method [``int``, input]

   ``P`` -- global order of accuracy of the embedding [``int``, input]

   ``C`` -- array of length ``S`` containing the stage times
   [``realtype``, input] 

   ``AI`` -- array of length ``S*S`` containing the IRK coefficients
   (stored in row-major, "C", order) [``realtype``, input]

   ``AE`` -- array of length ``S*S`` containing the ERK coefficients
   (stored in row-major, "C", order) [``realtype``, input]

   ``B`` -- array of length ``S`` containing the solution coefficients 
   [``realtype``, input]

   ``BEMBED`` -- array of length ``S`` containing the embedding
   coefficients [``realtype``, input]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]



.. _Fcn.FARKSETADAPTIVITYMETHOD:

``FARKSETADAPTIVITYMETHOD``
"""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSETADAPTIVITYMETHOD(METHOD, PARAMS, IER)``

:Description:  Interface to the routine :ref:`Fcn.ARKodeSetAdaptivityMethod`.

:Arguments: ``METHOD`` -- flag specifying the method [``int``, input]

   ``PARAMS`` -- array of length 9 containing the adaptivity parameters 
   [``realtype``, input]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]



.. _Fcn.FARKSETADAPTIVITYCONSTANTS:

``FARKSETADAPTIVITYCONSTANTS``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSETADAPTIVITYCONSTANTS(ETAMX1, ETAMXF, ETACF, SMALLNEF, IER)``

:Description:  Interface to the routine :ref:`Fcn.ARKodeSetAdaptivityConstants`.

:Arguments: ``ETAMX1`` -- max change for the first step [``realtype``, input]

   ``ETAMXF`` -- step change on error failure [``realtype``, input]

   ``ETACF`` -- step change on a convergence failure [``realtype``, input]

   ``SMALLNEF`` -- No. of error failures before enforcing ``ETAMXF`` [``int``, input]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]



.. _Fcn.FARKSETNEWTONCONSTANTS:

``FARKSETNEWTONCONSTANTS``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSETNEWTONCONSTANTS(CRDOWN, RDIV, IER)``

:Description:  Interface to the routine :ref:`Fcn.ARKodeSetNewtonConstants`.

:Arguments: ``CRDOWN`` -- convergence rate estimation constant [``realtype``, input]

   ``RDIV`` -- divergence bound [``realtype``, input]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]



.. _Fcn.FARKSETLSETUPCONSTANTS:

``FARKSETLSETUPCONSTANTS``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSETLSETUPCONSTANTS(DGMAX, MSBP, IER)``

:Description:  Interface to the routine :ref:`Fcn.ARKodeSetLSetupConstants`.

:Arguments: ``DGMAX`` -- maximum allowable gamma ratio [``realtype``, input]

   ``MSBP`` -- maximum number of time steps between ``lsetup`` calls [``int``, input]

   ``IER`` -- return flag (0 success, :math:`\ne 0` failure) [``int``, output]



.. _FInterface.LinearSolver:

Linear solver specification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the case of using either an implicit or ImEx method, the solution
of each Runge-Kutta stage may involve the solution of linear systems
related to the Jacobian :math:`J = \frac{\partial f_I}{\partial y}` of
the implicit portion of the ODE system. ``ARKode`` presently includes
seven choices for the treatment of these systems, and the user of
``FARKODE`` must call a routine with a specific name to make the
desired choice. 


[**S**] Dense treatment of the linear system
"""""""""""""""""""""""""""""""""""""""""""""""

To use the direct dense linear solver based on the internal ``ARKode``
implementation, the user must call the ``FARKDENSE`` routine.

.. _Fcn.FARKDENSE:

``FARKDENSE``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKDENSE(NEQ, IER)``

:Description:  Interfaces with the :ref:`Fcn.ARKDense` function to
   specify use of the dense direct linear solver.

:Arguments: ``NEQ`` -- size of the ODE system [``long int``, input]

   ``IER`` -- return flag (0 if success, -1 if a memory allocation
   error occurred, -2 for an illegal input) [``int``, output]


Alteratively, to use the LAPACK-based direct dense linear solver, a
user must call the similar ``FARKLAPACKDENSE`` routine.

.. _Fcn.FARKLAPACKDENSE:

``FARKLAPACKDENSE``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKLAPACKDENSE(NEQ, IER)``

:Description:  Interfaces with the :ref:`Fcn.ARKLapackDense` function
   to specify use of the LAPACK the dense direct linear solver.

:Arguments: ``NEQ`` -- size of the ODE system [``int``, input]

   ``IER`` -- return flag (0 if success, -1 if a memory allocation
   error occurred, -2 for an illegal input) [``int``, output]


As an option when using either of these dense linear solvers, the user
may supply a routine that computes a dense approximation of the system
Jacobian :math:`J = \frac{\partial f_I}{\partial y}`. If supplied, it
must have one of the following forms:

.. _Fcn.FARKDJAC:

``FARKDJAC``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKDJAC(NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)``

:Description:  Interface to provide a user-supplied dense Jacobian
   approximation function (of type :ref:`Fcn.ARKDenseJacFn`), to be
   used by the :ref:`Fcn.FARKDENSE` solver.

:Arguments: ``NEQ`` -- size of the ODE system [``long int``, input]

   ``T`` -- current value of the independent variable [``realtype``, input]

   ``Y`` -- array containing values of the dependent state variables [``realtype``, input]

   ``FY`` -- array containing values of the dependent state derivatives [``realtype``, input]

   ``DJAC`` -- 2D array containing the Jacobian entries [``realtype`` of size ``(NEQ,NEQ)``, output]

   ``H`` -- current step size [``realtype``, input]

   ``IPAR`` -- array containing integer user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``long int``, input]

   ``RPAR`` -- array containing real user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``realtype``, input]

   ``WK1``, ``WK2``, ``WK3``  -- array containing temporary workspace
   of same size as ``Y`` [``realtype``, input]

   ``IER`` -- return flag (0 if success, >0 if a recoverable error
   occurred, <0 if an unrecoverable error occurred) [``int``, output]

:Notes: Typically this routine will use only ``NEQ``, ``T``, ``Y``,
   and ``DJAC``. It must compute the Jacobian and store it column-wise in ``DJAC``.


.. _Fcn.FARKLDJAC:

``FARKLDJAC``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKLDJAC(NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)``

:Description:  Interface to provide a user-supplied dense Jacobian
   approximation function (of type :ref:`Fcn.ARKLapackJacFn`), to be
   used by the :ref:`Fcn.FARKLAPACKDENSE` solver.

:Arguments: these all match those for :ref:`Fcn.FARKDJAC`.


If either of the above routines (:ref:`Fcn.FARKDJAC` or
:ref:`Fcn.FARKLDJAC`) uses difference quotient approximations, it may
need to use the error weight array ``EWT`` and current stepsize ``H``
in the calculation of suitable increments. The array ``EWT`` can be
obtained by calling :ref:`Fcn.FARKGETERRWEIGHTS` using one of the work
arrays as temporary storage for ``EWT``. It may also need the unit
roundoff, which can be obtained as the optional output ``ROUT(6)``,
passed from the calling program to this routine using either ``RPAR``
or a common block. 

If the :ref:`Fcn.FARKDJAC` routine is provided, then, following the
call to :ref:`Fcn.FARKDENSE`, the user must call the routine
``FARKDENSESETJAC``. 

.. _Fcn.FARKDENSESETJAC:

``FARKDENSESETJAC``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKDENSESETJAC(FLAG, IER)``

:Description:  Interface to the :ref:`Fcn.ARKDenseSetJacFn` function,
   specifying to use the user-supplied routine :ref:`Fcn.FARKDJAC` for
   the Jacobian approximation.

:Arguments: ``FLAG`` -- any nonzero value specifies to use :ref:`Fcn.FARKDJAC` [``int``, input]

   ``IER`` -- return flag (0 if success, :math:`\ne 0` if an error
   occurred) [``int``, output]


Similarly, if the :ref:`Fcn.FARKLDJAC` routine is provided, then,
following the call to :ref:`Fcn.FARKLAPACKDENSE`, the user must call
the routine ``FARKLAPACKDENSESETJAC``. 

.. _Fcn.FARKLAPACKDENSESETJAC:

``FARKLAPACKDENSESETJAC``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKLAPACKDENSESETJAC(FLAG, IER)``

:Description:  Interface to the :ref:`Fcn.ARKLapackSetJacFn` function,
   specifying to use the user-supplied routine :ref:`Fcn.FARKLDJAC` for
   the Jacobian approximation.

:Arguments: ``FLAG`` -- any nonzero value specifies to use
   :ref:`Fcn.FARKLDJAC` [``int``, input]

   ``IER`` -- return flag (0 if success, :math:`\ne 0` if an error
   occurred) [``int``, output]





[**S**] Band treatment of the linear system
"""""""""""""""""""""""""""""""""""""""""""""""

To use the direct band linear solver based on the internal ``ARKode``
implementation, the user must call the ``FARKBAND`` routine.

.. _Fcn.FARKBAND:

``FARKBAND``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKBAND(NEQ, MU, ML, IER)``

:Description:  Interfaces with the :ref:`Fcn.ARKBand` function to
   specify use of the dense banded linear solver.

:Arguments: ``NEQ`` -- size of the ODE system [``long int``, input]

   ``MU`` -- upper half-bandwidth [``long int``, input]

   ``ML`` -- lower half-bandwidth [``long int``, input]

   ``IER`` -- return flag (0 if success, -1 if a memory allocation
   error occurred, -2 for an illegal input) [``int``, output]


Alteratively, to use the LAPACK-based direct banded linear solver, a
user must call the similar ``FARKLAPACKBAND`` routine.


.. _Fcn.FARKLAPACKBAND:

``FARKLAPACKBAND``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKLAPACKBAND(NEQ, MU, ML, IER)``

:Description:  Interfaces with the :ref:`Fcn.ARKLapackBand` function
   to specify use of the dense banded linear solver.

:Arguments: ``NEQ`` -- size of the ODE system [``int``, input]

   ``MU`` -- upper half-bandwidth [``int``, input]

   ``ML`` -- lower half-bandwidth [``int``, input]

   ``IER`` -- return flag (0 if success, -1 if a memory allocation
   error occurred, -2 for an illegal input) [``int``, output]


As an option when using either of these banded linear solvers, the user
may supply a routine that computes a banded approximation of the
linear system Jacobian :math:`J = \frac{\partial f_I}{\partial y}`. If
supplied, it must have one of the following forms:

.. _Fcn.FARKBJAC:

``FARKBJAC``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKBJAC(NEQ, MU, ML, MDIM, T, Y, FY, BJAC,
   H, IPAR, RPAR, WK1, WK2, WK3, IER)``

:Description:  Interface to provide a user-supplied band Jacobian
   approximation function (of type :ref:`Fcn.ARKBandJacFn`), to be
   used by the :ref:`Fcn.FARKBAND` solver.

:Arguments: ``NEQ`` -- size of the ODE system [``long int``, input]

   ``MU``   -- upper half-bandwidth [``long int``, input]

   ``ML``   -- lower half-bandwidth [``long int``, input]

   ``MDIM`` -- leading dimension of ``BJAC`` array [``long int``, input]

   ``T``    -- current value of the independent variable [``realtype``, input]

   ``Y``    -- array containing dependent state variables [``realtype``, input]

   ``FY``   -- array containing dependent state derivatives [``realtype``, input]

   ``BJAC`` -- 2D array containing the Jacobian entries [``realtype`` of size
   ``(MDIM,NEQ)``, output]

   ``H``    -- current step size [``realtype``, input]

   ``IPAR`` -- array containing integer user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``long int``, input]

   ``RPAR`` -- array containing real user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``realtype``, input]

   ``WK1``, ``WK2``, ``WK3``  -- array containing temporary workspace
   of same size as ``Y`` [``realtype``, input]

   ``IER`` -- return flag (0 if success, >0 if a recoverable error
   occurred, <0 if an unrecoverable error occurred) [``int``, output]

:Notes: Typically this routine will use only ``NEQ``, ``MU``, ``ML``,
   ``T``, ``Y``, and ``BJAC``. It must load the ``MDIM`` by ``N``
   array ``BJAC`` with the Jacobian matrix at the current
   :math:`(t,y)` in band form.  Store in ``BJAC(k,j)`` the Jacobian
   element :math:`J_{i,j}` with ``k = i - j + MU + 1`` (or ``k = 1
   ... ML+MU+1``) and ``j = 1 ... N``.


.. _Fcn.FARKLBJAC:

``FARKLBJAC``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKLBJAC(NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)``

:Description:  Interface to provide a user-supplied banded Jacobian
   approximation function (of type :ref:`Fcn.ARKLapackJacFn`), to be
   used by the :ref:`Fcn.FARKLAPACKBAND` solver.

:Arguments: these all match those for :ref:`Fcn.FARKBJAC`.


If either of the above routines (:ref:`Fcn.FARKBJAC` or
:ref:`Fcn.FARKLBJAC`) uses difference quotient approximations, it may
need to use the error weight array ``EWT`` and current stepsize ``H``
in the calculation of suitable increments. The array ``EWT`` can be
obtained by calling :ref:`Fcn.FARKGETERRWEIGHTS` using one of the work
arrays as temporary storage for ``EWT``. It may also need the unit
roundoff, which can be obtained as the optional output ``ROUT(6)``,
passed from the calling program to this routine using either ``RPAR``
or a common block. 

If the :ref:`Fcn.FARKBJAC` routine is provided, then, following the
call to :ref:`Fcn.FARKBAND`, the user must call the routine
``FARKBANDSETJAC``. 

.. _Fcn.FARKBANDSETJAC:

``FARKBANDSETJAC``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKBANDSETJAC(FLAG, IER)``

:Description:  Interface to the :ref:`Fcn.ARKBandSetJacFn` function,
   specifying to use the user-supplied routine :ref:`Fcn.FARKBJAC` for
   the Jacobian approximation.

:Arguments: ``FLAG`` -- any nonzero value specifies to use
   :ref:`Fcn.FARKBJAC` [``int``, input] 

   ``IER`` -- return flag (0 if success, :math:`\ne 0` if an error
   occurred) [``int``, output]


Similarly, if the :ref:`Fcn.FARKLBJAC` routine is provided, then,
following the call to :ref:`Fcn.FARKLAPACKBAND`, the user must call
the routine ``FARKLAPACKBANDSETJAC``. 

.. _Fcn.FARKLAPACKBANDSETJAC:

``FARKLAPACKBANDSETJAC``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKLAPACKBANDSETJAC(FLAG, IER)``

:Description:  Interface to the :ref:`Fcn.ARKLapackSetJacFn` function,
   specifying to use the user-supplied routine :ref:`Fcn.FARKLBJAC` for
   the Jacobian approximation.

:Arguments: ``FLAG`` -- any nonzero value specifies to use
   :ref:`Fcn.FARKLBJAC` [``int``, input]

   ``IER`` -- return flag (0 if success, :math:`\ne 0` if an error
   occurred) [``int``, output]




[**S**][**P**] SPGMR treatment of the linear systems
"""""""""""""""""""""""""""""""""""""""""""""""""""""

For the Scaled Preconditioned GMRES solution of the linear systems,
the user must call the ``FARKSPGMR`` routine.

.. _Fcn.FARKSPGMR:

``FARKSPGMR``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSPGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)``

:Description:  Interfaces with the :ref:`Fcn.ARKSpgmr` and
   ``ARKSpilsSet*`` routines to specify use of the ``SPGMR`` iterative
   linear solver.

:Arguments: 
   ``IPRETYPE`` -- preconditioner type [``int``, input]: 

      0 = none 

      1 = left only

      2 = right only

      3 = both sides

   ``IGSTYPE`` -- Gram-schmidt process type [``int``, input]: 

      1 = modified G-S

      2 = classical G-S

   ``MAXL`` -- maximum Krylov subspace dimension (0 for default) [``int``; input].

   ``DELT`` -- linear convergence tolerance factor (0.0 for default) [``realtype``, input].

   ``IER`` -- return flag (0 if success, -1 if a memory allocation
   error occurred, -2 for an illegal input) [``int``, output]


For descriptions of the optional user-supplied routines for use with
:ref:`Fcn.FARKSPGMR` see the section :ref:`FInterface.SpilsUserSupplied`.





[**S**][**P**] SPBCG treatment of the linear systems
"""""""""""""""""""""""""""""""""""""""""""""""""""""

For the Scaled Preconditioned Bi-CGStab solution of the linear systems,
the user must call the ``FARKSPBCG`` routine.

.. _Fcn.FARKSPBCG:

``FARKSPBCG``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSPBCG(IPRETYPE, MAXL, DELT, IER)``

:Description:  Interfaces with the :ref:`Fcn.ARKSpbcg` and
   ``ARKSpilsSet*`` routines to specify use of the ``SPBCG`` iterative
   linear solver.

:Arguments: The arguments are the same as those with the same names
   for :ref:`Fcn.FARKSPGMR`.

For descriptions of the optional user-supplied routines for use with
:ref:`Fcn.FARKSPBCG` see the section :ref:`FInterface.SpilsUserSupplied`.





[**S**][**P**] SPTFQMR treatment of the linear systems
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For the Scaled Preconditioned TFQMR solution of the linear systems,
the user must call the ``FARKSPTFQMR`` routine.

.. _Fcn.FARKSPTFQMR:

``FARKSPTFQMR``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSPTFQMR(IPRETYPE, MAXL, DELT, IER)``

:Description:  Interfaces with the :ref:`Fcn.ARKSptfqmr` and
   ``ARKSpilsSet*`` routines to specify use of the ``SPTFQMR`` iterative
   linear solver.

:Arguments: The arguments are the same as those with the same names
   for :ref:`Fcn.FARKSPGMR`.

For descriptions of the optional user-supplied routines for use with
:ref:`Fcn.FARKSPTFQMR` see the next section.



.. _FInterface.SpilsUserSupplied:

[**S**][**P**] User-supplied routines for SPGMR/SPBCG/SPTFQMR
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

With treatment of the linear systems by any of the Krylov iterative
solvers, there are three optional user-supplied routines --
:ref:`Fcn.FARKJTIMES`, :ref:`Fcn.FARKPSET` and :ref:`Fcn.FARKPSOL`.
The specifications of these functions are given below.

As an option when using the ``SPGMR``, ``SPBCG`` or ``SPTFQMR`` linear
solvers, the user may supply a routine that computes the product of
the system Jacobian :math:`J = \frac{\partial f_I}{\partial y}` and a
given vector :math:`v`.  If supplied, it must have the following form:

.. _Fcn.FARKJTIMES:

``FARKJTIMES``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKJTIMES(V, FJV, T, Y, FY, H, IPAR, RPAR, WORK, IER)``

:Description:  Interface to provide a user-supplied
   Jacobian-times-vector product approximation function (of type
   :ref:`Fcn.ARKSpilsJacTimesVecFn`), to be used by one of the Krylov
   iterative linear solvers.

:Arguments: ``V`` -- array containing the vector to multiply
   [``realtype``, input]

   ``FJV``  -- array containing resulting product vector
   [``realtype``, output] 

   ``T``    -- current value of the independent variable
   [``realtype``, input] 

   ``Y``    -- array containing dependent state variables
   [``realtype``, input] 

   ``FY``   -- array containing dependent state derivatives
   [``realtype``, input] 

   ``H``    -- current step size [``realtype``, input]

   ``IPAR`` -- array containing integer user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``long int``, input]

   ``RPAR`` -- array containing real user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``realtype``, input]

   ``WORK`` -- array containing temporary workspace of same size as
   ``Y`` [``realtype``, input]

   ``IER``  -- return flag  (0 if success, :math:`\ne 0` if an error)
   [``int``, output]

:Notes: Typically this routine will use only ``NEQ``, ``T``, ``Y``,
   ``V``, and ``FJV``.  It must compute the product vector :math:`Jv`,
   where :math:`v` is given in ``V``, and the product is stored in
   ``FJV``.


If this routine has been supplied by the user, then, following the
call to :ref:`Fcn.FARKSPGMR`, :ref:`Fcn.FARKSPBCG` or
:ref:`Fcn.FARKSPTFQMR`, the user must call the routine
``FARKSPILSSETJAC`` with ``FLAG`` :math:`\ne 0` to specify use of the
user-supplied Jacobian-times-vector function.


.. _Fcn.FARKSPILSSETJAC:

``FARKSPILSSETJAC``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSPILSSETJAC(FLAG, IER)``

:Description:  Interface to the function 
   :ref:`Fcn.ARKSpilsSetJacTimesVecFn` to specify use of the
   user-supplied Jacobian-times-vector function :ref:`Fcn.FARKJTIMES`.

:Arguments: ``FLAG`` -- flag denoting to use ``FARKJTIMES`` routine [``int``, input]

   ``IER``  -- return flag  (0 if success, :math:`\ne 0` if an error)
   [``int``, output]


If preconditioning is to be performed during the Krylov solver
(i.e. the solver was set up with ``IPRETYPE`` :math:`\ne 0`), then the
user must also call the routine ``FARKSPILSSETPREC`` with ``FLAG``
:math:`\ne 0`. 

.. _Fcn.FARKSPILSSETPREC:

``FARKSPILSSETPREC``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSPILSSETJAC(FLAG, IER)``

:Description:  Interface to the function 
   :ref:`Fcn.ARKSpilsSetPreconditioner` to specify use of the
   user-supplied preconditioner setup and solve functions,
   :ref:`Fcn.FARKPSET` and :ref:`Fcn.FARKPSOL`, respectively.

:Arguments: ``FLAG`` -- flag denoting use of user-supplied
   preconditioning routines [``int``, input] 

   ``IER``  -- return flag  (0 if success, :math:`\ne 0` if an error)
   [``int``, output]


In addition, the user must provide the following two routines to
implement the preconditioner setup and solve functions to be used
within the solve.

.. _Fcn.FARKPSET:

``FARKPSET``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKPSET(T,Y,FY,JOK,JCUR,GAMMA,H,IPAR,RPAR,V1,V2,V3,IER)``

:Description:  User-supplied preconditioner setup routine (of type
   :ref:`Fcn.ARKSpilsPrecSetupFn`). 

:Arguments: ``T`` -- current value of the independent variable
   [``realtype``, input] 

   ``Y`` -- current dependent state variable array [``realtype``, input]

   ``FY`` -- current dependent state variable derivative array [``realtype``, input]

   ``JOK`` -- flag indicating whether Jacobian-related data needs to be 
   recomputed [int, input]:
  
      0 = recompute, 

      1 = reuse with the current value of ``GAMMA``.

   ``JCUR`` -- return flag to denote if Jacobian data was recomputed
   (1=yes, 0=no)  [``realtype``, output]

   ``GAMMA`` -- Jacobian scaling factor [``realtype``, input]

   ``H`` -- current step size [``realtype``, input]

   ``IPAR`` -- array containing integer user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``long int``, input/output]

   ``RPAR`` -- array containing real user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``realtype``, input/output]

   ``V1``, ``V2``, ``V3`` -- arrays containing temporary workspace of
   same size as ``Y`` [``realtype``, input]

   ``IER``  -- return flag  (0 if success, >0 if a recoverable
   failure, <0 if a non-recoverable failure) [``int``, output]

:Notes: This routine must set up the preconditioner ``P`` to be used
   in the subsequent call to :ref:`Fcn.FARKPSOL`.  The preconditioner
   (or the product of the left and right preconditioners if using
   both) should be an approximation to the matrix  :math:`M - \gamma
   J`, where :math:`M` is the system mass matrix, :math:`\gamma` is
   the input ``GAMMA``, and :math:`J = \frac{\partial f_I}{\partial y}`.


.. _Fcn.FARKPSOL:

``FARKPSOL``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKPSOL(T,Y,FY,R,Z,GAMMA,DELTA,LR,IPAR,RPAR,VT,IER)``

:Description:  User-supplied preconditioner solve routine (of type
   :ref:`Fcn.ARKSpilsPrecSolveFn`). 

:Arguments: ``T`` -- current value of the independent variable
   [``realtype``, input] 

   ``Y`` -- current dependent state variable array [``realtype``, input]

   ``FY`` -- current dependent state variable derivative array [``realtype``, input]

   ``R`` -- right-hand side array [``realtype``, input]

   ``Z`` -- solution array [``realtype``, output]

   ``GAMMA`` -- Jacobian scaling factor [``realtype``, input]

   ``DELTA`` -- desired residual tolerance [``realtype``, input]

   ``LR`` -- flag denoting to solve the right or left preconditioner
   system:

      1 = left preconditioner

      2 = right preconditioner

   ``IPAR`` -- array containing integer user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``long int``, input/output]

   ``RPAR`` -- array containing real user data that was passed to
   :ref:`Fcn.FARKMALLOC` [``realtype``, input/output]

   ``VT`` -- array containing temporary workspace of same size as ``Y``  
   [``realtype``, input]

   ``IER``  -- return flag  (0 if success, >0 if a recoverable
   failure, <0 if a non-recoverable failure) [``int``, output]

:Notes: Typically this routine will use only ``NEQ``, ``T``, ``Y``,
   ``GAMMA``, ``R``, ``LR``, and ``Z``.  It must solve the
   preconditioner linear system :math:`Pz = r`.  The preconditioner
   (or the product of the left and right preconditioners if both are
   nontrivial) should be an approximation to the matrix  :math:`M - \gamma
   J`, where :math:`M` is the system mass matrix, :math:`\gamma` is
   the input ``GAMMA``, and :math:`J = \frac{\partial f_I}{\partial y}`.


Notes:

(a) If the user's :ref:`Fcn.FARKJTIMES` or :ref:`Fcn.FARKPSET` routine
    uses difference quotient approximations, it may need to use the
    error weight array ``EWT``, the current stepsize ``H``, and/or the
    unit roundoff, in the calculation of suitable increments. Also, If
    :ref:`Fcn.FARKPSOL` uses an iterative method in its solution, the
    residual vector :math:`\rho = r - Pz` of the system should be made
    less than :math:`\delta =` ``DELTA`` in the weighted l2 norm, i.e.
    
    .. math::
       \left(\sum_i \left(\rho_i * EWT_i\right)^2 \right)^{1/2} < \delta.

(b) If needed in :ref:`Fcn.FARKJTIMES`, :ref:`Fcn.FARKPSOL`, or
    :ref:`Fcn.FARKPSET`, the error weight array ``EWT`` can be
    obtained by calling :ref:`Fcn.FARKGETERRWEIGHTS` using one of the
    work arrays as temporary storage for ``EWT``. 

(c) If needed in :ref:`Fcn.FARKJTIMES`, :ref:`Fcn.FARKPSOL`, or
    :ref:`Fcn.FARKPSET`, the unit roundoff can be obtained as the
    optional output ``ROUT(6)`` (available after the call to
    :ref:`Fcn.FARKMALLOC`) and can be passed using either the ``RPAR``
    user data array or a common block. 




.. _FInterface.Solution:

Problem solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Carrying out the integration is accomplished by making calls to
:ref:`Fcn.FARKODE`.

.. _Fcn.FARKODE:

``FARKODE``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKODE(TOUT, T, Y, ITASK, IER)``

:Description:  Fortran interface to the C routine :ref:`Fcn.ARKode`
   for performing the solve, along with many of the ``ARK*Get*``
   routines for reporting on solver statistics.

:Arguments: ``TOUT`` -- next value of :math:`t` at which a solution is
   desired [``realtype``, input]

   ``T`` -- current value of independent variable reached by the solver
   [``realtype``, output] 

   ``Y`` -- array containing dependent state variables on output
   [``realtype``, output] 

   ``ITASK`` -- task indicator [``int``, input]:

      1 = normal mode (overshoot ``TOUT`` and interpolate)

      2 = one-step mode (return after each internal step taken)

      3 = normal ``tstop`` mode (like 1, but integration never
      proceeds past ``TSTOP``, which must be specified through a
      preceding call to :ref:`Fcn.FARKSETRIN` using the key
      ``STOP_TIME``)

      4 = one step ``tstop`` (like 2, but integration never goes past
      ``TSTOP``) 

   ``IER`` -- completion flag [int, output]: 

      0 = success, 

      1 = tstop return, 

      2 = root return, 

      values -1 ... -10 are failure modes (see :ref:`Fcn.ARKode` and
      :ref:`Constants`).

:Notes: The current values of the optional outputs are immediately
   available in ``IOUT`` and ``ROUT`` upon return from this function
   (see :ref:`FInterface.OptionalIntegerOutputTable` and
   :ref:`FInterface.OptionalRealOutputTable`).



.. _FInterface.AdditionalOutput:

Additional solution output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After a successful return from :ref:`Fcn.FARKODE`, the routine
:ref:`Fcn.FARKDKY` may be used to obtain a derivative of the solution,
of order up to 3, at any :math:`t` within the last step taken. 

.. _Fcn.FARKDKY:

``FARKDKY``
"""""""""""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKDKY(T, K, DKY, IER)``

:Description:  Fortran interface to the C routine :ref:`Fcn.ARKDKY`
   for interpolating output of the solution or its derivatives at any
   point within the last step taken.

:Arguments: ``T`` -- time at which solution derivative is desired,
   within the interval :math:`[t_n-h,t_n]`, [``realtype``, input].

   ``K`` -- derivative order :math:`(0 \le k \le 3)` [``int``, input]

   ``DKY`` -- array containing the computed K-th derivative of
   :math:`y` [``realtype``, output] 

   ``IER`` -- return flag (0 if success, <0 if an illegal argument)
   [``int``, output]




.. _FInterface.ReInit:

Problem reinitialization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To re-initialize the ``ARKode`` solver for the solution of a new
problem of the same size as one already solved, the user must call
:ref:`Fcn.FARKREINIT`. 

.. _Fcn.FARKREINIT:

``FARKREINIT``
"""""""""""""""

:Definition: ``SUBROUTINE FARKREINIT(T0, Y0, IMEX, IATOL, RTOL, ATOL, IER)`` 

:Description:  Re-initializes the Fortran interface to the ``ARKode``
   solver.

:Arguments: The arguments have the same names and meanings as those of
   :ref:`Fcn.FARKMALLOC`.

:Notes: This routine performs no memory allocation, instead using the
   existing memory created by the previous :ref:`Fcn.FARKMALLOC`
   call.  The call to specify the linear system solution method may or
   may not be needed.


Following a call to :ref:`Fcn.FARKREINIT`, a call to specify the
linear system solver must be made if the choice of linear solver is
being changed. Otherwise, a call to reinitialize the linear solver
last used may or may not be needed, depending on changes in the inputs
to it. 

In the case of the ``BAND`` solver, for any change in the
half-bandwidth parameters, call :ref:`Fcn.FARKBAND` (or
:ref:`Fcn.FARKLAPACKBAND`) again described above.

In the case of ``spgmr``, for a change of inputs other than ``MAXL``,
the user may call the routine :ref:`Fcn.FCVSPGMRREINIT` to
reinitialize ``SPGMR`` without reallocating its memory, as follows: 

.. _Fcn.FARKSPGMRREINIT:

``FARKSPGMRREINIT``
""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSPGMRREINIT(IPRETYPE, IGSTYPE, DELT, IER)`` 

:Description:  Re-initializes the Fortran interface to the ``SPGMR``
   linear solver.

:Arguments: The arguments have the same names and meanings as those of
   :ref:`Fcn.FARKSPGMR`.

However, if ``MAXL`` is being changed, then the user should call
:ref:`Fcn.FARKSPGMR` instead.

In the case of ``SPBCG``, for a change in any inputs, the user can
reinitialize ``SPBCG`` without reallocating its memory by calling
:ref:`Fcn.FARKSPBCGREINIT`, as follows:

.. _Fcn.FARKSPBCGREINIT:

``FARKSPBCGREINIT``
""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSPBCGREINIT(IPRETYPE, MAXL, DELT, IER)`` 

:Description:  Re-initializes the Fortran interface to the ``SPBCG``
   linear solver.

:Arguments: The arguments have the same names and meanings as those of
   :ref:`Fcn.FARKSPBCG`.

In the case of ``SPTFQMR``, for a change in any inputs, the user can
reinitialize ``SPTFQMR`` without reallocating its memory by calling
:ref:`Fcn.FARKSPTFQMRREINIT`, as follows:

.. _Fcn.FARKSPTFQMRREINIT:

``FARKSPTFQMRREINIT``
""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKSPTFQMRREINIT(IPRETYPE, MAXL, DELT, IER)`` 

:Description:  Re-initializes the Fortran interface to the ``SPBTFQMR``
   linear solver.

:Arguments: The arguments have the same names and meanings as those of
   :ref:`Fcn.FARKSPTFQMR`.




.. _FInterface.Deallocation:

Memory deallocation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To free the internal memory created by :ref:`Fcn.FARKMALLOC`, the user
may call :ref:`Fcn.FARKFREE`, as follows:

.. _Fcn.FARKFREE:

``FARKFREE``
""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKFREE()`` 

:Description:  Frees the internal memory created by :ref:`Fcn.FARKMALLOC`.

:Arguments: None.




.. _FInterface.OptionalOutputs:

FARKODE optional output
-----------------------------

The optional inputs to ``FARKODE`` have already been described in the
section :ref:`FInterface.OptionalInputs`.  


``IOUT`` and ``ROUT`` arrays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The optional outputs from the :ref:`Fcn.ARKode` solver are accessed
not through individual functions, but rather through a pair of arrays,
``IOUT`` (``long int`` type) of dimension at least 22, and ``ROUT``
(``realtype`` type) of dimension at least 6. These arrays are owned
(and allocated) by the user and are passed as arguments to
:ref:`Fcn.FARKMALLOC`. 

:ref:`FInterface.IOUTTable` and
:ref:`FInterface.ROUTTable` list the entries in these
arrays associated with the main ``ARKode`` solver, along with the
relevant ``ARKode`` function that is actually called to extract the
optional output.  Similarly,
:ref:`FInterface.DlsIOUTTable` lists the ``IOUT``
entries associated with the main ``ARKDENSE`` and ``ARKBAND`` direct
linear solvers, and :ref:`FInterface.SpilsIOUTTable`
lists the ``IOUT`` entries associated with the main ``ARKSPGMR``,
``ARKSPBCG`` and ``ARKSPTFQMR`` iterative linear solvers.

For more details on the optional inputs and outputs to ``ARKode``, see
the sections :ref:`CInterface.OptionalInputs` and
:ref:`CInterface.OptionalOutputs`.



.. _FInterface.IOUTTable:

Table: Optional FARKODE integer outputs
""""""""""""""""""""""""""""""""""""""""

   ==============  ===============  ===================================================
   ``IOUT`` Index  Optional output  ``ARKode`` function
   ==============  ===============  ===================================================
   1               LENRW            :ref:`Fcn.ARKodeGetWorkSpace`
   2               LENIW            :ref:`Fcn.ARKodeGetWorkSpace`
   3               NST              :ref:`Fcn.ARKodeGetNumSteps`
   4               NST_STB          :ref:`Fcn.ARKodeGetNumExpSteps`
   5               NST_ACC          :ref:`Fcn.ARKodeGetNumAccSteps`
   6               NST_CNV          :ref:`Fcn.ARKodeGetNumConvSteps`
   7               NFE              :ref:`Fcn.ARKodeGetNumRhsEvals` (:math:`f_E` calls)
   8               NFI              :ref:`Fcn.ARKodeGetNumRhsEvals` (:math:`f_I` calls)
   9               NSETUPS          :ref:`Fcn.ARKodeGetNumLinSolvSetups`
   10              NETF             :ref:`Fcn.ARKodeGetNumErrTestFails`
   11              NNI              :ref:`Fcn.ARKodeGetNumNonlinSolvIters`
   12              NCFN             :ref:`Fcn.ARKodeGetNumNonlinSolvConvFails`
   13              NGE              :ref:`Fcn.ARKodeGetNumGEvals`
   ==============  ===============  ===================================================



.. _FInterface.ROUTTable:

Table: Optional FARKODE real outputs 
""""""""""""""""""""""""""""""""""""""""

   ==============  ===============  ===================================================
   ``ROUT`` Index  Optional output  ``ARKode`` function
   ==============  ===============  ===================================================
   1               H0U              :ref:`Fcn.ARKodeGetActualInitStep`
   2               HU               :ref:`Fcn.ARKodeGetLastStep`
   3               HCUR             :ref:`Fcn.ARKodeGetCurrentStep`
   4               TCUR             :ref:`Fcn.ARKodeGetCurrentTime`
   5               TOLSF            :ref:`Fcn.ARKodeGetTolScaleFactor`
   6               UROUND           ``UNIT_ROUNDOFF`` (see :ref:`CInterface.DataTypes`)
   ==============  ===============  ===================================================



.. _FInterface.DlsIOUTTable:

Table: Optional ARKDENSE and ARKBAND outputs
""""""""""""""""""""""""""""""""""""""""""""""

   ==============  ===============  ===================================================
   ``IOUT`` Index  Optional output  ``ARKode`` function
   ==============  ===============  ===================================================
   14              LENRWLS          :ref:`Fcn.ARKDlsGetWorkSpace`
   15              LENIWLS          :ref:`Fcn.ARKDlsGetWorkSpace`
   16              LSTF             :ref:`Fcn.ARKDlsGetLastFlag`
   17              NFELS            :ref:`Fcn.ARKDlsGetNumRhsEvals`
   18              NJE              :ref:`Fcn.ARKDlsGetNumJacEvals`
   ==============  ===============  ===================================================



.. _FInterface.SpilsIOUTTable:

Table: Optional ARKSPGMR, ARKSPBCG and ARKSPTFQMR outputs 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

   ==============  ===============  ===================================================
   ``IOUT`` Index  Optional output  ``ARKode`` function
   ==============  ===============  ===================================================
   14              LENRWLS          :ref:`Fcn.ARKSpilsGetWorkSpace`
   15              LENIWLS          :ref:`Fcn.ARKSpilsGetWorkSpace`
   16              LSTF             :ref:`Fcn.ARKSpilsGetLastFlag`
   17              NFELS            :ref:`Fcn.ARKSpilsGetNumRhsEvals`
   18              NJTV             :ref:`Fcn.ARKSpilsGetNumJtimesEvals`
   19              NPE              :ref:`Fcn.ARKSpilsGetNumPrecEvals`
   20              NPS              :ref:`Fcn.ARKSpilsGetNumPrecSolves`
   21              NLI              :ref:`Fcn.ARKSpilsGetNumLinIters`
   22              NCFL             :ref:`Fcn.ARKSpilsGetNumConvFails`
   ==============  ===============  ===================================================



Additional optional output routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


In addition to the optional inputs communicated through ``FARKSET*``
calls and the optional outputs extracted from ``IOUT`` and ``ROUT``,
the following user-callable routines are available: 

To obtain the error weight array ``EWT``, containing the
multiplicative error weights used the WRMS norms, the user may call
the routine :ref:`Fcn.FARKGETERRWEIGHTS` as follows:


.. _Fcn.FARKGETERRWEIGHTS:

``FARKGETERRWEIGHTS``
""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKGETERRWEIGHTS(EWT, IER)`` 

:Description:  Retrieves the current error weight vector (interfaces
   with :ref:`Fcn.ARKodeGetErrWeights`).

:Arguments: ``EWT`` -- array containing the error weight vector
   [``realtype``, output] 

   ``IER``  -- return flag  (0 if success, :math:`\ne 0` if an error)
   [``int``, output]

:Notes: The array ``EWT``, of length ``NEQ`` if using
   ``NVECTOR_SERIAL`` or ``NLOCAL`` if using ``NVECTOR_PARALLEL``,
   must already have been declared by the user. 

Similarly, to obtain the estimated local errors, following a
successful call to :ref:`Fcn.FARKSOLVE`, the user may call the routine
:ref:`Fcn.FARKGETESTLOCALERR` as follows:


.. _Fcn.FARKGETESTLOCALERR:

``FARKGETESTLOCALERR``
""""""""""""""""""""""

:Definition: ``SUBROUTINE FARKGETESTLOCALERR(ELE, IER)`` 

:Description:  Retrieves the current local truncation error estimate
   vector (interfaces with :ref:`Fcn.ARKodeGetEstLocalErrors`).

:Arguments: ``ELE`` -- array with the estimated local error vector
   [``realtype``, output] 

   ``IER``  -- return flag  (0 if success, :math:`\ne 0` if an error)
   [``int``, output]

:Notes: The array ``ELE``, of length ``NEQ`` if using
   ``NVECTOR_SERIAL`` or ``NLOCAL`` if using ``NVECTOR_PARALLEL``,
   must already have been declared by the user. 




.. _FInterface.RootFinding:

Usage of the FARKROOT interface to rootfinding
-----------------------------------------------

(to be added)




.. _FInterface.BandPre:

Usage of the FARKBP interface to ARKBANDPRE
-----------------------------------------------

(to be added)



.. _FInterface.BBDPre:

Usage of the FARKBBD interface to ARKBBDPRE
-----------------------------------------------

(to be added)



