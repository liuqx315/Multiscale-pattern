..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _FInterface.OptionalOutputs:

FARKODE optional output
==============================

We note that the optional inputs to FARKODE have already been
described in the section :ref:`FInterface.OptionalInputs`.



IOUT and ROUT arrays
----------------------------

In the Fortran interface, the optional outputs from the
:f:func:`FARKODE()` solver are accessed not through individual
functions, but rather through a pair of user-allocated arrays, *IOUT*
(having ``long int`` type) of dimension at least 22, and *ROUT*
(having ``realtype`` type) of dimension at least 6.  These arrays must
be allocated by the user program that calls :f:func:`FARKODE()`, that
passes them through the Fortran interface as arguments to
:f:func:`FARKMALLOC()`.  Following this call, :f:func:`FARKODE()` will
modify the entries of these arrays to contain all optional output
values provided to a Fortran user.

In the following tables, :ref:`FInterface.IOUTTable` and
:ref:`FInterface.ROUTTable`, we list the entries in these
arrays by index, naming them according to their role with the main
ARKode solver, and list the relevant ARKode C/C++ function that is
actually called to extract the output value.  Similarly, optional
integer output values that are specific to the ARKDENSE and ARKBAND
linear solvers are listed in :ref:`FInterface.DlsIOUTTable`, optional
integer output values that are specific to the ARKKLU and ARKSUPERLUMT
linear solvers are listed in :ref:`FInterface.SlsIOUTTable`, while
integer optional output values specific to the ARKSPGMR,
ARKSPBCG, ARKSPTFQMR, ARKSPFGMR and ARKPCG iterative linear solvers
are listed in :ref:`FInterface.SpilsIOUTTable`.

For more details on the optional inputs and outputs to ARKode, see
the sections :ref:`CInterface.OptionalInputs` and
:ref:`CInterface.OptionalOutputs`.



.. _FInterface.IOUTTable:

Table: Optional FARKODE integer outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

==============  ===============  =========================================================
*IOUT* Index    Optional output  ARKode function
==============  ===============  =========================================================
1               LENRW            :c:func:`ARKodeGetWorkSpace()`
2               LENIW            :c:func:`ARKodeGetWorkSpace()`
3               NST              :c:func:`ARKodeGetNumSteps()`
4               NST_STB          :c:func:`ARKodeGetNumExpSteps()`
5               NST_ACC          :c:func:`ARKodeGetNumAccSteps()`
6               NST_ATT          :c:func:`ARKodeGetNumStepAttempts()`
7               NFE              :c:func:`ARKodeGetNumRhsEvals()` (num :math:`f_E` calls)
8               NFI              :c:func:`ARKodeGetNumRhsEvals()` (num :math:`f_I` calls)
9               NSETUPS          :c:func:`ARKodeGetNumLinSolvSetups()`
10              NETF             :c:func:`ARKodeGetNumErrTestFails()`
11              NNI              :c:func:`ARKodeGetNumNonlinSolvIters()`
12              NCFN             :c:func:`ARKodeGetNumNonlinSolvConvFails()`
13              NGE              :c:func:`ARKodeGetNumGEvals()`
==============  ===============  =========================================================



.. _FInterface.ROUTTable:

Table: Optional FARKODE real outputs 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

==============  ===============  ===============================================================
*ROUT* Index    Optional output  ARKode function
==============  ===============  ===============================================================
1               H0U              :c:func:`ARKodeGetActualInitStep()`
2               HU               :c:func:`ARKodeGetLastStep()`
3               HCUR             :c:func:`ARKodeGetCurrentStep()`
4               TCUR             :c:func:`ARKodeGetCurrentTime()`
5               TOLSF            :c:func:`ARKodeGetTolScaleFactor()`
6               UROUND           ``UNIT_ROUNDOFF`` (see the section :ref:`CInterface.DataTypes`)
==============  ===============  ===============================================================



.. _FInterface.DlsIOUTTable:

Table: Optional ARKDENSE and ARKBAND outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

==============  ===============  ===================================================
*IOUT* Index    Optional output  ARKode function
==============  ===============  ===================================================
14              LENRWLS          :c:func:`ARKDlsGetWorkSpace()`
15              LENIWLS          :c:func:`ARKDlsGetWorkSpace()`
16              LSTF             :c:func:`ARKDlsGetLastFlag()`
17              NFELS            :c:func:`ARKDlsGetNumRhsEvals()`
18              NJE              :c:func:`ARKDlsGetNumJacEvals()`
==============  ===============  ===================================================



.. _FInterface.SlsIOUTTable:

Table: Optional ARKKLU and ARKSUPERLUMT outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

==============  ===============  ===================================================
*IOUT* Index    Optional output  ARKode function
==============  ===============  ===================================================
16              LSTF             :c:func:`ARKSlsGetLastFlag()`
18              NJE              :c:func:`ARKSlsGetNumJacEvals()`
==============  ===============  ===================================================



.. _FInterface.SpilsIOUTTable:

Table: Optional ARKSPGMR, ARKSPBCG, ARKSPTFQMR, ARKSPFGMR and ARKPCG outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

==============  ===============  ===================================================
*IOUT* Index    Optional output  ARKode function
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
---------------------------------------------

In addition to the optional inputs communicated through FARKSET*
calls and the optional outputs extracted from *IOUT* and *ROUT*,
the following user-callable routines are available.


To obtain the error weight array *EWT*, containing the
multiplicative error weights used in the WRMS norms, the user may call 
the routine :f:func:`FARKGETERRWEIGHTS()` as follows:


.. f:subroutine:: FARKGETERRWEIGHTS(EWT, IER)
   
   Retrieves the current error weight vector (interfaces
   with :c:func:`ARKodeGetErrWeights()`).
      
   **Arguments:** 
      * *EWT* (``realtype``, output) -- array containing the error
	weight vector. 
      * *IER*  (``int``, output) -- return flag  (0 if success,
	:math:`\ne 0` if an error). 
      
   **Notes:**
   The array *EWT* must have already been allocated by the user, of
   the same size as the solution array *Y*.



Similarly, to obtain the estimated local truncation errors, following
a successful call to :f:func:`FARKODE()`, the user may call the
routine :f:func:`FARKGETESTLOCALERR()` as follows:


.. f:subroutine:: FARKGETESTLOCALERR(ELE, IER)
   
   Retrieves the current local truncation error estimate
   vector (interfaces with :c:func:`ARKodeGetEstLocalErrors()`).
      
   **Arguments:** 
      * *ELE* (``realtype``, output) -- array with the estimated local
	truncation error vector. 
      * *IER*  (``int``, output) -- return flag  (0 if success,
	:math:`\ne 0` if an error).
      
   **Notes:**
   The array *ELE* must have already been allocated by the user, of
   the same size as the solution array *Y*.

