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

The optional inputs to FARKODE have already been described in the
section :ref:`FInterface.OptionalInputs`.  


IOUT and ROUT arrays
----------------------------

The optional outputs from the :c:func:`ARKode()` solver are accessed
not through individual functions, but rather through a pair of arrays,
IOUT (``long int`` type) of dimension at least 22, and ROUT
(``realtype`` type) of dimension at least 6. These arrays are owned
(and allocated) by the user and are passed as arguments to
:f:func:`FARKMALLOC()`. 

:ref:`FInterface.IOUTTable` and
:ref:`FInterface.ROUTTable` list the entries in these
arrays associated with the main ARKode solver, along with the
relevant ARKode function that is actually called to extract the
optional output.  Similarly,
:ref:`FInterface.DlsIOUTTable` lists the IOUT
entries associated with the main ARKDENSE and ARKBAND direct
linear solvers, and :ref:`FInterface.SpilsIOUTTable`
lists the IOUT entries associated with the main ARKSPGMR,
ARKSPBCG, ARKSPTFQMR, ARKSPFGMR and ARKPCG iterative linear solvers.

For more details on the optional inputs and outputs to ARKode, see
the sections :ref:`CInterface.OptionalInputs` and
:ref:`CInterface.OptionalOutputs`.



.. _FInterface.IOUTTable:

Table: Optional FARKODE integer outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

==============  ===============  =====================================================
IOUT Index      Optional output  ARKode function
==============  ===============  =====================================================
1               LENRW            :c:func:`ARKodeGetWorkSpace()`
2               LENIW            :c:func:`ARKodeGetWorkSpace()`
3               NST              :c:func:`ARKodeGetNumSteps()`
4               NST_STB          :c:func:`ARKodeGetNumExpSteps()`
5               NST_ACC          :c:func:`ARKodeGetNumAccSteps()`
6               NST_ATT          :c:func:`ARKodeGetNumStepAttempts()`
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Table: Optional ARKSPGMR, ARKSPBCG, ARKSPTFQMR, ARKSPFGMR and ARKPCG outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
---------------------------------------------


In addition to the optional inputs communicated through FARKSET*
calls and the optional outputs extracted from IOUT and ROUT,
the following user-callable routines are available: 

To obtain the error weight array EWT, containing the
multiplicative error weights used the WRMS norms, the user may call
the routine :f:func:`FARKGETERRWEIGHTS()` as follows:



.. f:subroutine:: FARKGETERRWEIGHTS(EWT, IER)
   
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
successful call to :f:func:`FARKODE()`, the user may call the routine
:f:func:`FARKGETESTLOCALERR()` as follows:



.. f:subroutine:: FARKGETESTLOCALERR(ELE, IER)
   
   Retrieves the current local truncation error estimate
   vector (interfaces with :c:func:`ARKodeGetEstLocalErrors()`).
      
   **Arguments:** 
      * ELE (``realtype``, output) -- array with the estimated local error vector
      * IER  (``int``, output) -- return flag  (0 if success, :math:`\ne 0` if an error)
      
   **Notes:**
   The array ELE, of length NEQ if using NVECTOR_SERIAL or NLOCAL
   if using NVECTOR_PARALLEL, must already have been declared by
   the user.  

