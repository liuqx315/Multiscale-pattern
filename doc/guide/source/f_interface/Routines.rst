:tocdepth: 3


.. _FInterface.Routines:

FARKODE routines
===========================

The user-callable functions comprising the FARKODE solver interface,
with the corresponding ARKODE functions, are as follows:

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

  - :c:func:`FARKADAPTSET()` interfaces to :c:func:`ARKodeSetAdaptivityFn()`.

  - :c:func:`FARKEXPSTABSET()` interfaces to :c:func:`ARKodeSetStabilityFn()`.

  ..
     - :c:func:`FARKSETDIAGNOSTICS()` interfaces to :c:func:`ARKodeSetDiagnostics()`.

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

  - :c:func:`FARKLAPACKDENSE()` interfaces to :c:func:`ARKLapackDense()`.

  - :c:func:`FARKDENSESETJAC()` interfaces to :c:func:`ARKDlsSetDenseJacFn()`.

  - :c:func:`FARKBAND()` interfaces to :c:func:`ARKBand()`.

  - :c:func:`FARKLAPACKBAND()` interfaces to :c:func:`ARKLapackBand()`.

  - :c:func:`FARKBANDSETJAC()` interfaces to :c:func:`ARKDlsSetBandJacFn()`.

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

  - :c:func:`FARKPCG()` interfaces to :c:func:`ARKPcg()` and the PCG optional input
    functions (see :ref:`CInterface.ARKSpilsInputTable`).

  - :c:func:`FARKPCGREINIT()` interfaces to the PCG optional input
    functions.

  - :c:func:`FARKSPILSSETJAC()` interfaces to :c:func:`ARKSpilsSetJacTimesVecFn()`.

  - :c:func:`FARKSPILSSETPREC()` interfaces to :c:func:`ARKSpilsSetPreconditioner()`.


As with the native C interface, the FARKode solver interface requires
user-supplied functions to specify the ODE problem to be solved.  In
contrast to the case of direct use of ARKode, and of most Fortran ODE
solvers, the names of all user-supplied routines here are fixed, in
order to maximize portability for the resulting mixed-language program. 
As a result, whether using a purely implicit, purely explicit, or
mixed implicit-explicit solver, two routines must be provided by the
user (though one of which may do nothing):

.. cssclass:: table-bordered

+--------------------------+-----------------------------------+
| FARKODE routine          | ARKode interface                  |
| (FORTRAN, user-supplied) | function type                     |
+==========================+===================================+
| :c:func:`FARKIFUN()`     | :c:func:`ARKRhsFn()`              |
+--------------------------+-----------------------------------+
| :c:func:`FARKEFUN()`     | :c:func:`ARKRhsFn()`              |
+--------------------------+-----------------------------------+

In addition, as with the native C interface a user may provide
additional routines to assist in the solution process.  Each of the
following user-supplied routines is activated by calling the specified
"activation" routine: 

.. cssclass:: table-bordered

+--------------------------+-----------------------------------+-------------------------------+
| FARKODE routine          | ARKode interface                  | FARKODE "activation" routine  |
| (FORTRAN, user-supplied) | function type                     |                               |
+==========================+===================================+===============================+
| :c:func:`FARKDJAC()`     | :c:func:`ARKDlsDenseJacFn()`      | :c:func:`FARKDENSESETJAC()`   |
+--------------------------+-----------------------------------+-------------------------------+
| :c:func:`FARKBJAC()`     | :c:func:`ARKDlsBandJacFn()`       |  :c:func:`FARKBANDSETJAC()`   |
+--------------------------+-----------------------------------+-------------------------------+
| :c:func:`FARKPSET()`     | :c:func:`ARKSpilsPrecSetupFn()`   |  :c:func:`FARKSPILSSETPREC()` |
+--------------------------+-----------------------------------+-------------------------------+
| :c:func:`FARKPSOL()`     | :c:func:`ARKSpilsPrecSolveFn()`   |  :c:func:`FARKSPILSSETPREC()` |
+--------------------------+-----------------------------------+-------------------------------+
| :c:func:`FARKJTIMES()`   | :c:func:`ARKSpilsJacTimesVecFn()` |  :c:func:`FARKSPILSSETJAC()`  |
+--------------------------+-----------------------------------+-------------------------------+
| :c:func:`FARKEWT()`      | :c:func:`ARKEwtFn()`              |  :c:func:`FARKEWTSET()`       |
+--------------------------+-----------------------------------+-------------------------------+
| :c:func:`FARKADAPT()`    | :c:func:`ARKAdaptFn()`            |  :c:func:`FARKADAPTSET()`     |
+--------------------------+-----------------------------------+-------------------------------+
| :c:func:`FARKEXPSTAB()`  | :c:func:`ARKExpStabFn()`          |  :c:func:`FARKEXPSTABSET()`   |
+--------------------------+-----------------------------------+-------------------------------+
