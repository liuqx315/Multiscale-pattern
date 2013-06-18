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
the following user-supplied functions.  Each of which is listed with
the corresponding internal interface function which calls it (and its
type within ARKode):

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
| :c:func:`FARKBJAC()`     | FARKBandJac            | :c:func:`ARKDlsBandJacFn()`       |
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


