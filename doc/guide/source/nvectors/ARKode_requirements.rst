..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3



.. _NVectors.ARKode:

NVECTOR functions required by ARKode
==========================================

In the table below, we list the vector functions in the ``N_Vector``
module that are called within the ARKode package.  The table also
shows, for each function, which ARKode module uses the function:
the ARKde column shows function usage within the main integrator
module,  while the remaining columns show function usage within 
the ARKode linear solver modules (ARKSPILS stands for any of
ARKSPGMR, ARKSPBCG, ARKSPTFQMR, ARKSPFGMR or ARKPCG), the ARKBANDPRE
and ARKBBDPRE preconditioner modules, and the FARKODE module.


.. cssclass:: table-bordered

==================  ======  ========  =======  ========  ==========  =========  =======
Routine             ARKode  ARKDENSE  ARKBAND  ARKSPILS  ARKBANDPRE  ARKBBDPRE  FARKODE
==================  ======  ========  =======  ========  ==========  =========  =======
N_VAbs              X
N_VAddConst         X
N_VClone            X                          X
N_VCloneEmpty                                                                   X
N_VConst            X       X         X        X
N_VDestroy          X                          X                                X
N_VDiv              X                          X
N_VDotProd          X                          X
N_VGetArrayPointer          X         X                  X           X          X
N_VInv              X
N_VLinearSum        X       X                  X
N_VMaxNorm          X
N_VMin              X                                                           X
N_VProd                                        X
N_VScale            X       X         X        X         X           X
N_VSetArrayPointer          X                                                   X
N_VSpace            X
N_VWrmsNorm         X       X         X        X         X           X
==================  ======  ========  =======  ========  ==========  =========  =======


At this point, we should emphasize that the ARKode user does not need
to know anything about the usage of vector functions by the ARKode
code modules in order to use ARKode.  Instead, this information is
provided primarily for users interested in constructing a custom
``N_Vector`` module.  We note that a number of ``N_Vector`` functions
from the section :ref:`NVectors.Description` are not listed in the
above table.  Therefore a user-supplied ``N_Vector`` module for ARKode
could safely omit these functions from their implementation.
Furthermore, since the :c:func:`N_VSpace()` function is only
informational, while it should be supplied, the return values may be
arbitrary. 

