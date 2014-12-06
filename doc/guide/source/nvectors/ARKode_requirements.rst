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
shows, for each function, which ARKode module uses the function.  
The ARKode column shows function usage within the main integrator
module,  while the remaining columns show function usage within 
the ARKode linear solvers, the ARKBANDPRE and ARKBBDPRE
preconditioner modules, and the FARKODE module.  Here ARKDLS stands
for ARKDENSE and ARKBAND; ARKSPILS stands for ARKSPGMR, ARKSPBCG,
ARKSPTFQMR, ARKSPFGMR and ARKPCG; and ARKSLS stands for ARKKLU and
ARKSUPERLUMT.

At this point, we should emphasize that the user does not need to know
anything about ARKode's usage of vector functions in order to use
ARKode.  Instead, this information is provided primarily for users
interested in constructing a custom ``N_Vector`` module.  We note that
a number of ``N_Vector`` functions from the section
:ref:`NVectors.Description` are not listed in the above table.
Therefore a user-supplied ``N_Vector`` module for ARKode could safely
omit these functions from their implementation. 



.. cssclass:: table-bordered

==================  =============  ======  ======  ========  ==========  =========  =============
Routine             ARKode         ARKDLS  ARKSLS  ARKSPILS  ARKBANDPRE  ARKBBDPRE  FARKODE
==================  =============  ======  ======  ========  ==========  =========  =============
N_VAbs              X                                                               X
N_VAddConst         X                                                               X
N_VClone            X                              X                                X
N_VCloneEmpty                                                                       X
N_VConst            X              X       X       X                                X
N_VDestroy          X                              X                                X
N_VDiv              X                              X                                X
N_VDotProd          X\ :sup:`(a)`                  X                                X\ :sup:`(a)`
N_VGetArrayPointer                 X       X                 X           X          X
N_VInv              X                                                               X
N_VLinearSum        X              X               X                                X
N_VMaxNorm          X                                                               X
N_VMin              X                                                               X
N_VProd                                            X
N_VScale            X              X       X       X         X           X          X
N_VSetArrayPointer                 X                                                X
N_VSpace            X\ :sup:`(b)`                                                   X\ :sup:`(b)`
N_VWrmsNorm         X              X               X         X           X          X
==================  =============  ======  ======  ========  ==========  =========  =============

(a) The :c:func:`N_VDotProd()` function is only used by the main
    ARKode integrator module when the fixed-point nonlinear solver is
    specified; when solving an explicit problem or when using a Newton
    solver with direct or sparse linear solver, it need not be
    supplied by the ``N_Vector`` implementation.

(b) The :c:func:`N_VSpace()` function is only informational, and need
    not be supplied by the ``N_Vector`` implementation.



