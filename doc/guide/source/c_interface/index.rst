..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _CInterface:

=======================================
Using ARKode for C and C++ Applications
=======================================

This chapter is concerned with the use of ARKode for the solution
of initial value problems (IVPs) in a C or C++ language setting.  The
following sections treat the header files and the layout of the user's
main program, and provide descriptions of the ARKode user-callable 
functions and user-supplied functions. 

The example programs described in the companion document [R2013]_ may
be helpful. Those codes may be used as templates for new codes and are
included in the ARKode package ``examples`` subdirectory.

Users with applications written in Fortran should see the chapter
:ref:`FortranInterface`, that describes the Fortran/C interface
module, and may look to the Fortran example programs also described in
the companion document [R2013]_.  These codes are also located in the
ARKode package ``examples`` directory.

The user should be aware that not all linear solver and
preconditioning modules are compatible with all NVECTOR
implementations.  For example, NVECTOR_PARALLEL is not compatible with
the direct dense or direct band linear solvers since these linear
solver modules need to form the complete system Jacobian on a single
processor.  Specifically, the following ARKode modules can only be
used with NVECTOR_SERIAL, NVECTOR_OPENMP or NVECTOR_PTHREADS: ARKDENSE,
ARKBAND (using either the internal or the LAPACK implementation),
ARKKLU, ARKSUPERLUMT and ARKBANDPRE. Also, the preconditioner
module ARKBBDPRE can only be used with NVECTOR_PARALLEL. 

ARKode uses various constants for both input and output. These are
defined as needed in this chapter, but for convenience the full list
is provided separately in the section :ref:`Constants`. 

The relevant information on using ARKode's C and C++ interfaces is
detailed in the following sub-sections:

.. toctree::
   :maxdepth: 1

   General
   Skeleton
   User_callable
   User_supplied
   Preconditioners
