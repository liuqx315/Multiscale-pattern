..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _FortranInterface:

=====================================================
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

.. _FInterface.Portability:

Important notes on portability
--------------------------------------

In this package, the names of the interface functions, and the names
of the Fortran user routines called by them, appear as dummy names
which are mapped to actual values by a series of definitions in the
header files ``farkode.h``, ``farkroot.h``, ``farkbp.h``, and
``farkbbd.h``. By default, those mapping definitions depend in turn on
the C macro ``F77_FUNC`` defined in the header file
``sundials_config.h`` and decided upon at configuration time (see
:ref:`Installation`). 


.. _FInterface.DataTypes:

Fortran Data Types
-------------------------

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



Details on the FARKode interface are provided in the following sub-sections:

.. toctree::
   :maxdepth: 1

   Routines
   Usage
   Optional_output
   Rootfinding
   Preconditioning
