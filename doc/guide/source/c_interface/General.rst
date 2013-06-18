:tocdepth: 3

.. _CInterface.Headers:

Access to library and header files
===========================================

At this point, it is assumed that the installation of ARKode,
following the procedure described in the section :ref:`Installation`,
has been completed successfully. 

Regardless of where the user's application program resides, its
associated compilation and load commands must make reference to the
appropriate locations for the library and header files required by
ARKode. The relevant library files are 

- ``libdir/libsundialsarkode.lib``,
- ``libdir/libsundials_nvec*.lib`` (one or two files), 

where the file extension ``.lib`` is typically ``.so`` for shared
libraries and ``.a`` for static libraries. The relevant header files are
located in the subdirectories 

- ``incdir/include/arkode``
- ``incdir/include/sundials``
- ``incdir/include/nvector``

The directories ``libdir`` and ``incdir`` are the install library and
include directories, respectively.  For a default installation, these
are ``instdir/lib`` and ``instdir/include``, respectively, where ``instdir``
is the directory where SUNDIALS was installed (see the section
:ref:`Installation`).


.. _CInterface.DataTypes:

Data Types
===========================================

The ``sundials_types.h`` file contains the definition of the type
``realtype``, which is used by the SUNDIALS solvers for all
floating-point data.  The type ":index:`realtype`" can be ``float``,
``double``, or ``long double``, depending on how SUNDIALS was
installed (with the default being ``double``). The user can change the
precision of the SUNDIALS solvers floating-point arithmetic at the
configuration stage (see the section :ref:`Installation`). 

Additionally, based on the current precision, ``sundials_types.h``
defines the values :index:`BIG_REAL` to be the largest value
representable as a ``realtype``, :index:`SMALL_REAL` to be the
smallest positive value representable as a ``realtype``, and
:index:`UNIT_ROUNDOFF` to be the difference between 1.0 and the
minimum ``realtype`` greater than 1.0.  

Within SUNDIALS, real constants may be set to have the appropriate
precision by way of a macro called :index:`RCONST`.  It is this macro
that needs the ability to branch on the definition ``realtype``.  In
ANSI C, a floating-point constant with no suffix is stored as a
``double``. Placing the suffix "F" at the end of a floating point
constant makes it a ``float``, whereas using the suffix "L" makes it a
``long double``. For example,

.. code-block:: c

   #define A 1.0 
   #define B 1.0F 
   #define C 1.0L

defines ``A`` to be a ``double`` constant equal to 1.0, ``B`` to be a
``float`` constant equal to 1.0, and ``C`` to be a ``long double`` constant
equal to 1.0.  The macro call ``RCONST(1.0)`` automatically expands to
1.0 if ``realtype`` is ``double``, to ``1.0F`` if ``realtype`` is ``float``, or
to ``1.0L`` if ``realtype`` is ``long double``. SUNDIALS uses the ``RCONST``
macro internally to declare all of its floating-point constants. 

A user program which uses the type ``realtype`` and the ``RCONST`` macro
to handle floating-point constants is precision-independent except for
any calls to precision-specific standard math library functions (Our
example programs use both ``realtype`` and ``RCONST``).  Users can,
however, use the types ``double``, ``float``, or ``long double`` in their
code (assuming that this usage is consistent with the ``typedef`` for
``realtype``).  Thus, a previously existing piece of ANSI C code can use
SUNDIALS without modifying the code to use ``realtype``, so long as
the SUNDIALS libraries have been compiled using the same precision
(for details see the section :ref:`Installation`). 

SUNDIALS also defines a type ":index:`booleantype`", that can have
values ``TRUE`` and ``FALSE``, which is used for logic arguments
within the library.



Header Files
===========================================

The calling program must include several header files so that various
macros and data types can be used. The header file that is always
required is: 

- ``arkode.h``, the main header file for ARKode, which defines the
  several types and various constants, and includes function
  prototypes. 

Note that ``arkode.h`` includes ``sundials_types.h`` directly, which
defines the types ``realtype`` and ``booleantype`` and the
constants ``FALSE`` and ``TRUE``, so a user program does not need to
include ``sundials_types.h`` directly. 

The calling program must also include an NVECTOR implementation
header file (see the chapter :ref:`NVectors` for details).  For the two
NVECTOR implementations that are included in the ARKode package, the
corresponding header files are: 

* ``nvector_serial.h``, which defines the serial implementation
  NVECTOR_SERIAL; 
* ``nvector_parallel.h``, which defines the parallel (MPI)
  implementation, NVECTOR_PARALLEL.

Note that both these files in turn include the header file
``sundials_nvector.h`` which defines the abstract ``N_Vector`` data
type.

Finally, if the user includes a non-trivial implicit component to their
ODE system (and hence requires a Newton solver for the resulting
nonlinear systems of equations), then a linear solver module header
file will be required. The header files corresponding to the various
linear solvers availble for use with ARKode are: 

- ``arkode_dense.h``, which is used with the dense direct linear solver; 
- ``arkode_band.h``, which is used with the band direct linear solver;
- ``arkode_lapack.h``, which is used with LAPACK implementations of dense
  or band direct linear solvers; 
- ``arkode_spgmr.h``, which is used with the scaled, preconditioned GMRES
  Krylov linear solver SPGMR;
- ``arkode_spbcgs.h``, which is used with the scaled, preconditioned
  Bi-CGStab Krylov linear solver SPBCG;
- ``arkode_sptfqmr.h``, which is used with the scaled, preconditioned
  TFQMR Krylov solver SPTFQMR.
- ``arkode_pcg.h``, which is used with the preconditioned
  conjugate gradient linear solver PCG;

The header files for the dense and banded linear solvers (both
internal and LAPACK) include the file ``arkode_direct.h``, which defines
common functions.  This in turn includes a file (``sundials_direct.h``)
which defines the matrix type for these direct linear solvers
(``DlsMat``), as well as various functions and macros for acting on and
accessing entries of such matrices. 

The header files for the Krylov iterative solvers each include
``arkode_spils.h`` which defines common functions and which in turn
includes a header file (``sundials_iterative.h``) which enumerates the
preconditioning type and the choices for the Gram-Schmidt process (for
the SPGMR solver). 

Other headers may be needed, according to the choice of
preconditioner, etc.  For example, if preconditioning for an iterative
linear solver were performed using a block-diagonal
matrix, the header ``sundials_dense.h`` would need to be included for
access to the underlying generic dense linear solver to be used for
preconditioning. 
