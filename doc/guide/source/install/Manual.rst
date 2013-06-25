:tocdepth: 3

.. _Installation.Manual:

Manually building SUNDIALS
================================

With the addition of CMake support, the installation of the SUNDIALS
package on almost any platform was greatly simplified. However, if for
whatever reason, neither of the two procedures described above is
convenient (for example for users who prefer to own the build process
or otherwise incorporate SUNDIALS or one of its solvers in a larger
project with its own build system), we provide here a few directions
for a completely manual installation. 

The following files are required to compile a SUNDIALS solver module:

* public header files are located under ``SRCDIR/include/SOLVER``
* implementation header files and source files are located under
  ``SRCDIR/src/SOLVER``
* (optional) Fortran/C interface files are located under
  ``SRCDIR/src/SOLVER/fcmix`` 
* shared public header files are located under 
  ``SRCDIR/include/sundials``
* shared source files are located under ``SRCDIR/src/sundials``
* (optional) NVECTOR_SERIAL header and source files are located under 
  ``SRCDIR/include/nvector`` and ``SRCDIR/src/nvec_ser``
* (optional) NVECTOR_PARALLEL header and source are files located
  under ``SRCDIR/include/nvector`` and ``SRCDIR/src/nvec_par``
* the configuration header file, ``sundials_config.h`` (see below)

A sample header file that, appropriately modified, can be used as
``sundials_config.h`` (otherwise created automatically by the
configure or CMake scripts), is provided below. 

.. code-block:: c

   /* SUNDIALS configuration header file */
   #define SUNDIALS_PACKAGE_VERSION "2.5.0"

   #define SUNDIALS_F77_FUNC(name,NAME) name ## _
   #define SUNDIALS_F77_FUNC_(name,NAME) name ## _

   #define SUNDIALS_DOUBLE_PRECISION 1

   #define SUNDIALS_USE_GENERIC_MATH

   #define SUNDIALS_BLAS_LAPACK 1

   #define SUNDIALS_MPI_COMM_F2C 1

   #define SUNDIALS_EXPORT

The various preprocessor macros defined within ``sundials_config.h``
have the following uses: 

* Precision of the SUNDIALS ``realtype`` type

  Only one of the macros :index:`SUNDIALS_SINGLE_PRECISION`,
  :index:`SUNDIALS_DOUBLE_PRECISION` and
  :index:`SUNDIALS_EXTENDED_PRECISION` should be defined to indicate
  if the SUNDIALS ``realtype`` type is   an alias for ``float``,
  ``double``, or ``long double``, respectively. 

* Use of generic math functions

  If :index:`SUNDIALS_USE_GENERIC_MATH` is defined, then the functions
  in ``sundials_math.h`` and ``sundials_math.c`` will use the ``pow``,
  ``sqrt``, ``fabs``, and ``exp`` functions from the standard math
  library (see ``math.h``), regardless of the definition of
  ``realtype``. Otherwise, if ``realtype`` is defined to be an alias
  for the ``float`` C-type, then SUNDIALS will use ``powf``,
  ``sqrtf``, ``fabsf``, and ``expf``. If ``realtype`` is instead
  defined to be a synonym for the ``long double`` C-type, then
  ``powl``, ``sqrtl``, ``fabsl``, and ``expl`` will be used. 

  Note: Although the ``powf/powl``, ``sqrtf/sqrtl``,
  ``fabsf/fabsl``, and ``expf/expl`` routines are not
  specified in the ANSI C standard, they are ISO C99
  requirements. Consequently, these routines will only be used if
  available. 

* Fortran name-mangling scheme

  The macros given below are used to transform the C-language function
  names defined in the Fortran-C interface modules in a manner
  consistent with the preferred Fortran compiler, thus allowing native
  C functions to be called from within a Fortran subroutine. The
  name-mangling scheme is specified by appropriately defining the
  following parameterized macros (using the stringization operator,
  ``##``, if necessary): 

  * :index:`SUNDIALS_F77_FUNC(name,NAME)`
  * :index:`SUNDIALS_F77_FUNC_(name,NAME)`

  For example, to specify that mangled C-language function names
  should be lowercase with one underscore appended, include

  .. code-block:: c

     #define SUNDIALS_F77_FUNC(name,NAME) name ## _
     #define SUNDIALS_F77_FUNC_(name,NAME) name ## _

  in the ``sundials_config.h`` header file.

* Availability of BLAS/LAPACK libraries

  If working libraries for BLAS and LAPACK are available, then the
  macro :index:`SUNDIALS_BLAS_LAPACK` should be set to 1; otherwise it 
  should have the value 0.

* Use of an MPI communicator other than ``MPI_COMM_WORLD`` in Fortran 

  If the macro :index:`SUNDIALS_MPI_COMM_F2C` is defined, then the MPI
  implementation used to build SUNDIALS defines the type ``MPI_Fint``
  and the function ``MPI_Comm_f2c``, and it is possible to use MPI
  communicators other than ``MPI_COMM_WORLD`` with the Fortran-C
  interface modules. 

* The macro :index:`SUNDIALS_EXPORT` is used when marking SUNDIALS API
  functions for export/import. When building shared SUNDIALS libraries
  under Windows, use 

  .. code-block:: c

     #define SUNDIALS_EXPORT __declspec(dllexport)

  When linking to shared SUNDIALS libraries under Windows, use

  .. code-block:: c

     #define SUNDIALS_EXPORT __declspec(dllimport)

  In all other cases (other platforms or static libraries under
  Windows), the ``SUNDIALS_EXPORT`` macro is empty.

