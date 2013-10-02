..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _Installation.Autotools:

Autotools-based installation
=========================================

The installation procedure outlined below will work on commodity
LINUX/UNIX systems without modification.  However, users are still
encouraged to carefully read this entire section before attempting to
install the SUNDIALS suite, in case non-default choices are desired
for compilers, compilation options, installation location, etc. The
user may invoke the configuration script with the help flag to view a
complete listing of available options, by issuing the command 

.. code-block:: bash

   % ./configure --help

from within ``SRCDIR``.

The installation steps for SUNDIALS can be as simple as the following:

.. code-block:: bash

   % cd SRCDIR
   % ./configure
   % make
   % make install

in which case the SUNDIALS header files and libraries are installed
under ``/usr/local/include`` and ``/usr/local/lib``,
respectively. Note that, by default, the example programs are not
built and installed.  To delete all temporary files created by
building SUNDIALS, issue 

.. code-block:: bash

   % make clean

To prepare the SUNDIALS distribution for a new install (using, for
example, different options and/or installation destinations), issue 

.. code-block:: bash

   % make distclean

The above steps are for an "in-source" build. For an "out-of-source"
build (recommended), the procedure is simply:

.. code-block:: bash

   % cd BUILDDIR
   % SRCDIR/configure
   % make
   % make install

Note that, in this case, ``make clean`` and ``make distclean`` are
irrelevant. Indeed, if disk space is a priority, the entire ``BUILDDIR``
can be purged after the installation completes. For a new install, a
new ``BUILDDIR`` directory can be created and used.




Configuration options
------------------------------

The installation procedure given above will generally work without
modification; however, if the system includes multiple MPI
implementations, then certain configure script-related options may be
used to indicate which MPI implementation should be used. Also, if the
user wants to use non-default language compilers, then, again, the
necessary shell environment variables must be appropriately
redefined. The remainder of this section provides explanations of
available configure script options. 


General options
^^^^^^^^^^^^^^^^^^^^^^^^^

:index:`--prefix=PREFIX <--prefix=PREFIX (autotools option)>`

   Location for architecture-independent files.

   Default: ``PREFIX=/usr/local``

:index:`--exec-prefix=EPREFIX <--exec-prefix=EPREFIX (autotools option)>`

   Location for architecture-dependent files.

   Default: ``EPREFIX=/usr/local``

:index:`--includedir=DIR <--includedir=DIR (autotools option)>`

   Alternate location for installation of header files. 

   Default: ``DIR=PREFIX/include``

:index:`--libdir=DIR <--libdir=DIR (autotools option)>`

   Alternate location for installation of libraries.

   Default: ``DIR=EPREFIX/lib``

:index:`--disable-solver <--disable-solver (autotools option)>`

   Although each existing solver module is built 
   by default, support for a given solver can be explicitly disabled
   using this option. The valid values for solver are: arkode, cvode,
   cvodes, ida, idas, and kinsol.

:index:`--enable-examples <--enable-examples (autotools option)>`
 
   Available example programs are not built by 
   default. Use this option to enable compilation of all pertinent
   example programs. Upon completion of the ``make`` command, the
   example executables will be created under solver-specific
   subdirectories of ``BUILDDIR/examples``: 

   ``BUILDDIR/examples/SOLVER/serial``: serial C examples

   ``BUILDDIR/examples/SOLVER/parallel``: parallel C examples

   ``BUILDDIR/examples/SOLVER/fcmix_serial``: serial Fortran examples

   ``BUILDDIR/examples/SOLVER/fcmix_parallel``: parallel Fortran
   examples

   `Note`: Some of these subdirectories may not exist depending upon
   the solver and/or the configuration options given. 

:index:`--with-examples-instdir=EXINSTDIR <--with-examples-instdir=EXINSTDIR (autotools option)>`
 
   Alternate location for example
   executables and sample output files (valid only if examples are
   enabled). Note that installation of example files can be completely
   disabled by issuing ``EXINSTDIR=no`` (in case building the examples
   is desired only as a test of the SUNDIALS libraries). 

   Default: ``DIR=EPREFIX/examples``

:index:`--with-cppflags=ARG <--with-cppflags=ARG (autotools option)>`

   Specify additional C preprocessor flags (e.g.,
   ``--with-cppflags=-I<include_dir``> if necessary header files are
   located in nonstandard locations). 

:index:`--with-cflags=ARG <--with-cflags=ARG (autotools option)>`

   Specify additional C compilation flags.

:index:`--with-ldflags=ARG <--with-ldflags=ARG (autotools option)>`

   Specify additional linker flags (e.g., 
   ``--with-ldflags=-L<lib_dir>`` if required libraries are located in
   nonstandard locations). 

:index:`--with-libs=ARG <--with-libs=ARG (autotools option)>`

   Specify additional libraries to be used (e.g.,
   ``--with-libs=-lfoo`` to link with the library named ``libfoo.a`` or
   ``libfoo.so``). 

:index:`--with-precision=ARG <--with-precision=ARG (autotools option)>`

   By default, SUNDIALS will define a real number
   (internally referred to as ``realtype``) to be a double-precision
   floating-point numeric data type (``double`` C-type); however, this
   option may be used to build SUNDIALS with ``realtype`` defined
   instead as a single-precision floating-point numeric data type
   (``float`` C-type) if ``--with-precision=single``, or as a ``long
   double`` C-type if ``--with-precision=extended``. 

   Default ``double``:

   Users should not build SUNDIALS with support for single-precision
   floating-point arithmetic on 32- or 64-bit systems.  This will
   almost certainly result in unreliable numerical solutions. The
   configuration option ``--with-precision=single`` is intended for
   systems on which single-precision arithmetic involves at least 14
   decimal digits. 


Options for Fortran support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:index:`--disable-fcmix <--disable-fcmix (autotools option)>`

   Using this option will disable all Fortran
   support. The FARKODE, FCVODE, FKINSOL, FIDA and FNVECTOR modules
   will not be built, regardless of availability. 

:index:`--with-fflags=ARG <--with-fflags=ARG (autotools option)>`

   Specify additional Fortran compilation flags.


Options for MPI support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following configuration options are only applicable to the
parallel SUNDIALS packages: 

:index:`--disable-mpi <--disable-mpi (autotools option)>`

   Using this option will completely disable MPI support.

:index:`--with-mpicc=ARG <--with-mpicc=ARG (autotools option)>`

:index:`--with-mpif77=ARG <--with-mpif77=ARG (autotools option)>`

   By default, the configuration utility script will
   use the MPI compiler scripts named ``mpicc`` and ``mpif77`` to
   compile the parallelized SUNDIALS subroutines; however, for reasons
   of compatibility, different executable names may be specified via
   the above options. Also, ``--with-mpif77=no`` can be used to
   disable the use of MPI compiler scripts, thus causing the serial C
   and Fortran compilers to be used to compile the parallelized
   SUNDIALS functions and examples. 

:index:`--with-mpi-root=MPIDIR <--with-mpi-root=MPIDIR (autotools option)>`

   This option may be used to specify which MPI
   implementation should be used. The SUNDIALS configuration script
   will automatically check under the subdirectories ``MPIDIR/include``
   and ``MPIDIR/lib`` for the necessary header files and
   libraries. The subdirectory ``MPIDIR/bin`` will also be searched
   for the C and Fortran MPI compiler scripts, unless the user
   uses ``--with-mpicc=no`` or ``--with-mpif77=no``.

:index:`--with-mpi-incdir=INCDIR <--with-mpi-incdir=INCDIR (autotools option)>`

:index:`--with-mpi-libdir=LIBDIR <--with-mpi-libdir=LIBDIR (autotools option)>`

:index:`--with-mpi-libs=LIBS <--with-mpi-libs=LIBS (autotools option)>`

   These options may be used if the user would
   prefer not to use a preexisting MPI compiler script, but instead
   would rather use a serial complier and provide the flags necessary
   to compile the MPI-aware subroutines in SUNDIALS.

   Often an MPI implementation will have unique library names and so
   it may be necessary to specify the appropriate libraries to use
   (e.g., ``--with-mpi-libs=-lmpich``). 

   Default: ``INCDIR=MPIDIR/include`` and ``LIBDIR=MPIDIR/lib``

:index:`--with-mpi-flags=ARG <--with-mpi-flags=ARG (autotools option)>`

   Specify additional MPI-specific flags.


Options for library support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, only static libraries are built, but the following option
may be used to build shared libraries on supported platforms.

:index:`--enable-shared <--enable-shared (autotools option)>`

   Using this particular option will result in both
   static and shared versions of the available SUNDIALS libraries
   being built if the system supports shared libraries. To build only
   shared libraries also specify ``--disable-static``.

   Note: The FARKODE, FCVODE, FKINSOL and FIDA libraries can only be
   built as static libraries because they contain references to
   externally defined symbols, namely user-supplied Fortran
   subroutines.  Although the Fortran interfaces to the serial and
   parallel implementations of the supplied NVECTOR module do not
   contain any unresolvable external symbols, the libraries are still
   built as static libraries for the purpose of consistency.


Options for BLAS/LAPACK support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``configure`` script will attempt to automatically determine the
proper libraries to be linked for support of the BLAS/LAPACK linear
solver module. If these are not found, or if BLAS and/or LAPACK
libraries are installed in a non-standard location, the following
options can be used: 

:index:`--with-blas=BLASDIR <--with-blas=BLASDIR (autotools option)>`

   Specify the BLAS library.

   Default: none

:index:`--with-lapack=LAPACKDIR <--with-lapack=LAPACKDIR (autotools option)>`

   Specify the LAPACK library.

   Default: none


Environment variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following environment variables can be locally (re)defined for use
during the configuration of SUNDIALS. See the next section for
illustrations of these. 

:index:`CC <CC (env. variable)>`

:index:`F77 <F77 (env. variable)>`

   Since the configuration script uses the first C and Fortran
   compilers found in the current executable search path, then each
   relevant shell variable (CC and F77) must be locally (re)defined in
   order to use a different compiler. For example, to use ``xcc``
   (executable name of chosen compiler) as the C language compiler,
   use ``CC=xcc`` in the ``configure`` step. 

:index:`CFLAGS <CFLAGS (env. variable)>`

:index:`FFLAGS <FFLAGS (env. variable)>`

   Use these environment variables to override the default C
   and Fortran compilation flags. 




Configuration examples
--------------------------------------

The following examples are meant to help demonstrate proper usage of
the configure options. 

To build SUNDIALS using the default C and Fortran compilers, and
default ``mpicc`` and ``mpif77`` parallel compilers, enable
compilation of examples, and install libraries, headers, and example
sources under appropriate subdirectories of
``/home/myname/sundials/``, use 

.. code-block:: bash

   % configure --prefix=/home/myname/sundials --enable-examples

To disable installation of the examples, use:

.. code-block::  bash

   % configure --prefix=/home/myname/sundials \
               --enable-examples --with-examples-instdir=no

The following example builds SUNDIALS using ``gcc`` as the serial C
compiler, ``gfortran`` as the serial Fortran compiler, the default
``mpicc`` as the parallel C compiler, the default ``mpif77`` as the
parallel Fortran compiler, and appends the ``-O3`` compilaton flag to
the list of default flags: 

.. code-block:: bash

   % configure CC=gcc F77=gfortran --with-cflags=-O3 --with-fflags=-O3 \
               --with-mpicc=mpicc --with-mpif77=mpif77

The next example again builds SUNDIALS using ``gcc`` as the serial C
compiler, but the ``--with-mpicc=no`` option explicitly disables the
use of the corresponding MPI compiler script. In addition, since the 
``--with-mpi-root`` option is given, the compilation flags 
``-I/usr/apps/mpich/1.2.4/include`` and
``-L/usr/apps/mpich/1.2.4/lib`` are passed to ``gcc`` when compiling
the MPI-enabled functions. The ``--with-mpi-libs`` option is required
so that the configure script can check if ``gcc`` can link with the 
appropriate MPI library. The ``--disable-lapack`` option explicitly
disables support for BLAS/LAPACK, while the ``--disable-fcmix``
explicitly disables building the FCMIX interfaces. Note that, because
of the last two options, no Fortran-related settings are checked for.

.. code-block:: bash

   % configure CC=gcc --with-mpicc=no \
               --with-mpi-root=/usr/apps/mpich/1.2.4 \
               --with-mpi-libs=-lmpich \
               --disable-lapack --disable-fcmix

Finally, a minimal configuration and installation of SUNDIALS in
``/home/myname/sundials/`` (serial only, no Fortran support, no
examples) can be obtained with: 

.. code-block:: bash

   % configure --prefix=/home/myname/sundials \
               --disable-mpi --disable-lapack --disable-fcmix
