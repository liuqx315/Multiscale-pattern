.. _Installation:

ARKode Installation Procedure
=============================

The installation of ARKode is accomplished by installing the SUNDIALS
suite as a whole, according to the instructions that follow. The
installation procedure remains the same, whether or not the downloaded
file contains solvers other than ARKode. [#f1]_ 

The SUNDIALS suite (or individual solvers) are distributed as
compressed archives (``.tar.gz``). The name of the distribution
archive is of the form ``SOLVER-X.Y.Z.tar.gz``, where ``SOLVER`` is
one of: ``sundials``, ``arkode``, ``cvode``, ``cvodes``, ``ida``,
``idas``, or ``kinsol``, and ``X.Y.Z`` represents the version number
(of the SUNDIALS suite or of the individual solver). To begin the
installation, first uncompress and expand the sources, by issuing

.. code-block:: bash

   % tar -zxf SOLVER-X.Y.Z.tar.gz

This will extract source files under a directory ``SOLVER-X.Y.Z.``

Starting with version 2.4.0 of SUNDIALS, two installation methods are
provided: a standard LINUX/UNIX autotools-based method, and a new method
based on :index:`CMake`. Before providing detailed explanations on the
installation procedure for the two approaches, we begin with a few
common observations:

* In the remainder of this chapter, we make the following
  distinctions:

  ``SRCDIR`` 
     is the directory ``SOLVER-X.Y.Z`` created above; i.e. the
     directory containing the SUNDIALS sources.

  ``BUILDDIR`` 
     is the (temporary) directory under which SUNDIALS is built.

  ``INSTDIR`` 
     is the directory under which the SUNDIALS exported
     header files and libraries will be installed. Typically, header
     files are exported under a directory ``INSTDIR/include`` while
     libraries are installed under ``INSTDIR/lib``, with ``INSTDIR``
     specified at configuration time. 

* For the CMake-based installation, in-source builds are prohibited;
  in other words, the build directory ``BUILDDIR`` can **not** be the
  same as ``SRCDIR`` and such an attempt will lead to an error.  For
  autotools-based installation, in-source builds are allowed, although
  even in that case we recommend using a separate ``BUILDDIR``. Indeed,
  this prevents "polluting" the source tree and allows efficient
  builds for different configurations and/or options. 

* The installation directory ``INSTDIR`` can not be the same as the
  source directory ``SRCDIR``. 

* By default, only the libraries and header files are exported to the
  installation directory ``INSTDIR``.  If enabled by the user (with the
  appropriate option to ``configure`` or toggle for CMake), the
  examples distributed with SUNDIALS will be built together with the
  solver libraries but the installation step will result in exporting
  (by default in a subdirectory of the installation directory) the
  example sources and sample outputs together with automatically
  generated configuration files that reference the installed SUNDIALS
  headers and libraries. As such, these configuration files for the
  SUNDIALS examples can be used as "templates" for your own
  problems. The ``configure`` script will install makefiles. CMake
  installs ``CMakeLists.txt`` files and also (as an option available
  only under Unix/Linux) makefiles. Note that both installation
  approaches also allow the option of building the SUNDIALS examples
  without having to install them. (This can be used as a sanity check
  for the freshly built libraries.) 

* Even if generation of shared libraries is enabled, only static
  libraries are created for the FCMIX modules.  Because of the use of
  fixed names for the Fortran user-provided subroutines, FCMIX shared
  libraries would result in "undefined symbol" errors at link time.



.. rubric:: Footnotes

.. [#f1] Files for both the serial and parallel versions of ARKode are
	 included in the distribution. For users in a serial computing
	 environment, the files specific to parallel environments
	 (which may be deleted) are as follows: all files in
	 ``src/nvec_par/``;  ``nvector parallel.h`` (in
	 ``include/nvector/``); ``arkode_bbdpre.c``,
	 ``arkode_bbdpre_impl.h`` (in ``src/arkode/``);
	 ``arkode_bbdpre.h`` (in ``include/arkode/``); ``farkbbd.c``,
	 ``farkbbd.h`` (in ``src/arkode/fcmix/``); all files in
	 ``examples/arkode/parallel/``; all files in
	 ``examples/arkode/fcmix_parallel/``. (By "serial version" of
	 ARKode we mean the ARKode solver with the serial NVECTOR
	 module attached, and similarly for “parallel version”.) 



.. _Installation.Autotools:


Autotools-based installation
------------------------------

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
^^^^^^^^^^^^^^^^^^^^^^^^^

The installation procedure given above will generally work without
modification; however, if the system includes multiple MPI
implementations, then certain configure script-related options may be
used to indicate which MPI implementation should be used. Also, if the
user wants to use non-default language compilers, then, again, the
necessary shell environment variables must be appropriately
redefined. The remainder of this section provides explanations of
available configure script options. 


General options
"""""""""""""""""

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
""""""""""""""""""""""""""""""

:index:`--disable-fcmix <--disable-fcmix (autotools option)>`

   Using this option will disable all Fortran
   support. The FARKODE, FCVODE, FKINSOL, FIDA and FNVECTOR modules
   will not be built, regardless of availability. 

:index:`--with-fflags=ARG <--with-fflags=ARG (autotools option)>`

   Specify additional Fortran compilation flags.


Options for MPI support
""""""""""""""""""""""""""

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
""""""""""""""""""""""""""""""

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
"""""""""""""""""""""""""""""""""""

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
""""""""""""""""""""""""

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^

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



.. _Installation.CMake:

CMake-based installation
-------------------------

Support for CMake-based installation has been added to SUNDIALS
primarily to provide a platform-independent build system. Like
autotools, CMake can generate a Unix Makefile. Unlike autotools, CMake
can also create KDevelop, Visual Studio, and (Apple) XCode project
files from the same configuration file. In addition, CMake provides a
GUI front end and therefore the installation process is more
interactive than when using autotools. 

The installation options are very similar to the options mentioned
above (although their default values may differ
slightly). Practically, all configurations supported by the
autotools-based installation approach are also possible with CMake,
the only notable exception being cross-compilation, which is currently
not implemented in the CMake approach. 

The SUNDIALS build process requires CMake version 2.4.x or higher and
a working compiler. On Unix-like operating systems, it also requires
Make (and ``curses``, including its development libraries, for the GUI
front end to CMake, ``ccmake``), while on Windows it requires Visual
Studio. While many Linux distributions offer CMake, the version
included is probably out of date. Many new CMake features have been
added recently, and you should download the latest version from 
http://www.cmake.org/HTML/Download.html. Build instructions for
Cmake (only necessary for Unix-like systems) can be found on the CMake
website. Once CMake is installed, Linux/Unix users will be able to use
``ccmake``, while Windows users will be able to use ``CMakeSetup``. 

As noted above, when using CMake to configure, build and install
SUNDIALS, it is always required to use a separate build
directory. While in-source builds are possible, they are explicitly
prohibited by the SUNDIALS CMake scripts (one of the reasons being
that, unlike autotools, CMake does not provide a ``make distclean``
procedure and it is therefore difficult to clean-up the source tree
after an in-source build).



Configuring, building, and installing on Unix-like systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These instructions use :index:`ccmake` from the CMake installed
location. ``ccmake`` is a Curses based GUI for CMake. To run it, go to
the build directory and specify as an argument the source directory: 

.. code-block:: bash

   % mkdir BUILDDIR
   % cd BUILDDIR
   % ccmake SRCDIR

About ``ccmake``:

* Iterative process

  * Select values, run configure (press the ``<c>`` key)
  * Set the settings, run configure, set the settings, run configure,
    etc. 

* Repeat until all values are set and the generate option is available
  (press the ``<g>`` key) 
* Some variables (advanced variables) are not visible right away
* To see advanced varables, toggle to advanced mode (press the ``<t>``
  key) 
* To set a variable, move the cursor to the variable and press
  ``<enter>`` 

  * If it is a boolean (``ON/OFF``) it will flip the value
  * If it is string or file, it will allow editing of the string
  * For file and directories, the ``<tab>`` key can be used to
    complete 

* To search for a variable press the ``</>`` key, and to repeat the
  search, press the ``<n>`` key 

CMake will now generate makefiles including all dependencies and all
rules to build SUNDIALS on this system.  You should not, however, try
to move the build directory to another location on this system or to
another system. Once you have makefiles you should be able to just
type: 

.. code-block:: bash

   % make

To install SUNDIALS in the installation directory specified at
configuration time, simply run 

.. code-block:: bash

   % make install



Configuring, building, and installing on Windows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These instructions use :index:`CMakeSetup` from the CMake install
location. Make sure to select the appropriate source and the build
directory.  Also, make sure to pick the appropriate generator (e.g. on
Visual Studio 6, pick the Visual Studio 6 generator). Some CMake
versions will ask you to select the generator the first time you press
Configure instead of having a drop-down menu in the main dialog.

About ``CMakeSetup``:

* Iterative process

  * Select values, press the Configure button
  * Set the settings, run configure, set the settings, run configure,
    etc. 

* Repeat until all values are set and the ``OK`` button becomes available. 
* Some variables (advanced variables) are not visible right away
* To see advanced varables, toggle to advanced mode ("Show Advanced
  Values" toggle).  
* To set the value of a variable, click on that value.

  * If it is boolean (``ON/OFF``), a drop-down menu will appear for
    changing the value.  
  * If it is file or directory, an ellipsis button will appear ("...")
    on the far right of the entry.  Clicking this button will bring up
    the file or directory selection dialog.  
  * If it is a string, it will become an editable string.

CMake will now create Visual Studio project files. You should now be
able to open the SUNDIALS project (or workspace) file. Make sure to
select the appropriate build type (Debug, Release, ...). To build
SUNDIALS, simply build the ``ALL_BUILD`` target. To install SUNDIALS,
simply run the ``INSTALL`` target within the build system.



Configuration options
^^^^^^^^^^^^^^^^^^^^^^^^^

A complete list of all available options for a CMake-based SUNDIALS
configuration is provide below.  Note that the default values shown
are for a typical configuration on a Linux system and are provided as
illustration only. Some of them will be different on different
systems. 

:index:`BUILD_ARKODE <BUILD_ARKODE (CMake option)>` 

   Build the ARKODE library 

   Default: ``ON``

:index:`BUILD_CVODE <BUILD_CVODE (CMake option)>`

   Build the CVODE library

   Default: ``ON``

:index:`BUILD_CVODES <BUILD_CVODES (CMake option)>` 

   Build the CVODES library

   Default: ``ON``

:index:`BUILD_IDA <BUILD_IDA (CMake option)>` 

   Build the IDA library

   Default: ``ON``

:index:`BUILD_IDAS <BUILD_IDAS (CMake option)>` 

   Build the IDAS library

   Default: ``ON``

:index:`BUILD_KINSOL <BUILD_KINSOL (CMake option)>` 

   Build the KINSOL library

   Default: ``ON``

:index:`BUILD_SHARED_LIBS <BUILD_SHARED_LIBS (CMake option)>` 

   Build shared libraries

   Default: ``OFF``

:index:`BUILD_STATIC_LIBS <BUILD_STATIC_LIBS (CMake option)>` 

   Build static libraries

   Default: ``ON``

:index:`CMAKE_BUILD_TYPE <CMAKE_BUILD_TYPE (CMake option)>` 

   Choose the type of build, options are: 
   ``None`` (``CMAKE_C_FLAGS`` used), ``Debug``, ``Release``,
   ``RelWithDebInfo``, and ``MinSizeRel``

   Default:

:index:`CMAKE_C_COMPILER <CMAKE_C_COMPILER (CMake option)>` 

   C compiler

   Default: ``/usr/bin/gcc``

:index:`CMAKE_C_FLAGS <CMAKE_C_FLAGS (CMake option)>` 

   Flags for C compiler

   Default:

:index:`CMAKE_C_FLAGS_DEBUG <CMAKE_C_FLAGS_DEBUG (CMake option)>` 

   Flags used by the compiler during debug
   builds

   Default: ``-g``

:index:`CMAKE_C_FLAGS_MINSIZEREL <CMAKE_C_FLAGS_MINSIZEREL (CMake option)>` 

   Flags used by the compiler during
   release minsize builds

   Default: ``-Os -DNDEBUG``

:index:`CMAKE_C_FLAGS_RELEASE <CMAKE_C_FLAGS_RELEASE (CMake option)>` 

   Flags used by the compiler during release
   builds

   Default: ``-O3 -DNDEBUG``

:index:`CMAKE_BACKWARDS_COMPATIBILITY <CMAKE_BACKWARDS_COMPATIBILITY (CMake option)>` 

   For backwards compatibility, what
   version of CMake commands and syntax should this version of CMake
   allow. 

   Default: ``2.4``

:index:`CMAKE_Fortran_COMPILER <CMAKE_Fortran_COMPILER (CMake option)>` 

   Fortran compiler

   Default: ``/usr/bin/g77``

   Note: Fortran support (and all related options) are triggered only
   if either Fortran-C support is enabled (``FCMIX_ENABLE`` is ``ON``) or
   BLAS/LAPACK support is enabled (``LAPACK_ENABLE`` is ``ON``). 

:index:`CMAKE_Fortran_FLAGS <CMAKE_Fortran_FLAGS (CMake option)>` 

   Flags for Fortran compiler

   Default:

:index:`CMAKE_Fortran_FLAGS_DEBUG <CMAKE_Fortran_FLAGS_DEBUG (CMake option)>` 

   Flags used by the compiler during debug
   builds

   Default:

:index:`CMAKE_Fortran_FLAGS_MINSIZEREL <CMAKE_Fortran_FLAGS_MINSIZEREL (CMake option)>` 

   Flags used by the compiler during
   release minsize builds 

   Default:

:index:`CMAKE_Fortran_FLAGS_RELEASE <CMAKE_Fortran_FLAGS_RELEASE (CMake option)>` 

   Flags used by the compiler during
   release builds

   Default:

:index:`CMAKE_INSTALL_PREFIX <CMAKE_INSTALL_PREFIX (CMake option)>` 

   Install path prefix, prepended onto install
   directories

   Default: ``/usr/local``

   Note: The user must have write access to the location specified
   through this option. Exported SUNDIALS header files and libraries
   will be installed under subdirectories ``include`` and ``lib`` of
   ``CMAKE_INSTALL_PREFIX``, respectively. 

:index:`EXAMPLES_ENABLE <EXAMPLES_ENABLE (CMake option)>` 

   Build the SUNDIALS examples

   Default: ``OFF``

   Note: setting this option to ``ON`` will trigger additional options
   related to how and where example programs will be installed.

:index:`EXAMPLES_GENERATE_MAKEFILES <EXAMPLES_GENERATE_MAKEFILES (CMake option)>` 

   Create Makefiles for building the
   examples

   Default: ``ON``

   Note: This option is triggered only if enabling the building and
   installing of the example programs (i.e., both ``EXAMPLES_ENABLE``
   and ``EXAMPLEs_INSTALL`` are set to ``ON``) and if configuration is
   done on a Unix-like system. If enabled, makefiles for the
   compilation of the example programs (using the installed SUNDIALS
   libraries) will be automatically generated and exported to the
   directory specified by ``EXAMPLES_INSTALL_PATH``. 

:index:`EXAMPLES_INSTALL <EXAMPLES_INSTALL (CMake option)>` 

   Install example files

   Default: ``ON``

   Note: This option is triggered only if building example programs is
   enabled (``EXAMPLES_ENABLE`` is set to ``ON``). If the user
   requires installation of example programs then the sources and
   sample output files for all SUNDIALS modules that are currently
   enabled will be exported to the directory specified by
   ``EXAMPLES_INSTALL_PATH``. A CMake configuration script will also
   be automatically generated and exported to the same
   directory. Additionally, if the configuration is done under a
   Unix-like system, an additional option
   (``EXAMPLES_GENERATE_MAKEFILES``) will be triggered.  

:index:`EXAMPLES_INSTALL_PATH <EXAMPLES_INSTALL_PATH (CMake option)>` 

   Output directory for installing example
   files

   Default: ``/usr/local/examples``

   Note: The actual default value for this option will be an
   ``examples`` subdirectory created under ``CMAKE_INSTALL_PREFIX``.

:index:`EXAMPLES_USE_STATIC_LIBS <EXAMPLES_USE_STATIC_LIBS (CMake option)>` 

   Link examples using the static libraries 

   Default: ``OFF``

   Note: This option is triggered only if building shared libraries is
   enabled (``BUILD_SHARED_LIBS`` is ``ON``).

:index:`FCMIX_ENABLE <FCMIX_ENABLE (CMake option)>` 

   Enable Fortran-C support

   Default: ``OFF``

:index:`LAPACK_ENABLE <LAPACK_ENABLE (CMake option)>` 

   Enable LAPACK support

   Default: ``OFF``

   Note: Setting this option to ``ON`` will trigger the two additional
   options see below. 

:index:`LAPACK_LIBRARIES <LAPACK_LIBRARIES (CMake option)>` 

   LAPACK (and BLAS) libraries

   Default: ``/usr/lib/liblapack.so;/usr/lib/libblas.so``

:index:`LAPACK_LINKER_FLAGS <LAPACK_LINKER_FLAGS (CMake option)>` 

   LAPACK (and BLAS) required linker flags

   Default: ``-lg2c``

:index:`MPI_ENABLE <MPI_ENABLE (CMake option)>` 

   Enable MPI support

   Default: ``OFF``

   Note: Setting this option to ``ON`` will trigger several additional
   options related to MPI. 

:index:`MPI_MPICC <MPI_MPICC (CMake option)>` 

   ``mpicc`` program

   Default: ``/home/radu/apps/mpich1/gcc/bin/mpicc``

   Note: This option is triggered only if using MPI compiler scripts
   (``MPI_USE_MPISCRIPTS`` is ``ON``). 

:index:`MPI_MPIF77 <MPI_MPIF77 (CMake option)>` 

   ``mpif77`` program

   Default: ``/home/radu/apps/mpich1/gcc/bin/mpif77``

   Note: This option is triggered only if using MPI compiler scripts
   (``MPI_USE_MPISCRIPTS`` is ``ON``) and Fortran-C support is enabled
   (``FCMIX_ENABLE`` is ``ON``). 

:index:`MPI_INCLUDE_PATH <MPI_INCLUDE_PATH (CMake option)>` 

   Path to MPI header files

   Default: ``/home/radu/apps/mpich1/gcc/include``

   Note: This option is triggered only if not using MPI compiler
   scripts (``MPI_USE_MPISCRIPTS`` is ``OFF``).

:index:`MPI_LIBRARIES <MPI_LIBRARIES (CMake option)>` 

   MPI libraries

   Default: ``/home/radu/apps/mpich1/gcc/lib/libmpich.a``

   Note: This option is triggered only if not using MPI compiler
   scripts (``MPI_USE_MPISCRIPTS`` is ``OFF``).

:index:`MPI_USE_MPISCRIPTS <MPI_USE_MPISCRIPTS (CMake option)>` 

   Use MPI compiler scripts

   Default: ``ON``

:index:`SUNDIALS_PRECISION <SUNDIALS_PRECISION (CMake option)>` 

   Precision used in SUNDIALS, options are: ``double``, ``single`` or
   ``extended``

   Default: ``double``

:index:`USE_GENERIC_MATH <USE_GENERIC_MATH (CMake option)>` 

   Use generic (``stdc``) math libraries

   Default: ``ON``




Manually building SUNDIALS
-----------------------------

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



Installed libraries and exported header files
-----------------------------------------------

Using the standard SUNDIALS build system, the command

.. code-block:: bash

   % make install

will install the libraries under ``LIBDIR`` and the public header
files under ``INCLUDEDIR``. The default values for these directories
are ``INSTDIR/lib`` and ``INSTDIR/include``, respectively, but can
be changed using the configure script options ``--prefix``,
``--exec-prefix``, ``--includedir`` and ``--libdir`` (see the section
:ref:`Installation.Autotools`) or the appropriate CMake options (see
the section :ref:`Installation.CMake`). For example, a global
installation of SUNDIALS on a LINUX/UNIX system could be accomplished
using

.. code-block:: bash

   % configure --prefix=/opt/sundials-2.5.0

Although all installed libraries reside under ``LIBDIR``, the public
header files are further organized into subdirectories under
``INCLUDEDIR``. 

The installed libraries and exported header files are listed for
reference in the :ref:`Table: SUNDIALS libraries and header files
<Installation.Table>`. The file extension ``.LIB`` is typically ``.so``
for shared libraries and ``.a`` for static libraries. Note that, in
this table names are relative to ``LIBDIR`` for libraries and to
``INCLUDEDIR`` for header files.  

A typical user program need not explicitly include any of the shared
SUNDIALS header files from under the ``INCLUDEDIR/sundials``
directory since they are explicitly included by the appropriate solver
header files (e.g., ``arkode_dense.h`` includes
``sundials_dense.h``). However, it is both legal and safe to do so
(e.g., the functions declared in ``sundials_dense.h`` could be used in
building a preconditioner).



.. _Installation.Table:

Table: SUNDIALS libraries and header files

+--------------------------------+---------------------------------+
| Shared            Libraries    | n/a                             |
+--------------------------------+---------------------------------+
| Shared            Header files | sundials/sundials_config.h,     |
|                                | sundials/sundials_types.h,      |
|                                | sundials/sundials_math.h,       |
|                                | sundials/sundials_nvector.h,    |
|                                | sundials/sundials_fnvector.h,   |
|                                | sundials/sundials_direct.h,     |
|                                | sundials/sundials_lapack.h,     |
|                                | sundials/sundials_dense.h,      |
|                                | sundials/sundials_band.h,       |
|                                | sundials/sundials_iterative.h,  |
|                                | sundials/sundials_spgmr.h,      |
|                                | sundials/sundials_spbcgs.h,     |
|                                | sundials/sundials_sptfqmr.h     |
+--------------------------------+---------------------------------+
| Serial NVECTOR    Libraries    | libsundials_nvecserial.LIB,     |
|                                | libsundials_fnvecserial.a       |
+--------------------------------+---------------------------------+
| Serial NVECTOR    Header files | nvector/nvector_serial.h        |
+--------------------------------+---------------------------------+
| Parallel NVECTOR  Libraries    | libsundials_nvecparallel.LIB,   |
|                                | libsundials_fnvecparallel.a     |
+--------------------------------+---------------------------------+
| Parallel NVECTOR  Header files | nvector/nvector_parallel.h      |
+--------------------------------+---------------------------------+
| ARKODE            Libraries    | libsundials_arkode.LIB,         |
|                                | libsundials_farkode.a           |
+--------------------------------+---------------------------------+
| ARKODE            Header files | arkode/arkode.h,                |
|                                | arkode/arkode_impl.h,           |
|                                | arkode/arkode_direct.h,         |
|                                | arkode/arkode_lapack.h,         |
|                                | arkode/arkode_dense.h,          |
|                                | arkode/arkode_band.h,           |
|                                | arkode/arkode_spils.h,          |
|                                | arkode/arkode_spgmr.h,          |
|                                | arkode/arkode_spbcgs.h,         |
|                                | arkode/arkode_sptfqmr.h,        |
|                                | arkode/arkode_bandpre.h,        |
|                                | arkode/arkode_bbdpre.h          |
+--------------------------------+---------------------------------+
| CVODE             Libraries    | libsundials_cvode.LIB,          |
|                                | libsundials_fcvoce.a            |
+--------------------------------+---------------------------------+
| CVODE             Header files | cvode/cvode.h,                  |
|                                | cvode/cvode_impl.h,             |
|                                | cvode/cvode_direct.h,           |
|                                | cvode/cvode_lapack.h,           |
|                                | cvode/cvode_dense.h,            |
|                                | cvode/cvode_band.h,             |
|                                | cvode/cvode_diag.h,             |
|                                | cvode/cvode_spils.h,            |
|                                | cvode/cvode_spgmr.h,            |
|                                | cvode/cvode_spbcgs.h,           |
|                                | cvode/cvode_sptfqmr.h,          |
|                                | cvode/cvode_bandpre.h,          |
|                                | cvode/cvode_bbdpre.h            |
+--------------------------------+---------------------------------+
| CVODES            Libraries    | libsundials_cvodes.LIB          |
+--------------------------------+---------------------------------+
| CVODES            Header files | cvodes/cvodes.h,                |
|                                | cvodes/cvodes_impl.h,           |
|                                | cvodes/cvodes_direct.h,         |
|                                | cvodes/cvodes_lapack.h,         |
|                                | cvodes/cvodes_dense.h,          |
|                                | cvodes/cvodes_band.h,           |
|                                | cvodes/cvodes_diag.h,           |
|                                | cvodes/cvodes_spils.h,          |
|                                | cvodes/cvodes_spgmr.h,          |
|                                | cvodes/cvodes_spbcgs.h,         |
|                                | cvodes/cvodes_sptfqmr.h,        |
|                                | cvodes/cvodes_bandpre.h,        |
|                                | cvodes/cvodes_bbdpre.h          |
+--------------------------------+---------------------------------+
| IDA               Libraries    | libsundials_ida.LIB,            |
|                                | libsundials_fida.a              |
+--------------------------------+---------------------------------+
| IDA               Header files | ida/ida.h,                      |
|                                | ida/ida_impl.h,                 |
|                                | ida/ida_direct.h,               |
|                                | ida/ida_lapack.h,               |
|                                | ida/ida_dense.h,                |
|                                | ida/ida_band.h,                 |
|                                | ida/ida_spils.h,                |
|                                | ida/ida_spgmr.h,                |
|                                | ida/ida_spbcgs.h,               |
|                                | ida/ida_sptfqmr.h,              |
|                                | ida/ida_bbdpre.h                |
+--------------------------------+---------------------------------+
| IDAS              Libraries    | libsundials_idas.LIB            |
+--------------------------------+---------------------------------+
| IDAS              Header files | idas/idas.h,                    |
|                                | idas/idas_impl.h,               |
|                                | idas/idas_direct.h,             |
|                                | idas/idas_lapack.h,             |
|                                | idas/idas_dense.h,              |
|                                | idas/idas_band.h,               |
|                                | idas/idas_spils.h,              |
|                                | idas/idas_spgmr.h,              |
|                                | idas/idas_spbcgs.h,             |
|                                | idas/idas_sptfqmr.h,            |
|                                | idas/idas_bbdpre.h              |
+--------------------------------+---------------------------------+
| KINSOL            Libraries    | libsundials_kinsol.LIB,         |
|                                | libsundials_fkinsol.a           |
+--------------------------------+---------------------------------+
| KINSOL            Header files | kinsol/kinsol.h,                |
|                                | kinsol/kinsol_impl.h,           |
|                                | kinsol/kinsol_direct.h,         |
|                                | kinsol/kinsol_lapack.h,         |
|                                | kinsol/kinsol_dense.h,          |
|                                | kinsol/kinsol_band.h,           |
|                                | kinsol/kinsol_spils.h,          |
|                                | kinsol/kinsol_spgmr.h,          |
|                                | kinsol/kinsol_spbcgs.h,         |
|                                | kinsol/kinsol_sptfqmr.h,        |
|                                | kinsol/kinsol_bbdpre.h          |
+--------------------------------+---------------------------------+
