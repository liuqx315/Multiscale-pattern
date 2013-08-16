:tocdepth: 3

.. _Installation.CMake:

CMake-based installation
======================================

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
----------------------------------------------------------------

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
----------------------------------------------------------------

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
----------------------------------------------------------------

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

