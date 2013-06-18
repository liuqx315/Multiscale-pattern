:tocdepth: 3

.. _Installation:

=============================
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


Further details on the autotools- and CMake-based installation
procedures, instructions for manual compilation, and a roadmap of the
resulting installed libraries and exported header files, are provided
in the following subsections: 

.. toctree::
   :maxdepth: 1

   Autotools
   CMake
   Manual
   Installed




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


