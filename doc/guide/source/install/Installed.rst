..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _Installation.Results:

Installed libraries and exported header files
====================================================

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
---------------------------------------------

.. cssclass:: table-bordered

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
|                                | sundials/sundials_spfgmr.h,     |
|                                | sundials/sundials_pcg.h,        |
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
|                                | arkode/arkode_spfgmr.h,         |
|                                | arkode/arkode_pcg.h,            |
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
