:tocdepth: 3

.. _NVectors:

=======================
Vector Data Structures
=======================

The SUNDIALS library comes packaged with two NVECTOR implementations,
one designed for serial simulations and the second for
distributed-memory parallel simulations.  Both implementations assume
that the process-local data is stored contiguously, and they in turn
provide a variety of standard vector algebra operations that may be
performed on the data.  

In addition, SUNDIALS provides a simple interface for generic vectors
(akin to a C++ *abstract base class*).  All of the major SUNDIALS
solvers (CVODE, IDA, KINSOL, ARKODE) in turn are constructed to only
depend on these generic vector operations, making them immediately
extensible to new user-defined vector objects.  The only exceptions to
this rule relate to the dense and banded linear system solvers, and
band linear system preconditioners, that will be discussed in the next
section, :ref:`LinearSolvers`.

Details on the generic NVECTOR module, as well as the serial and
parallel vector modules provided by SUNDIALS, are provided in the
following sub-sections:

.. toctree::
   :maxdepth: 1

   Description
   NVector_Serial
   NVector_Parallel
   ARKode_requirements
