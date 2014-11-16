..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------
   ARKode documentation master file, created by
   sphinx-quickstart on Sat Dec 22 20:38:03 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

====================================
ARKode Documentation
====================================

This is the documentation for ARKode, an adaptive step time
integration package for stiff, nonstiff and multi-rate systems of
ordinary differential equations (ODEs).  
The ARKode solver is a component of the `SUNDIALS
<https://computation.llnl.gov/casc/sundials/main.html>`_ suite of
nonlinear and differential/algebraic equation solvers. It is designed
to have a similar user experience to the `CVODE
<https://computation.llnl.gov/casc/sundials/description/description.html#descr_cvode>`_
solver, including user modes to allow adaptive integration to specified
output times, return after each internal step and root-finding
capabilities, and for calculations in serial and using either
shared-memory parallelism (via OpenMP or Pthreads) or
distributed-memory parallelism (via MPI). The default integration and
solver options should apply to most users, though complete control
over all internal parameters and time adaptivity algorithms is enabled
through optional interface routines.

ARKode is written in C, with C++ and Fortran interfaces.  

Due to its similarities in both function and design with the CVODE
package, a significant portion of this documentation has been directly
adapted from the CVODE documentation [HS2012]_. 

ARKode is developed by `Southern Methodist University
<http://www.smu.edu>`_, with support by the `US Department of Energy
<http://www.doe.gov>`_ through the `FASTMath
<http://www.fastmath-scidac.org/>`_ SciDAC Institute, under subcontract
B598130 from `Lawrence Livermore National Laboratory
<http://www.llnl.gov>`_. 



.. only:: html

   Documentation sections:

.. toctree::
   :maxdepth: 1

   Introduction
   Mathematics
   Organization
   c_interface/index.rst
   f_interface/index.rst
   nvectors/index.rst
   linear_solvers/index.rst
   Install
   Constants
   Butcher
   References

.. only:: html

   * :ref:`genindex`
   * :ref:`search`


