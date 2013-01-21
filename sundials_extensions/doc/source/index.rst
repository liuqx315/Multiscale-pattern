.. ARKode documentation master file, created by
   sphinx-quickstart on Sat Dec 22 20:38:03 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the ARKode documentation!
====================================

This is the documentation for ARKode, an adaptive step time
integration package for stiff, nonstiff and multi-rate systems of
ordinary differential equations (ODEs).  
The ARKode solver is a component of the `SUNDIALS
<https://computation.llnl.gov/casc/sundials/main.html>`_ suite of
nonlinear and differential/algebraic equation solvers. It is designed
to have a similar user experience to the `CVODE
<https://computation.llnl.gov/casc/sundials/description/description.html#descr_cvode>`_
solver, with user modes to allow adaptive integration to specified
output times, return after each internal step and root-finding
capabilities, for calculations both in serial and parallel (via
MPI). The default integration and solver options should apply to most
users, though complete control over all internal parameters and time
adaptivity algorithms is enabled through optional interface routines.  

Due to its similarities in both function and design with the CVODE
package, a significant portion of this documentation has been directly
adapted from the CVODE documentation [HindmarshSerban2012]_.

ARKode is written in C, with C++ and Fortran interfaces.

ARKode is developed by `Southern Methodist University
<http://www.smu.edu>`_, with support by the `US Department of Energy
<http://www.doe.gov>`_ through the `FASTMath
<http://www.fastmath-scidac.org/>`_ SciDAC Institute, under subcontract
B598130 from `Lawrence Livermore National Laboratory
<http://www.llnl.gov>`_. 


Table of Contents
-------------------

.. toctree::
   :maxdepth: 2

   Introduction
   Mathematics
   Organization
   CInterface
   FortranInterface
   NVectors
   LinearSolvers
   Installation
   Constants
   Examples
   References


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

