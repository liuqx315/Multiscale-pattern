..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _LinearSolvers.custom:

Providing Alternate Linear Solver Modules
==================================================


The central ARKode module interfaces with the linear solver module
using calls to one of four routines. These are denoted here by
:c:func:`linit()`, :c:func:`lsetup()`, :c:func:`lsolve()`, and
:c:func:`lfree()`. Briefly, their purposes are as follows:

* :c:func:`linit()`: initializes and allocate memory specific to the
  linear solver; 
* :c:func:`lsetup()`: evaluates and preprocesses the Jacobian or
  preconditioner; 
* :c:func:`lsolve()`: solves the linear system;
* :c:func:`lfree()`: frees the linear solver memory.

A linear solver module must also provide a user-callable specification
routine (like those described in the section
:ref:`CInterface.LinearSolvers`) which will attach the above four
routines to the main ARKode memory block. The ARKode memory block is a
structure defined in the header file ``arkode_impl.h``. A pointer to
such a structure is defined as the type ``ARKodeMem``. The four
fields in a ``ARKodeMem`` structure that must point to the linear
solver's functions are ``ark_linit``, ``ark_lsetup``, ``ark_lsolve``,
and ``ark_lfree``, respectively. Note that of the four interface
routines, only the :c:func:`lsolve()` routine is required. The
:c:func:`lfree()` routine must be provided only if the solver
specification routine makes any memory allocation. The linear
solver specification function must also set the value of the field
``ark_setupNonNull`` in the ARKode memory block -- to ``TRUE`` if
:c:func:`lsetup()` is used, or ``FALSE`` otherwise. 

For consistency with the existing ARKode linear solver modules, we
recommend that the return value of the specification function be 0 for
a successful return or a negative value if an error occurs (e.g. if
the pointer to the main ARKode memory block is ``NULL``, an input is
illegal, the NVECTOR implementation is not compatible, a memory
allocation fails, etc). 

To facilitate data exchange between the four interface functions, the
field ``ark_lmem`` in the ARKode memory block can be used to attach a
linear solver-specific memory block. That memory should be allocated
in the linear solver specification function. 

These four routines that interface between ARKode and the linear
solver module necessarily have fixed call sequences.  Thus, a user
wishing to implement another linear solver within the ARKode package
must adhere to this set of interfaces. The following is a complete
description of the call list for each of these routines. Note that the
call list of each routine includes a pointer to the main ARKode memory
block, by which the routine can access various data related to the
ARKode solution. The contents of this memory block are given in the
file ``arkode_impl.h`` (but not reproduced here, for the sake of
space).


Initialization function
-----------------------------------

The type definition of :c:func:`linit()` is

.. c:function:: typedef int (*linit)(ARKodeMem ark_mem)

   Completes initializations for the specific linear solver, such as
   counters and statistics. 

   **Arguments:**
      * `ark_mem` -- pointer to the ARKode memory block.
   
   **Return value:**  Should return 0 if it has successfully
   initialized the ARKode linear solver and -1 otherwise.



Setup function
-----------------------------------

   
The type definition of :c:func:`lsetup()` is

.. c:function:: typedef int (*lsetup)(ARKodeMem ark_mem, int convfail, N_Vector ypred, N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)

   Prepares the linear solver for subsequent calls to
   :c:func:`lsolve()`. It may re-compute Jacobian-related data is it
   deems necessary.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `convfail` -- an input flag used to indicate any problem that
	occurred during the solution of the nonlinear equation on the
	current time step for which the linear solver is being
	used. This flag can be used to help decide whether the
	Jacobian data kept by a linear solver needs to be
	updated or not. Its possible values are:

        - ARK_NO_FAILURES: this value is passed if either this is the
	  first call for this step, or the local error test failed on
	  the previous attempt at this step (but the Newton iteration
	  converged).
        - ARK_FAIL_BAD_J: this value is passed if (a) the previous
	  Newton corrector iteration did not converge and the linear
	  solver's setup routine indicated that its Jacobian-related
	  data is not current, or (b) during the previous Newton
	  corrector iteration, the linear solver's solve routine
	  failed in a recoverable manner and the linear solver's setup
	  routine indicated that its Jacobian-related data is not
	  current. 
        - ARK_FAIL_OTHER: this value is passed if during the current
	  internal step try, the previous Newton iteration failed to
	  converge even though the linear solver was using current
	  Jacobian-related data.

      * `ypred` -- is the predicted :math:`y` vector for the current
	ARKode internal step. 
      * `fpred` -- is the value of the implicit right-hand side at
	`ypred`, i.e. :math:`f_I(t_n,ypred)`. 
      * `jcurPtr` -- is a pointer to a boolean to be filled in by
	:c:func:`lsetup()`. The function should set ``*jcurPtr = TRUE``
        if its Jacobian data is current after the call and should set
	``*jcurPtr = FALSE`` if its Jacobian data is not current. If
	:c:func:`lsetup()` calls for re-evaluation of Jacobian data
	(based on `convfail` and ARKode state data), it should return
	``*jcurPtr = TRUE`` unconditionally; otherwise an infinite
	loop can result.
      * `vtemp1`, `vtemp2`, `vtemp3` -- are temporary variables of
	type ``N_Vector`` provided for use by :c:func:`lsetup()`. 
   
   **Return value:** 
   Should return 0 if successful, a positive value
   for a recoverable error, and a negative value for an unrecoverable
   error.





Solve function
-----------------------------------

The type definition of :c:func:`lsolve()` is

.. c:function:: typedef int (*lsolve)(ARKodeMem ark_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector fcur)

   Solves the linear equation :math:`A x = b`, where :math:`A` arises
   in the Newton iteration :eq:`Newton_system` and gives
   some approximation to :math:`M - \gamma J`, :math:`J = \frac{\partial
   f_I}{\partial y}(t_n, ycur)`.  Note, the right-hand side vector
   :math:`b` is input, and :math:`\gamma` is available as
   ``ark_mem->ark_gamma``. 

   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `b` -- is the right-hand side vector :math:`b`. The solution
	is also to be returned in the vector :math:`b`. 
      * `weight` -- is a vector that contains the error weights. These
	are the :math:`w_i` of :ref:`CInterface.ErrorWeight`.
      * `ycur` -- is a vector that contains the solver's current
	approximation to :math:`y(t_n)`. 
      * `fcur` -- is a vector that contains :math:`f_I(t_n, ycur)`.

   **Return value:**  Should return 0 if successful, a positive value
   for a recoverable error, and a negative value for an unrecoverable
   error. 



Memory deallocation function
-----------------------------------

The type definition of :c:func:`lfree()` is

.. c:function:: typedef void (*lfree)(ARKodeMem ark_mem)

   free up any memory allocated by the linear solver.

   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.

   **Return value:**  None

   **Notes:**  This routine is called once a problem has been
   completed and the linear solver is no longer needed.
