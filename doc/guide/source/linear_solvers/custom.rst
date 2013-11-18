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


Newton system linear solver
------------------------------

The central ARKode module interfaces with the Newton system linear
solver module using calls to one of four routines. These are denoted
here by :c:func:`linit()`, :c:func:`lsetup()`, :c:func:`lsolve()`, and
:c:func:`lfree()`. Briefly, their purposes are as follows:

* :c:func:`linit()`: initializes and allocates memory specific to the
  linear solver; 
* :c:func:`lsetup()`: evaluates and preprocesses the Jacobian or
  preconditioner in preparation for solves; 
* :c:func:`lsolve()`: solves the linear system;
* :c:func:`lfree()`: frees the linear solver memory.

A linear solver module must also provide a user-callable specification
routine (like those described in the section
:ref:`CInterface.LinearSolvers`) which will attach the above four
routines to the main ARKode memory block. The ARKode memory block is a
structure defined in the header file ``arkode_impl.h``. A pointer to
such a structure is defined as the type ``ARKodeMem``. The four
fields in the ``ARKodeMem`` structure that refer to the Newton system
linear solver's functions are ``ark_linit``, ``ark_lsetup``,
``ark_lsolve``, and ``ark_lfree``, respectively.  Note that of these
interface routines, only the :c:func:`lsolve()` routine is
required. The :c:func:`lfree()` routine must be provided only if the
solver specification routine makes any memory allocation.  The linear
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




Mass matrix linear solver
------------------------------

Similarly, for problems involving a non-identity mass matrix
:math:`M\ne I`, the main ARKode module interfaces with the mass matrix
linear solver module using calls to one of four routines:
:c:func:`minit()`, :c:func:`msetup()`, :c:func:`msolve()`, and
:c:func:`mfree()`. Briefly, their purposes are as follows: 

* :c:func:`minit()`: initializes and allocates memory specific to the
  mass matrix linear solver; 
* :c:func:`msetup()`: evaluates and preprocesses the mass matrix or
  associated preconditioner in preparation for solves; 
* :c:func:`msolve()`: solves the mass matrix system;
* :c:func:`mfree()`: frees the mass matrix linear solver memory.

As with the Newton system linear solver, a mass matrix linear solver
module must also provide a user-callable specification routine (like
those described in the section :ref:`CInterface.LinearSolvers`) which
will attach the above four routines to the main ARKode memory
block.  The four fields in the ``ARKodeMem`` structure that refer to
the mass matrix system linear solver's functions are ``ark_minit``,
``ark_msetup``, ``ark_msolve``, and ``ark_mfree``, respectively.  As
with the Newton system solver, only :c:func:`msolve()` is required,
and :c:func:`mfree()` must be provided only if the solver
specification routine makes any memory allocation.  The mass matrix
linear solver specification function must also set the value of the
field ``ark_MassSetupNonNull`` in the ARKode memory block -- to
``TRUE`` if :c:func:`msetup()` is used, or ``FALSE`` otherwise. 

For consistency with the existing ARKode linear solver modules, we
recommend that the return value of the specification function be 0 for
a successful return or a negative value if an error occurs (e.g. if
the pointer to the main ARKode memory block is ``NULL``, an input is
illegal, the NVECTOR implementation is not compatible, a memory
allocation fails, etc). 

To facilitate data exchange between the four interface functions, the
field ``ark_mass_mem`` in the ARKode memory block can be used to
attach a linear solver-specific memory block.  That memory should be
allocated in the linear solver specification function. 



These above routines that interface between ARKode and the Newton
system or mass matrix linear solver module necessarily have fixed call
sequences.  Thus, a user wishing to implement another linear solver
within the ARKode package must adhere to this set of interfaces.  The
following is a complete description of the call list for each of these
routines.  Note that the call list of each routine includes a pointer
to the main ARKode memory block, by which the routine can access
various data related to the ARKode solution. The contents of this
memory block are given in the file ``arkode_impl.h`` (but not
reproduced here, for the sake of space).





Initialization function
-----------------------------------

The type definition of :c:func:`linit()` is

.. c:function:: typedef int (*linit)(ARKodeMem ark_mem)

   Completes initializations for the specific linear solver, such as
   counters and statistics. 

   **Arguments:**
      * *ark_mem* -- pointer to the ARKode memory block.
   
   **Return value:**  Should return 0 if it has successfully
   initialized the ARKode linear solver and -1 otherwise.


Similarly, the type definition of :c:func:`minit()` is

.. c:function:: typedef int (*minit)(ARKodeMem ark_mem)

   Completes initializations for the specific linear solver, such as
   counters and statistics. 

   **Arguments:**
      * *ark_mem* -- pointer to the ARKode memory block.
   
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
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *convfail* -- an input flag used to indicate any problem that
	occurred during the solution of the nonlinear equation on the
	current time step for which the linear solver is being
	used. This flag can be used to help decide whether the
	Jacobian data kept by a linear solver needs to be
	updated or not. Its possible values are:

        - *ARK_NO_FAILURES*: this value is passed if either this is the
	  first call for this step, or the local error test failed on
	  the previous attempt at this step (but the Newton iteration
	  converged).
        - *ARK_FAIL_BAD_J*: this value is passed if (a) the previous
	  Newton corrector iteration did not converge and the linear
	  solver's setup routine indicated that its Jacobian-related
	  data is not current, or (b) during the previous Newton
	  corrector iteration, the linear solver's solve routine
	  failed in a recoverable manner and the linear solver's setup
	  routine indicated that its Jacobian-related data is not
	  current. 
        - *ARK_FAIL_OTHER*: this value is passed if during the current
	  internal step try, the previous Newton iteration failed to
	  converge even though the linear solver was using current
	  Jacobian-related data.

      * *ypred* -- is the predicted :math:`y` vector for the current
	ARKode internal step. 
      * *fpred* -- is the value of the implicit right-hand side at
	*ypred*, i.e. :math:`f_I(t_n,ypred)`. 
      * *jcurPtr* -- is a pointer to a boolean to be filled in by
	:c:func:`lsetup()`. The function should set ``*jcurPtr = TRUE``
        if its Jacobian data is current after the call and should set
	``*jcurPtr = FALSE`` if its Jacobian data is not current. If
	:c:func:`lsetup()` calls for re-evaluation of Jacobian data
	(based on *convfail* and ARKode state data), it should return
	``*jcurPtr = TRUE`` unconditionally; otherwise an infinite
	loop can result.
      * *vtemp1*, *vtemp2*, *vtemp3* -- are temporary variables of
	type ``N_Vector`` provided for use by :c:func:`lsetup()`. 
   
   **Return value:** 
   Should return 0 if successful, a positive value
   for a recoverable error, and a negative value for an unrecoverable
   error.


Similarly, the type definition of :c:func:`msetup()` is

.. c:function:: typedef int (*msetup)(ARKodeMem ark_mem, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)

   Prepares the mass matrix linear solver for subsequent calls to
   :c:func:`msolve()`. It may re-compute mass-matrix-related data is
   it deems necessary.
   
   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *vtemp1*, *vtemp2*, *vtemp3* -- are temporary variables of
	type ``N_Vector`` provided for use by :c:func:`msetup()`. 
   
   **Return value:** 
   Should return 0 if successful, a positive value
   for a recoverable error, and a negative value for an unrecoverable
   error.





Solve function
-----------------------------------

The type definition of :c:func:`lsolve()` is

.. c:function:: typedef int (*lsolve)(ARKodeMem ark_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector fcur)

   Solves the linear equation :math:`{\mathcal A} x = b`, where
   :math:`{\mathcal A}` arises  in the Newton iteration (see the
   section :ref:`Mathematics.Linear`) and gives some approximation to
   :math:`M - \gamma J`, :math:`J = \frac{\partial}{\partial y} f_I(t_n, ycur)`.  
   Note, the right-hand side vector 
   :math:`b` is input, and :math:`\gamma` is available as
   ``ark_mem->ark_gamma``. 

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *b* -- is the right-hand side vector :math:`b`. The solution
	is also to be returned in the vector :math:`b`. 
      * *weight* -- is a vector that contains the residual weights. These
	are the :math:`rwt_i` of :ref:`CInterface.ResidualWeight`.
      * *ycur* -- is a vector that contains the solver's current
	approximation to :math:`y(t_n)`. 
      * *fcur* -- is a vector that contains :math:`f_I(t_n, ycur)`.

   **Return value:**  Should return 0 if successful, a positive value
   for a recoverable error, and a negative value for an unrecoverable
   error. 


Similarly, the type definition of :c:func:`msolve()` is

.. c:function:: typedef int (*msolve)(ARKodeMem ark_mem, N_Vector b, N_Vector weight)

   Solves the linear equation :math:`M x = b`, where :math:`M` is the
   system mass matrix.  Note, the right-hand side vector :math:`b` is
   input, and holds the solution :math:`x` on output.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.
      * *b* -- is the right-hand side vector :math:`b`. The solution
	is also to be returned in the vector :math:`b`. 
      * *weight* -- is a vector that contains the error weights. These
	are the :math:`rwt_i` of :ref:`CInterface.ResidualWeight`.

   **Return value:**  Should return 0 if successful, a positive value
   for a recoverable error, and a negative value for an unrecoverable
   error. 



Memory deallocation function
-----------------------------------

The type definition of :c:func:`lfree()` is

.. c:function:: typedef void (*lfree)(ARKodeMem ark_mem)

   free up any memory allocated by the linear solver.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.

   **Return value:**  None

   **Notes:**  This routine is called once a problem has been
   completed and the linear solver is no longer needed.


Similarly, the type definition of :c:func:`mfree()` is

.. c:function:: typedef void (*mfree)(ARKodeMem ark_mem)

   free up any memory allocated by the mass matrix linear solver.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKode memory block.

   **Return value:**  None

   **Notes:**  This routine is called once a problem has been
   completed and the mass matrix solver is no longer needed.
