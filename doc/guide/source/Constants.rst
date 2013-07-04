:tocdepth: 3

.. _Constants:

================
ARKode Constants
================

Below we list all input and output constants used by the main solver
and linear solver modules, together with their numerical values and a
short description of their meaning. 


ARKode input constants
==========================

ARKode main solver module
^^^^^^^^^^^^^^^^^^^^^^^^^^

  :index:`ARK_NORMAL` (1): 
     Solver returns at a specified output time.

  :index:`ARK_ONE_STEP`  (2): 
     Solver returns after each successful step.


Iterative linear solver module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  :index:`PREC_NONE`  (0): 
     No preconditioning.

  :index:`PREC_LEFT`  (1): 
     Preconditioning on the left only.

  :index:`PREC_RIGHT`  (2): 
     Preconditioning on the right only.

  :index:`PREC_BOTH`  (3): 
     Preconditioning on both the left and the right.

  :index:`MODIFIED_GS`  (1): 
     Use modified Gram-Schmidt procedure.

  :index:`CLASSICAL_GS`  (2): 
     Use classical Gram-Schmidt procedure.




ARKode output constants
==========================

ARKode main solver module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  :index:`ARK_SUCCESS`  (0): 
     Successful function return.

  :index:`ARK_TSTOP_RETURN`  (1): 
     ARKode succeeded by reachign the specified
     stopping point.

  :index:`ARK_ROOT_RETURN`  (2): 
     ARKode succeeded and found one more more roots.

  :index:`ARK_WARNING`  (99): 
     ARKode succeeded but an unusual situation occurred.

  :index:`ARK_TOO_MUCH_WORK`  (-1): 
     The solver took ``mxstep`` internal steps
     but could not reach ``tout``.

  :index:`ARK_TOO_MUCH_ACC`  (-2): 
     The solver could not satisfy the accuracy
     demanded by the user for some internal step.

  :index:`ARK_ERR_FAILURE`  (-3): 
     Error test failures occurred too many times
     during one internal time step, or the minimum step size was
     reached. 

  :index:`ARK_CONV_FAILURE`  (-4): 
     Convergence test failures occurred too many
     times during one internal time step, or the minimum step size was
     reached. 

  :index:`ARK_LINIT_FAIL`  (-5): 
     The linear solver's initialization function failed.

  :index:`ARK_LSETUP_FAIL`  (-6): 
     The linear solver's setup function failed in
     an unrecoverable manner.

  :index:`ARK_LSOLVE_FAIL`  (-7): 
     The linear solver's solve function failed in 
     an unrecoverable manner.

  :index:`ARK_RHSFUNC_FAIL`  (-8): 
     The right-hand side function failed in an
     unrecoverable manner.

  :index:`ARK_FIRST_RHSFUNC_ERR`  (-9): 
     The right-hand side function failed 
     at the first call.

  :index:`ARK_REPTD_RHSFUNC_ERR`  (-10): 
     The right-hand side function had 
     repeated recoverable errors.

  :index:`ARK_UNREC_RHSFUNC_ERR`  (-11): 
     The right-hand side function had a
     recoverable error, but no recovery is possible.

  :index:`ARK_RTFUNC_FAIL`  (-12): 
     The rootfinding function failed in an
     unrecoverable manner.

  :index:`ARK_MEM_FAIL`  (-20): 
     A memory allocation failed.

  :index:`ARK_MEM_NULL`  (-21): 
     The ``arkode_mem`` argument was ``NULL``.

  :index:`ARK_ILL_INPUT`  (-22): 
     One of the function inputs is illegal.

  :index:`ARK_NO_MALLOC`  (-23): 
     The ARKode memory block was not allocated by 
     a call to :c:func:`ARKodeMalloc()`.

  :index:`ARK_BAD_K`  (-24): 
     The derivative order :math:`k` is larger than allowed.

  :index:`ARK_BAD_T`  (-25): 
     The time :math:`t` is outside the last step taken.

  :index:`ARK_BAD_DKY`  (-26): 
     The output derivative vector is ``NULL``.

  :index:`ARK_TOO_CLOSE`  (-27): 
     The output and initial times are too close to 
     each other.


ARKDLS linear solver modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  :index:`ARKDLS_SUCCESS`  (0): 
     Successful function return.

  :index:`ARKDLS_MEM_NULL`  (-1): 
     The ``arkode_mem`` argument was ``NULL``.

  :index:`ARKDLS_LMEM_NULL`  (-2): 
     The ARKDLS linear solver has not been initialized.

  :index:`AKRDLS_ILL_INPUT`  (-3): 
     The ARKDLS solver is not compatible with
     the current NVECTOR module.

  :index:`ARKDLS_MEM_FAIL`  (-4): 
     A memory allocation request failed.

  :index:`ARKDLS_JACFUNC_UNRECVR`  (-5): 
     The Jacobian function failed in an
     unrecoverable manner.

  :index:`ARKDLS_JACFUNC_RECVR`  (-6): 
     The Jacobian function had a recoverable error.



ARKSPILS linear solver modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  :index:`ARKSPILS_SUCCESS`  (0): 
     Successful function return.

  :index:`ARKSPILS_MEM_NULL`  (-1): 
     The ``arkode_mem`` argument was ``NULL``.

  :index:`ARKSPILS_LMEM_NULL`  (-2): 
     The ARKSPILS linear solver has not been initialized.

  :index:`AKRSPILS_ILL_INPUT`  (-3): 
     The ARKSPILS solver is not compatible with
     the current NVECTOR module, or an input value was illegal.

  :index:`ARKSPILS_MEM_FAIL`  (-4): 
     A memory allocation request failed.

  :index:`ARKSPILS_PMEM_FAIL`  (-5): 
     The preconditioner module has not been initialized.



ARKSPGMR generic linear solver module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


  :index:`SPGMR_SUCCESS`  (0): 
     Converged.

  :index:`SPGMR_RES_REDUCED`  (1): 
     No convergence, but the residual norm was
     reduced. 

  :index:`SPGMR_CONV_FAIL`  (2): 
     Failure to converge.

  :index:`SPGMR_QRFACT_FAIL`  (3): 
     A singular matrix was found during the 
     QR factorization.

  :index:`SPGMR_PSOLVE_FAIL_REC`  (4): 
     The preconditioner solve function 
     failed recoverably.

  :index:`SPGMR_ATIMES_FAIL_REC`  (5): 
     The Jacobian-times-vector function 
     failed recoverably.

  :index:`SPGMR_PSET_FAIL_REC`  (6): 
     The preconditioner setup function failed 
     recoverably.

  :index:`SPGMR_MEM_NULL`  (-1): 
     The SPGMR memory is ``NULL``

  :index:`SPGMR_ATIMES_FAIL_UNREC`  (-2): 
     The Jacobian-times-vector function
     failed unrecoverably.

  :index:`SPGMR_PSOLVE_FAIL_UNREC`  (-3): 
     The preconditioner solve function 
     failed unrecoverably.

  :index:`SPGMR_GS_FAIL`  (-4): 
     Failure in the Gram-Schmidt procedure.

  :index:`SPGMR_QRSOL_FAIL`  (-5): 
     The matrix :MATH:`R` was found to be
     singular during the QR solve phase.

  :index:`SPGMR_PSET_FAIL_UNREC`  (-6): 
     The preconditioner setup function 
     failed unrecoverably.



ARKSPBCG generic linear solver module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  :index:`SPBCG_SUCCESS`  (0): 
     Converged.

  :index:`SPBCG_RES_REDUCED`  (1): 
     No convergence, but the residual norm 
     was reduced.

  :index:`SPBCG_CONV_FAIL`  (2): 
     Failure to converge.

  :index:`SPBCG_PSOLVE_FAIL_REC`  (3): 
     The preconditioner solve function 
     failed recoverably.

  :index:`SPBCG_ATIMES_FAIL_REC`  (4): 
     The Jacobian-times-vector function 
     failed recoverably.

  :index:`SPBCG_PSET_FAIL_REC`  (5): 
     The preconditioner setup function 
     failed recoverably.

  :index:`SPBCG_MEM_NULL`  (-1): 
     The SPBCG memory is ``NULL``

  :index:`SPBCG_ATIMES_FAIL_UNREC`  (-2): 
     The Jacobian-times-vector function 
     failed unrecoverably.

  :index:`SPBCG_PSOLVE_FAIL_UNREC`  (-3): 
     The preconditioner solve function 
     failed unrecoverably.

  :index:`SPBCG_PSET_FAIL_UNREC`  (-4): 
     The preconditioner setup function 
     failed unrecoverably.



ARKSPTFQMR generic linear solver module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  :index:`SPTFQMR_SUCCESS`  (0): 
     Converged.

  :index:`SPTFQMR_RES_REDUCED`  (1): 
     No convergence, but the residual norm 
     was reduced.

  :index:`SPTFQMR_CONV_FAIL`  (2): 
     Failure to converge.

  :index:`SPTFQMR_PSOLVE_FAIL_REC`  (3): 
     The preconditioner solve function 
     failed recoverably.

  :index:`SPTFQMR_ATIMES_FAIL_REC`  (4): 
     The Jacobian-times-vector function 
     failed recoverably.

  :index:`SPTFQMR_PSET_FAIL_REC`  (5): 
     The preconditioner setup function 
     failed recoverably.

  :index:`SPTFQMR_MEM_NULL`  (-1): 
     The SPTFQMR memory is ``NULL``

  :index:`SPTFQMR_ATIMES_FAIL_UNREC`  (-2): 
     The Jacobian-times-vector 
     function failed.

  :index:`SPTFQMR_PSOLVE_FAIL_UNREC`  (-3): 
     The preconditioner solve function 
     failed unrecoverably.

  :index:`SPTFQMR_PSET_FAIL_UNREC`  (-4): 
     The preconditioner setup function 
     failed unrecoverably.


ARKSPFGMR generic linear solver module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


  :index:`SPFGMR_SUCCESS`  (0): 
     Converged.

  :index:`SPFGMR_RES_REDUCED`  (1): 
     No convergence, but the residual norm was
     reduced. 

  :index:`SPFGMR_CONV_FAIL`  (2): 
     Failure to converge.

  :index:`SPFGMR_QRFACT_FAIL`  (3): 
     A singular matrix was found during the 
     QR factorization.

  :index:`SPFGMR_PSOLVE_FAIL_REC`  (4): 
     The preconditioner solve function 
     failed recoverably.

  :index:`SPFGMR_ATIMES_FAIL_REC`  (5): 
     The Jacobian-times-vector function 
     failed recoverably.

  :index:`SPFGMR_PSET_FAIL_REC`  (6): 
     The preconditioner setup function failed 
     recoverably.

  :index:`SPFGMR_MEM_NULL`  (-1): 
     The SPFGMR memory is ``NULL``

  :index:`SPFGMR_ATIMES_FAIL_UNREC`  (-2): 
     The Jacobian-times-vector function
     failed unrecoverably.

  :index:`SPFGMR_PSOLVE_FAIL_UNREC`  (-3): 
     The preconditioner solve function 
     failed unrecoverably.

  :index:`SPFGMR_GS_FAIL`  (-4): 
     Failure in the Gram-Schmidt procedure.

  :index:`SPFGMR_QRSOL_FAIL`  (-5): 
     The matrix :MATH:`R` was found to be
     singular during the QR solve phase.

  :index:`SPFGMR_PSET_FAIL_UNREC`  (-6): 
     The preconditioner setup function 
     failed unrecoverably.



ARKPCG generic linear solver module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  :index:`PCG_SUCCESS`  (0): 
     Converged.

  :index:`PCG_RES_REDUCED`  (1): 
     No convergence, but the residual norm 
     was reduced.

  :index:`PCG_CONV_FAIL`  (2): 
     Failure to converge.

  :index:`PCG_PSOLVE_FAIL_REC`  (3): 
     The preconditioner solve function 
     failed recoverably.

  :index:`PCG_ATIMES_FAIL_REC`  (4): 
     The Jacobian-times-vector function 
     failed recoverably.

  :index:`PCG_PSET_FAIL_REC`  (5): 
     The preconditioner setup function 
     failed recoverably.

  :index:`PCG_MEM_NULL`  (-1): 
     The PCG memory is ``NULL``

  :index:`PCG_ATIMES_FAIL_UNREC`  (-2): 
     The Jacobian-times-vector function 
     failed unrecoverably.

  :index:`PCG_PSOLVE_FAIL_UNREC`  (-3): 
     The preconditioner solve function 
     failed unrecoverably.

  :index:`PCG_PSET_FAIL_UNREC`  (-4): 
     The preconditioner setup function 
     failed unrecoverably.


