/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * This is the header file for the ARKODE diagonal linear solver, ARKDIAG.
 * -----------------------------------------------------------------
 ***** UNTOUCHED *****
 * -----------------------------------------------------------------
 */

#ifndef _ARKDIAG_H
#define _ARKDIAG_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_nvector.h>

/*
 * -----------------------------------------------------------------
 * Function : ARKDiag
 * -----------------------------------------------------------------
 * A call to the ARKDiag function links the main integrator with
 * the ARKDIAG linear solver.
 *
 * arkode_mem is the pointer to the integrator memory returned by
 *           ARKodeCreate.
 *
 * The return value of ARKDiag is one of:
 *    ARKDIAG_SUCCESS   if successful
 *    ARKDIAG_MEM_NULL  if the arkode memory was NULL
 *    ARKDIAG_MEM_FAIL  if there was a memory allocation failure
 *    ARKDIAG_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int ARKDiag(void *arkode_mem);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the ARKDIAG linear solver
 * -----------------------------------------------------------------
 *
 * ARKDiagGetWorkSpace returns the real and integer workspace used
 *                    by ARKDIAG.
 * ARKDiagGetNumRhsEvals returns the number of calls to the user
 *                      f routine due to finite difference Jacobian
 *                      evaluation.
 *                      Note: The number of diagonal approximate
 *                      Jacobians formed is equal to the number of
 *                      ARKDiagSetup calls. This number is available
 *                      through ARKodeGetNumLinSolvSetups.
 * ARKDiagGetLastFlag returns the last error flag set by any of
 *                   the ARKDIAG interface functions.
 *
 * The return value of ARKDiagGet* is one of:
 *    ARKDIAG_SUCCESS   if successful
 *    ARKDIAG_MEM_NULL  if the arkode memory was NULL
 *    ARKDIAG_LMEM_NULL if the arkdiag memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int ARKDiagGetWorkSpace(void *arkode_mem, long int *lenrwLS, long int *leniwLS);
SUNDIALS_EXPORT int ARKDiagGetNumRhsEvals(void *arkode_mem, long int *nfevalsLS);
SUNDIALS_EXPORT int ARKDiagGetLastFlag(void *arkode_mem, long int *flag);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a ARKDIAG return flag
 * -----------------------------------------------------------------
 */
  
SUNDIALS_EXPORT char *ARKDiagGetReturnFlagName(long int flag);

/*
 * -----------------------------------------------------------------
 * ARKDIAG return values
 * -----------------------------------------------------------------
 */

#define ARKDIAG_SUCCESS          0
#define ARKDIAG_MEM_NULL        -1
#define ARKDIAG_LMEM_NULL       -2
#define ARKDIAG_ILL_INPUT       -3
#define ARKDIAG_MEM_FAIL        -4

/* Additional last_flag values */

#define ARKDIAG_INV_FAIL        -5
#define ARKDIAG_RHSFUNC_UNRECVR -6
#define ARKDIAG_RHSFUNC_RECVR   -7

#ifdef __cplusplus
}
#endif

#endif
