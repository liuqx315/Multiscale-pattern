/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * This is the header file for the ARKODE dense linear solver, ARKDENSE.
 * -----------------------------------------------------------------
 ***** UNTOUCHED *****
 * -----------------------------------------------------------------
 */

#ifndef _ARKDENSE_H
#define _ARKDENSE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <arkode/arkode_direct.h>
#include <sundials/sundials_dense.h>

/*
 * -----------------------------------------------------------------
 * Function: ARKDense
 * -----------------------------------------------------------------
 * A call to the ARKDense function links the main integrator with
 * the ARKDENSE linear solver.
 *
 * arkode_mem is the pointer to the integrator memory returned by
 *           ARKodeCreate.
 *
 * N is the size of the ODE system.
 *
 * The return value of ARKDense is one of:
 *    ARKDLS_SUCCESS   if successful
 *    ARKDLS_MEM_NULL  if the arkode memory was NULL
 *    ARKDLS_MEM_FAIL  if there was a memory allocation failure
 *    ARKDLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int ARKDense(void *arkode_mem, long int N);

#ifdef __cplusplus
}
#endif

#endif
