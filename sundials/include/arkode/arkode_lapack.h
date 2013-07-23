/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
----------------------------------------------------------------- 
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Header file for the ARKODE dense linear solver ARKLAPACK.
---------------------------------------------------------------*/

#ifndef _ARKLAPACK_H
#define _ARKLAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <arkode/arkode_direct.h>
#include <sundials/sundials_lapack.h>


/*===============================================================
  EXPORTED FUNCTIONS
===============================================================*/

/*---------------------------------------------------------------
 ARKLapackDense:

 A call to the ARKLapackDense function links the main integrator
 with the ARKLAPACK linear solver using dense Jacobians.

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 N is the size of the ODE system.

 The return value of ARKLapackDense is one of:
    ARKLAPACK_SUCCESS   if successful
    ARKLAPACK_MEM_NULL  if the ARKODE memory was NULL
    ARKLAPACK_MEM_FAIL  if there was a memory allocation failure
    ARKLAPACK_ILL_INPUT if a required vector operation is missing
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKLapackDense(void *arkode_mem, int N);


/*---------------------------------------------------------------
 ARKLapackBand:

 A call to the ARKLapackBand function links the main integrator
 with the ARKLAPACK linear solver using banded Jacobians. 

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 N is the size of the ODE system.

 mupper is the upper bandwidth of the band Jacobian approximation.

 mlower is the lower bandwidth of the band Jacobian approximation.

 The return value of ARKLapackBand is one of:
    ARKLAPACK_SUCCESS   if successful
    ARKLAPACK_MEM_NULL  if the ARKODE memory was NULL
    ARKLAPACK_MEM_FAIL  if there was a memory allocation failure
    ARKLAPACK_ILL_INPUT if a required vector operation is missing 
                        or if a bandwidth has an illegal value.
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKLapackBand(void *arkode_mem, int N, int mupper, int mlower);

#ifdef __cplusplus
}
#endif

#endif
