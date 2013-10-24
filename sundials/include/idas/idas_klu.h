/*
 * -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2013, Lawrence Livermore National Security
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the IDAKLU linear solver module.
 * -----------------------------------------------------------------
 */

#ifndef _IDASKLU_H
#define _IDASKLU_H

#include <idas/idas_sparse.h>
#include <sundials/sundials_sparse.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * Function : IDAKLU
 * -----------------------------------------------------------------
 * A call to the IDAKLU function links the main integrator      
 * with the IDAKLU linear solver module.                        
 *                                                                
 * ida_mem is the pointer to integrator memory returned by        
 *     IDACreate.             
 *
 *                                                                
 * IDAKLU returns:                                              
 *     IDASLU_SUCCESS   = 0  if successful                              
 *     IDASLU_LMEM_FAIL = -1 if there was a memory allocation failure   
 *     IDASLU_ILL_INPUT = -2 if NVECTOR found incompatible           
 *                                                                
 * NOTE: The KLU linear solver assumes a serial implementation  
 *       of the NVECTOR package. Therefore, IDAKLU will first
 *       test for a compatible N_Vector internal representation
 *       by checking that the functions N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int IDAKLU(void *ida_mem, int n, int nnz); 

/*
 * -----------------------------------------------------------------
 * Function: IDAKLUB
 * -----------------------------------------------------------------
 * IDAKLUB links the main IDAS integrator with the IDAKLU
 * linear solver for the backward integration.
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int IDAKLUB(void *ida_mem, int which, int nB, int nnzB);



  
#ifdef __cplusplus
}
#endif

#endif
