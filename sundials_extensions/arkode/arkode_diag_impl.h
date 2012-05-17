/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * Implementation header file for the diagonal linear solver, ARKDIAG.
 * -----------------------------------------------------------------
 ***** UNTOUCHED *****
 * -----------------------------------------------------------------
 */

#ifndef _ARKDIAG_IMPL_H
#define _ARKDIAG_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <arkode/arkode_diag.h>

/*
 * -----------------------------------------------------------------
 * Types: ARKDiagMemRec, ARKDiagMem
 * -----------------------------------------------------------------
 * The type ARKDiagMem is pointer to a ARKDiagMemRec.
 * This structure contains ARKDiag solver-specific data.
 * -----------------------------------------------------------------
 */

typedef struct {

  realtype di_gammasv; /* gammasv = gamma at the last call to setup */
                       /* or solve                                  */

  N_Vector di_M;       /* M = (I - gamma J)^{-1} , gamma = h / l1   */

  N_Vector di_bit;     /* temporary storage vector                  */

  N_Vector di_bitcomp; /* temporary storage vector                  */

  long int di_nfeDI;   /* no. of calls to f due to difference 
                          quotient diagonal Jacobian approximation  */

  long int di_last_flag;    /* last error return flag                    */

} ARKDiagMemRec, *ARKDiagMem;

/* Error Messages */

#define MSGDG_ARKMEM_NULL "Integrator memory is NULL."
#define MSGDG_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGDG_MEM_FAIL "A memory request failed."
#define MSGDG_LMEM_NULL "ARKDIAG memory is NULL."
#define MSGDG_RHSFUNC_FAILED "The right-hand side routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
