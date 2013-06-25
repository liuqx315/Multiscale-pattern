/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
----------------------------------------------------------------- 
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 This is the header file for the ARKODE scaled preconditioned 
 CG linear solver, ARKPCG.
---------------------------------------------------------------*/

#ifndef _ARKPCG_H
#define _ARKPCG_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <arkode/arkode_spils.h>
#include <sundials/sundials_pcg.h>

/*---------------------------------------------------------------
 ARKPcg:

 A call to the ARKPcg function links the main ARKODE integrator
 with the ARKPCG linear solver.

 arkode_mem is the pointer to the integrator memory returned by
           ARKodeCreate.

 pretype   is the type of user preconditioning to be done.
           This must be one of the four enumeration constants
           PREC_NONE, PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined
           in iterative.h. These correspond to no preconditioning,
           left preconditioning only, right preconditioning
           only, and both left and right preconditioning,
           respectively.  However, since PCG requires a symmetric 
	   linear operator, this flag is checked for any one of 
	   PREC_LEFT, PREC_RIGHT and PREC_BOTH -- if any are 
	   found then preconditioning is used.  It is assumed 
	   that the preconditioner implements a symmetric linear 
	   operator.

 maxl      is the maximum Krylov dimension. This is an
           optional input to the ARKPCG solver. Pass 0 to
           use the default value ARKPCG_MAXL=5.

 The return value of ARKPcg is one of:
    ARKSPILS_SUCCESS   if successful
    ARKSPILS_MEM_NULL  if the arkode memory was NULL
    ARKSPILS_MEM_FAIL  if there was a memory allocation failure
    ARKSPILS_ILL_INPUT if a required vector operation is missing
 The above constants are defined in arkode_spils.h

---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKPcg(void *arkode_mem, int pretype, int maxl);


#ifdef __cplusplus
}
#endif

#endif
