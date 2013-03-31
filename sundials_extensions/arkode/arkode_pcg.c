/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
----------------------------------------------------------------- 
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 This is the implementation file for the ARKPCG linear solver.
---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_spils_impl.h"
#include "arkode_impl.h"
#include <arkode/pcg/sundials_pcg.h>
#include <sundials/sundials_math.h>

#define NO_DEBUG_OUTPUT
#ifdef DEBUG_OUTPUT
#include <nvector/nvector_serial.h>
#endif


/* Constants */
#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* ARKPCG linit, lsetup, lsolve, and lfree routines */
static int ARKPcgInit(ARKodeMem ark_mem);
static int ARKPcgSetup(ARKodeMem ark_mem, int convfail, 
		       N_Vector ypred, N_Vector fpred, 
		       booleantype *jcurPtr, N_Vector vtemp1,
		       N_Vector vtemp2, N_Vector vtemp3);
static int ARKPcgSolve(ARKodeMem ark_mem, N_Vector b, 
		       N_Vector weight, N_Vector ynow, 
		       N_Vector fnow);
static void ARKPcgFree(ARKodeMem ark_mem);


/*---------------------------------------------------------------
 ARKPcg:

 This routine initializes the memory record and sets various 
 function fields specific to the Pcg linear solver module. 
 ARKPcg first calls the existing lfree routine if this is not 
 NULL. It then sets the ark_linit, ark_lsetup, ark_lsolve, 
 ark_lfree fields in (*arkode_mem) to be ARKPcgInit, 
 ARKPcgSetup, ARKPcgSolve, and ARKPcgFree, respectively. 
 It allocates memory for a structure of type ARKSpilsMemRec and 
 sets the ark_lmem field in (*arkode_mem) to the address of 
 this structure. It sets setupNonNull in (*arkode_mem),
 and sets various fields in the ARKSpilsMemRec structure.
 Finally, ARKPcg allocates memory for ytemp and x, and calls
 PcgMalloc to allocate memory for the Pcg solver.
---------------------------------------------------------------*/
int ARKPcg(void *arkode_mem, int pretype, int maxl)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  PcgMem pcg_mem;
  int mxl;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKPCG", 
		    "ARKPcg", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check if N_VDotProd and N_VWrmsNorm is present */
  if ((ark_mem->ark_tempv->ops->nvdotprod == NULL) ||
      (ark_mem->ark_tempv->ops->nvwrmsnorm == NULL)) {
    ARKProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKPCG", 
		    "ARKPcg", MSGS_BAD_NVECTOR);
    return(ARKSPILS_ILL_INPUT);
  }

  if (ark_mem->ark_lfree != NULL) ark_mem->ark_lfree(ark_mem);

  /* Set four main function fields in ark_mem */
  ark_mem->ark_linit  = ARKPcgInit;
  ark_mem->ark_lsetup = ARKPcgSetup;
  ark_mem->ark_lsolve = ARKPcgSolve;
  ark_mem->ark_lfree  = ARKPcgFree;

  /* Get memory for ARKSpilsMemRec */
  arkspils_mem = NULL;
  arkspils_mem = (ARKSpilsMem) malloc(sizeof(struct ARKSpilsMemRec));
  if (arkspils_mem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKPCG", 
		    "ARKPcg", MSGS_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* Set ILS type */
  arkspils_mem->s_type = SPILS_PCG;

  /* Set Pcg parameters that have been passed in call sequence */
  arkspils_mem->s_pretype = pretype;
  mxl = arkspils_mem->s_maxl = (maxl <= 0) ? ARKSPILS_MAXL : maxl;

  /* Set defaults for Jacobian-related fileds */
  arkspils_mem->s_jtimesDQ = TRUE;
  arkspils_mem->s_jtimes   = NULL;
  arkspils_mem->s_j_data   = NULL;

  /* Set defaults for preconditioner-related fields */
  arkspils_mem->s_pset   = NULL;
  arkspils_mem->s_psolve = NULL;
  arkspils_mem->s_pfree  = NULL;
  arkspils_mem->s_P_data = ark_mem->ark_user_data;

  /* Set default values for the rest of the Pcg parameters */
  arkspils_mem->s_eplifac   = ARKSPILS_EPLIN;
  arkspils_mem->s_last_flag = ARKSPILS_SUCCESS;
  ark_mem->ark_setupNonNull = FALSE;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    ARKProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKPCG", 
		    "ARKPcg", MSGS_BAD_PRETYPE);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_ILL_INPUT);
  }

  /* Allocate memory for ytemp and x */
  arkspils_mem->s_ytemp = N_VClone(ark_mem->ark_tempv);
  if (arkspils_mem->s_ytemp == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKPCG", 
		    "ARKPcg", MSGS_MEM_FAIL);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }

  arkspils_mem->s_x = N_VClone(ark_mem->ark_tempv);
  if (arkspils_mem->s_x == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKPCG", 
		    "ARKPcg", MSGS_MEM_FAIL);
    N_VDestroy(arkspils_mem->s_ytemp);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, arkspils_mem->s_ytemp);
  arkspils_mem->s_sqrtN = RSqrt(N_VDotProd(arkspils_mem->s_ytemp, 
					   arkspils_mem->s_ytemp));

  /* Call PcgMalloc to allocate workspace for Pcg */
  pcg_mem = NULL;
  pcg_mem = PcgMalloc(mxl, ark_mem->ark_tempv);
  if (pcg_mem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKPCG", 
		    "ARKPcg", MSGS_MEM_FAIL);
    N_VDestroy(arkspils_mem->s_ytemp);
    N_VDestroy(arkspils_mem->s_x);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }
  
  /* Attach PCG memory to spils memory structure */
  arkspils_mem->s_spils_mem = (void *) pcg_mem;

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_lmem = arkspils_mem;

  return(ARKSPILS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKPcgInit:

 This routine does remaining initializations specific to the 
 Pcg linear solver.
---------------------------------------------------------------*/
static int ARKPcgInit(ARKodeMem ark_mem)
{
  ARKSpilsMem arkspils_mem;
  PcgMem pcg_mem;

  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;
  pcg_mem = (PcgMem) arkspils_mem->s_spils_mem;

  /* Initialize counters */
  arkspils_mem->s_npe = arkspils_mem->s_nli = 0;
  arkspils_mem->s_nps = arkspils_mem->s_ncfl = 0;
  arkspils_mem->s_nstlpre = arkspils_mem->s_njtimes = 0;
  arkspils_mem->s_nfes = 0;

  /* Check for legal combination pretype - psolve */
  if ((arkspils_mem->s_pretype != PREC_NONE) && 
      (arkspils_mem->s_psolve == NULL)) {
    ARKProcessError(ark_mem, -1, "ARKPCG", 
		    "ARKPcgInit", MSGS_PSOLVE_REQ);
    arkspils_mem->s_last_flag = ARKSPILS_ILL_INPUT;
    return(-1);
  }

  /* Set setupNonNull = TRUE iff there is preconditioning
     (pretype != PREC_NONE)  and there is a preconditioning
     setup phase (pset != NULL) */
  ark_mem->ark_setupNonNull = (arkspils_mem->s_pretype != PREC_NONE) 
    && (arkspils_mem->s_pset != NULL);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (arkspils_mem->s_jtimesDQ) {
    arkspils_mem->s_jtimes = ARKSpilsDQJtimes;
    arkspils_mem->s_j_data = ark_mem;
  } else {
    arkspils_mem->s_j_data = ark_mem->ark_user_data;
  }

  /*  Set maxl in the PCG memory in case it was changed by the user */
  pcg_mem->l_max = arkspils_mem->s_maxl;

  arkspils_mem->s_last_flag = ARKSPILS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 ARKPcgSetup:

 This routine does the setup operations for the Pcg linear 
 solver. It makes a decision as to whether or not to signal for 
 reevaluation of Jacobian data in the pset routine, based on 
 various state variables, then it calls pset. If we signal for 
 reevaluation, then we reset jcur = *jcurPtr to TRUE, regardless 
 of the pset output. In any case, if jcur == TRUE, we increment 
 npe and save nst in nstlpre.
---------------------------------------------------------------*/
static int ARKPcgSetup(ARKodeMem ark_mem, int convfail, 
		       N_Vector ypred, N_Vector fpred, 
		       booleantype *jcurPtr, N_Vector vtemp1,
		       N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  int  retval;
  ARKSpilsMem arkspils_mem;

  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((ark_mem->ark_gamma/ark_mem->ark_gammap) - ONE);
  jbad = (ark_mem->ark_nst == 0) || 
    (ark_mem->ark_nst > arkspils_mem->s_nstlpre + ARKSPILS_MSBPRE) ||
    ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKSPILS_DGMAX)) ||
    (convfail == ARK_FAIL_OTHER);
  *jcurPtr = jbad;
  jok = !jbad;

  /* Call pset routine and possibly reset jcur */
  retval = arkspils_mem->s_pset(ark_mem->ark_tn, ypred, fpred, jok, 
				jcurPtr, ark_mem->ark_gamma, 
				arkspils_mem->s_P_data, vtemp1, 
				vtemp2, vtemp3);
  if (retval < 0) {
    ARKProcessError(ark_mem, PCG_PSET_FAIL_UNREC, "ARKPCG", 
		    "ARKPcgSetup", MSGS_PSET_FAILED);
    arkspils_mem->s_last_flag = PCG_PSET_FAIL_UNREC;
  }
  if (retval > 0) {
    arkspils_mem->s_last_flag = PCG_PSET_FAIL_REC;
  }

  if (jbad) *jcurPtr = TRUE;

  /* If jcur = TRUE, increment npe and save nst value */
  if (*jcurPtr) {
    arkspils_mem->s_npe++;
    arkspils_mem->s_nstlpre = ark_mem->ark_nst;
  }

  arkspils_mem->s_last_flag = PCG_SUCCESS;

  /* Return the same value that pset returned */
  return(retval);
}


/*---------------------------------------------------------------
 ARKPcgSolve:

 This routine handles the call to the generic solver PcgSolve
 for the solution of the linear system Ax = b with the PCG 
 method. The solution x is returned in the vector b.

 If the WRMS norm of b is small, we return x = b (if this is the 
 first Newton iteration) or x = 0 (if a later Newton iteration).

 Otherwise, we set the tolerance parameter and initial guess 
 (x = 0), call PcgSolve, and copy the solution x into b. The 
 x-scaling and b-scaling arrays are both equal to weight.

 The counters nli, nps, and ncfl are incremented, and the return 
 value is set according to the success of PcgSolve. The 
 success flag is returned if PcgSolve converged, or if this is 
 the first Newton iteration and the residual norm was reduced 
 below its initial value.
---------------------------------------------------------------*/
static int ARKPcgSolve(ARKodeMem ark_mem, N_Vector b, 
		       N_Vector weight, N_Vector ynow, 
		       N_Vector fnow)
{
  realtype bnorm, res_norm;
  ARKSpilsMem arkspils_mem;
  PcgMem pcg_mem;
  int nli_inc, nps_inc, retval;
  
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;
  pcg_mem = (PcgMem) arkspils_mem->s_spils_mem;

  /* Test norm(b); if small, return x = 0 or x = b */
  arkspils_mem->s_deltar = arkspils_mem->s_eplifac * ark_mem->ark_eLTE; 
  bnorm = N_VWrmsNorm(b, weight);
  if (bnorm <= arkspils_mem->s_deltar) {
    if (ark_mem->ark_mnewt > 0) N_VConst(ZERO, b); 
    return(0);
  }

#ifdef DEBUG_OUTPUT
  printf("ARKPcgSolve: bnorm = %19.16g, s_deltar = %19.16g,  s_sqrtN = %19.16g\n", 
	 bnorm, arkspils_mem->s_deltar, arkspils_mem->s_sqrtN);
#endif
      
  /* Set vectors ycur and fcur for use by the Atimes and Psolve routines */
  arkspils_mem->s_ycur = ynow;
  arkspils_mem->s_fcur = fnow;

  /* Set inputs delta and initial guess x = 0 to PcgSolve */
  arkspils_mem->s_delta = arkspils_mem->s_deltar * arkspils_mem->s_sqrtN;
  N_VConst(ZERO, arkspils_mem->s_x);

#ifdef DEBUG_OUTPUT
 printf("ARKPcgSolve: delta = %19.16g\n", arkspils_mem->s_delta);
 printf("  ||guess|| = %19.16g\n", N_VWrmsNorm(arkspils_mem->s_x, weight));
#endif
  
  /* Call PcgSolve and copy x to b */
  retval = PcgSolve(pcg_mem, ark_mem, arkspils_mem->s_x, b, 
		    arkspils_mem->s_pretype, arkspils_mem->s_delta,
		    ark_mem, weight, ARKSpilsAtimes, 
		    ARKSpilsPSolve, &res_norm, &nli_inc, &nps_inc);
  N_VScale(ONE, arkspils_mem->s_x, b);
  
#ifdef DEBUG_OUTPUT
 printf("ARKPcgSolve: retval = %i,  res_norm = %19.16g,  nli_inc = %i,  nps_inc = %i\n", 
	retval, res_norm, nli_inc, nps_inc);
 printf("  ||x|| = %19.16g\n", N_VWrmsNorm(arkspils_mem->s_x, weight));
#endif
  
  /* Increment counters nli, nps, and ncfl */
  arkspils_mem->s_nli += nli_inc;
  arkspils_mem->s_nps += nps_inc;
  if (retval != PCG_SUCCESS) arkspils_mem->s_ncfl++;

  /* Interpret return value from SpbcgSolve */
  arkspils_mem->s_last_flag = retval;

  switch(retval) {

  case PCG_SUCCESS:
    return 0;
    break;
  case PCG_RES_REDUCED:
    /* allow reduction but not solution on first Newton iteration, 
       otherwise return with a recoverable failure */
    if (ark_mem->ark_mnewt == 0) return(0);
    else                         return(1);
    break;
  case PCG_CONV_FAIL:
    return 1;
    break;
  case PCG_PSOLVE_FAIL_REC:
    return 1;
    break;
  case PCG_ATIMES_FAIL_REC:
    return 1;
    break;
  case PCG_MEM_NULL:
    return -1;
    break;
  case PCG_ATIMES_FAIL_UNREC:
    ARKProcessError(ark_mem, PCG_ATIMES_FAIL_UNREC, "ARKPCG", 
		    "ARKPcgSolve", MSGS_JTIMES_FAILED);    
    return -1;
    break;
  case PCG_PSOLVE_FAIL_UNREC:
    ARKProcessError(ark_mem, PCG_PSOLVE_FAIL_UNREC, "ARKPCG", 
		    "ARKPcgSolve", MSGS_PSOLVE_FAILED);
    return -1;
    break;
  }

  return 0;
}


/*---------------------------------------------------------------
 ARKPcgFree:

 This routine frees memory specific to the Pcg linear solver.
---------------------------------------------------------------*/
static void ARKPcgFree(ARKodeMem ark_mem)
{
  ARKSpilsMem arkspils_mem;
  PcgMem pcg_mem;

  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  N_VDestroy(arkspils_mem->s_ytemp);
  N_VDestroy(arkspils_mem->s_x);

  pcg_mem = (PcgMem) arkspils_mem->s_spils_mem;
  PcgFree(pcg_mem);

  if (arkspils_mem->s_pfree != NULL) (arkspils_mem->s_pfree)(ark_mem);

  free(arkspils_mem);
  ark_mem->ark_lmem = NULL;
}


/*---------------------------------------------------------------
     EOF
---------------------------------------------------------------*/
