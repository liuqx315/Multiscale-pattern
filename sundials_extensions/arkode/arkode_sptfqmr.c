/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * This is the implementation file for the ARKSPTFQMR linear solver.
 * -----------------------------------------------------------------
 ***** UNTOUCHED *****
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <arkode/arkode_sptfqmr.h>
#include "arkode_spils_impl.h"
#include "arkode_impl.h"

#include <sundials/sundials_sptfqmr.h>
#include <sundials/sundials_math.h>

/* Other Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* ARKSPTFQMR linit, lsetup, lsolve, and lfree routines */

static int ARKSptfqmrInit(ARKodeMem ark_mem);

static int ARKSptfqmrSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
			  N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
			  N_Vector vtemp2, N_Vector vtemp3);

static int ARKSptfqmrSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
			  N_Vector ynow, N_Vector fnow);

static void ARKSptfqmrFree(ARKodeMem ark_mem);


/* Readability Replacements */

#define tq           (ark_mem->ark_tq)
#define nst          (ark_mem->ark_nst)
#define tn           (ark_mem->ark_tn)
#define gamma        (ark_mem->ark_gamma)
#define gammap       (ark_mem->ark_gammap)
#define f            (ark_mem->ark_f)
#define user_data    (ark_mem->ark_user_data)
#define ewt          (ark_mem->ark_ewt)
#define errfp        (ark_mem->ark_errfp)
#define mnewt        (ark_mem->ark_mnewt)
#define linit        (ark_mem->ark_linit)
#define lsetup       (ark_mem->ark_lsetup)
#define lsolve       (ark_mem->ark_lsolve)
#define lfree        (ark_mem->ark_lfree)
#define lmem         (ark_mem->ark_lmem)
#define vec_tmpl     (ark_mem->ark_tempv)
#define setupNonNull (ark_mem->ark_setupNonNull)

#define sqrtN       (arkspils_mem->s_sqrtN)   
#define ytemp       (arkspils_mem->s_ytemp)
#define x           (arkspils_mem->s_x)
#define ycur        (arkspils_mem->s_ycur)
#define fcur        (arkspils_mem->s_fcur)
#define delta       (arkspils_mem->s_delta)
#define deltar      (arkspils_mem->s_deltar)
#define npe         (arkspils_mem->s_npe)
#define nli         (arkspils_mem->s_nli)
#define nps         (arkspils_mem->s_nps)
#define ncfl        (arkspils_mem->s_ncfl)
#define nstlpre     (arkspils_mem->s_nstlpre)
#define njtimes     (arkspils_mem->s_njtimes)
#define nfes        (arkspils_mem->s_nfes)
#define spils_mem   (arkspils_mem->s_spils_mem)

#define jtimesDQ (arkspils_mem->s_jtimesDQ)
#define jtimes  (arkspils_mem->s_jtimes)
#define j_data  (arkspils_mem->s_j_data)

#define last_flag   (arkspils_mem->s_last_flag)

/*
 * -----------------------------------------------------------------
 * Function : ARKSptfqmr
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the Sptfqmr linear solver module. ARKSptfqmr first
 * calls the existing lfree routine if this is not NULL. It then sets
 * the ark_linit, ark_lsetup, ark_lsolve, ark_lfree fields in (*arkode_mem)
 * to be ARKSptfqmrInit, ARKSptfqmrSetup, ARKSptfqmrSolve, and ARKSptfqmrFree,
 * respectively. It allocates memory for a structure of type
 * ARKSpilsMemRec and sets the ark_lmem field in (*arkode_mem) to the
 * address of this structure. It sets setupNonNull in (*arkode_mem),
 * and sets various fields in the ARKSpilsMemRec structure.
 * Finally, ARKSptfqmr allocates memory for ytemp and x, and calls
 * SptfqmrMalloc to allocate memory for the Sptfqmr solver.
 * -----------------------------------------------------------------
 */

int ARKSptfqmr(void *arkode_mem, int pretype, int maxl)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  SptfqmrMem sptfqmr_mem;
  int mxl;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPTFQMR", "ARKSptfqmr", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check if N_VDotProd is present */
  if (vec_tmpl->ops->nvdotprod == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPTFQMR", "ARKSptfqmr", MSGS_BAD_NVECTOR);
    return(ARKSPILS_ILL_INPUT);
  }

  if (lfree != NULL) lfree(ark_mem);

  /* Set four main function fields in ark_mem */
  linit  = ARKSptfqmrInit;
  lsetup = ARKSptfqmrSetup;
  lsolve = ARKSptfqmrSolve;
  lfree  = ARKSptfqmrFree;

  /* Get memory for ARKSpilsMemRec */
  arkspils_mem = NULL;
  arkspils_mem = (ARKSpilsMem) malloc(sizeof(struct ARKSpilsMemRec));
  if (arkspils_mem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPTFQMR", "ARKSptfqmr", MSGS_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* Set ILS type */
  arkspils_mem->s_type = SPILS_SPTFQMR;

  /* Set Sptfqmr parameters that have been passed in call sequence */
  arkspils_mem->s_pretype = pretype;
  mxl = arkspils_mem->s_maxl = (maxl <= 0) ? ARKSPILS_MAXL : maxl;

  /* Set defaults for Jacobian-related fileds */
  jtimesDQ = TRUE;
  jtimes   = NULL;
  j_data   = NULL;

  /* Set defaults for preconditioner-related fields */
  arkspils_mem->s_pset   = NULL;
  arkspils_mem->s_psolve = NULL;
  arkspils_mem->s_pfree  = NULL;
  arkspils_mem->s_P_data = ark_mem->ark_user_data;

  /* Set default values for the rest of the Sptfqmr parameters */
  arkspils_mem->s_eplifac = ARKSPILS_EPLIN;

  arkspils_mem->s_last_flag = ARKSPILS_SUCCESS;

  setupNonNull = FALSE;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    ARKProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPTFQMR", "ARKSptfqmr", MSGS_BAD_PRETYPE);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_ILL_INPUT);
  }

  /* Allocate memory for ytemp and x */

  ytemp = N_VClone(vec_tmpl);
  if (ytemp == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPTFQMR", "ARKSptfqmr", MSGS_MEM_FAIL);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }

  x = N_VClone(vec_tmpl);
  if (x == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPTFQMR", "ARKSptfqmr", MSGS_MEM_FAIL);
    N_VDestroy(ytemp);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }

  /* Compute sqrtN from a dot product */
  N_VConst(ONE, ytemp);
  sqrtN = RSqrt(N_VDotProd(ytemp, ytemp));

  /* Call SptfqmrMalloc to allocate workspace for Sptfqmr */
  sptfqmr_mem = NULL;
  sptfqmr_mem = SptfqmrMalloc(mxl, vec_tmpl);
  if (sptfqmr_mem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKSPTFQMR", "ARKSptfqmr", MSGS_MEM_FAIL);
    N_VDestroy(ytemp);
    N_VDestroy(x);
    free(arkspils_mem); arkspils_mem = NULL;
    return(ARKSPILS_MEM_FAIL);
  }
  
  /* Attach SPTFQMR memory to spils memory structure */
  spils_mem = (void *) sptfqmr_mem;

  /* Attach linear solver memory to integrator memory */
  lmem = arkspils_mem;

  return(ARKSPILS_SUCCESS);
}

/* Additional readability replacements */

#define pretype (arkspils_mem->s_pretype)
#define eplifac (arkspils_mem->s_eplifac)
#define maxl    (arkspils_mem->s_maxl)
#define psolve  (arkspils_mem->s_psolve)
#define pset    (arkspils_mem->s_pset)
#define P_data  (arkspils_mem->s_P_data)

/*
 * -----------------------------------------------------------------
 * Function : ARKSptfqmrInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the Sptfqmr
 * linear solver.
 * -----------------------------------------------------------------
 */

static int ARKSptfqmrInit(ARKodeMem ark_mem)
{
  ARKSpilsMem arkspils_mem;
  SptfqmrMem sptfqmr_mem;

  arkspils_mem = (ARKSpilsMem) lmem;
  sptfqmr_mem = (SptfqmrMem) spils_mem;

  /* Initialize counters */
  npe = nli = nps = ncfl = nstlpre = 0;
  njtimes = nfes = 0;

  /* Check for legal combination pretype - psolve */
  if ((pretype != PREC_NONE) && (psolve == NULL)) {
    ARKProcessError(ark_mem, -1, "ARKSPTFQMR", "ARKSptfqmrInit", MSGS_PSOLVE_REQ);
    last_flag = ARKSPILS_ILL_INPUT;
    return(-1);
  }

  /* Set setupNonNull = TRUE iff there is preconditioning
     (pretype != PREC_NONE)  and there is a preconditioning
     setup phase (pset != NULL) */
  setupNonNull = (pretype != PREC_NONE) && (pset != NULL);

  /* Set Jacobian-related fields, based on jtimesDQ */
  if (jtimesDQ) {
    jtimes = ARKSpilsDQJtimes;
    j_data = ark_mem;
  } else {
    j_data = user_data;
  }

  /*  Set maxl in the SPTFQMR memory in case it was changed by the user */
  sptfqmr_mem->l_max  = maxl;

  last_flag = ARKSPILS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : ARKSptfqmrSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the Sptfqmr linear solver.
 * It makes a decision as to whether or not to signal for reevaluation
 * of Jacobian data in the pset routine, based on various state
 * variables, then it calls pset. If we signal for reevaluation,
 * then we reset jcur = *jcurPtr to TRUE, regardless of the pset output.
 * In any case, if jcur == TRUE, we increment npe and save nst in nstlpre.
 * -----------------------------------------------------------------
 */

static int ARKSptfqmrSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
			  N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
			  N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  int  retval;
  ARKSpilsMem arkspils_mem;

  arkspils_mem = (ARKSpilsMem) lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlpre + ARKSPILS_MSBPRE) ||
      ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKSPILS_DGMAX)) ||
      (convfail == ARK_FAIL_OTHER);
  *jcurPtr = jbad;
  jok = !jbad;

  /* Call pset routine and possibly reset jcur */
  retval = pset(tn, ypred, fpred, jok, jcurPtr, gamma, P_data, 
                vtemp1, vtemp2, vtemp3);
  if (retval < 0) {
    ARKProcessError(ark_mem, SPTFQMR_PSET_FAIL_UNREC, "ARKSPTFQMR", "ARKSptfqmrSetup", MSGS_PSET_FAILED);
    last_flag = SPTFQMR_PSET_FAIL_UNREC;
  }
  if (retval > 0) {
    last_flag = SPTFQMR_PSET_FAIL_REC;
  }

  if (jbad) *jcurPtr = TRUE;

  /* If jcur = TRUE, increment npe and save nst value */
  if (*jcurPtr) {
    npe++;
    nstlpre = nst;
  }

  last_flag = SPTFQMR_SUCCESS;

  /* Return the same value that pset returned */
  return(retval);
}

/*
 * -----------------------------------------------------------------
 * Function : ARKSptfqmrSolve
 * -----------------------------------------------------------------
 * This routine handles the call to the generic solver SptfqmrSolve
 * for the solution of the linear system Ax = b with the SPTFQMR method.
 * The solution x is returned in the vector b.
 *
 * If the WRMS norm of b is small, we return x = b (if this is the first
 * Newton iteration) or x = 0 (if a later Newton iteration).
 *
 * Otherwise, we set the tolerance parameter and initial guess (x = 0),
 * call SptfqmrSolve, and copy the solution x into b. The x-scaling and
 * b-scaling arrays are both equal to weight.
 *
 * The counters nli, nps, and ncfl are incremented, and the return value
 * is set according to the success of SptfqmrSolve. The success flag is
 * returned if SptfqmrSolve converged, or if this is the first Newton
 * iteration and the residual norm was reduced below its initial value.
 * -----------------------------------------------------------------
 */

static int ARKSptfqmrSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
			  N_Vector ynow, N_Vector fnow)
{
  realtype bnorm, res_norm;
  ARKSpilsMem arkspils_mem;
  SptfqmrMem sptfqmr_mem;
  int nli_inc, nps_inc, retval;
  
  arkspils_mem = (ARKSpilsMem) lmem;

  sptfqmr_mem = (SptfqmrMem) spils_mem;

  /* Test norm(b); if small, return x = 0 or x = b */
  deltar = eplifac * tq[4]; 

  bnorm = N_VWrmsNorm(b, weight);
  if (bnorm <= deltar) {
    if (mnewt > 0) N_VConst(ZERO, b); 
    return(0);
  }

  /* Set vectors ycur and fcur for use by the Atimes and Psolve routines */
  ycur = ynow;
  fcur = fnow;

  /* Set inputs delta and initial guess x = 0 to SptfqmrSolve */  
  delta = deltar * sqrtN;
  N_VConst(ZERO, x);
  
  /* Call SptfqmrSolve and copy x to b */
  retval = SptfqmrSolve(sptfqmr_mem, ark_mem, x, b, pretype, delta,
                        ark_mem, weight, weight, ARKSpilsAtimes, ARKSpilsPSolve,
                        &res_norm, &nli_inc, &nps_inc);

  N_VScale(ONE, x, b);
  
  /* Increment counters nli, nps, and ncfl */
  nli += nli_inc;
  nps += nps_inc;
  if (retval != SPTFQMR_SUCCESS) ncfl++;

  /* Interpret return value from SpgmrSolve */

  last_flag = retval;

  switch(retval) {

  case SPTFQMR_SUCCESS:
    return(0);
    break;
  case SPTFQMR_RES_REDUCED:
    if (mnewt == 0) return(0);
    else            return(1);
    break;
  case SPTFQMR_CONV_FAIL:
    return(1);
    break;
  case SPTFQMR_PSOLVE_FAIL_REC:
    return(1);
    break;
  case SPTFQMR_ATIMES_FAIL_REC:
    return(1);
    break;
  case SPTFQMR_MEM_NULL:
    return(-1);
    break;
  case SPTFQMR_ATIMES_FAIL_UNREC:
    ARKProcessError(ark_mem, SPTFQMR_ATIMES_FAIL_UNREC, "ARKSPTFQMR", "ARKSptfqmrSolve", MSGS_JTIMES_FAILED);    
    return(-1);
    break;
  case SPTFQMR_PSOLVE_FAIL_UNREC:
    ARKProcessError(ark_mem, SPTFQMR_PSOLVE_FAIL_UNREC, "ARKSPTFQMR", "ARKSptfqmrSolve", MSGS_PSOLVE_FAILED);
    return(-1);
    break;
  }

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : ARKSptfqmrFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the Sptfqmr linear solver.
 * -----------------------------------------------------------------
 */

static void ARKSptfqmrFree(ARKodeMem ark_mem)
{
  ARKSpilsMem arkspils_mem;
  SptfqmrMem sptfqmr_mem;
    
  arkspils_mem = (ARKSpilsMem) lmem;

  N_VDestroy(ytemp);
  N_VDestroy(x);

  sptfqmr_mem = (SptfqmrMem) spils_mem;
  SptfqmrFree(sptfqmr_mem);

  if (arkspils_mem->s_pfree != NULL) (arkspils_mem->s_pfree)(ark_mem);

  free(arkspils_mem);
  ark_mem->ark_lmem = NULL;

  return;
}

