/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * This is the implementation file for the ARKSPILS linear solvers.
 * -----------------------------------------------------------------
 ***** UNTOUCHED *****
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_spils_impl.h"

/* Private constants */

#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define ONE    RCONST(1.0)

/* Algorithmic constants */

#define MAX_ITERS  3  /* max. number of attempts to recover in DQ J*v */

/* Readability Replacements */

#define lrw1      (ark_mem->ark_lrw1)
#define liw1      (ark_mem->ark_liw1)
#define tq        (ark_mem->ark_tq)
#define tn        (ark_mem->ark_tn)
#define h         (ark_mem->ark_h)
#define gamma     (ark_mem->ark_gamma)
#define nfe       (ark_mem->ark_nfe)
#define f         (ark_mem->ark_f)
#define user_data (ark_mem->ark_user_data)
#define ewt       (ark_mem->ark_ewt)
#define lmem      (ark_mem->ark_lmem)

#define ils_type  (arkspils_mem->s_type)
#define sqrtN     (arkspils_mem->s_sqrtN)   
#define ytemp     (arkspils_mem->s_ytemp)
#define x         (arkspils_mem->s_x)
#define ycur      (arkspils_mem->s_ycur)
#define fcur      (arkspils_mem->s_fcur)
#define delta     (arkspils_mem->s_delta)
#define npe       (arkspils_mem->s_npe)
#define nli       (arkspils_mem->s_nli)
#define nps       (arkspils_mem->s_nps)
#define ncfl      (arkspils_mem->s_ncfl)
#define njtimes   (arkspils_mem->s_njtimes)
#define nfes      (arkspils_mem->s_nfes)

#define jtimesDQ  (arkspils_mem->s_jtimesDQ)
#define jtimes    (arkspils_mem->s_jtimes)
#define j_data    (arkspils_mem->s_j_data)

#define last_flag (arkspils_mem->s_last_flag)

/*
 * -----------------------------------------------------------------
 * OPTIONAL INPUT and OUTPUT
 * -----------------------------------------------------------------
 */


/*
 * -----------------------------------------------------------------
 * ARKSpilsSetPrecType
 * -----------------------------------------------------------------
 */

int ARKSpilsSetPrecType(void *arkode_mem, int pretype)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsSetPrecType", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsSetPrecType", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  /* Check for legal pretype */ 
  if ((pretype != PREC_NONE) && (pretype != PREC_LEFT) &&
      (pretype != PREC_RIGHT) && (pretype != PREC_BOTH)) {
    ARKProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", "ARKSpilsSetPrecType", MSGS_BAD_PRETYPE);
    return(ARKSPILS_ILL_INPUT);
  }

  arkspils_mem->s_pretype = pretype;

  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsSetGSType
 * -----------------------------------------------------------------
 */

int ARKSpilsSetGSType(void *arkode_mem, int gstype)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsSetGSType", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsSetGSType", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  if (ils_type != SPILS_SPGMR) {
    ARKProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", "ARKSpilsSetGSType", MSGS_BAD_LSTYPE);
    return(ARKSPILS_ILL_INPUT);
  }

  /* Check for legal gstype */
  if ((gstype != MODIFIED_GS) && (gstype != CLASSICAL_GS)) {
    ARKProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", "ARKSpilsSetGSType", MSGS_BAD_GSTYPE);
    return(ARKSPILS_ILL_INPUT);
  }

  arkspils_mem->s_gstype = gstype;

  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * Function : ARKSpilsSetMaxl
 * -----------------------------------------------------------------
 */

int ARKSpilsSetMaxl(void *arkode_mem, int maxl)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  int mxl;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsSetMaxl", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(NULL, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsSetMaxl", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  if (ils_type == SPILS_SPGMR) {
    ARKProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", "ARKSpilsSetMaxl", MSGS_BAD_LSTYPE);
    return(ARKSPILS_ILL_INPUT);
  }

  mxl = (maxl <= 0) ? ARKSPILS_MAXL : maxl;
  arkspils_mem->s_maxl = mxl;

  return(ARKSPILS_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * ARKSpilsSetEpsLin
 * -----------------------------------------------------------------
 */

int ARKSpilsSetEpsLin(void *arkode_mem, realtype eplifac)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsSetEpsLin", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsSetEpsLin", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  /* Check for legal eplifac */
  if(eplifac < ZERO) {
    ARKProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKSPILS", "ARKSpilsSetEpsLin", MSGS_BAD_EPLIN);
    return(ARKSPILS_ILL_INPUT);
  }

  arkspils_mem->s_eplifac = (eplifac == ZERO) ? ARKSPILS_EPLIN : eplifac;

  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsSetPrecSetupFn
 * -----------------------------------------------------------------
 */

int ARKSpilsSetPreconditioner(void *arkode_mem, 
                             ARKSpilsPrecSetupFn pset, ARKSpilsPrecSolveFn psolve)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsSetPreconditioner", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsSetPreconditioner", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  arkspils_mem->s_pset = pset;
  arkspils_mem->s_psolve = psolve;

  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsSetJacTimesVecFn
 * -----------------------------------------------------------------
 */

int ARKSpilsSetJacTimesVecFn(void *arkode_mem, ARKSpilsJacTimesVecFn jtv)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsSetJacTimesVecFn", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsSetJacTimesVecFn", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  if (jtv != NULL) {
    jtimesDQ = FALSE;
    jtimes = jtv;
  } else {
    jtimesDQ = TRUE;
  }

  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsGetWorkSpace
 * -----------------------------------------------------------------
 */

int ARKSpilsGetWorkSpace(void *arkode_mem, long int *lenrwLS, long int *leniwLS)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  int maxl;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsGetWorkSpace", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsGetWorkSpace", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  
  switch(ils_type) {
  case SPILS_SPGMR:
    maxl = arkspils_mem->s_maxl;
    *lenrwLS = lrw1*(maxl + 5) + maxl*(maxl + 4) + 1;
    *leniwLS = liw1*(maxl + 5);
    break;
  case SPILS_SPBCG:
    *lenrwLS = lrw1 * 9;
    *leniwLS = liw1 * 9;
    break;
  case SPILS_SPTFQMR:
    *lenrwLS = lrw1*11;
    *leniwLS = liw1*11;
    break;
  }


  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsGetNumPrecEvals
 * -----------------------------------------------------------------
 */

int ARKSpilsGetNumPrecEvals(void *arkode_mem, long int *npevals)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsGetNumPrecEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsGetNumPrecEvals", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  *npevals = npe;

  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsGetNumPrecSolves
 * -----------------------------------------------------------------
 */

int ARKSpilsGetNumPrecSolves(void *arkode_mem, long int *npsolves)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsGetNumPrecSolves", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsGetNumPrecSolves", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  *npsolves = nps;

  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsGetNumLinIters
 * -----------------------------------------------------------------
 */

int ARKSpilsGetNumLinIters(void *arkode_mem, long int *nliters)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsGetNumLinIters", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsGetNumLinIters", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  *nliters = nli;

  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsGetNumConvFails
 * -----------------------------------------------------------------
 */

int ARKSpilsGetNumConvFails(void *arkode_mem, long int *nlcfails)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsGetNumConvFails", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsGetNumConvFails", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  *nlcfails = ncfl;

  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsGetNumJtimesEvals
 * -----------------------------------------------------------------
 */

int ARKSpilsGetNumJtimesEvals(void *arkode_mem, long int *njvevals)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsGetNumJtimesEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsGetNumJtimesEvals", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  *njvevals = njtimes;

  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsGetNumRhsEvals
 * -----------------------------------------------------------------
 */

int ARKSpilsGetNumRhsEvals(void *arkode_mem, long int *nfevalsLS)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsGetNumRhsEvals", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsGetNumRhsEvals", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  *nfevalsLS = nfes;

  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsGetLastFlag
 * -----------------------------------------------------------------
 */

int ARKSpilsGetLastFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKSPILS", "ARKSpilsGetLastFlag", MSGS_ARKMEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKSPILS", "ARKSpilsGetLastFlag", MSGS_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) lmem;

  *flag = last_flag;

  return(ARKSPILS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *ARKSpilsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case ARKSPILS_SUCCESS:
    sprintf(name,"ARKSPILS_SUCCESS");
    break; 
  case ARKSPILS_MEM_NULL:
    sprintf(name,"ARKSPILS_MEM_NULL");
    break;
  case ARKSPILS_LMEM_NULL:
    sprintf(name,"ARKSPILS_LMEM_NULL");
    break;
  case ARKSPILS_ILL_INPUT:
    sprintf(name,"ARKSPILS_ILL_INPUT");
    break;
  case ARKSPILS_MEM_FAIL:
    sprintf(name,"ARKSPILS_MEM_FAIL");
    break;
  case ARKSPILS_PMEM_NULL:
    sprintf(name,"ARKSPILS_PMEM_NULL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * -----------------------------------------------------------------
 * ARKSPILS private functions
 * -----------------------------------------------------------------
 */


/* Additional readability Replacements */

#define pretype (arkspils_mem->s_pretype)
#define eplifac (arkspils_mem->s_eplifac)
#define maxl    (arkspils_mem->s_maxl)
#define psolve  (arkspils_mem->s_psolve)
#define P_data  (arkspils_mem->s_P_data)

/*
 * -----------------------------------------------------------------
 * ARKSpilsAtimes
 * -----------------------------------------------------------------
 * This routine generates the matrix-vector product z = Mv, where
 * M = I - gamma*J. The product J*v is obtained by calling the jtimes 
 * routine. It is then scaled by -gamma and added to v to obtain M*v.
 * The return value is the same as the value returned by jtimes --
 * 0 if successful, nonzero otherwise.
 * -----------------------------------------------------------------
 */

int ARKSpilsAtimes(void *arkode_mem, N_Vector v, N_Vector z)
{
  ARKodeMem   ark_mem;
  ARKSpilsMem arkspils_mem;
  int jtflag;

  ark_mem = (ARKodeMem) arkode_mem;
  arkspils_mem = (ARKSpilsMem) lmem;

  jtflag = jtimes(v, z, tn, ycur, fcur, j_data, ytemp);
  njtimes++;
  if (jtflag != 0) return(jtflag);

  N_VLinearSum(ONE, v, -gamma, z, z);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsPSolve
 * -----------------------------------------------------------------
 * This routine interfaces between the generic Sp***Solve routine
 * (within the SPGMR, SPBCG, or SPTFQMR solver) and the
 * user's psolve routine.  It passes to psolve all required state 
 * information from arkode_mem.  Its return value is the same as that
 * returned by psolve. Note that the generic SP*** solver guarantees
 * that ARKSpilsPSolve will not be called in the case in which
 * preconditioning is not done. This is the only case in which the
 * user's psolve routine is allowed to be NULL.
 * -----------------------------------------------------------------
 */

int ARKSpilsPSolve(void *arkode_mem, N_Vector r, N_Vector z, int lr)
{
  ARKodeMem   ark_mem;
  ARKSpilsMem arkspils_mem;
  int retval;

  ark_mem = (ARKodeMem) arkode_mem;
  arkspils_mem = (ARKSpilsMem)lmem;

  /* This call is counted in nps within the ARKSp***Solve routine */
  retval = psolve(tn, ycur, fcur, r, z, gamma, delta, lr, P_data, ytemp);

  return(retval);     
}

/*
 * -----------------------------------------------------------------
 * ARKSpilsDQJtimes
 * -----------------------------------------------------------------
 * This routine generates a difference quotient approximation to
 * the Jacobian times vector f_y(t,y) * v. The approximation is 
 * Jv = vnrm[f(y + v/vnrm) - f(y)], where vnrm = (WRMS norm of v) is
 * input, i.e. the WRMS norm of v/vnrm is 1.
 * -----------------------------------------------------------------
 */

int ARKSpilsDQJtimes(N_Vector v, N_Vector Jv, realtype t, 
                    N_Vector y, N_Vector fy,
                    void *data, N_Vector work)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  realtype sig, siginv;
  int iter, retval;

  /* data is arkode_mem */
  ark_mem = (ARKodeMem) data;
  arkspils_mem = (ARKSpilsMem) lmem;

  /* Initialize perturbation to 1/||v|| */
  sig = ONE/N_VWrmsNorm(v, ewt);

  for (iter=0; iter<MAX_ITERS; iter++) {

    /* Set work = y + sig*v */
    N_VLinearSum(sig, v, ONE, y, work);

    /* Set Jv = f(tn, y+sig*v) */
    retval = f(t, work, Jv, user_data); 
    nfes++;
    if (retval == 0) break;
    if (retval < 0)  return(-1);

    sig *= PT25;
  }

  if (retval > 0) return(+1);

  /* Replace Jv by (Jv - fy)/sig */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, fy, Jv);

  return(0);
}
