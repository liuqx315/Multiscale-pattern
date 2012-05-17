/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * This is the implementation file for the ARKDIAG linear solver.
 * -----------------------------------------------------------------
 ***** UNTOUCHED *****
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "arkode_diag_impl.h"
#include "arkode_impl.h"

/* Other Constants */
  
#define FRACT RCONST(0.1)
#define ONE   RCONST(1.0)

/* ARKDIAG linit, lsetup, lsolve, and lfree routines */

static int ARKDiagInit(ARKodeMem ark_mem);

static int ARKDiagSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
                       N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                       N_Vector vtemp2, N_Vector vtemp3);

static int ARKDiagSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur);

static void ARKDiagFree(ARKodeMem ark_mem);

/* Readability Replacements */

#define lrw1      (ark_mem->ark_lrw1)
#define liw1      (ark_mem->ark_liw1)
#define f         (ark_mem->ark_f)
#define uround    (ark_mem->ark_uround)
#define tn        (ark_mem->ark_tn)
#define h         (ark_mem->ark_h)
#define rl1       (ark_mem->ark_rl1)
#define gamma     (ark_mem->ark_gamma)
#define ewt       (ark_mem->ark_ewt)
#define nfe       (ark_mem->ark_nfe)
#define zn        (ark_mem->ark_zn)
#define linit     (ark_mem->ark_linit)
#define lsetup    (ark_mem->ark_lsetup)
#define lsolve    (ark_mem->ark_lsolve)
#define lfree     (ark_mem->ark_lfree)
#define lmem      (ark_mem->ark_lmem)
#define vec_tmpl  (ark_mem->ark_tempv)
#define setupNonNull   (ark_mem->ark_setupNonNull)

#define gammasv   (arkdiag_mem->di_gammasv)
#define M         (arkdiag_mem->di_M)
#define bit       (arkdiag_mem->di_bit)
#define bitcomp   (arkdiag_mem->di_bitcomp)
#define nfeDI     (arkdiag_mem->di_nfeDI)
#define last_flag (arkdiag_mem->di_last_flag)

/*
 * -----------------------------------------------------------------
 * ARKDiag 
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the diagonal linear solver module.  ARKDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the ark_linit, ark_lsetup, ark_lsolve, ark_lfree fields in (*arkode_mem)
 * to be ARKDiagInit, ARKDiagSetup, ARKDiagSolve, and ARKDiagFree,
 * respectively.  It allocates memory for a structure of type
 * ARKDiagMemRec and sets the ark_lmem field in (*arkode_mem) to the
 * address of this structure.  It sets setupNonNull in (*arkode_mem) to
 * TRUE.  Finally, it allocates memory for M, bit, and bitcomp.
 * The ARKDiag return value is SUCCESS = 0, LMEM_FAIL = -1, or 
 * LIN_ILL_INPUT=-2.
 * -----------------------------------------------------------------
 */
  
int ARKDiag(void *arkode_mem)
{
  ARKodeMem ark_mem;
  ARKDiagMem arkdiag_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKDIAG_MEM_NULL, "ARKDIAG", "ARKDiag", MSGDG_ARKMEM_NULL);
    return(ARKDIAG_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check if N_VCompare and N_VInvTest are present */
  if(vec_tmpl->ops->nvcompare == NULL ||
     vec_tmpl->ops->nvinvtest == NULL) {
    ARKProcessError(ark_mem, ARKDIAG_ILL_INPUT, "ARKDIAG", "ARKDiag", MSGDG_BAD_NVECTOR);
    return(ARKDIAG_ILL_INPUT);
  }

  if (lfree != NULL) lfree(ark_mem);
  
  /* Set four main function fields in ark_mem */
  linit  = ARKDiagInit;
  lsetup = ARKDiagSetup;
  lsolve = ARKDiagSolve;
  lfree  = ARKDiagFree;

  /* Get memory for ARKDiagMemRec */
  arkdiag_mem = NULL;
  arkdiag_mem = (ARKDiagMem) malloc(sizeof(ARKDiagMemRec));
  if (arkdiag_mem == NULL) {
    ARKProcessError(ark_mem, ARKDIAG_MEM_FAIL, "ARKDIAG", "ARKDiag", MSGDG_MEM_FAIL);
    return(ARKDIAG_MEM_FAIL);
  }

  last_flag = ARKDIAG_SUCCESS;

  /* Set flag setupNonNull = TRUE */
  setupNonNull = TRUE;

  /* Allocate memory for M, bit, and bitcomp */
    
  M = N_VClone(vec_tmpl);
  if (M == NULL) {
    ARKProcessError(ark_mem, ARKDIAG_MEM_FAIL, "ARKDIAG", "ARKDiag", MSGDG_MEM_FAIL);
    free(arkdiag_mem); arkdiag_mem = NULL;
    return(ARKDIAG_MEM_FAIL);
  }

  bit = N_VClone(vec_tmpl);
  if (bit == NULL) {
    ARKProcessError(ark_mem, ARKDIAG_MEM_FAIL, "ARKDIAG", "ARKDiag", MSGDG_MEM_FAIL);
    N_VDestroy(M);
    free(arkdiag_mem); arkdiag_mem = NULL;
    return(ARKDIAG_MEM_FAIL);
  }

  bitcomp = N_VClone(vec_tmpl);
  if (bitcomp == NULL) {
    ARKProcessError(ark_mem, ARKDIAG_MEM_FAIL, "ARKDIAG", "ARKDiag", MSGDG_MEM_FAIL);
    N_VDestroy(M);
    N_VDestroy(bit);
    free(arkdiag_mem); arkdiag_mem = NULL;
    return(ARKDIAG_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = arkdiag_mem;

  return(ARKDIAG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKDiagGetWorkSpace
 * -----------------------------------------------------------------
 */

int ARKDiagGetWorkSpace(void *arkode_mem, long int *lenrwLS, long int *leniwLS)
{
  ARKodeMem ark_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKDIAG_MEM_NULL, "ARKDIAG", "ARKDiagGetWorkSpace", MSGDG_ARKMEM_NULL);
    return(ARKDIAG_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *lenrwLS = 3*lrw1;
  *leniwLS = 3*liw1;

  return(ARKDIAG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKDiagGetNumRhsEvals
 * -----------------------------------------------------------------
 */

int ARKDiagGetNumRhsEvals(void *arkode_mem, long int *nfevalsLS)
{
  ARKodeMem ark_mem;
  ARKDiagMem arkdiag_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKDIAG_MEM_NULL, "ARKDIAG", "ARKDiagGetNumRhsEvals", MSGDG_ARKMEM_NULL);
    return(ARKDIAG_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKDIAG_LMEM_NULL, "ARKDIAG", "ARKDiagGetNumRhsEvals", MSGDG_LMEM_NULL);
    return(ARKDIAG_LMEM_NULL);
  }
  arkdiag_mem = (ARKDiagMem) lmem;

  *nfevalsLS = nfeDI;

  return(ARKDIAG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKDiagGetLastFlag
 * -----------------------------------------------------------------
 */

int ARKDiagGetLastFlag(void *arkode_mem, long int *flag)
{
  ARKodeMem ark_mem;
  ARKDiagMem arkdiag_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKDIAG_MEM_NULL, "ARKDIAG", "ARKDiagGetLastFlag", MSGDG_ARKMEM_NULL);
    return(ARKDIAG_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (lmem == NULL) {
    ARKProcessError(ark_mem, ARKDIAG_LMEM_NULL, "ARKDIAG", "ARKDiagGetLastFlag", MSGDG_LMEM_NULL);
    return(ARKDIAG_LMEM_NULL);
  }
  arkdiag_mem = (ARKDiagMem) lmem;

  *flag = last_flag;

  return(ARKDIAG_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKDiagGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *ARKDiagGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case ARKDIAG_SUCCESS:
    sprintf(name,"ARKDIAG_SUCCESS");
    break;   
  case ARKDIAG_MEM_NULL:
    sprintf(name,"ARKDIAG_MEM_NULL");
    break;
  case ARKDIAG_LMEM_NULL:
    sprintf(name,"ARKDIAG_LMEM_NULL");
    break;
  case ARKDIAG_ILL_INPUT:
    sprintf(name,"ARKDIAG_ILL_INPUT");
    break;
  case ARKDIAG_MEM_FAIL:
    sprintf(name,"ARKDIAG_MEM_FAIL");
    break;
  case ARKDIAG_INV_FAIL:
    sprintf(name,"ARKDIAG_INV_FAIL");
    break;
  case ARKDIAG_RHSFUNC_UNRECVR:
    sprintf(name,"ARKDIAG_RHSFUNC_UNRECVR");
    break;
  case ARKDIAG_RHSFUNC_RECVR:
    sprintf(name,"ARKDIAG_RHSFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * -----------------------------------------------------------------
 * ARKDiagInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the diagonal
 * linear solver.
 * -----------------------------------------------------------------
 */

static int ARKDiagInit(ARKodeMem ark_mem)
{
  ARKDiagMem arkdiag_mem;

  arkdiag_mem = (ARKDiagMem) lmem;

  nfeDI = 0;

  last_flag = ARKDIAG_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * ARKDiagSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the diagonal linear 
 * solver.  It constructs a diagonal approximation to the Newton matrix 
 * M = I - gamma*J, updates counters, and inverts M.
 * -----------------------------------------------------------------
 */

static int ARKDiagSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
                       N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                       N_Vector vtemp2, N_Vector vtemp3)
{
  realtype r;
  N_Vector ftemp, y;
  booleantype invOK;
  ARKDiagMem arkdiag_mem;
  int retval;

  arkdiag_mem = (ARKDiagMem) lmem;

  /* Rename work vectors for use as temporary values of y and f */
  ftemp = vtemp1;
  y     = vtemp2;

  /* Form y with perturbation = FRACT*(func. iter. correction) */
  r = FRACT * rl1;
  N_VLinearSum(h, fpred, -ONE, zn[1], ftemp);
  N_VLinearSum(r, ftemp, ONE, ypred, y);

  /* Evaluate f at perturbed y */
  retval = f(tn, y, M, ark_mem->ark_user_data);
  nfeDI++;
  if (retval < 0) {
    ARKProcessError(ark_mem, ARKDIAG_RHSFUNC_UNRECVR, "ARKDIAG", "ARKDiagSetup", MSGDG_RHSFUNC_FAILED);
    last_flag = ARKDIAG_RHSFUNC_UNRECVR;
    return(-1);
  }
  if (retval > 0) {
    last_flag = ARKDIAG_RHSFUNC_RECVR;
    return(1);
  }

  /* Construct M = I - gamma*J with J = diag(deltaf_i/deltay_i) */
  N_VLinearSum(ONE, M, -ONE, fpred, M);
  N_VLinearSum(FRACT, ftemp, -h, M, M);
  N_VProd(ftemp, ewt, y);
  /* Protect against deltay_i being at roundoff level */
  N_VCompare(uround, y, bit);
  N_VAddConst(bit, -ONE, bitcomp);
  N_VProd(ftemp, bit, y);
  N_VLinearSum(FRACT, y, -ONE, bitcomp, y);
  N_VDiv(M, y, M);
  N_VProd(M, bit, M);
  N_VLinearSum(ONE, M, -ONE, bitcomp, M);

  /* Invert M with test for zero components */
  invOK = N_VInvTest(M, M);
  if (!invOK) {
    last_flag = ARKDIAG_INV_FAIL;
    return(1);
  }

  /* Set jcur = TRUE, save gamma in gammasv, and return */
  *jcurPtr = TRUE;
  gammasv = gamma;
  last_flag = ARKDIAG_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * ARKDiagSolve
 * -----------------------------------------------------------------
 * This routine performs the solve operation for the diagonal linear
 * solver.  If necessary it first updates gamma in M = I - gamma*J.
 * -----------------------------------------------------------------
 */

static int ARKDiagSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur)
{
  booleantype invOK;
  realtype r;
  ARKDiagMem arkdiag_mem;

  arkdiag_mem = (ARKDiagMem) lmem;
  
  /* If gamma has changed, update factor in M, and save gamma value */

  if (gammasv != gamma) {
    r = gamma / gammasv;
    N_VInv(M, M);
    N_VAddConst(M, -ONE, M);
    N_VScale(r, M, M);
    N_VAddConst(M, ONE, M);
    invOK = N_VInvTest(M, M);
    if (!invOK) {
      last_flag = ARKDIAG_INV_FAIL;
      return (1);
    }
    gammasv = gamma;
  }

  /* Apply M-inverse to b */
  N_VProd(b, M, b);

  last_flag = ARKDIAG_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * ARKDiagFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the diagonal linear solver.
 * -----------------------------------------------------------------
 */

static void ARKDiagFree(ARKodeMem ark_mem)
{
  ARKDiagMem arkdiag_mem;
  
  arkdiag_mem = (ARKDiagMem) lmem;

  N_VDestroy(M);
  N_VDestroy(bit);
  N_VDestroy(bitcomp);
  free(arkdiag_mem);
  ark_mem->ark_lmem = NULL;
}
