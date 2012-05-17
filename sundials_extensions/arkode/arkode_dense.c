/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * This is the impleentation file for the ARKDENSE linear solver.
 * -----------------------------------------------------------------
 ***** UNTOUCHED *****
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <arkode/arkode_dense.h>
#include "arkode_direct_impl.h"
#include "arkode_impl.h"

#include <sundials/sundials_math.h>

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* ARKDENSE linit, lsetup, lsolve, and lfree routines */
 
static int arkDenseInit(ARKodeMem ark_mem);

static int arkDenseSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, 
                        N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

static int arkDenseSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                        N_Vector ycur, N_Vector fcur);

static void arkDenseFree(ARKodeMem ark_mem);

/* Readability Replacements */

#define lmm       (ark_mem->ark_lmm)
#define f         (ark_mem->ark_f)
#define nst       (ark_mem->ark_nst)
#define tn        (ark_mem->ark_tn)
#define h         (ark_mem->ark_h)
#define gamma     (ark_mem->ark_gamma)
#define gammap    (ark_mem->ark_gammap)
#define gamrat    (ark_mem->ark_gamrat)
#define ewt       (ark_mem->ark_ewt)
#define linit     (ark_mem->ark_linit)
#define lsetup    (ark_mem->ark_lsetup)
#define lsolve    (ark_mem->ark_lsolve)
#define lfree     (ark_mem->ark_lfree)
#define lmem      (ark_mem->ark_lmem)
#define vec_tmpl     (ark_mem->ark_tempv)
#define setupNonNull (ark_mem->ark_setupNonNull)

#define mtype     (arkdls_mem->d_type)
#define n         (arkdls_mem->d_n)
#define jacDQ     (arkdls_mem->d_jacDQ)
#define jac       (arkdls_mem->d_djac)
#define M         (arkdls_mem->d_M)
#define lpivots   (arkdls_mem->d_lpivots)
#define savedJ    (arkdls_mem->d_savedJ)
#define nstlj     (arkdls_mem->d_nstlj)
#define nje       (arkdls_mem->d_nje)
#define nfeDQ     (arkdls_mem->d_nfeDQ)
#define J_data    (arkdls_mem->d_J_data)
#define last_flag (arkdls_mem->d_last_flag)
                  
/*
 * -----------------------------------------------------------------
 * ARKDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense linear solver module.  ARKDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the ark_linit, ark_lsetup, ark_lsolve, ark_lfree fields in (*arkode_mem)
 * to be arkDenseInit, arkDenseSetup, arkDenseSolve, and arkDenseFree,
 * respectively.  It allocates memory for a structure of type
 * ARKDlsMemRec and sets the ark_lmem field in (*arkode_mem) to the
 * address of this structure.  It sets setupNonNull in (*arkode_mem) to
 * TRUE, and the d_jac field to the default arkDlsDenseDQJac.
 * Finally, it allocates memory for M, savedJ, and lpivots.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, ARKDense will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int ARKDense(void *arkode_mem, long int N)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKDLS_MEM_NULL, "ARKDENSE", "ARKDense", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL ||
      vec_tmpl->ops->nvsetarraypointer == NULL) {
    ARKProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDENSE", "ARKDense", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(ark_mem);

  /* Set four main function fields in ark_mem */
  linit  = arkDenseInit;
  lsetup = arkDenseSetup;
  lsolve = arkDenseSolve;
  lfree  = arkDenseFree;

  /* Get memory for ARKDlsMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMem) malloc(sizeof(struct ARKDlsMemRec));
  if (arkdls_mem == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDENSE", "ARKDense", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;

  /* Initialize Jacobian-related data */
  jacDQ = TRUE;
  jac = NULL;
  J_data = NULL;

  last_flag = ARKDLS_SUCCESS;

  setupNonNull = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for M, savedJ, and pivot array */

  M = NULL;
  M = NewDenseMat(N, N);
  if (M == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDENSE", "ARKDense", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  savedJ = NULL;
  savedJ = NewDenseMat(N, N);
  if (savedJ == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDENSE", "ARKDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  lpivots = NULL;
  lpivots = NewLintArray(N);
  if (lpivots == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDENSE", "ARKDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    DestroyMat(savedJ);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = arkdls_mem;

  return(ARKDLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * arkDenseInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int arkDenseInit(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;

  arkdls_mem = (ARKDlsMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;

  /* Set Jacobian function and data, depending on jacDQ */
  if (jacDQ) {
    jac = arkDlsDenseDQJac;
    J_data = ark_mem;
  } else {
    J_data = ark_mem->ark_user_data;
  }

  last_flag = ARKDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * arkDenseSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the dense linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy.  In any case, it constructs the Newton matrix 
 * M = I - gamma*J, updates counters, and calls the dense LU 
 * factorization routine.
 * -----------------------------------------------------------------
 */

static int arkDenseSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, 
                        N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  long int ier;
  ARKDlsMem arkdls_mem;
  int retval;

  arkdls_mem = (ARKDlsMem) lmem;
 
  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
 
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlj + ARKD_MSBJ) ||
         ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKD_DGMAX)) ||
         (convfail == ARK_FAIL_OTHER);
  jok = !jbad;
 
  if (jok) {

    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    DenseCopy(savedJ, M);

  } else {

    /* If jok = FALSE, call jac routine for new J value */
    nje++;
    nstlj = nst;
    *jcurPtr = TRUE;
    SetToZero(M);

    retval = jac(n, tn, ypred, fpred, M, J_data, vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      ARKProcessError(ark_mem, ARKDLS_JACFUNC_UNRECVR, "ARKDENSE", "arkDenseSetup", MSGD_JACFUNC_FAILED);
      last_flag = ARKDLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      last_flag = ARKDLS_JACFUNC_RECVR;
      return(1);
    }

    DenseCopy(M, savedJ);

  }
  
  /* Scale and add I to get M = I - gamma*J */
  DenseScale(-gamma, M);
  AddIdentity(M);

  /* Do LU factorization of M */
  ier = DenseGETRF(M, lpivots); 

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * arkDenseSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.  The returned value is 0.
 * -----------------------------------------------------------------
 */

static int arkDenseSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                        N_Vector ycur, N_Vector fcur)
{
  ARKDlsMem arkdls_mem;
  realtype *bd;

  arkdls_mem = (ARKDlsMem) lmem;
  
  bd = N_VGetArrayPointer(b);

  DenseGETRS(M, lpivots, bd);

  /* If ARK_BDF, scale the correction to account for change in gamma */
  if ((lmm == ARK_BDF) && (gamrat != ONE)) {
    N_VScale(TWO/(ONE + gamrat), b, b);
  }
  
  last_flag = ARKDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * arkDenseFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void arkDenseFree(ARKodeMem ark_mem)
{
  ARKDlsMem  arkdls_mem;

  arkdls_mem = (ARKDlsMem) lmem;
  
  DestroyMat(M);
  DestroyMat(savedJ);
  DestroyArray(lpivots);
  free(arkdls_mem);
  ark_mem->ark_lmem = NULL;
}

