/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * This is the implementation file for the ARKBAND linear solver.
 * -----------------------------------------------------------------
 ***** UNTOUCHED *****
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <arkode/arkode_band.h>
#include "arkode_direct_impl.h"
#include "arkode_impl.h"

#include <sundials/sundials_math.h>

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* ARKBAND linit, lsetup, lsolve, and lfree routines */

static int arkBandInit(ARKodeMem ark_mem);

static int arkBandSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
                       N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                       N_Vector vtemp2, N_Vector vtemp3);

static int arkBandSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur);

static void arkBandFree(ARKodeMem ark_mem);

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
#define nfe       (ark_mem->ark_nfe)
#define linit     (ark_mem->ark_linit)
#define lsetup    (ark_mem->ark_lsetup)
#define lsolve    (ark_mem->ark_lsolve)
#define lfree     (ark_mem->ark_lfree)
#define lmem      (ark_mem->ark_lmem)
#define vec_tmpl      (ark_mem->ark_tempv)
#define setupNonNull  (ark_mem->ark_setupNonNull)

#define mtype      (arkdls_mem->d_type)
#define n          (arkdls_mem->d_n)
#define jacDQ      (arkdls_mem->d_jacDQ)
#define jac        (arkdls_mem->d_bjac)
#define M          (arkdls_mem->d_M)
#define mu         (arkdls_mem->d_mu)
#define ml         (arkdls_mem->d_ml)
#define smu        (arkdls_mem->d_smu)
#define lpivots    (arkdls_mem->d_lpivots)
#define savedJ     (arkdls_mem->d_savedJ)
#define nstlj      (arkdls_mem->d_nstlj)
#define nje        (arkdls_mem->d_nje)
#define nfeDQ       (arkdls_mem->d_nfeDQ)
#define J_data     (arkdls_mem->d_J_data)
#define last_flag  (arkdls_mem->d_last_flag)

/*
 * -----------------------------------------------------------------
 * ARKBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module.  ARKBand first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * ark_linit, ark_lsetup, ark_lsolve, and ark_lfree fields in (*arkode_mem)
 * to be arkBandInit, arkBandSetup, arkBandSolve, and arkBandFree,
 * respectively.  It allocates memory for a structure of type
 * ARKDlsMemRec and sets the ark_lmem field in (*arkode_mem) to the
 * address of this structure.  It sets setupNonNull in (*arkode_mem) to be
 * TRUE, d_mu to be mupper, d_ml to be mlower, and the d_jac field to be 
 * arkDlsBandDQJac.
 * Finally, it allocates memory for M, savedJ, and pivot.  The ARKBand
 * return value is SUCCESS = 0, LMEM_FAIL = -1, or LIN_ILL_INPUT = -2.
 *
 * NOTE: The band linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, ARKBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */
                  
int ARKBand(void *arkode_mem, long int N, long int mupper, long int mlower)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKDLS_MEM_NULL, "ARKBAND", "ARKBand", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL) {
    ARKProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKBAND", "ARKBand", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  if (lfree != NULL) lfree(ark_mem);

  /* Set four main function fields in ark_mem */  
  linit  = arkBandInit;
  lsetup = arkBandSetup;
  lsolve = arkBandSolve;
  lfree  = arkBandFree;
  
  /* Get memory for ARKDlsMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMem) malloc(sizeof(struct ARKDlsMemRec));
  if (arkdls_mem == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKBand", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_BAND;
  
  /* Initialize Jacobian-related data */
  jacDQ = TRUE;
  jac = NULL;
  J_data = NULL;

  last_flag = ARKDLS_SUCCESS;

  setupNonNull = TRUE;
  
  /* Load problem dimension */
  n = N;

  /* Load half-bandwiths in arkdls_mem */
  ml = mlower;
  mu = mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= N) || (mu >= N)) {
    ARKProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKBAND", "ARKBand", MSGD_BAD_SIZES);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  smu = MIN(N-1, mu + ml);

  /* Allocate memory for M, savedJ, and pivot arrays */
  M = NULL;
  M = NewBandMat(N, mu, ml, smu);
  if (M == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKBand", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  savedJ = NULL;
  savedJ = NewBandMat(N, mu, ml, mu);
  if (savedJ == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKBand", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  lpivots = NULL;
  lpivots = NewLintArray(N);
  if (lpivots == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKBand", MSGD_MEM_FAIL);
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
 * arkBandInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the band
 * linear solver.
 * -----------------------------------------------------------------
 */

static int arkBandInit(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;

  arkdls_mem = (ARKDlsMem) lmem;

  nje   = 0;
  nfeDQ  = 0;
  nstlj = 0;

  /* Set Jacobian function and data, depending on jacDQ */
  if (jacDQ) {
    jac = arkDlsBandDQJac;
    J_data = ark_mem;
  } else {
    J_data = ark_mem->ark_user_data;
  }

  last_flag = ARKDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * arkBandSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy.  In any case, it constructs the Newton matrix 
 * M = I - gamma*J, updates counters, and calls the band LU 
 * factorization routine.
 * -----------------------------------------------------------------
 */

static int arkBandSetup(ARKodeMem ark_mem, int convfail, N_Vector ypred,
                       N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1,
                       N_Vector vtemp2, N_Vector vtemp3)
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
    BandCopy(savedJ, M, mu, ml);

  } else {

    /* If jok = FALSE, call jac routine for new J value */
    nje++;
    nstlj = nst;
    *jcurPtr = TRUE;
    SetToZero(M); 

    retval = jac(n, mu, ml, tn, ypred, fpred, M, J_data, vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      ARKProcessError(ark_mem, ARKDLS_JACFUNC_UNRECVR, "ARKBAND", "arkBandSetup", MSGD_JACFUNC_FAILED);
      last_flag = ARKDLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      last_flag = ARKDLS_JACFUNC_RECVR;
      return(1);
    }

    BandCopy(M, savedJ, mu, ml);

  }
  
  /* Scale and add I to get M = I - gamma*J */
  BandScale(-gamma, M);
  AddIdentity(M);

  /* Do LU factorization of M */
  ier = BandGBTRF(M, lpivots);

  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) {
    last_flag = ier;
    return(1);
  }
  last_flag = ARKDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * arkBandSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the band linear solver
 * by calling the band backsolve routine.  The return value is 0.
 * -----------------------------------------------------------------
 */

static int arkBandSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur)
{
  ARKDlsMem arkdls_mem;
  realtype *bd;

  arkdls_mem = (ARKDlsMem) lmem;

  bd = N_VGetArrayPointer(b);

  BandGBTRS(M, lpivots, bd);

  /* If ARK_BDF, scale the correction to account for change in gamma */
  if ((lmm == ARK_BDF) && (gamrat != ONE)) {
    N_VScale(TWO/(ONE + gamrat), b, b);
  }

  last_flag = ARKDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * arkBandFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the band linear solver.
 * -----------------------------------------------------------------
 */

static void arkBandFree(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;

  arkdls_mem = (ARKDlsMem) lmem;

  DestroyMat(M);
  DestroyMat(savedJ);
  DestroyArray(lpivots);
  free(arkdls_mem);
  ark_mem->ark_lmem = NULL;
}

