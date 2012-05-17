/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * This is the implementation file for a ARKODE dense linear solver
 * using BLAS and LAPACK functions.
 * -----------------------------------------------------------------
 ***** UNTOUCHED *****
 * -----------------------------------------------------------------
 */

/*
 * NOTE: the only operation that does not use Blas/Lapack functions
 *       is matrix plus mass matrix (in calculating M-gamma*J in lsetup)
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <arkode/arkode_lapack.h>
#include "arkode_direct_impl.h"
#include "arkode_impl.h"

#include <sundials/sundials_math.h>

/* Constant */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* ARKLAPACK DENSE linit, lsetup, lsolve, and lfree routines */ 
static int arkLapackDenseInit(ARKodeMem ark_mem);
static int arkLapackDenseSetup(ARKodeMem ark_mem, int convfail, 
                              N_Vector yP, N_Vector fctP, 
                              booleantype *jcurPtr,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int arkLapackDenseSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector fctC);
static void arkLapackDenseFree(ARKodeMem ark_mem);

/* ARKLAPACK BAND linit, lsetup, lsolve, and lfree routines */ 
static int arkLapackBandInit(ARKodeMem ark_mem);
static int arkLapackBandSetup(ARKodeMem ark_mem, int convfail, 
                             N_Vector yP, N_Vector fctP, 
                             booleantype *jcurPtr,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int arkLapackBandSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                             N_Vector yC, N_Vector fctC);
static void arkLapackBandFree(ARKodeMem ark_mem);

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define lmm            (ark_mem->ark_lmm)
#define f              (ark_mem->ark_f)
#define nst            (ark_mem->ark_nst)
#define tn             (ark_mem->ark_tn)
#define h              (ark_mem->ark_h)
#define gamma          (ark_mem->ark_gamma)
#define gammap         (ark_mem->ark_gammap)
#define gamrat         (ark_mem->ark_gamrat)
#define ewt            (ark_mem->ark_ewt)

#define linit          (ark_mem->ark_linit)
#define lsetup         (ark_mem->ark_lsetup)
#define lsolve         (ark_mem->ark_lsolve)
#define lfree          (ark_mem->ark_lfree)
#define lmem           (ark_mem->ark_lmem)
#define tempv          (ark_mem->ark_tempv)
#define setupNonNull   (ark_mem->ark_setupNonNull)

#define mtype          (arkdls_mem->d_type)
#define n              (arkdls_mem->d_n)
#define ml             (arkdls_mem->d_ml)
#define mu             (arkdls_mem->d_mu)
#define smu            (arkdls_mem->d_smu)
#define jacDQ          (arkdls_mem->d_jacDQ)
#define djac           (arkdls_mem->d_djac)
#define bjac           (arkdls_mem->d_bjac)
#define M              (arkdls_mem->d_M)
#define savedJ         (arkdls_mem->d_savedJ)
#define pivots         (arkdls_mem->d_pivots)
#define nstlj          (arkdls_mem->d_nstlj)
#define nje            (arkdls_mem->d_nje)
#define nfeDQ          (arkdls_mem->d_nfeDQ)
#define J_data         (arkdls_mem->d_J_data)
#define last_flag      (arkdls_mem->d_last_flag)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * -----------------------------------------------------------------
 * ARKLapackDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the linear solver module.  ARKLapackDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the ark_linit, ark_lsetup, ark_lsolve, ark_lfree fields in (*arkode_mem)
 * to be arkLapackDenseInit, arkLapackDenseSetup, arkLapackDenseSolve, 
 * and arkLapackDenseFree, respectively.  It allocates memory for a 
 * structure of type ARKDlsMemRec and sets the ark_lmem field in 
 * (*arkode_mem) to the address of this structure.  It sets setupNonNull 
 * in (*arkode_mem) to TRUE, and the d_jac field to the default 
 * arkDlsDenseDQJac. Finally, it allocates memory for M, pivots, and 
 * savedJ.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, ARKLapackDense will first 
 *       test for a compatible N_Vector internal representation 
 *       by checking that N_VGetArrayPointer and N_VSetArrayPointer 
 *       exist.
 * -----------------------------------------------------------------
 */
int ARKLapackDense(void *arkode_mem, int N)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKDLS_MEM_NULL, "ARKLAPACK", "ARKLapackDense", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the LAPACK solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    ARKProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKLAPACK", "ARKLapackDense", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(ark_mem);

  /* Set four main function fields in ark_mem */
  linit  = arkLapackDenseInit;
  lsetup = arkLapackDenseSetup;
  lsolve = arkLapackDenseSolve;
  lfree  = arkLapackDenseFree;

  /* Get memory for ARKDlsMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMem) malloc(sizeof(struct ARKDlsMemRec));
  if (arkdls_mem == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", "ARKLapackDense", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;

  /* Initialize Jacobian-related data */
  jacDQ  = TRUE;
  djac   = NULL;
  J_data = NULL;

  last_flag = ARKDLS_SUCCESS;
  setupNonNull = TRUE;

  /* Set problem dimension */
  n = (long int) N;

  /* Allocate memory for M, pivot array, and savedJ */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = NewDenseMat(n, n);
  if (M == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", "ARKLapackDense", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", "ARKLapackDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  savedJ = NewDenseMat(n, n);
  if (savedJ == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", "ARKLapackDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    DestroyArray(pivots);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = arkdls_mem;

  return(ARKDLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * ARKLapackBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module. It first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * ark_linit, ark_lsetup, ark_lsolve, and ark_lfree fields in (*arkode_mem)
 * to be arkLapackBandInit, arkLapackBandSetup, arkLapackBandSolve, 
 * and arkLapackBandFree, respectively.  It allocates memory for a 
 * structure of type ARKLapackBandMemRec and sets the ark_lmem field in 
 * (*arkode_mem) to the address of this structure.  It sets setupNonNull 
 * in (*arkode_mem) to be TRUE, mu to be mupper, ml to be mlower, and 
 * the jacE and jacI field to NULL.
 * Finally, it allocates memory for M, pivots, and savedJ.  
 * The ARKLapackBand return value is ARKDLS_SUCCESS = 0, 
 * ARKDLS_MEM_FAIL = -1, or ARKDLS_ILL_INPUT = -2.
 *
 * NOTE: The ARKLAPACK linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, ARKLapackBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */                  
int ARKLapackBand(void *arkode_mem, int N, int mupper, int mlower)
{
  ARKodeMem ark_mem;
  ARKDlsMem arkdls_mem;

  /* Return immediately if arkode_mem is NULL */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKDLS_MEM_NULL, "ARKLAPACK", "ARKLapackBand", MSGD_ARKMEM_NULL);
    return(ARKDLS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (tempv->ops->nvgetarraypointer == NULL) {
    ARKProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKLAPACK", "ARKLapackBand", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  if (lfree != NULL) lfree(ark_mem);

  /* Set four main function fields in ark_mem */  
  linit  = arkLapackBandInit;
  lsetup = arkLapackBandSetup;
  lsolve = arkLapackBandSolve;
  lfree  = arkLapackBandFree;
  
  /* Get memory for ARKDlsMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMem) malloc(sizeof(struct ARKDlsMemRec));
  if (arkdls_mem == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", "ARKLapackBand", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_BAND;

  /* Initialize Jacobian-related data */
  jacDQ  = TRUE;
  bjac   = NULL;
  J_data = NULL;

  last_flag = ARKDLS_SUCCESS;
  setupNonNull = TRUE;
  
  /* Load problem dimension */
  n = (long int) N;

  /* Load half-bandwiths in arkdls_mem */
  ml = (long int) mlower;
  mu = (long int) mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= n) || (mu >= n)) {
    ARKProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKLAPACK", "ARKLapackBand", MSGD_BAD_SIZES);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  smu = MIN(n-1, mu + ml);

  /* Allocate memory for M, pivot array, and savedJ */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = NewBandMat(n, mu, ml, smu);
  if (M == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", "ARKLapackBand", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }  
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", "ARKLapackBand", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  savedJ = NewBandMat(n, mu, ml, smu);
  if (savedJ == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKLAPACK", "ARKLapackBand", MSGD_MEM_FAIL);
    DestroyMat(M);
    DestroyArray(pivots);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = arkdls_mem;

  return(ARKDLS_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH DENSE JACOBIANS
 * =================================================================
 */

/*
 * arkLapackDenseInit does remaining initializations specific to the dense
 * linear solver.
 */
static int arkLapackDenseInit(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;

  arkdls_mem = (ARKDlsMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;

  /* Set Jacobian function and data, depending on jacDQ */
  if (jacDQ) {
    djac = arkDlsDenseDQJac;
    J_data = ark_mem;
  } else {
    J_data = ark_mem->ark_user_data;
  }

  last_flag = ARKDLS_SUCCESS;
  return(0);
}

/*
 * arkLapackDenseSetup does the setup operations for the dense linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy. In any case, it constructs the Newton matrix M = I - gamma*J
 * updates counters, and calls the dense LU factorization routine.
 */
static int arkLapackDenseSetup(ARKodeMem ark_mem, int convfail,
                              N_Vector yP, N_Vector fctP,
                              booleantype *jcurPtr,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  ARKDlsMem arkdls_mem;
  realtype dgamma, fact;
  booleantype jbad, jok;
  int ier, retval, one = 1;
  int intn, lenmat;

  arkdls_mem = (ARKDlsMem) lmem;
  intn = (int) n;
  lenmat = M->ldata ;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlj + ARKD_MSBJ) ||
    ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKD_DGMAX)) ||
    (convfail == ARK_FAIL_OTHER);
  jok = !jbad;
  
  if (jok) {
    
    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    dcopy_f77(&lenmat, savedJ->data, &one, M->data, &one);
    
  } else {
    
    /* If jok = FALSE, call jac routine for new J value */
    nje++;
    nstlj = nst;
    *jcurPtr = TRUE;
    SetToZero(M);

    retval = djac(n, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);

    if (retval == 0) {
      dcopy_f77(&lenmat, M->data, &one, savedJ->data, &one);
    } else if (retval < 0) {
      ARKProcessError(ark_mem, ARKDLS_JACFUNC_UNRECVR, "ARKLAPACK", "arkLapackDenseSetup", MSGD_JACFUNC_FAILED);
      last_flag = ARKDLS_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      last_flag = ARKDLS_JACFUNC_RECVR;
      return(1);
    }
    
  }

  /* Scale J by - gamma */
  fact = -gamma;
  dscal_f77(&lenmat, &fact, M->data, &one);
  
  /* Add identity to get M = I - gamma*J*/
  AddIdentity(M);

  /* Do LU factorization of M */
  dgetrf_f77(&intn, &intn, M->data, &intn, pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = (long int) ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * arkLapackDenseSolve handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.
 */
static int arkLapackDenseSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector fctC)
{
  ARKDlsMem arkdls_mem;
  realtype *bd, fact;
  int ier, one = 1;
  int intn;

  arkdls_mem = (ARKDlsMem) lmem;
  intn = (int) n;

  bd = N_VGetArrayPointer(b);

  dgetrs_f77("N", &intn, &one, M->data, &intn, pivots, bd, &intn, &ier, 1); 

  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm == ARK_BDF) && (gamrat != ONE)) {
    fact = TWO/(ONE + gamrat);
    dscal_f77(&intn, &fact, bd, &one); 
  }
  
  last_flag = ARKDLS_SUCCESS;
  return(0);
}

/*
 * arkLapackDenseFree frees memory specific to the dense linear solver.
 */
static void arkLapackDenseFree(ARKodeMem ark_mem)
{
  ARKDlsMem  arkdls_mem;

  arkdls_mem = (ARKDlsMem) lmem;
  
  DestroyMat(M);
  DestroyArray(pivots);
  DestroyMat(savedJ);
  free(arkdls_mem); 
  arkdls_mem = NULL;
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH BAND JACOBIANS
 * =================================================================
 */

/*
 * arkLapackBandInit does remaining initializations specific to the band
 * linear solver.
 */
static int arkLapackBandInit(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;

  arkdls_mem = (ARKDlsMem) lmem;

  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;

  /* Set Jacobian function and data, depending on jacDQ */
  if (jacDQ) {
    bjac = arkDlsBandDQJac;
    J_data = ark_mem;
  } else {
    J_data = ark_mem->ark_user_data;
  }
  
  last_flag = ARKDLS_SUCCESS;
  return(0);
}

/*
 * arkLapackBandSetup does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy. In any case, it constructs the Newton matrix M = I - gamma*J, 
 * updates counters, and calls the band LU factorization routine.
 */
static int arkLapackBandSetup(ARKodeMem ark_mem, int convfail, 
                             N_Vector yP, N_Vector fctP, 
                             booleantype *jcurPtr,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  ARKDlsMem arkdls_mem;
  realtype dgamma, fact;
  booleantype jbad, jok;
  int ier, retval, one = 1;
  int intn, iml, imu, lenmat, ldmat;

  arkdls_mem = (ARKDlsMem) lmem;
  intn = (int) n;
  iml = (int) ml;
  imu = (int) mu;
  lenmat = M->ldata;
  ldmat = M->ldim;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlj + ARKD_MSBJ) ||
    ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKD_DGMAX)) ||
    (convfail == ARK_FAIL_OTHER);
  jok = !jbad;
  
  if (jok) {
    
    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    dcopy_f77(&lenmat, savedJ->data, &one, M->data, &one);
    
  } else {
    
    /* If jok = FALSE, call jac routine for new J value */
    nje++;
    nstlj = nst;
    *jcurPtr = TRUE;
    SetToZero(M);

    retval = bjac(n, mu, ml, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
    if (retval == 0) {
      dcopy_f77(&lenmat, M->data, &one, savedJ->data, &one);
    } else if (retval < 0) {
      ARKProcessError(ark_mem, ARKDLS_JACFUNC_UNRECVR, "ARKLAPACK", "arkLapackBandSetup", MSGD_JACFUNC_FAILED);
      last_flag = ARKDLS_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      last_flag = ARKDLS_JACFUNC_RECVR;
      return(1);
    }
    
  }
  
  /* Scale J by - gamma */
  fact = -gamma;
  dscal_f77(&lenmat, &fact, M->data, &one);
  
  /* Add identity to get M = I - gamma*J*/
  AddIdentity(M);
  
  /* Do LU factorization of M */
  dgbtrf_f77(&intn, &intn, &iml, &imu, M->data, &ldmat, pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = (long int) ier;
  if (ier > 0) return(1);
  return(0);

}

/*
 * arkLapackBandSolve handles the solve operation for the band linear solver
 * by calling the band backsolve routine.
 */
static int arkLapackBandSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                             N_Vector yC, N_Vector fctC)
{
  ARKDlsMem arkdls_mem;
  realtype *bd, fact;
  int ier, one = 1;
  int intn, iml, imu, ldmat;

  arkdls_mem = (ARKDlsMem) lmem;
  intn = (int) n;
  iml = (int) ml;
  imu = (int) mu;
  ldmat = M->ldim;

  bd = N_VGetArrayPointer(b);

  dgbtrs_f77("N", &intn, &iml, &imu, &one, M->data, &ldmat, pivots, bd, &intn, &ier, 1);
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm == ARK_BDF) && (gamrat != ONE)) {
    fact = TWO/(ONE + gamrat);
    dscal_f77(&intn, &fact, bd, &one); 
  }

  last_flag = ARKDLS_SUCCESS;
  return(0);
}

/*
 * arkLapackBandFree frees memory specific to the band linear solver.
 */
static void arkLapackBandFree(ARKodeMem ark_mem)
{
  ARKDlsMem  arkdls_mem;

  arkdls_mem = (ARKDlsMem) lmem;
  
  DestroyMat(M);
  DestroyArray(pivots);
  DestroyMat(savedJ);
  free(arkdls_mem); 
  arkdls_mem = NULL;
}

