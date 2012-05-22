/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
----------------------------------------------------------------- 
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 This is the implementation file for the ARKBAND linear solver.
---------------------------------------------------------------*/

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


/*---------------------------------------------------------------
 ARKBand:

 This routine initializes the memory record and sets various function
 fields specific to the band linear solver module.  ARKBand first calls
 the existing lfree routine if this is not NULL.  It then sets the
 ark_linit, ark_lsetup, ark_lsolve, and ark_lfree fields in (*arkode_mem)
 to be arkBandInit, arkBandSetup, arkBandSolve, and arkBandFree,
 respectively.  It allocates memory for a structure of type
 ARKDlsMemRec and sets the ark_lmem field in (*arkode_mem) to the
 address of this structure.  It sets setupNonNull in (*arkode_mem) to be
 TRUE, d_mu to be mupper, d_ml to be mlower, and the d_jac field to be 
 arkDlsBandDQJac.
 Finally, it allocates memory for M, savedJ, and pivot.  The ARKBand
 return value is SUCCESS = 0, LMEM_FAIL = -1, or LIN_ILL_INPUT = -2.

 NOTE: The band linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, ARKBand will first 
       test for compatible a compatible N_Vector internal
       representation by checking that the function 
       N_VGetArrayPointer exists.
---------------------------------------------------------------*/
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
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL) {
    ARKProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKBAND", "ARKBand", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  if (ark_mem->ark_lfree != NULL) ark_mem->ark_lfree(ark_mem);

  /* Set four main function fields in ark_mem */  
  ark_mem->ark_linit  = arkBandInit;
  ark_mem->ark_lsetup = arkBandSetup;
  ark_mem->ark_lsolve = arkBandSolve;
  ark_mem->ark_lfree  = arkBandFree;
  
  /* Get memory for ARKDlsMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMem) malloc(sizeof(struct ARKDlsMemRec));
  if (arkdls_mem == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKBand", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* Set matrix type */
  arkdls_mem->d_type = SUNDIALS_BAND;
  
  /* Initialize Jacobian-related data */
  arkdls_mem->d_jacDQ = TRUE;
  arkdls_mem->d_bjac = NULL;
  arkdls_mem->d_J_data = NULL;

  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;

  ark_mem->ark_setupNonNull = TRUE;
  
  /* Load problem dimension */
  arkdls_mem->d_n = N;

  /* Load half-bandwiths in arkdls_mem */
  arkdls_mem->d_ml = mlower;
  arkdls_mem->d_mu = mupper;

  /* Test ml and mu for legality */
  if ((arkdls_mem->d_ml < 0) || (arkdls_mem->d_mu < 0) || (arkdls_mem->d_ml >= N) || (arkdls_mem->d_mu >= N)) {
    ARKProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKBAND", "ARKBand", MSGD_BAD_SIZES);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  arkdls_mem->d_smu = MIN(N-1, arkdls_mem->d_mu + arkdls_mem->d_ml);

  /* Allocate memory for M, savedJ, and pivot arrays */
  arkdls_mem->d_M = NULL;
  arkdls_mem->d_M = NewBandMat(N, arkdls_mem->d_mu, arkdls_mem->d_ml, arkdls_mem->d_smu);
  if (arkdls_mem->d_M == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKBand", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  arkdls_mem->d_savedJ = NULL;
  arkdls_mem->d_savedJ = NewBandMat(N, arkdls_mem->d_mu, arkdls_mem->d_ml, arkdls_mem->d_mu);
  if (arkdls_mem->d_savedJ == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKBand", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  arkdls_mem->d_lpivots = NULL;
  arkdls_mem->d_lpivots = NewLintArray(N);
  if (arkdls_mem->d_lpivots == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKBAND", "ARKBand", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    DestroyMat(arkdls_mem->d_savedJ);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  ark_mem->ark_lmem = arkdls_mem;

  return(ARKDLS_SUCCESS);
}


/*---------------------------------------------------------------
 arkBandInit:

 This routine does remaining initializations specific to the band
 linear solver.
---------------------------------------------------------------*/
static int arkBandInit(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  arkdls_mem->d_nje   = 0;
  arkdls_mem->d_nfeDQ  = 0;
  arkdls_mem->d_nstlj = 0;

  /* Set Jacobian function and data, depending on jacDQ */
  if (arkdls_mem->d_jacDQ) {
    arkdls_mem->d_bjac = arkDlsBandDQJac;
    arkdls_mem->d_J_data = ark_mem;
  } else {
    arkdls_mem->d_J_data = ark_mem->ark_user_data;
  }

  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkBandSetup:

 This routine does the setup operations for the band linear 
 solver. It makes a decision whether or not to call the Jacobian 
 evaluation routine based on various state variables, and if not 
 it uses the saved copy.  In any case, it constructs the Newton 
 matrix  M = I - gamma*J, updates counters, and calls the band 
 LU factorization routine.
---------------------------------------------------------------*/
static int arkBandSetup(ARKodeMem ark_mem, int convfail, 
			N_Vector ypred, N_Vector fpred, 
			booleantype *jcurPtr, N_Vector vtemp1,
			N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  long int ier;
  ARKDlsMem arkdls_mem;
  int retval;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((ark_mem->ark_gamma/ark_mem->ark_gammap) - ONE);
  jbad = (ark_mem->ark_nst == 0) || (ark_mem->ark_nst > arkdls_mem->d_nstlj + ARKD_MSBJ) ||
         ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKD_DGMAX)) ||
         (convfail == ARK_FAIL_OTHER);
  jok = !jbad;
  
  if (jok) {

    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    BandCopy(arkdls_mem->d_savedJ, arkdls_mem->d_M, arkdls_mem->d_mu, arkdls_mem->d_ml);

  } else {

    /* If jok = FALSE, call jac routine for new J value */
    arkdls_mem->d_nje++;
    arkdls_mem->d_nstlj = ark_mem->ark_nst;
    *jcurPtr = TRUE;
    SetToZero(arkdls_mem->d_M); 

    retval = arkdls_mem->d_bjac(arkdls_mem->d_n, arkdls_mem->d_mu, arkdls_mem->d_ml, ark_mem->ark_tn, ypred, fpred, arkdls_mem->d_M, arkdls_mem->d_J_data, vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      ARKProcessError(ark_mem, ARKDLS_JACFUNC_UNRECVR, "ARKBAND", "arkBandSetup", MSGD_JACFUNC_FAILED);
      arkdls_mem->d_last_flag = ARKDLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      arkdls_mem->d_last_flag = ARKDLS_JACFUNC_RECVR;
      return(1);
    }

    BandCopy(arkdls_mem->d_M, arkdls_mem->d_savedJ, arkdls_mem->d_mu, arkdls_mem->d_ml);

  }
  
  /* Scale and add I to get M = I - gamma*J */
  BandScale(-ark_mem->ark_gamma, arkdls_mem->d_M);
  AddIdentity(arkdls_mem->d_M);

  /* Do LU factorization of M */
  ier = BandGBTRF(arkdls_mem->d_M, arkdls_mem->d_lpivots);

  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) {
    arkdls_mem->d_last_flag = ier;
    return(1);
  }
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkBandSolve:

 This routine handles the solve operation for the band linear solver
 by calling the band backsolve routine.  The return value is 0.
---------------------------------------------------------------*/
static int arkBandSolve(ARKodeMem ark_mem, N_Vector b, N_Vector weight,
                       N_Vector ycur, N_Vector fcur)
{
  ARKDlsMem arkdls_mem;
  realtype *bd;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  bd = N_VGetArrayPointer(b);

  BandGBTRS(arkdls_mem->d_M, arkdls_mem->d_lpivots, bd);

  /* scale the correction to account for change in gamma */
  if (ark_mem->ark_gamrat != ONE) 
    N_VScale(TWO/(ONE + ark_mem->ark_gamrat), b, b);

  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkBandFree:

 This routine frees memory specific to the band linear solver.
---------------------------------------------------------------*/
static void arkBandFree(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;

  DestroyMat(arkdls_mem->d_M);
  DestroyMat(arkdls_mem->d_savedJ);
  DestroyArray(arkdls_mem->d_lpivots);
  free(arkdls_mem);
  ark_mem->ark_lmem = NULL;
}


/*---------------------------------------------------------------
    EOF
---------------------------------------------------------------*/
