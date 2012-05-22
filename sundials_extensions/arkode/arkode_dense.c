/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
----------------------------------------------------------------- 
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 This is the impleentation file for the ARKDENSE linear solver.
---------------------------------------------------------------*/

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

                  
/*---------------------------------------------------------------
 ARKDense:

 This routine initializes the memory record and sets various 
 function fields specific to the dense linear solver module.  
 ARKDense first calls the existing lfree routine if this is not 
 NULL.  Then it sets the ark_linit, ark_lsetup, ark_lsolve, 
 ark_lfree fields in (*arkode_mem) to be arkDenseInit, 
 arkDenseSetup, arkDenseSolve, and arkDenseFree, respectively.  
 It allocates memory for a structure of type ARKDlsMemRec and 
 sets the ark_lmem field in (*arkode_mem) to the address of this 
 structure.  It sets setupNonNull in (*arkode_mem) to TRUE, and 
 the d_jac field to the default arkDlsDenseDQJac.  Finally, it
 allocates memory for M, savedJ, and lpivots. The return value 
 is SUCCESS = 0, or LMEM_FAIL = -1.

 NOTE: The dense linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, ARKDense will first 
       test for compatible a compatible N_Vector internal
       representation by checking that N_VGetArrayPointer and
       N_VSetArrayPointer exist.
---------------------------------------------------------------*/
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
  if (ark_mem->ark_tempv->ops->nvgetarraypointer == NULL ||
      ark_mem->ark_tempv->ops->nvsetarraypointer == NULL) {
    ARKProcessError(ark_mem, ARKDLS_ILL_INPUT, "ARKDENSE", "ARKDense", MSGD_BAD_NVECTOR);
    return(ARKDLS_ILL_INPUT);
  }

  if (ark_mem->ark_lfree !=NULL) ark_mem->ark_lfree(ark_mem);

  /* Set four main function fields in ark_mem */
  ark_mem->ark_linit  = arkDenseInit;
  ark_mem->ark_lsetup = arkDenseSetup;
  ark_mem->ark_lsolve = arkDenseSolve;
  ark_mem->ark_lfree  = arkDenseFree;

  /* Get memory for ARKDlsMemRec */
  arkdls_mem = NULL;
  arkdls_mem = (ARKDlsMem) malloc(sizeof(struct ARKDlsMemRec));
  if (arkdls_mem == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDENSE", "ARKDense", MSGD_MEM_FAIL);
    return(ARKDLS_MEM_FAIL);
  }

  /* Set matrix type */
  arkdls_mem->d_type = SUNDIALS_DENSE;

  /* Initialize Jacobian-related data */
  arkdls_mem->d_jacDQ = TRUE;
  arkdls_mem->d_djac = NULL;
  arkdls_mem->d_J_data = NULL;

  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;

  ark_mem->ark_setupNonNull = TRUE;

  /* Set problem dimension */
  arkdls_mem->d_n = N;

  /* Allocate memory for M, savedJ, and pivot array */

  arkdls_mem->d_M = NULL;
  arkdls_mem->d_M = NewDenseMat(N, N);
  if (arkdls_mem->d_M == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDENSE", "ARKDense", MSGD_MEM_FAIL);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  arkdls_mem->d_savedJ = NULL;
  arkdls_mem->d_savedJ = NewDenseMat(N, N);
  if (arkdls_mem->d_savedJ == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDENSE", "ARKDense", MSGD_MEM_FAIL);
    DestroyMat(arkdls_mem->d_M);
    free(arkdls_mem); arkdls_mem = NULL;
    return(ARKDLS_MEM_FAIL);
  }
  arkdls_mem->d_lpivots = NULL;
  arkdls_mem->d_lpivots = NewLintArray(N);
  if (arkdls_mem->d_lpivots == NULL) {
    ARKProcessError(ark_mem, ARKDLS_MEM_FAIL, "ARKDENSE", "ARKDense", MSGD_MEM_FAIL);
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
 arkDenseInit:

 This routine does remaining initializations specific to the 
 dense linear solver.
---------------------------------------------------------------*/
static int arkDenseInit(ARKodeMem ark_mem)
{
  ARKDlsMem arkdls_mem;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;
  
  arkdls_mem->d_nje   = 0;
  arkdls_mem->d_nfeDQ = 0;
  arkdls_mem->d_nstlj = 0;

  /* Set Jacobian function and data, depending on jacDQ */
  if (arkdls_mem->d_jacDQ) {
    arkdls_mem->d_djac = arkDlsDenseDQJac;
    arkdls_mem->d_J_data = ark_mem;
  } else {
    arkdls_mem->d_J_data = ark_mem->ark_user_data;
  }

  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkDenseSetup:

 This routine does the setup operations for the dense linear solver.
 It makes a decision whether or not to call the Jacobian evaluation
 routine based on various state variables, and if not it uses the 
 saved copy.  In any case, it constructs the Newton matrix 
 M = I - gamma*J, updates counters, and calls the dense LU 
 factorization routine.
---------------------------------------------------------------*/
static int arkDenseSetup(ARKodeMem ark_mem, int convfail, 
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
  jbad = (ark_mem->ark_nst == 0) || 
    (ark_mem->ark_nst > arkdls_mem->d_nstlj + ARKD_MSBJ) ||
    ((convfail == ARK_FAIL_BAD_J) && (dgamma < ARKD_DGMAX)) ||
    (convfail == ARK_FAIL_OTHER);
  jok = !jbad;
 
  if (jok) {

    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    DenseCopy(arkdls_mem->d_savedJ, arkdls_mem->d_M);

  } else {

    /* If jok = FALSE, call jac routine for new J value */
    arkdls_mem->d_nje++;
    arkdls_mem->d_nstlj = ark_mem->ark_nst;
    *jcurPtr = TRUE;
    SetToZero(arkdls_mem->d_M);

    retval = arkdls_mem->d_djac(arkdls_mem->d_n, ark_mem->ark_tn, ypred, 
				fpred, arkdls_mem->d_M, arkdls_mem->d_J_data, 
				vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      ARKProcessError(ark_mem, ARKDLS_JACFUNC_UNRECVR, "ARKDENSE", 
		      "arkDenseSetup",  MSGD_JACFUNC_FAILED);
      arkdls_mem->d_last_flag = ARKDLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      arkdls_mem->d_last_flag = ARKDLS_JACFUNC_RECVR;
      return(1);
    }

    DenseCopy(arkdls_mem->d_M, arkdls_mem->d_savedJ);

  }
  
  /* Scale and add I to get M = I - gamma*J */
  DenseScale(-ark_mem->ark_gamma, arkdls_mem->d_M);
  AddIdentity(arkdls_mem->d_M);

  /* Do LU factorization of M */
  ier = DenseGETRF(arkdls_mem->d_M, arkdls_mem->d_lpivots); 

  /* Return 0 if the LU was complete; otherwise return 1 */
  arkdls_mem->d_last_flag = ier;
  if (ier > 0) return(1);
  return(0);
}


/*---------------------------------------------------------------
 arkDenseSolve:

 This routine handles the solve operation for the dense linear 
 solver by calling the dense backsolve routine.  The returned 
 value is 0.
---------------------------------------------------------------*/
static int arkDenseSolve(ARKodeMem ark_mem, N_Vector b, 
			 N_Vector weight, N_Vector ycur, 
			 N_Vector fcur)
{
  ARKDlsMem arkdls_mem;
  realtype *bd;

  arkdls_mem = (ARKDlsMem) ark_mem->ark_lmem;
  
  bd = N_VGetArrayPointer(b);

  DenseGETRS(arkdls_mem->d_M, arkdls_mem->d_lpivots, bd);

  /* If ARK_BDF, scale the correction to account for change in gamma */
  if ((ark_mem->ark_lmm == ARK_BDF) && (ark_mem->ark_gamrat != ONE)) {
    N_VScale(TWO/(ONE + ark_mem->ark_gamrat), b, b);
  }
  
  arkdls_mem->d_last_flag = ARKDLS_SUCCESS;
  return(0);
}


/*---------------------------------------------------------------
 arkDenseFree:

 This routine frees memory specific to the dense linear solver.
---------------------------------------------------------------*/
static void arkDenseFree(ARKodeMem ark_mem)
{
  ARKDlsMem  arkdls_mem;

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
