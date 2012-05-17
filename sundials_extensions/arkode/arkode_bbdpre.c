/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * This file contains implementations of routines for a
 * band-block-diagonal preconditioner, i.e. a block-diagonal
 * matrix with banded blocks, for use with ARKODE, a ARKSPILS linear
 * solver, and the parallel implementation of NVECTOR.
 * -----------------------------------------------------------------
 ***** UNTOUCHED *****
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_bbdpre_impl.h"
#include "arkode_spils_impl.h"

#include <arkode/arkode_sptfqmr.h>
#include <arkode/arkode_spbcgs.h>
#include <arkode/arkode_spgmr.h>

#include <sundials/sundials_math.h>

#define MIN_INC_MULT RCONST(1000.0)

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)

/* Prototypes of functions ARKBBDPrecSetup and ARKBBDPrecSolve */

static int ARKBBDPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                          booleantype jok, booleantype *jcurPtr, 
                          realtype gamma, void *bbd_data, 
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int ARKBBDPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                          N_Vector r, N_Vector z, 
                          realtype gamma, realtype delta,
                          int lr, void *bbd_data, N_Vector tmp);

/* Prototype for ARKBBDPrecFree */
static void ARKBBDPrecFree(ARKodeMem ark_mem);


/* Prototype for difference quotient Jacobian calculation routine */

static int ARKBBDDQJac(ARKBBDPrecData pdata, realtype t, 
                      N_Vector y, N_Vector gy, 
                      N_Vector ytemp, N_Vector gtemp);

/* Redability replacements */

#define uround   (ark_mem->ark_uround)
#define vec_tmpl (ark_mem->ark_tempv)

/*
 * -----------------------------------------------------------------
 * User-Callable Functions: initialization, reinit and free
 * -----------------------------------------------------------------
 */

int ARKBBDPrecInit(void *arkode_mem, long int Nlocal, 
                   long int mudq, long int mldq,
                   long int mukeep, long int mlkeep, 
                   realtype dqrely, 
                   ARKLocalFn gloc, ARKCommFn cfn)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  ARKBBDPrecData pdata;
  long int muk, mlk, storage_mu;
  int flag;

  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKBBDPRE", "ARKBBDPrecInit", MSGBBD_MEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if one of the SPILS linear solvers has been attached */
  if (ark_mem->ark_lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKBBDPRE", "ARKBBDPrecInit", MSGBBD_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Test if the NVECTOR package is compatible with the BLOCK BAND preconditioner */
  if(vec_tmpl->ops->nvgetarraypointer == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_ILL_INPUT, "ARKBBDPRE", "ARKBBDPrecInit", MSGBBD_BAD_NVECTOR);
    return(ARKSPILS_ILL_INPUT);
  }

  /* Allocate data memory */
  pdata = NULL;
  pdata = (ARKBBDPrecData) malloc(sizeof *pdata);  
  if (pdata == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBBDPRE", "ARKBBDPrecInit", MSGBBD_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* Set pointers to gloc and cfn; load half-bandwidths */
  pdata->arkode_mem = arkode_mem;
  pdata->gloc = gloc;
  pdata->cfn = cfn;
  pdata->mudq = MIN(Nlocal-1, MAX(0,mudq));
  pdata->mldq = MIN(Nlocal-1, MAX(0,mldq));
  muk = MIN(Nlocal-1, MAX(0,mukeep));
  mlk = MIN(Nlocal-1, MAX(0,mlkeep));
  pdata->mukeep = muk;
  pdata->mlkeep = mlk;

  /* Allocate memory for saved Jacobian */
  pdata->savedJ = NewBandMat(Nlocal, muk, mlk, muk);
  if (pdata->savedJ == NULL) { 
    free(pdata); pdata = NULL; 
    ARKProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBBDPRE", "ARKBBDPrecInit", MSGBBD_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL); 
  }

  /* Allocate memory for preconditioner matrix */
  storage_mu = MIN(Nlocal-1, muk + mlk);
  pdata->savedP = NULL;
  pdata->savedP = NewBandMat(Nlocal, muk, mlk, storage_mu);
  if (pdata->savedP == NULL) {
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    ARKProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBBDPRE", "ARKBBDPrecInit", MSGBBD_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }
  /* Allocate memory for lpivots */
  pdata->lpivots = NULL;
  pdata->lpivots = NewLintArray(Nlocal);
  if (pdata->lpivots == NULL) {
    DestroyMat(pdata->savedP);
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    ARKProcessError(ark_mem, ARKSPILS_MEM_FAIL, "ARKBBDPRE", "ARKBBDPrecInit", MSGBBD_MEM_FAIL);
    return(ARKSPILS_MEM_FAIL);
  }

  /* Set pdata->dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? dqrely : RSqrt(uround);

  /* Store Nlocal to be used in ARKBBDPrecSetup */
  pdata->n_local = Nlocal;

  /* Set work space sizes and initialize nge */
  pdata->rpwsize = Nlocal*(muk + 2*mlk + storage_mu + 2);
  pdata->ipwsize = Nlocal;
  pdata->nge = 0;

  /* Overwrite the P_data field in the SPILS memory */
  arkspils_mem->s_P_data = pdata;

  /* Attach the pfree function */
  arkspils_mem->s_pfree = ARKBBDPrecFree;

  /* Attach preconditioner solve and setup functions */
  flag = ARKSpilsSetPreconditioner(arkode_mem, ARKBBDPrecSetup, ARKBBDPrecSolve);

  return(flag);
}


int ARKBBDPrecReInit(void *arkode_mem, 
                    long int mudq, long int mldq, 
                    realtype dqrely)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  ARKBBDPrecData pdata;
  long int Nlocal;

  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKBBDPRE", "ARKBBDPrecReInit", MSGBBD_MEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Test if one of the SPILS linear solvers has been attached */
  if (ark_mem->ark_lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKBBDPRE", "ARKBBDPrecReInit", MSGBBD_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  /* Test if the preconditioner data is non-NULL */
  if (arkspils_mem->s_P_data == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_PMEM_NULL, "ARKBBDPRE", "ARKBBDPrecReInit", MSGBBD_PMEM_NULL);
    return(ARKSPILS_PMEM_NULL);
  } 
  pdata = (ARKBBDPrecData) arkspils_mem->s_P_data;

  /* Load half-bandwidths */
  Nlocal = pdata->n_local;
  pdata->mudq = MIN(Nlocal-1, MAX(0,mudq));
  pdata->mldq = MIN(Nlocal-1, MAX(0,mldq));

  /* Set pdata->dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? dqrely : RSqrt(uround);

  /* Re-initialize nge */
  pdata->nge = 0;

  return(ARKSPILS_SUCCESS);
}

int ARKBBDPrecGetWorkSpace(void *arkode_mem, long int *lenrwBBDP, long int *leniwBBDP)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  ARKBBDPrecData pdata;

  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKBBDPRE", "ARKBBDPrecGetWorkSpace", MSGBBD_MEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKBBDPRE", "ARKBBDPrecGetWorkSpace", MSGBBD_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  if (arkspils_mem->s_P_data == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_PMEM_NULL, "ARKBBDPRE", "ARKBBDPrecGetWorkSpace", MSGBBD_PMEM_NULL);
    return(ARKSPILS_PMEM_NULL);
  } 
  pdata = (ARKBBDPrecData) arkspils_mem->s_P_data;

  *lenrwBBDP = pdata->rpwsize;
  *leniwBBDP = pdata->ipwsize;

  return(ARKSPILS_SUCCESS);
}

int ARKBBDPrecGetNumGfnEvals(void *arkode_mem, long int *ngevalsBBDP)
{
  ARKodeMem ark_mem;
  ARKSpilsMem arkspils_mem;
  ARKBBDPrecData pdata;

  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARKSPILS_MEM_NULL, "ARKBBDPRE", "ARKBBDPrecGetNumGfnEvals", MSGBBD_MEM_NULL);
    return(ARKSPILS_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_lmem == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_LMEM_NULL, "ARKBBDPRE", "ARKBBDPrecGetNumGfnEvals", MSGBBD_LMEM_NULL);
    return(ARKSPILS_LMEM_NULL);
  }
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;

  if (arkspils_mem->s_P_data == NULL) {
    ARKProcessError(ark_mem, ARKSPILS_PMEM_NULL, "ARKBBDPRE", "ARKBBDPrecGetNumGfnEvals", MSGBBD_PMEM_NULL);
    return(ARKSPILS_PMEM_NULL);
  } 
  pdata = (ARKBBDPrecData) arkspils_mem->s_P_data;

  *ngevalsBBDP = pdata->nge;

  return(ARKSPILS_SUCCESS);
}

/* Readability Replacements */

#define Nlocal  (pdata->n_local)
#define mudq    (pdata->mudq)
#define mldq    (pdata->mldq)
#define mukeep  (pdata->mukeep)
#define mlkeep  (pdata->mlkeep)
#define dqrely  (pdata->dqrely)
#define gloc    (pdata->gloc)
#define cfn     (pdata->cfn)
#define savedJ  (pdata->savedJ)
#define savedP  (pdata->savedP)
#define lpivots (pdata->lpivots)
#define nge     (pdata->nge)

/*
 * -----------------------------------------------------------------
 * Function : ARKBBDPrecSetup                                      
 * -----------------------------------------------------------------
 * ARKBBDPrecSetup generates and factors a banded block of the
 * preconditioner matrix on each processor, via calls to the
 * user-supplied gloc and cfn functions. It uses difference
 * quotient approximations to the Jacobian elements.
 *
 * ARKBBDPrecSetup calculates a new J,if necessary, then calculates
 * P = I - gamma*J, and does an LU factorization of P.
 *
 * The parameters of ARKBBDPrecSetup used here are as follows:
 *
 * t       is the current value of the independent variable.
 *
 * y       is the current value of the dependent variable vector,
 *         namely the predicted value of y(t).
 *
 * fy      is the vector f(t,y).
 *
 * jok     is an input flag indicating whether Jacobian-related
 *         data needs to be recomputed, as follows:
 *           jok == FALSE means recompute Jacobian-related data
 *                  from scratch.
 *           jok == TRUE  means that Jacobian data from the
 *                  previous ARKBBDPrecon call can be reused
 *                  (with the current value of gamma).
 *         A ARKBBDPrecon call with jok == TRUE should only occur
 *         after a call with jok == FALSE.
 *
 * jcurPtr is a pointer to an output integer flag which is
 *         set by ARKBBDPrecon as follows:
 *           *jcurPtr = TRUE if Jacobian data was recomputed.
 *           *jcurPtr = FALSE if Jacobian data was not recomputed,
 *                      but saved data was reused.
 *
 * gamma   is the scalar appearing in the Newton matrix.
 *
 * bbd_data is a pointer to the preconditioner data set by
 *          ARKBBDPrecInit
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated
 *           for NVectors which are be used by ARKBBDPrecSetup
 *           as temporary storage or work space.
 *
 * Return value:
 * The value returned by this ARKBBDPrecSetup function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried).
 * -----------------------------------------------------------------
 */

static int ARKBBDPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                          booleantype jok, booleantype *jcurPtr, 
                          realtype gamma, void *bbd_data, 
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  long int ier;
  ARKBBDPrecData pdata;
  ARKodeMem ark_mem;
  int retval;

  pdata = (ARKBBDPrecData) bbd_data;

  ark_mem = (ARKodeMem) pdata->arkode_mem;

  if (jok) {

    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    BandCopy(savedJ, savedP, mukeep, mlkeep);

  } else {

    /* Otherwise call ARKBBDDQJac for new J value */
    *jcurPtr = TRUE;
    SetToZero(savedJ);

    retval = ARKBBDDQJac(pdata, t, y, tmp1, tmp2, tmp3);
    if (retval < 0) {
      ARKProcessError(ark_mem, -1, "ARKBBDPRE", "ARKBBDPrecSetup", MSGBBD_FUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    BandCopy(savedJ, savedP, mukeep, mlkeep);

  }
  
  /* Scale and add I to get P = I - gamma*J */
  BandScale(-gamma, savedP);
  AddIdentity(savedP);
 
  /* Do LU factorization of P in place */
  ier = BandGBTRF(savedP, lpivots);
 
  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : ARKBBDPrecSolve
 * -----------------------------------------------------------------
 * ARKBBDPrecSolve solves a linear system P z = r, with the
 * band-block-diagonal preconditioner matrix P generated and
 * factored by ARKBBDPrecSetup.
 *
 * The parameters of ARKBBDPrecSolve used here are as follows:
 *
 * r is the right-hand side vector of the linear system.
 *
 * bbd_data is a pointer to the preconditioner data set by
 *   ARKBBDPrecInit.
 *
 * z is the output vector computed by ARKBBDPrecSolve.
 *
 * The value returned by the ARKBBDPrecSolve function is always 0,
 * indicating success.
 * -----------------------------------------------------------------
 */

static int ARKBBDPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                          N_Vector r, N_Vector z, 
                          realtype gamma, realtype delta,
                          int lr, void *bbd_data, N_Vector tmp)
{
  ARKBBDPrecData pdata;
  realtype *zd;

  pdata = (ARKBBDPrecData) bbd_data;

  /* Copy r to z, then do backsolve and return */
  N_VScale(ONE, r, z);
  
  zd = N_VGetArrayPointer(z);

  BandGBTRS(savedP, lpivots, zd);

  return(0);
}


static void ARKBBDPrecFree(ARKodeMem ark_mem)
{
  ARKSpilsMem arkspils_mem;
  ARKBBDPrecData pdata;
  
  if (ark_mem->ark_lmem == NULL) return;
  arkspils_mem = (ARKSpilsMem) ark_mem->ark_lmem;
  
  if (arkspils_mem->s_P_data == NULL) return;
  pdata = (ARKBBDPrecData) arkspils_mem->s_P_data;

  DestroyMat(savedJ);
  DestroyMat(savedP);
  DestroyArray(lpivots);

  free(pdata);
  pdata = NULL;
}


#define ewt       (ark_mem->ark_ewt)
#define h         (ark_mem->ark_h)
#define user_data (ark_mem->ark_user_data)

/*
 * -----------------------------------------------------------------
 * Function : ARKBBDDQJac
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation
 * to the local block of the Jacobian of g(t,y). It assumes that a
 * band matrix of type DlsMat is stored columnwise, and that elements
 * within each column are contiguous. All matrix elements are generated
 * as difference quotients, by way of calls to the user routine gloc.
 * By virtue of the band structure, the number of these calls is
 * bandwidth + 1, where bandwidth = mldq + mudq + 1.
 * But the band matrix kept has bandwidth = mlkeep + mukeep + 1.
 * This routine also assumes that the local elements of a vector are
 * stored contiguously.
 * -----------------------------------------------------------------
 */

static int ARKBBDDQJac(ARKBBDPrecData pdata, realtype t, 
                      N_Vector y, N_Vector gy, 
                      N_Vector ytemp, N_Vector gtemp)
{
  ARKodeMem ark_mem;
  realtype gnorm, minInc, inc, inc_inv;
  long int group, i, j, width, ngroups, i1, i2;
  realtype *y_data, *ewt_data, *gy_data, *gtemp_data, *ytemp_data, *col_j;
  int retval;

  ark_mem = (ARKodeMem) pdata->arkode_mem;

  /* Load ytemp with y = predicted solution vector */
  N_VScale(ONE, y, ytemp);

  /* Call cfn and gloc to get base value of g(t,y) */
  if (cfn != NULL) {
    retval = cfn(Nlocal, t, y, user_data);
    if (retval != 0) return(retval);
  }

  retval = gloc(Nlocal, t, ytemp, gy, user_data);
  nge++;
  if (retval != 0) return(retval);

  /* Obtain pointers to the data for various vectors */
  y_data     =  N_VGetArrayPointer(y);
  gy_data    =  N_VGetArrayPointer(gy);
  ewt_data   =  N_VGetArrayPointer(ewt);
  ytemp_data =  N_VGetArrayPointer(ytemp);
  gtemp_data =  N_VGetArrayPointer(gtemp);

  /* Set minimum increment based on uround and norm of g */
  gnorm = N_VWrmsNorm(gy, ewt);
  minInc = (gnorm != ZERO) ?
           (MIN_INC_MULT * ABS(h) * uround * Nlocal * gnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mldq + mudq + 1;
  ngroups = MIN(width, Nlocal);

  /* Loop over groups */  
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < Nlocal; j+=width) {
      inc = MAX(dqrely*ABS(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate g with incremented y */
    retval = gloc(Nlocal, t, ytemp, gtemp, user_data);
    nge++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < Nlocal; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(savedJ,j);
      inc = MAX(dqrely*ABS(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = MAX(0, j-mukeep);
      i2 = MIN(j+mlkeep, Nlocal-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) =
          inc_inv * (gtemp_data[i] - gy_data[i]);
    }
  }

  return(0);
}
