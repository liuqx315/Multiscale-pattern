/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
----------------------------------------------------------------- 
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 This is the implementation file for the main ARKODE integrator.
 It is independent of the ARKODE linear solver in use.
---------------------------------------------------------------*/

/*===============================================================
             Import Header Files                                 
===============================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>


/*===============================================================
             Private Helper Functions Prototypes
===============================================================*/
static booleantype ARKCheckNvector(N_Vector tmpl);

static int ARKInitialSetup(ARKodeMem ark_mem);
static int ARKSetButcherTables(ARKodeMem ark_mem);
static int ARKCheckButcherTables(ARKodeMem ark_mem);

static booleantype ARKAllocVectors(ARKodeMem ark_mem, 
				   N_Vector tmpl);
static int ARKAllocRKVectors(ARKodeMem ark_mem);
static void ARKFreeVectors(ARKodeMem ark_mem);

static int ARKEwtSetSS(ARKodeMem ark_mem, N_Vector ycur, 
		       N_Vector weight);
static int ARKEwtSetSV(ARKodeMem ark_mem, N_Vector ycur, 
		       N_Vector weight);

static int ARKHin(ARKodeMem ark_mem, realtype tout);
static realtype ARKUpperBoundH0(ARKodeMem ark_mem, 
				realtype tdist);
static int ARKYddNorm(ARKodeMem ark_mem, realtype hg, 
		      realtype *yddnrm);

static int ARKStep(ARKodeMem ark_mem);
static int ARKStep2(ARKodeMem ark_mem);

static int ARKAdapt(ARKodeMem ark_mem, realtype *hnew);
static int ARKAdaptPID(ARKodeMem ark_mem, realtype *hnew);
static int ARKAdaptPI(ARKodeMem ark_mem, realtype *hnew);
static int ARKAdaptI(ARKodeMem ark_mem, realtype *hnew);
static int ARKAdaptExpGus(ARKodeMem ark_mem, realtype *hnew);
static int ARKAdaptImpGus(ARKodeMem ark_mem, realtype *hnew);
static int ARKAdaptImExGus(ARKodeMem ark_mem, realtype *hnew);

static void ARKAdjustParams(ARKodeMem ark_mem);
static void ARKAdjustOrder(ARKodeMem ark_mem, int deltaq);
static void ARKAdjustBDF(ARKodeMem ark_mem, int deltaq);
static void ARKIncreaseBDF(ARKodeMem ark_mem);
static void ARKDecreaseBDF(ARKodeMem ark_mem);

static void ARKAdjustParams2(ARKodeMem ark_mem);

static void ARKRescale(ARKodeMem ark_mem);

static int ARKDenseEval(ARKodeMem ark_mem, realtype tau,
			int d, int order, N_Vector yout);
static void ARKPredict(ARKodeMem ark_mem);
static void ARKPredict2(ARKodeMem ark_mem, int istage);

static void ARKSet2(ARKodeMem ark_mem);

static void ARKSet(ARKodeMem ark_mem);
static void ARKSetBDF(ARKodeMem ark_mem);
static void ARKSetTqBDF(ARKodeMem ark_mem, realtype hsum, 
			realtype alpha0, realtype alpha0_hat, 
			realtype xi_inv, realtype xistar_inv);

static int ARKNlsResid(ARKodeMem ark_mem, N_Vector y, 
		       N_Vector fy, N_Vector r);
static int ARKNls(ARKodeMem ark_mem, int nflag);
static int ARKNlsNewton(ARKodeMem ark_mem, int nflag);
static int ARKLs(ARKodeMem ark_mem, int nflag);

static int ARKHandleNFlag(ARKodeMem ark_mem, int *nflagPtr, 
			  realtype saved_t, int *ncfPtr);

static void ARKRestore(ARKodeMem ark_mem);

static int ARKDoErrorTest(ARKodeMem ark_mem, int *nflagPtr,
			  realtype saved_t, int *nefPtr, 
			  realtype *dsmPtr);
static int ARKDoErrorTest2(ARKodeMem ark_mem, int *nflagPtr,
			   realtype saved_t, int *nefPtr, 
			   realtype dsm);

static realtype ARKComputeSolutions(ARKodeMem ark_mem);

static void ARKCompleteStep(ARKodeMem ark_mem);
static int ARKCompleteStage(ARKodeMem ark_mem);

static void ARKPrepareNextStep(ARKodeMem ark_mem, realtype dsm);
static void ARKSetEta(ARKodeMem ark_mem);
static realtype ARKComputeEtaqm1(ARKodeMem ark_mem);
static realtype ARKComputeEtaqp1(ARKodeMem ark_mem);
static void ARKChooseEta(ARKodeMem ark_mem);

static int ARKHandleFailure(ARKodeMem ark_mem,int flag);

static int ARKRootCheck1(ARKodeMem ark_mem);
static int ARKRootCheck2(ARKodeMem ark_mem);
static int ARKRootCheck3(ARKodeMem ark_mem);
static int ARKRootfind(ARKodeMem ark_mem);

static int ARKFullRHS(ARKodeMem ark_mem, realtype t, 
		      N_Vector y, N_Vector tmp, N_Vector f);


/*===============================================================
  EXPORTED FUNCTIONS IMPLEMENTATION
===============================================================*/

/*---------------------------------------------------------------
 ARKodeCreate:

 ARKodeCreate creates an internal memory block for a problem to 
 be solved by ARKODE. If successful, ARKodeCreate returns a 
 pointer to the problem memory. This pointer should be passed to
 ARKodeInit. If an initialization error occurs, ARKodeCreate 
 prints an error message to standard err and returns NULL. 
---------------------------------------------------------------*/
void *ARKodeCreate()
{
  int i, j, iret;
  ARKodeMem ark_mem;

  ark_mem = NULL;
  ark_mem = (ARKodeMem) malloc(sizeof(struct ARKodeMemRec));
  if (ark_mem == NULL) {
    ARKProcessError(NULL, 0, "ARKODE", "ARKodeCreate", 
		    MSGARK_ARKMEM_FAIL);
    return(NULL);
  }

  /* Zero out ark_mem */
  memset(ark_mem, 0, sizeof(struct ARKodeMemRec));

  /* Set uround */
  ark_mem->ark_uround = UNIT_ROUNDOFF;

  /* Set default values for integrator optional inputs */
  iret = ARKodeSetDefaults((void *)ark_mem);
  if (iret != ARK_SUCCESS) {
    ARKProcessError(NULL, 0, "ARKODE", "ARKodeCreate", 
		    "Error setting default solver options");
    return(NULL);
  }

  /* Initialize internal RK parameters and coefficients */
  ark_mem->ark_stages = 0;
  ark_mem->ark_istage = 0;
  ark_mem->ark_q = 0;
  ark_mem->ark_p = 0;
  ark_mem->ark_stagesE = 0;
  ark_mem->ark_qE = 0;
  ark_mem->ark_pE = 0;
  for (i=0; i<ARK_S_MAX; i++) {
    for (j=0; j<ARK_S_MAX; j++) {
      ark_mem->ark_Ae[i][j] = ZERO;
      ark_mem->ark_Ai[i][j] = ZERO;
    }
    ark_mem->ark_c[i]   = ZERO;
    ark_mem->ark_b[i]   = ZERO;
    ark_mem->ark_b2[i]  = ZERO;
    ark_mem->ark_Fi[i]  = NULL;
    ark_mem->ark_cE[i]  = ZERO;
    ark_mem->ark_bE[i]  = ZERO;
    ark_mem->ark_b2E[i] = ZERO;
    ark_mem->ark_Fe[i]  = NULL;
  }

  /* Initialize root finding variables */
  ark_mem->ark_glo     = NULL;
  ark_mem->ark_ghi     = NULL;
  ark_mem->ark_grout   = NULL;
  ark_mem->ark_iroots  = NULL;
  ark_mem->ark_rootdir = NULL;
  ark_mem->ark_gfun    = NULL;
  ark_mem->ark_nrtfn   = 0;
  ark_mem->ark_gactive = NULL;
  ark_mem->ark_mxgnull = 1;

  /* Set the saved values qmax_alloc and smax_alloc
     (no vectors allocated yet, so initialize to 0) */
  ark_mem->ark_qmax_alloc = 0;
  
  /* Initialize lrw and liw */
  /***** UPDATE THESE!!! *****/
  ark_mem->ark_lrw = 58 + 2*L_MAX + NUM_TESTS;
  ark_mem->ark_liw = 40;

  /* No mallocs have been done yet */
  ark_mem->ark_VabstolMallocDone = FALSE;
  ark_mem->ark_MallocDone        = FALSE;

  /* Return pointer to ARKODE memory block */
  return((void *)ark_mem);
}


/*---------------------------------------------------------------
 ARKodeInit:
 
 ARKodeInit allocates and initializes memory for a problem. All 
 problem inputs are checked for errors. If any error occurs during 
 initialization, it is reported to the file whose file pointer is 
 errfp and an error flag is returned. Otherwise, it returns 
 ARK_SUCCESS.
---------------------------------------------------------------*/
int ARKodeInit(void *arkode_mem, ARKRhsFn fe, ARKRhsFn fi, 
	       realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  booleantype nvectorOK, allocOK;
  long int lrw1, liw1;

  /* Check arkode_mem */
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKodeInit", 
		    MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check for legal input parameters */
  if (y0==NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeInit", MSGARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  /* Set implicit/explicit problem based on function pointers */
  if (fe == NULL) ark_mem->ark_implicit = TRUE;
  if (fi == NULL) ark_mem->ark_explicit = TRUE;

  /* Check that at least one of fe,fi is supplied and is to be used */
  if (ark_mem->ark_implicit && ark_mem->ark_explicit) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE",
  		    "ARKodeInit", MSGARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Test if all required vector operations are implemented */
  nvectorOK = ARKCheckNvector(y0);
  if (!nvectorOK) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeInit", MSGARK_BAD_NVECTOR);
    return(ARK_ILL_INPUT);
  }

  /* Set space requirements for one N_Vector */
  if (y0->ops->nvspace != NULL) {
    N_VSpace(y0, &lrw1, &liw1);
  } else {
    lrw1 = 0;
    liw1 = 0;
  }
  ark_mem->ark_lrw1 = lrw1;
  ark_mem->ark_liw1 = liw1;


  /* Allocate the solver vectors (using y0 as a template) */
  allocOK = ARKAllocVectors(ark_mem, y0);
  if (!allocOK) {
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKodeInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  /* All error checking is complete at this point */

  /* Copy the input parameters into ARKODE state */
  ark_mem->ark_fe = fe;
  ark_mem->ark_fi = fi;
  ark_mem->ark_tn = t0;

  /* Set step parameters */
  ark_mem->ark_q      = 1;
  ark_mem->ark_L      = 2;
  ark_mem->ark_qwait  = ark_mem->ark_L;
  ark_mem->ark_etamax = ETAMX1;
  ark_mem->ark_qold   = 0;
  ark_mem->ark_hold   = ZERO;
  ark_mem->ark_tolsf  = ONE;

  /* Set the linear solver addresses to NULL.
     (We check != NULL later, in ARKode.) */
  ark_mem->ark_linit  = NULL;
  ark_mem->ark_lsetup = NULL;
  ark_mem->ark_lsolve = NULL;
  ark_mem->ark_lfree  = NULL;
  ark_mem->ark_lmem   = NULL;

  /* Initialize zn[0] in the history array */
  N_VScale(ONE, y0, ark_mem->ark_zn[0]);

  /* Initialize error history (1.0 => met target accuracy) */
  ark_mem->ark_hadapt_ehist[0] = ONE;
  ark_mem->ark_hadapt_ehist[1] = ONE;
  ark_mem->ark_hadapt_ehist[2] = ONE;

  /* Initialize all the counters */
  ark_mem->ark_nst     = 0;
  ark_mem->ark_nst_acc = 0;
  ark_mem->ark_nst_exp = 0;
  ark_mem->ark_nst_con = 0;
  ark_mem->ark_nfe     = 0;
  ark_mem->ark_nfi     = 0;
  ark_mem->ark_ncfn    = 0;
  ark_mem->ark_netf    = 0;
  ark_mem->ark_nni     = 0;
  ark_mem->ark_nsetups = 0;
  ark_mem->ark_nhnil   = 0;
  ark_mem->ark_nstlp   = 0;
  ark_mem->ark_nge     = 0;
  ark_mem->ark_irfnd   = 0;

  /* Initialize other integrator optional outputs */
  ark_mem->ark_h0u      = ZERO;
  ark_mem->ark_next_h   = ZERO;
  ark_mem->ark_next_q   = 0;

  /* Problem has been successfully initialized */
  ark_mem->ark_MallocDone = TRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeReInit:

 ARKodeReInit re-initializes ARKODE's memory for a problem, 
 assuming it has already been allocated in a prior ARKodeInit 
 call.  All problem specification inputs are checked for errors.
 If any error occurs during initialization, it is reported to 
 the file whose file pointer is errfp.

 The return value is ARK_SUCCESS = 0 if no errors occurred, or
 a negative value otherwise.
---------------------------------------------------------------*/
int ARKodeReInit(void *arkode_mem, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
 
  /* Check arkode_mem */
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeReInit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check if arkode_mem was allocated */
  if (ark_mem->ark_MallocDone == FALSE) {
    ARKProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", 
		    "ARKodeReInit", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check for legal input parameters */
  if (y0 == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeReInit", MSGARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }
  
  /* Copy the input parameters into ARKODE state */
  ark_mem->ark_tn = t0;
  
  /* Set step parameters */
  ark_mem->ark_q      = 1;
  ark_mem->ark_L      = 2;
  ark_mem->ark_qwait  = ark_mem->ark_L;
  ark_mem->ark_etamax = ETAMX1;
  ark_mem->ark_qold   = 0;
  ark_mem->ark_hold   = ZERO;
  ark_mem->ark_tolsf  = ONE;

  /* Initialize zn[0] in the history array */
  N_VScale(ONE, y0, ark_mem->ark_zn[0]);
 
  /* Initialize error history (1.0 => met target accuracy) */
  ark_mem->ark_hadapt_ehist[0] = ONE;
  ark_mem->ark_hadapt_ehist[1] = ONE;
  ark_mem->ark_hadapt_ehist[2] = ONE;

  /* Initialize all the counters */
  ark_mem->ark_nst     = 0;
  ark_mem->ark_nst_acc = 0;
  ark_mem->ark_nst_exp = 0;
  ark_mem->ark_nst_con = 0;
  ark_mem->ark_nfe     = 0;
  ark_mem->ark_nfi     = 0;
  ark_mem->ark_ncfn    = 0;
  ark_mem->ark_netf    = 0;
  ark_mem->ark_nni     = 0;
  ark_mem->ark_nsetups = 0;
  ark_mem->ark_nhnil   = 0;
  ark_mem->ark_nstlp   = 0;
  ark_mem->ark_nge     = 0;
  ark_mem->ark_irfnd   = 0;

  /* Initialize other integrator optional outputs */
  ark_mem->ark_h0u    = ZERO;
  ark_mem->ark_next_h = ZERO;
  ark_mem->ark_next_q = 0;

  /* Problem has been successfully re-initialized */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSStolerances:
 ARKodeSVtolerances:
 ARKodeWFtolerances:

 These functions specify the integration tolerances. One of them
 SHOULD be called before the first call to ARKode; otherwise 
 default values of reltol=1e-4 and abstol=1e-9 will be used, 
 which may be entirely incorrect for a specific problem.

 ARKodeSStolerances specifies scalar relative and absolute 
   tolerances.
 ARKodeSVtolerances specifies scalar relative tolerance and a 
   vector absolute tolerance (a potentially different absolute 
   tolerance for each vector component).
 ARKodeWFtolerances specifies a user-provides function (of type
   ARKEwtFn) which will be called to set the error weight vector.
---------------------------------------------------------------*/
int ARKodeSStolerances(void *arkode_mem, realtype reltol, 
		       realtype abstol)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSStolerances", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_MallocDone == FALSE) {
    ARKProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", 
		    "ARKodeSStolerances", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check inputs */
  if (reltol < ZERO) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSStolerances", MSGARK_BAD_RELTOL);
    return(ARK_ILL_INPUT);
  }
  if (abstol < ZERO) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSStolerances", MSGARK_BAD_ABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Copy tolerances into memory */
  ark_mem->ark_reltol    = reltol;
  ark_mem->ark_Sabstol   = abstol;
  ark_mem->ark_itol      = ARK_SS;
  ark_mem->ark_user_efun = FALSE;
  ark_mem->ark_efun      = ARKEwtSet;
  ark_mem->ark_e_data    = NULL; /* set to arkode_mem in InitialSetup */

  return(ARK_SUCCESS);
}

int ARKodeSVtolerances(void *arkode_mem, realtype reltol, 
		       N_Vector abstol)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSVtolerances", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_MallocDone == FALSE) {
    ARKProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", 
		    "ARKodeSVtolerances", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check inputs */
  if (reltol < ZERO) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSVtolerances", MSGARK_BAD_RELTOL);
    return(ARK_ILL_INPUT);
  }
  if (N_VMin(abstol) < ZERO) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSVtolerances", MSGARK_BAD_ABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Copy tolerances into memory */
  if ( !(ark_mem->ark_VabstolMallocDone) ) {
    ark_mem->ark_Vabstol = N_VClone(ark_mem->ark_ewt);
    ark_mem->ark_lrw += ark_mem->ark_lrw1;
    ark_mem->ark_liw += ark_mem->ark_liw1;
    ark_mem->ark_VabstolMallocDone = TRUE;
  }
  N_VScale(ONE, abstol, ark_mem->ark_Vabstol);
  ark_mem->ark_reltol    = reltol;
  ark_mem->ark_itol      = ARK_SV;
  ark_mem->ark_user_efun = FALSE;
  ark_mem->ark_efun      = ARKEwtSet;
  ark_mem->ark_e_data    = NULL; /* set to arkode_mem in InitialSetup */

  return(ARK_SUCCESS);
}

int ARKodeWFtolerances(void *arkode_mem, ARKEwtFn efun)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeWFtolerances", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_MallocDone == FALSE) {
    ARKProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", 
		    "ARKodeWFtolerances", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Copy tolerance data into memory */
  ark_mem->ark_itol      = ARK_WF;
  ark_mem->ark_user_efun = TRUE;
  ark_mem->ark_efun      = efun;
  ark_mem->ark_e_data    = NULL; /* set to user_data in InitialSetup */

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeRootInit:

 ARKodeRootInit initializes a rootfinding problem to be solved
 during the integration of the ODE system.  It loads the root
 function pointer and the number of root functions, and allocates
 workspace memory.  The return value is ARK_SUCCESS = 0 if no 
 errors occurred, or a negative value otherwise.
---------------------------------------------------------------*/
int ARKodeRootInit(void *arkode_mem, int nrtfn, ARKRootFn g)
{
  ARKodeMem ark_mem;
  int i, nrt;

  /* Check arkode_mem pointer */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeRootInit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  nrt = (nrtfn < 0) ? 0 : nrtfn;

  /* If rerunning ARKodeRootInit() with a different number of root
     functions (changing number of gfun components), then free
     currently held memory resources */
  if ((nrt != ark_mem->ark_nrtfn) && (ark_mem->ark_nrtfn > 0)) {
    free(ark_mem->ark_glo);     ark_mem->ark_glo     = NULL;
    free(ark_mem->ark_ghi);     ark_mem->ark_ghi     = NULL;
    free(ark_mem->ark_grout);   ark_mem->ark_grout   = NULL;
    free(ark_mem->ark_iroots);  ark_mem->ark_iroots  = NULL;
    free(ark_mem->ark_rootdir); ark_mem->ark_rootdir = NULL;
    free(ark_mem->ark_gactive); ark_mem->ark_gactive = NULL;

    ark_mem->ark_lrw -= 3 * (ark_mem->ark_nrtfn);
    ark_mem->ark_liw -= 3 * (ark_mem->ark_nrtfn);
  }

  /* If ARKodeRootInit() was called with nrtfn == 0, then set 
     ark_nrtfn to zero and ark_gfun to NULL before returning */
  if (nrt == 0) {
    ark_mem->ark_nrtfn = nrt;
    ark_mem->ark_gfun = NULL;
    return(ARK_SUCCESS);
  }

  /* If rerunning ARKodeRootInit() with the same number of root 
     functions (not changing number of gfun components), then 
     check if the root function argument has changed */
  /* If g != NULL then return as currently reserved memory 
     resources will suffice */
  if (nrt == ark_mem->ark_nrtfn) {
    if (g != ark_mem->ark_gfun) {
      if (g == NULL) {
        free(ark_mem->ark_glo);     ark_mem->ark_glo     = NULL;
        free(ark_mem->ark_ghi);     ark_mem->ark_ghi     = NULL;
        free(ark_mem->ark_grout);   ark_mem->ark_grout   = NULL;
        free(ark_mem->ark_iroots);  ark_mem->ark_iroots  = NULL;
        free(ark_mem->ark_rootdir); ark_mem->ark_rootdir = NULL;
        free(ark_mem->ark_gactive); ark_mem->ark_gactive = NULL;

        ark_mem->ark_lrw -= 3*nrt;
        ark_mem->ark_liw -= 3*nrt;

        ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKodeRootInit", MSGARK_NULL_G);
        return(ARK_ILL_INPUT);
      }
      else {
        ark_mem->ark_gfun = g;
        return(ARK_SUCCESS);
      }
    }
    else return(ARK_SUCCESS);
  }

  /* Set variable values in ARKode memory block */
  ark_mem->ark_nrtfn = nrt;
  if (g == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeRootInit", MSGARK_NULL_G);
    return(ARK_ILL_INPUT);
  }
  else ark_mem->ark_gfun = g;

  /* Allocate necessary memory and return */
  ark_mem->ark_glo = NULL;
  ark_mem->ark_glo = (realtype *) malloc(nrt*sizeof(realtype));
  if (ark_mem->ark_glo == NULL) {
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  ark_mem->ark_ghi = NULL;
  ark_mem->ark_ghi = (realtype *) malloc(nrt*sizeof(realtype));
  if (ark_mem->ark_ghi == NULL) {
    free(ark_mem->ark_glo); ark_mem->ark_glo = NULL;
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  ark_mem->ark_grout = NULL;
  ark_mem->ark_grout = (realtype *) malloc(nrt*sizeof(realtype));
  if (ark_mem->ark_grout == NULL) {
    free(ark_mem->ark_glo); ark_mem->ark_glo = NULL;
    free(ark_mem->ark_ghi); ark_mem->ark_ghi = NULL;
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  ark_mem->ark_iroots = NULL;
  ark_mem->ark_iroots = (int *) malloc(nrt*sizeof(int));
  if (ark_mem->ark_iroots == NULL) {
    free(ark_mem->ark_glo); ark_mem->ark_glo = NULL; 
    free(ark_mem->ark_ghi); ark_mem->ark_ghi = NULL;
    free(ark_mem->ark_grout); ark_mem->ark_grout = NULL;
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  ark_mem->ark_rootdir = NULL;
  ark_mem->ark_rootdir = (int *) malloc(nrt*sizeof(int));
  if (ark_mem->ark_rootdir == NULL) {
    free(ark_mem->ark_glo); ark_mem->ark_glo = NULL; 
    free(ark_mem->ark_ghi); ark_mem->ark_ghi = NULL;
    free(ark_mem->ark_grout); ark_mem->ark_grout = NULL;
    free(ark_mem->ark_iroots); ark_mem->ark_iroots = NULL;
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }
  ark_mem->ark_gactive = NULL;
  ark_mem->ark_gactive = (booleantype *) malloc(nrt*sizeof(booleantype));
  if (ark_mem->ark_gactive == NULL) {
    free(ark_mem->ark_glo); ark_mem->ark_glo = NULL; 
    free(ark_mem->ark_ghi); ark_mem->ark_ghi = NULL;
    free(ark_mem->ark_grout); ark_mem->ark_grout = NULL;
    free(ark_mem->ark_iroots); ark_mem->ark_iroots = NULL;
    free(ark_mem->ark_rootdir); ark_mem->ark_rootdir = NULL;
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODES", 
		    "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  /* Set default values for rootdir (both directions) */
  for(i=0; i<nrt; i++) ark_mem->ark_rootdir[i] = 0;

  /* Set default values for gactive (all active) */
  for(i=0; i<nrt; i++) ark_mem->ark_gactive[i] = TRUE;

  ark_mem->ark_lrw += 3*nrt;
  ark_mem->ark_liw += 3*nrt;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKode:

 This routine is the main driver of the ARKODE package. 

 It integrates over a time interval defined by the user, by 
 calling ARKStep to do internal time steps.

 The first time that ARKode is called for a successfully 
 initialized problem, it computes a tentative initial step size.

 ARKode supports two modes as specified by itask: ARK_NORMAL and
 ARK_ONE_STEP.  In the ARK_NORMAL mode, the solver steps until 
 it reaches or passes tout and then interpolates to obtain 
 y(tout).  In the ARK_ONE_STEP mode, it takes one internal step 
 and returns.
---------------------------------------------------------------*/
int ARKode(void *arkode_mem, realtype tout, N_Vector yout, 
	   realtype *tret, int itask)
{
  ARKodeMem ark_mem;
  long int nstloc;
  int retval, hflag, kflag, istate, ir, ier, irfndp;
  int ewtsetOK;
  realtype troundoff, tout_hin, rh, nrm;
  booleantype inactive_roots;

  /*-------------------------------------
    1. Check and process inputs
  -------------------------------------*/

  /* Check if arkode_mem exists */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKode", 
		    MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check if arkode_mem was allocated */
  if (ark_mem->ark_MallocDone == FALSE) {
    ARKProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", "ARKode", 
		    MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }
  
  /* Check for yout != NULL */
  if ((ark_mem->ark_y = yout) == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
		    MSGARK_YOUT_NULL);
    return(ARK_ILL_INPUT);
  }

  /* Check for tret != NULL */
  if (tret == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
		    MSGARK_TRET_NULL);
    return(ARK_ILL_INPUT);
  }

  /* Check for valid itask */
  if ( (itask != ARK_NORMAL) && (itask != ARK_ONE_STEP) ) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
		    MSGARK_BAD_ITASK);
    return(ARK_ILL_INPUT);
  }

  if (itask == ARK_NORMAL) ark_mem->ark_toutc = tout;
  ark_mem->ark_taskc = itask;


  /*----------------------------------------
    2. Initializations performed only at
       the first step (nst=0):
       - initial setup
       - initialize Nordsieck history array
       - compute initial step size
       - check for approach to tstop
       - check for approach to a root
  ----------------------------------------*/
  if (ark_mem->ark_nst == 0) {

    ier = ARKInitialSetup(ark_mem);
    if (ier!= ARK_SUCCESS) return(ier);
    
    /* Copy f(t0,y0) into ark_zn[1]. */
    if (!ark_mem->ark_explicit)
      N_VScale(ONE, ark_mem->ark_fnew, ark_mem->ark_zn[1]);
    
    /* Set initial h (from H0 or ARKHin). */
    ark_mem->ark_h = ark_mem->ark_hin;
    if ( (ark_mem->ark_h != ZERO) && 
	 ((tout-ark_mem->ark_tn)*ark_mem->ark_h < ZERO) ) {
      ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
		      MSGARK_BAD_H0);
      return(ARK_ILL_INPUT);
    }
    if (ark_mem->ark_h == ZERO) {
      tout_hin = tout;
      if ( ark_mem->ark_tstopset && 
	   (tout-ark_mem->ark_tn)*(tout-ark_mem->ark_tstop) > 0 ) 
	tout_hin = ark_mem->ark_tstop; 
      hflag = ARKHin(ark_mem, tout_hin);
      if (hflag != ARK_SUCCESS) {
        istate = ARKHandleFailure(ark_mem, hflag);
        return(istate);
      }
    }
    rh = ABS(ark_mem->ark_h)*ark_mem->ark_hmax_inv;
    if (rh > ONE) ark_mem->ark_h /= rh;
    if (ABS(ark_mem->ark_h) < ark_mem->ark_hmin) 
      ark_mem->ark_h *= ark_mem->ark_hmin/ABS(ark_mem->ark_h);

    /* Check for approach to tstop */
    if (ark_mem->ark_tstopset) {
      if ( (ark_mem->ark_tstop - ark_mem->ark_tn)*ark_mem->ark_h < ZERO ) {
        ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
			MSGARK_BAD_TSTOP, ark_mem->ark_tstop, ark_mem->ark_tn);
        return(ARK_ILL_INPUT);
      }
      if ( (ark_mem->ark_tn + ark_mem->ark_h - ark_mem->ark_tstop)*ark_mem->ark_h > ZERO ) 
        ark_mem->ark_h = (ark_mem->ark_tstop - ark_mem->ark_tn)*(ONE-FOUR*ark_mem->ark_uround);
    }

    /* Scale zn[1] by h.*/
    ark_mem->ark_hscale = ark_mem->ark_h; 
    ark_mem->ark_h0u    = ark_mem->ark_h;
    ark_mem->ark_hprime = ark_mem->ark_h;
    N_VScale(ark_mem->ark_h, ark_mem->ark_zn[1], ark_mem->ark_zn[1]);

    /* Check for zeros of root function g at and near t0. */
    if (ark_mem->ark_nrtfn > 0) {
      retval = ARKRootCheck1(ark_mem);

      if (retval == ARK_RTFUNC_FAIL) {
        ARKProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "ARKRootCheck1", 
			MSGARK_RTFUNC_FAILED, ark_mem->ark_tn);
        return(ARK_RTFUNC_FAIL);
      }
    }

  } /* end of first call block */

  /*------------------------------------------------------
    3. At following steps, perform stop tests:
       - check for root in last step
       - check if we passed tstop
       - check if we passed tout (NORMAL mode)
       - check if current tn was returned (ONE_STEP mode)
       - check if we are close to tstop
         (adjust step size if needed)
  -------------------------------------------------------*/

  if (ark_mem->ark_nst > 0) {

    /* Estimate an infinitesimal time interval to be used as
       a roundoff for time quantities (based on current time 
       and step size) */
    troundoff = FUZZ_FACTOR*ark_mem->ark_uround *
      (ABS(ark_mem->ark_tn) + ABS(ark_mem->ark_h));

    /* First, check for a root in the last step taken, other than the
       last root found, if any.  If itask = ARK_ONE_STEP and y(tn) was not
       returned because of an intervening root, return y(tn) now.     */
    if (ark_mem->ark_nrtfn > 0) {

      irfndp = ark_mem->ark_irfnd;
      
      retval = ARKRootCheck2(ark_mem);

      if (retval == CLOSERT) {
        ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKRootCheck2", 
			MSGARK_CLOSE_ROOTS, ark_mem->ark_tlo);
        return(ARK_ILL_INPUT);
      } else if (retval == ARK_RTFUNC_FAIL) {
        ARKProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "ARKRootCheck2", 
			MSGARK_RTFUNC_FAILED, ark_mem->ark_tlo);
        return(ARK_RTFUNC_FAIL);
      } else if (retval == RTFOUND) {
        ark_mem->ark_tretlast = *tret = ark_mem->ark_tlo;
        return(ARK_ROOT_RETURN);
      }

      /* If tn is distinct from tretlast (within roundoff),
         check remaining interval for roots */
      if ( ABS(ark_mem->ark_tn - ark_mem->ark_tretlast) > troundoff ) {

        retval = ARKRootCheck3(ark_mem);

        if (retval == ARK_SUCCESS) {     /* no root found */
          ark_mem->ark_irfnd = 0;
          if ((irfndp == 1) && (itask == ARK_ONE_STEP)) {
            ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
            N_VScale(ONE, ark_mem->ark_zn[0], yout);
            return(ARK_SUCCESS);
          }
        } else if (retval == RTFOUND) {  /* a new root was found */
          ark_mem->ark_irfnd = 1;
          ark_mem->ark_tretlast = *tret = ark_mem->ark_tlo;
          return(ARK_ROOT_RETURN);
        } else if (retval == ARK_RTFUNC_FAIL) {  /* g failed */
          ARKProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "ARKRootCheck3", 
			  MSGARK_RTFUNC_FAILED, ark_mem->ark_tlo);
          return(ARK_RTFUNC_FAIL);
        }

      }

    } /* end of root stop check */

    /* In ARK_NORMAL mode, test if tout was reached */
    if ( (itask == ARK_NORMAL) && 
	 ((ark_mem->ark_tn-tout)*ark_mem->ark_h >= ZERO) ) {
      ark_mem->ark_tretlast = *tret = tout;
      ier =  ARKodeGetDky(ark_mem, tout, 0, yout);
      if (ier != ARK_SUCCESS) {
        ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
			"ARKode", MSGARK_BAD_TOUT, tout);
        return(ARK_ILL_INPUT);
      }
      return(ARK_SUCCESS);
    }

    /* In ARK_ONE_STEP mode, test if tn was returned */
    if ( itask == ARK_ONE_STEP && 
	 ABS(ark_mem->ark_tn - ark_mem->ark_tretlast) > troundoff ) {
      ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
      N_VScale(ONE, ark_mem->ark_zn[0], yout);
      return(ARK_SUCCESS);
    }

    /* Test for tn at tstop or near tstop */
    if ( ark_mem->ark_tstopset ) {

      if ( ABS(ark_mem->ark_tn - ark_mem->ark_tstop) <= troundoff) {
        ier =  ARKodeGetDky(ark_mem, ark_mem->ark_tstop, 0, yout);
        if (ier != ARK_SUCCESS) {
          ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
			  MSGARK_BAD_TSTOP, ark_mem->ark_tstop, ark_mem->ark_tn);
          return(ARK_ILL_INPUT);
        }
        ark_mem->ark_tretlast = *tret = ark_mem->ark_tstop;
        ark_mem->ark_tstopset = FALSE;
        return(ARK_TSTOP_RETURN);
      }
      
      /* If next step would overtake tstop, adjust stepsize */
      if ( (ark_mem->ark_tn + ark_mem->ark_hprime - ark_mem->ark_tstop)*ark_mem->ark_h > ZERO ) {
        ark_mem->ark_hprime = (ark_mem->ark_tstop - ark_mem->ark_tn)*(ONE-FOUR*ark_mem->ark_uround);
        ark_mem->ark_eta = ark_mem->ark_hprime/ark_mem->ark_h;
      }
    }
  } /* end stopping tests block */  

  /*--------------------------------------------------
    4. Looping point for internal steps
   
       4.1. check for errors (too many steps, too much
            accuracy requested, step size too small)
       4.2. take a new step (call ARKStep)
       4.3. stop on error 
       4.4. perform stop tests:
            - check for root in last step
            - check if tout was passed
            - check if close to tstop
            - check if in ONE_STEP mode (must return)
  --------------------------------------------------*/
  nstloc = 0;
  for(;;) {
   
    ark_mem->ark_next_h = ark_mem->ark_h;
    ark_mem->ark_next_q = ark_mem->ark_q;
    
    /* Reset and check ewt */
    if (ark_mem->ark_nst > 0) {
      ewtsetOK = ark_mem->ark_efun(ark_mem->ark_zn[0], 
				   ark_mem->ark_ewt, 
				   ark_mem->ark_e_data);
      if (ewtsetOK != 0) {
        if (ark_mem->ark_itol == ARK_WF) 
          ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
			  MSGARK_EWT_NOW_FAIL, ark_mem->ark_tn);
        else 
          ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", 
			  MSGARK_EWT_NOW_BAD, ark_mem->ark_tn);
	
        istate = ARK_ILL_INPUT;
        ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
        N_VScale(ONE, ark_mem->ark_zn[0], yout);
        break;
      }
    }
    
    /* Check for too many steps */
    if ( (ark_mem->ark_mxstep>0) && (nstloc >= ark_mem->ark_mxstep) ) {
      ARKProcessError(ark_mem, ARK_TOO_MUCH_WORK, "ARKODE", "ARKode", 
		      MSGARK_MAX_STEPS, ark_mem->ark_tn);
      istate = ARK_TOO_MUCH_WORK;
      ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
      N_VScale(ONE, ark_mem->ark_zn[0], yout);
      break;
    }

    /* Check for too much accuracy requested */
    nrm = N_VWrmsNorm(ark_mem->ark_zn[0], ark_mem->ark_ewt);
    ark_mem->ark_tolsf = ark_mem->ark_uround * nrm;
    if (ark_mem->ark_tolsf > ONE) {
      ARKProcessError(ark_mem, ARK_TOO_MUCH_ACC, "ARKODE", "ARKode", 
		      MSGARK_TOO_MUCH_ACC, ark_mem->ark_tn);
      istate = ARK_TOO_MUCH_ACC;
      ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
      N_VScale(ONE, ark_mem->ark_zn[0], yout);
      ark_mem->ark_tolsf *= TWO;
      break;
    } else {
      ark_mem->ark_tolsf = ONE;
    }

    /* Check for h below roundoff level in tn */
    if (ark_mem->ark_tn + ark_mem->ark_h == ark_mem->ark_tn) {
      ark_mem->ark_nhnil++;
      if (ark_mem->ark_nhnil <= ark_mem->ark_mxhnil) 
        ARKProcessError(ark_mem, ARK_WARNING, "ARKODE", "ARKode", 
			MSGARK_HNIL, ark_mem->ark_tn, ark_mem->ark_h);
      if (ark_mem->ark_nhnil == ark_mem->ark_mxhnil) 
        ARKProcessError(ark_mem, ARK_WARNING, "ARKODE", "ARKode", 
			MSGARK_HNIL_DONE);
    }

    /* Call ARKStep to take a step */
    kflag = ARKStep(ark_mem);

    /* Process failed step cases, and exit loop */
    if (kflag != ARK_SUCCESS) {
      istate = ARKHandleFailure(ark_mem, kflag);
      ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
      N_VScale(ONE, ark_mem->ark_zn[0], yout);
      break;
    }
    
    nstloc++;

    /* Check for root in last step taken. */
    if (ark_mem->ark_nrtfn > 0) {

      retval = ARKRootCheck3(ark_mem);
      if (retval == RTFOUND) {  /* A new root was found */
        ark_mem->ark_irfnd = 1;
        istate = ARK_ROOT_RETURN;
        ark_mem->ark_tretlast = *tret = ark_mem->ark_tlo;
        break;
      } else if (retval == ARK_RTFUNC_FAIL) { /* g failed */
        ARKProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "ARKRootCheck3", 
			MSGARK_RTFUNC_FAILED, ark_mem->ark_tlo);
        istate = ARK_RTFUNC_FAIL;
        break;
      }

      /* If we are at the end of the first step and we still have
       * some event functions that are inactive, issue a warning
       * as this may indicate a user error in the implementation
       * of the root function. */
      if (ark_mem->ark_nst==1) {
        inactive_roots = FALSE;
        for (ir=0; ir<ark_mem->ark_nrtfn; ir++) { 
          if (!ark_mem->ark_gactive[ir]) {
            inactive_roots = TRUE;
            break;
          }
        }
        if ((ark_mem->ark_mxgnull > 0) && inactive_roots) {
          ARKProcessError(ark_mem, ARK_WARNING, "ARKODES", "ARKode", 
			  MSGARK_INACTIVE_ROOTS);
        }
      }
    }

    /* In NORMAL mode, check if tout reached */
    if ( (itask == ARK_NORMAL) &&
	 (ark_mem->ark_tn-tout)*ark_mem->ark_h >= ZERO ) {
      istate = ARK_SUCCESS;
      ark_mem->ark_tretlast = *tret = tout;
      (void) ARKodeGetDky(ark_mem, tout, 0, yout);
      ark_mem->ark_next_q = ark_mem->ark_qprime;
      ark_mem->ark_next_h = ark_mem->ark_hprime;
      break;
    }

    /* Check if tn is at tstop or near tstop */
    if ( ark_mem->ark_tstopset ) {
      troundoff = FUZZ_FACTOR*ark_mem->ark_uround * 
	(ABS(ark_mem->ark_tn) + ABS(ark_mem->ark_h));
      if ( ABS(ark_mem->ark_tn - ark_mem->ark_tstop) <= troundoff) {
        (void) ARKodeGetDky(ark_mem, ark_mem->ark_tstop, 0, yout);
        ark_mem->ark_tretlast = *tret = ark_mem->ark_tstop;
        ark_mem->ark_tstopset = FALSE;
        istate = ARK_TSTOP_RETURN;
        break;
      }
      if ( (ark_mem->ark_tn + ark_mem->ark_hprime - ark_mem->ark_tstop)*ark_mem->ark_h > ZERO ) {
        ark_mem->ark_hprime = (ark_mem->ark_tstop - ark_mem->ark_tn)*(ONE-FOUR*ark_mem->ark_uround);
        ark_mem->ark_eta = ark_mem->ark_hprime/ark_mem->ark_h;
      }
    }

    /* In ONE_STEP mode, copy y and exit loop */
    if (itask == ARK_ONE_STEP) {
      istate = ARK_SUCCESS;
      ark_mem->ark_tretlast = *tret = ark_mem->ark_tn;
      N_VScale(ONE, ark_mem->ark_zn[0], yout);
      ark_mem->ark_next_q = ark_mem->ark_qprime;
      ark_mem->ark_next_h = ark_mem->ark_hprime;
      break;
    }

  } /* end looping for internal steps */

  return(istate);
}


/*---------------------------------------------------------------
 ARKodeGetDky:

 This routine computes the k-th derivative of the interpolating
 polynomial at the time t and stores the result in the vector 
 dky. The formula is:
         q 
  dky = SUM c(j,k) * (t - tn)^(j-k) * h^(-j) * zn[j] , 
        j=k 
 where c(j,k) = j*(j-1)*...*(j-k+1), q is the current order, and
 zn[j] is the j-th column of the Nordsieck history array.

 This function is called by ARKode with k = 0 and t = tout, but
 may also be called directly by the user.
---------------------------------------------------------------*/
int ARKodeGetDky(void *arkode_mem, realtype t, int k, N_Vector dky)
{
  realtype s, c, r;
  realtype tfuzz, tp, tn1;
  int i, j;
  ARKodeMem ark_mem;
  
  /* Check all inputs for legality */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKodeGetDky", 
		    MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  if (dky == NULL) {
    ARKProcessError(ark_mem, ARK_BAD_DKY, "ARKODE", "ARKodeGetDky", 
		    MSGARK_NULL_DKY);
    return(ARK_BAD_DKY);
  }
  if ((k < 0) || (k > ark_mem->ark_q)) {
    ARKProcessError(ark_mem, ARK_BAD_K, "ARKODE", "ARKodeGetDky", 
		    MSGARK_BAD_K);
    return(ARK_BAD_K);
  }
  
  /* Allow for some slack */
  tfuzz = FUZZ_FACTOR * ark_mem->ark_uround * 
    (ABS(ark_mem->ark_tn) + ABS(ark_mem->ark_hold));
  if (ark_mem->ark_hold < ZERO) tfuzz = -tfuzz;
  tp = ark_mem->ark_tn - ark_mem->ark_hold - tfuzz;
  tn1 = ark_mem->ark_tn + tfuzz;
  if ((t-tp)*(t-tn1) > ZERO) {
    ARKProcessError(ark_mem, ARK_BAD_T, "ARKODE", "ARKodeGetDky", 
		    MSGARK_BAD_T, t, ark_mem->ark_tn-ark_mem->ark_hold, 
		    ark_mem->ark_tn);
    return(ARK_BAD_T);
  }

  /* Sum the differentiated interpolating polynomial */
  s = (t - ark_mem->ark_tn) / ark_mem->ark_h;
  for (j=ark_mem->ark_q; j >= k; j--) {
    c = ONE;
    for (i=j; i >= j-k+1; i--) c *= i;
    if (j == ark_mem->ark_q) {
      N_VScale(c, ark_mem->ark_zn[ark_mem->ark_q], dky);
    } else {
      N_VLinearSum(c, ark_mem->ark_zn[j], s, dky, dky);
    }
  }
  if (k == 0) return(ARK_SUCCESS);
  r = RPowerI(ark_mem->ark_h,-k);
  N_VScale(r, dky, dky);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeFree:

 This routine frees the problem memory allocated by ARKodeInit.
 Such memory includes all the vectors allocated by 
 ARKAllocVectors, and ARKAllocRKVectors, and the memory lmem for
 the linear solver (deallocated by a call to lfree).
---------------------------------------------------------------*/
void ARKodeFree(void **arkode_mem)
{
  ARKodeMem ark_mem;

  if (*arkode_mem == NULL) return;

  ark_mem = (ARKodeMem) (*arkode_mem);
  
  ARKFreeVectors(ark_mem);

  if (ark_mem->ark_lfree != NULL) 
    ark_mem->ark_lfree(ark_mem);

  if (ark_mem->ark_nrtfn > 0) {
    free(ark_mem->ark_glo);     ark_mem->ark_glo     = NULL;
    free(ark_mem->ark_ghi);     ark_mem->ark_ghi     = NULL;
    free(ark_mem->ark_grout);   ark_mem->ark_grout   = NULL;
    free(ark_mem->ark_iroots);  ark_mem->ark_iroots  = NULL;
    free(ark_mem->ark_rootdir); ark_mem->ark_rootdir = NULL;
    free(ark_mem->ark_gactive); ark_mem->ark_gactive = NULL;
  }

  free(*arkode_mem);
  *arkode_mem = NULL;
}



/*===============================================================
   Private Functions Implementation
===============================================================*/

/*---------------------------------------------------------------
 ARKCheckNvector:

 This routine checks if all required vector operations are 
 present.  If any of them is missing it returns FALSE.
---------------------------------------------------------------*/
static booleantype ARKCheckNvector(N_Vector tmpl)
{
  if ((tmpl->ops->nvclone     == NULL) ||
      (tmpl->ops->nvdestroy   == NULL) ||
      (tmpl->ops->nvlinearsum == NULL) ||
      (tmpl->ops->nvconst     == NULL) ||
      (tmpl->ops->nvprod      == NULL) ||
      (tmpl->ops->nvdiv       == NULL) ||
      (tmpl->ops->nvscale     == NULL) ||
      (tmpl->ops->nvabs       == NULL) ||
      (tmpl->ops->nvinv       == NULL) ||
      (tmpl->ops->nvaddconst  == NULL) ||
      (tmpl->ops->nvmaxnorm   == NULL) ||
      (tmpl->ops->nvwrmsnorm  == NULL) ||
      (tmpl->ops->nvmin       == NULL))
    return(FALSE);
  else
    return(TRUE);
}


/*---------------------------------------------------------------
 ARKAllocVectors

 This routine allocates the ARKODE vectors ewt, acor, tempv, 
 ftemp, zn[0], ..., zn[qmax], Fe[0], ..., Fe[stages-1], and 
 Fi[0], ..., Fi[stages-1].  If any of these vectors 
 already exist, they are left alone.  Otherwise, it will 
 allocate each vector by cloning the input vector. If all memory
 allocations are successful, ARKAllocVectors returns TRUE. 
 Otherwise all vector memory is freed and ARKAllocVectors 
 returns FALSE. This routine also updates the optional outputs 
 lrw and liw, which are (respectively) the lengths of the real 
 and integer work spaces.
---------------------------------------------------------------*/
static booleantype ARKAllocVectors(ARKodeMem ark_mem, N_Vector tmpl)
{
  int j;

  /* Allocate ewt if needed */
  if (ark_mem->ark_ewt == NULL) {
    ark_mem->ark_ewt = N_VClone(tmpl);
    if (ark_mem->ark_ewt == NULL) {
      ARKFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate acor if needed */
  if (ark_mem->ark_acor == NULL) {
    ark_mem->ark_acor = N_VClone(tmpl);
    if (ark_mem->ark_acor == NULL) {
      ARKFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate tempv if needed */
  if (ark_mem->ark_tempv == NULL) {
    ark_mem->ark_tempv = N_VClone(tmpl);
    if (ark_mem->ark_tempv == NULL) {
      ARKFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate ftemp if needed */
  if (ark_mem->ark_ftemp == NULL) {
    ark_mem->ark_ftemp = N_VClone(tmpl);
    if (ark_mem->ark_ftemp == NULL) {
      ARKFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate fold if needed */
  if (ark_mem->ark_fold == NULL) {
    ark_mem->ark_fold = N_VClone(tmpl);
    if (ark_mem->ark_fold == NULL) {
      ARKFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate fnew if needed */
  if (ark_mem->ark_fnew == NULL) {
    ark_mem->ark_fnew = N_VClone(tmpl);
    if (ark_mem->ark_fnew == NULL) {
      ARKFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate yold if needed */
  if (ark_mem->ark_yold == NULL) {
    ark_mem->ark_yold = N_VClone(tmpl);
    if (ark_mem->ark_yold == NULL) {
      ARKFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate ynew if needed */
  if (ark_mem->ark_ynew == NULL) {
    ark_mem->ark_ynew = N_VClone(tmpl);
    if (ark_mem->ark_ynew == NULL) {
      ARKFreeVectors(ark_mem);
      return(FALSE);
    } else {
      ark_mem->ark_lrw += ark_mem->ark_lrw1;
      ark_mem->ark_liw += ark_mem->ark_liw1;
    }
  }

  /* Allocate zn[0] ... zn[qmax] if needed */
  for (j=0; j <= ark_mem->ark_qmax; j++) {
    if (ark_mem->ark_zn[j] == NULL) {
      ark_mem->ark_zn[j] = N_VClone(tmpl);
      if (ark_mem->ark_zn[j] == NULL) {
	ARKFreeVectors(ark_mem);
	return(FALSE);
      } else {
	ark_mem->ark_lrw += ark_mem->ark_lrw1;
	ark_mem->ark_liw += ark_mem->ark_liw1;
      }
    }
  }

  /* Update the value of qmax_alloc */
  ark_mem->ark_qmax_alloc = 
    MAX(ark_mem->ark_qmax_alloc, ark_mem->ark_qmax);

  return(TRUE);
}


/*---------------------------------------------------------------
 ARKAllocRKVectors

 This routine allocates the ARKODE vectors 
 Fe[0], ..., Fe[stages-1], Fi[0], ..., Fi[stages-1], and fa and 
 fb (if necessary).

 If any of these vectors already exist, they are left alone.  
 Otherwise, it will allocate each vector by cloning the input 
 vector. If all memory allocations are successful, 
 ARKAllocRKVectors returns TRUE.  Otherwise all vector memory is 
 freed and ARKAllocRKVectors returns FALSE. This routine also 
 updates the optional outputs lrw and liw, which are 
 (respectively) the lengths of the real and integer work spaces.

 NOTE: this routine is separate from ARKAllocVectors, since that
 routine is called in ARKodeInit, which may be called by a user
 prior to setting their desired method order, 
 implicit/explicit/imex problem, and number of RK stages.  Hence
 this routine is called during ARKInitialStep, which is called 
 on the first call to the ARKode() integrator.
---------------------------------------------------------------*/
static int ARKAllocRKVectors(ARKodeMem ark_mem)
{
  int j;

  /* Allocate Fe[0] ... Fe[stages] if needed */
  if (!ark_mem->ark_implicit) {
    for (j=0; j<ark_mem->ark_stages; j++) {
      if (ark_mem->ark_Fe[j] == NULL) {
	ark_mem->ark_Fe[j] = N_VClone(ark_mem->ark_ewt);
	if (ark_mem->ark_Fe[j] == NULL) {
	  ARKFreeVectors(ark_mem);
	  return(ARK_MEM_FAIL);
	} else {
	  ark_mem->ark_lrw += ark_mem->ark_lrw1;
	  ark_mem->ark_liw += ark_mem->ark_liw1;
	}
      }
    }
  }

  /* Allocate Fi[0] ... Fi[stages] if needed */
  if (!ark_mem->ark_explicit) {
    for (j=0; j<ark_mem->ark_stages; j++) {
      if (ark_mem->ark_Fi[j] == NULL) {
	ark_mem->ark_Fi[j] = N_VClone(ark_mem->ark_ewt);
	if (ark_mem->ark_Fi[j] == NULL) {
	  ARKFreeVectors(ark_mem);
	  return(ARK_MEM_FAIL);
	} else {
	  ark_mem->ark_lrw += ark_mem->ark_lrw1;
	  ark_mem->ark_liw += ark_mem->ark_liw1;
	}
      }
    }
  }

  /* Allocate fa if needed */
  if (ark_mem->ark_dense_q > 3) {
    if (ark_mem->ark_fa == NULL) {
      ark_mem->ark_fa = N_VClone(ark_mem->ark_ewt);
      if (ark_mem->ark_fa == NULL) {
	ARKFreeVectors(ark_mem);
	return(FALSE);
      } else {
	ark_mem->ark_lrw += ark_mem->ark_lrw1;
	ark_mem->ark_liw += ark_mem->ark_liw1;
      }
    }
  }

  /* Allocate fb if needed */
  if (ark_mem->ark_dense_q == 5) {
    if (ark_mem->ark_fb == NULL) {
      ark_mem->ark_fb = N_VClone(ark_mem->ark_ewt);
      if (ark_mem->ark_fb == NULL) {
	ARKFreeVectors(ark_mem);
	return(FALSE);
      } else {
	ark_mem->ark_lrw += ark_mem->ark_lrw1;
	ark_mem->ark_liw += ark_mem->ark_liw1;
      }
    }
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKFreeVectors

 This routine frees the ARKODE vectors allocated in both
 ARKAllocVectors and ARKAllocRKVectors.
---------------------------------------------------------------*/
static void ARKFreeVectors(ARKodeMem ark_mem)
{
  int j;
  if (ark_mem->ark_ewt != NULL) {
    N_VDestroy(ark_mem->ark_ewt);
    ark_mem->ark_ewt = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_acor != NULL) {
    N_VDestroy(ark_mem->ark_acor);
    ark_mem->ark_acor = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_tempv != NULL) {
    N_VDestroy(ark_mem->ark_tempv);
    ark_mem->ark_tempv = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_ftemp != NULL) {
    N_VDestroy(ark_mem->ark_ftemp);
    ark_mem->ark_ftemp = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  for(j=0; j<=ark_mem->ark_qmax_alloc; j++) {
    if (ark_mem->ark_zn[j] != NULL) {
      N_VDestroy(ark_mem->ark_zn[j]);
      ark_mem->ark_zn[j] = NULL;
      ark_mem->ark_lrw -= ark_mem->ark_lrw1;
      ark_mem->ark_liw -= ark_mem->ark_liw1;
    }
  }
  if (ark_mem->ark_fold != NULL) {
    N_VDestroy(ark_mem->ark_fold);
    ark_mem->ark_fold = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_fnew != NULL) {
    N_VDestroy(ark_mem->ark_fnew);
    ark_mem->ark_fnew = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_yold != NULL) {
    N_VDestroy(ark_mem->ark_yold);
    ark_mem->ark_yold = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_ynew != NULL) {
    N_VDestroy(ark_mem->ark_ynew);
    ark_mem->ark_ynew = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_fa != NULL) {
    N_VDestroy(ark_mem->ark_fa);
    ark_mem->ark_fa = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  if (ark_mem->ark_fb != NULL) {
    N_VDestroy(ark_mem->ark_fb);
    ark_mem->ark_fb = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
  for(j=0; j<ARK_S_MAX; j++) {
    if (ark_mem->ark_Fe[j] != NULL) {
      N_VDestroy(ark_mem->ark_Fe[j]);
      ark_mem->ark_Fe[j] = NULL;
      ark_mem->ark_lrw -= ark_mem->ark_lrw1;
      ark_mem->ark_liw -= ark_mem->ark_liw1;
    }
    if (ark_mem->ark_Fi[j] != NULL) {
      N_VDestroy(ark_mem->ark_Fi[j]);
      ark_mem->ark_Fi[j] = NULL;
      ark_mem->ark_lrw -= ark_mem->ark_lrw1;
      ark_mem->ark_liw -= ark_mem->ark_liw1;
    }
  }
  if (ark_mem->ark_Vabstol != NULL) {
    N_VDestroy(ark_mem->ark_Vabstol);
    ark_mem->ark_Vabstol = NULL;
    ark_mem->ark_lrw -= ark_mem->ark_lrw1;
    ark_mem->ark_liw -= ark_mem->ark_liw1;
  }
}


/*---------------------------------------------------------------
 ARKInitialSetup

 This routine performs input consistency checks at the first step.
 If needed, it also checks the linear solver module and calls the
 linear solver initialization routine.
---------------------------------------------------------------*/
static int ARKInitialSetup(ARKodeMem ark_mem)
{
  int i, ier;

  /* Did the user specify tolerances? */
  if (ark_mem->ark_itol == ARK_NN) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKInitialSetup", MSGARK_NO_TOLS);
    return(ARK_ILL_INPUT);
  }

  /* Create Butcher tables */
  ier = ARKSetButcherTables(ark_mem);
  if (ier != ARK_SUCCESS) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKInitialStep", "Could not create Butcher table");
    return(ARK_ILL_INPUT);
  }

  /* Allocate ARK RHS vector memory */
  ier = ARKAllocRKVectors(ark_mem);
  if (ier != ARK_SUCCESS) {
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", 
		    "ARKInitialSetup", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  /* Fill in implicit guesses for first step with initial condition */
  if (!ark_mem->ark_explicit) 
    for (i=0; i<ark_mem->ark_stages; i++)
      N_VScale(ONE, ark_mem->ark_zn[0], ark_mem->ark_Fi[i]);

  /* Set data for efun */
  if (ark_mem->ark_user_efun) 
    ark_mem->ark_e_data = ark_mem->ark_user_data;
  else                        
    ark_mem->ark_e_data = ark_mem;

  /* Load initial error weights */
  ier = ark_mem->ark_efun(ark_mem->ark_zn[0], 
			  ark_mem->ark_ewt, 
			  ark_mem->ark_e_data);
  if (ier != 0) {
    if (ark_mem->ark_itol == ARK_WF) 
      ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "ARKInitialSetup", MSGARK_EWT_FAIL);
    else 
      ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "ARKInitialSetup", MSGARK_BAD_EWT);
    return(ARK_ILL_INPUT);
  }
  
  /* Check if lsolve function exists and call linit (if it exists) */
  if (ark_mem->ark_lsolve == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKInitialSetup", MSGARK_LSOLVE_NULL);
    return(ARK_ILL_INPUT);
  }
  if (ark_mem->ark_linit != NULL) {
    ier = ark_mem->ark_linit(ark_mem);
    if (ier != 0) {
      ARKProcessError(ark_mem, ARK_LINIT_FAIL, "ARKODE", 
		      "ARKInitialSetup", MSGARK_LINIT_FAIL);
      return(ARK_LINIT_FAIL);
    }
  }

  /* Fill initial ynew and fnew arrays */
  N_VScale(ONE, ark_mem->ark_zn[0], ark_mem->ark_ynew);
  ier = ARKFullRHS(ark_mem, ark_mem->ark_tn, ark_mem->ark_zn[0], 
		   ark_mem->ark_ftemp, ark_mem->ark_fnew);
  if (ier != 0)  return(ARK_RHSFUNC_FAIL);
  
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 PRIVATE FUNCTIONS FOR ARKODE
---------------------------------------------------------------*/

/*---------------------------------------------------------------
 ARKHin

 This routine computes a tentative initial step size h0. 
 If tout is too close to tn (= t0), then ARKHin returns 
 ARK_TOO_CLOSE and h remains uninitialized. Note that here tout 
 is either the value passed to ARKode at the first call or the 
 value of tstop (if tstop is enabled and it is closer to t0=tn 
 than tout). If the RHS function fails unrecoverably, ARKHin 
 returns ARK_RHSFUNC_FAIL. If the RHS function fails recoverably 
 too many times and recovery is not possible, ARKHin returns 
 ARK_REPTD_RHSFUNC_ERR. Otherwise, ARKHin sets h to the chosen 
 value h0 and returns ARK_SUCCESS.

 The algorithm used seeks to find h0 as a solution of
       (WRMS norm of (h0^2 ydd / 2)) = 1, 
 where ydd = estimated second derivative of y.  Although this 
 choice is based on an error expansion of the Backward Euler
 method, and hence results in an overly-conservative time step 
 for our higher-order ARK methods, it does find an order-of-
 magnitude estimate of the initial time scale of the solution.
 Since this method is only used on the first time step, the
 additional caution will not overly hinder solver efficiency.

 We start with an initial estimate equal to the geometric mean 
 of the lower and upper bounds on the step size.

 Loop up to MAX_ITERS times to find h0.
 Stop if new and previous values differ by a factor < 2.
 Stop if hnew/hg > 2 after one iteration, as this probably 
 means that the ydd value is bad because of cancellation error.        
  
 For each new proposed hg, we allow MAX_ITERS attempts to
 resolve a possible recoverable failure from f() by reducing
 the proposed stepsize by a factor of 0.2. If a legal stepsize
 still cannot be found, fall back on a previous value if 
 possible, or else return ARK_REPTD_RHSFUNC_ERR.

 Finally, we apply a bias (0.5) and verify that h0 is within 
 bounds.
---------------------------------------------------------------*/
static int ARKHin(ARKodeMem ark_mem, realtype tout)
{
  int retval, sign, count1, count2;
  realtype tdiff, tdist, tround, hlb, hub;
  realtype hg, hgs, hs, hnew, hrat, h0, yddnrm;
  booleantype hgOK, hnewOK;

  /* If tout is too close to tn, give up */
  if ((tdiff = tout-ark_mem->ark_tn) == ZERO) return(ARK_TOO_CLOSE);
  
  sign = (tdiff > ZERO) ? 1 : -1;
  tdist = ABS(tdiff);
  tround = ark_mem->ark_uround * MAX(ABS(ark_mem->ark_tn), ABS(tout));

  if (tdist < TWO*tround) return(ARK_TOO_CLOSE);
  
  /* Set lower and upper bounds on h0, and take geometric mean 
     as first trial value.
     Exit with this value if the bounds cross each other. */
  hlb = HLB_FACTOR * tround;
  hub = ARKUpperBoundH0(ark_mem, tdist);

  hg  = RSqrt(hlb*hub);

  if (hub < hlb) {
    if (sign == -1) ark_mem->ark_h = -hg;
    else            ark_mem->ark_h =  hg;
    return(ARK_SUCCESS);
  }
  
  /* Outer loop */
  hnewOK = FALSE;
  hs = hg;     /* safeguard against 'uninitialized variable' warning */
  for(count1 = 1; count1 <= MAX_ITERS; count1++) {

    /* Attempts to estimate ydd */
    hgOK = FALSE;

    for (count2 = 1; count2 <= MAX_ITERS; count2++) {
      hgs = hg*sign;
      retval = ARKYddNorm(ark_mem, hgs, &yddnrm);
      /* If f() failed unrecoverably, give up */
      if (retval < 0) return(ARK_RHSFUNC_FAIL);
      /* If successful, we can use ydd */
      if (retval == ARK_SUCCESS) {hgOK = TRUE; break;}
      /* f() failed recoverably; cut step size and test it again */
      hg *= POINT2;
    }

    /* If f() failed recoverably MAX_ITERS times */
    if (!hgOK) {
      /* Exit if this is the first or second pass. No recovery possible */
      if (count1 <= 2) return(ARK_REPTD_RHSFUNC_ERR);
      /* We have a fall-back option. The value hs is a previous hnew which
         passed through f(). Use it and break */
      hnew = hs;
      break;
    }

    /* The proposed step size is feasible. Save it. */
    hs = hg;

    /* If stopping criteria was met, or if this is the last pass, stop */
    if ( (hnewOK) || (count1 == MAX_ITERS))  {hnew = hg; break;}

    /* Propose new step size */
    hnew = (yddnrm*hub*hub > TWO) ? RSqrt(TWO/yddnrm) : RSqrt(hg*hub);
    hrat = hnew/hg;
    
    /* Accept hnew if it does not differ from hg by more than a factor of 2 */
    if ((hrat > HALF) && (hrat < TWO))  hnewOK = TRUE;

    /* After one pass, if ydd seems to be bad, use fall-back value. */
    if ((count1 > 1) && (hrat > TWO)) {
      hnew = hg;
      hnewOK = TRUE;
    }

    /* Send this value back through f() */
    hg = hnew;
  }

  /* Apply bounds, bias factor, and attach sign */
  h0 = H_BIAS*hnew;
  if (h0 < hlb) h0 = hlb;
  if (h0 > hub) h0 = hub;
  if (sign == -1) h0 = -h0;
  ark_mem->ark_h = h0;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKUpperBoundH0

 This routine sets an upper bound on abs(h0) based on
 tdist = tn - t0 and the values of y[i]/y'[i].
---------------------------------------------------------------*/
static realtype ARKUpperBoundH0(ARKodeMem ark_mem, realtype tdist)
{
  realtype hub_inv, hub;
  N_Vector temp1, temp2;

  /* Bound based on |y0|/|y0'| -- allow at most an increase of
   * HUB_FACTOR in y0 (based on a forward Euler step). The weight 
   * factor is used as a safeguard against zero components in y0. */
  temp1 = ark_mem->ark_tempv;
  temp2 = ark_mem->ark_acor;

  N_VAbs(ark_mem->ark_ynew, temp2);
  ark_mem->ark_efun(ark_mem->ark_ynew, temp1, ark_mem->ark_e_data);
  N_VInv(temp1, temp1);
  N_VLinearSum(HUB_FACTOR, temp2, ONE, temp1, temp1);

  N_VAbs(ark_mem->ark_fnew, temp2);

  N_VDiv(temp2, temp1, temp1);
  hub_inv = N_VMaxNorm(temp1);

  /* bound based on tdist -- allow at most a step of magnitude
   * HUB_FACTOR * tdist */
  hub = HUB_FACTOR*tdist;

  /* Use the smaller of the two */
  if (hub*hub_inv > ONE) hub = ONE/hub_inv;

  return(hub);
}


/*---------------------------------------------------------------
 ARKFullRHS

 This routine computes the full ODE RHS [fe(t,y) + fi(t,y)].

 Inputs (unchanged): ark_mem, t, y
 Output: f
 Temporary: tmp
---------------------------------------------------------------*/
static int ARKFullRHS(ARKodeMem ark_mem, realtype t, 
		      N_Vector y, N_Vector tmp, N_Vector f)
{
  int retval;

  /* explicit problem -- only call fe */
  if (ark_mem->ark_explicit) {

    retval = ark_mem->ark_fe(t, y, f, ark_mem->ark_user_data);
    ark_mem->ark_nfe++;
    if (retval != 0) {
      ARKProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", "ARKFullRHS", 
		      MSGARK_RHSFUNC_FAILED, ark_mem->ark_tn);
      return(ARK_RHSFUNC_FAIL);
    }


  /* implicit problem -- only call fi */
  } else if (ark_mem->ark_implicit) {

    retval = ark_mem->ark_fi(t, y, f, ark_mem->ark_user_data);
    ark_mem->ark_nfi++;
    if (retval != 0) {
      ARKProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", 
		      "ARKFullRHS", MSGARK_RHSFUNC_FAILED, ark_mem->ark_tn);
      return(ARK_RHSFUNC_FAIL);
    }


  /* imex problem */
  } else {
    
    /* explicit portion */
    retval = ark_mem->ark_fe(t, y, f, ark_mem->ark_user_data);
    ark_mem->ark_nfe++;
    if (retval != 0) {
      ARKProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", "ARKFullRHS", 
		      MSGARK_RHSFUNC_FAILED, ark_mem->ark_tn);
      return(ARK_RHSFUNC_FAIL);
    }

    /* implicit portion */
    retval = ark_mem->ark_fi(t, y, tmp, ark_mem->ark_user_data);
    ark_mem->ark_nfi++;
    if (retval != 0) {
      ARKProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", 
		      "ARKFullRHS", MSGARK_RHSFUNC_FAILED, ark_mem->ark_tn);
      return(ARK_RHSFUNC_FAIL);
    }
    N_VLinearSum(ONE, tmp, ONE, f, f);
  }
  
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKYddNorm

 This routine computes an estimate of the second derivative of y
 using a difference quotient, and returns its WRMS norm.
---------------------------------------------------------------*/
static int ARKYddNorm(ARKodeMem ark_mem, realtype hg, realtype *yddnrm)
{
  int retval;

  N_VLinearSum(hg, ark_mem->ark_fnew, ONE, 
	       ark_mem->ark_ynew, ark_mem->ark_y);
  retval = ARKFullRHS(ark_mem, ark_mem->ark_tn+hg, ark_mem->ark_y, 
		      ark_mem->ark_ftemp, ark_mem->ark_tempv);
  if (retval != 0) return(ARK_RHSFUNC_FAIL);
  N_VLinearSum(ONE, ark_mem->ark_tempv, -ONE, 
	       ark_mem->ark_fnew, ark_mem->ark_tempv);
  N_VScale(ONE/hg, ark_mem->ark_tempv, ark_mem->ark_tempv);

  *yddnrm = N_VWrmsNorm(ark_mem->ark_tempv, ark_mem->ark_ewt);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKStep

 This routine performs one internal arkode step, from tn to tn + h.
 It calls other routines to do all the work.

 The main operations done here are as follows:
 - preliminary adjustments if a new step size was chosen;
 - prediction of the Nordsieck history array zn at tn + h;
 - setting of multistep method coefficients and test quantities;
 - solution of the nonlinear system;
 - testing the local error;
 - updating zn and other state data if successful;
 - resetting stepsize and order for the next step.
 On a failure in the nonlinear system solution or error test, the
 step may be reattempted, depending on the nature of the failure.
---------------------------------------------------------------*/
static int ARKStep(ARKodeMem ark_mem)
{
  realtype saved_t, dsm;
  int ncf, nef;
  int nflag, kflag, eflag;
  
  saved_t = ark_mem->ark_tn;
  ncf = nef = 0;
  nflag = FIRST_CALL;

  if ((ark_mem->ark_nst > 0) && (ark_mem->ark_hprime != ark_mem->ark_h)) 
    ARKAdjustParams(ark_mem);
  
  /* Looping point for attempts to take a step */
  for(;;) {  

    /* Update current time with new time step.  If tstop is enabled, it is
       possible for tn + h to be past tstop by roundoff, and in that case, 
       we reset tn (after incrementing by h) to tstop. */
    ark_mem->ark_tn += ark_mem->ark_h;
    if (ark_mem->ark_tstopset) {
      if ((ark_mem->ark_tn - ark_mem->ark_tstop)*ark_mem->ark_h > ZERO) 
	ark_mem->ark_tn = ark_mem->ark_tstop;
    }

    /* Call predictor */
    ARKPredict(ark_mem);  

    /* Set up internal data structures for evaluation of BDF residual */
    ARKSet(ark_mem);

    nflag = ARKNls(ark_mem, nflag);
    kflag = ARKHandleNFlag(ark_mem, &nflag, saved_t, &ncf);

    /* Go back in loop if we need to predict again (nflag=PREV_CONV_FAIL)*/
    if (kflag == PREDICT_AGAIN) continue;

    /* Return if nonlinear solve failed and recovery not possible. */
    if (kflag != DO_ERROR_TEST) return(kflag);

    /* Perform error test (nflag=ARK_SUCCESS) */
    eflag = ARKDoErrorTest(ark_mem, &nflag, saved_t, &nef, &dsm);

    /* Go back in loop if we need to predict again (nflag=PREV_ERR_FAIL) */
    if (eflag == TRY_AGAIN)  continue;

    /* Return if error test failed and recovery not possible. */
    if (eflag != ARK_SUCCESS) return(eflag);

    /* Error test passed (eflag=ARK_SUCCESS), break from loop */
    break;

  }

  /* Nonlinear system solve and error test were both successful.
     Update data, and consider change of step and/or order.       */
  ARKCompleteStep(ark_mem); 
  ARKPrepareNextStep(ark_mem, dsm); 

  ark_mem->ark_etamax = ETAMX2;

  /*  Finally, we rescale the acor array to be the 
      estimated local error vector. */
  N_VScale(ark_mem->ark_tq[2], ark_mem->ark_acor, ark_mem->ark_acor);
  return(ARK_SUCCESS);
      
}


/*---------------------------------------------------------------
 ARKStep2

 This routine performs one internal arkode step, from tn to tn + h.
 It calls other routines to do all the work.

 The main operations done here are as follows:
 - preliminary adjustments if a new step size was chosen;
 - loop over internal time step stages:
   - predict stage solution (if implicit);
   - set up RHS data using old solutions;
   - solve the implicit system, or evaluate explicit update;
 - testing the local error;
 - reset stepsize for the next step.
 On a failure in the linear/nonlinear solvers or error test, the
 step may be reattempted, depending on the nature of the failure.
---------------------------------------------------------------*/
static int ARKStep2(ARKodeMem ark_mem)
{
  realtype saved_t, dsm;
  int retval, ncf, nef, is, nflag, kflag, eflag;
  
  saved_t = ark_mem->ark_tn;
  ncf = nef = 0;
  nflag = FIRST_CALL;

  if ((ark_mem->ark_nst > 0) && (ark_mem->ark_hprime != ark_mem->ark_h)) 
    ARKAdjustParams2(ark_mem);
  
  /* Looping point for attempts to take a step */
  for(;;) {  

    /* Loop over internal stages to the step */
    for (is=0; is<ark_mem->ark_stages; is++) {

      /* store current stage index */
      ark_mem->ark_istage = is;

      /* Looping point for attempts to solve a stage */
      for(;;) {  

	/* Update current time with new time step.  If tstop is enabled, it is
	   possible for tn + h to be past tstop by roundoff, and in that case, 
	   we reset tn (after incrementing by h) to tstop. */
	ark_mem->ark_tn = ark_mem->ark_told + ark_mem->ark_c[is]*ark_mem->ark_h;
	if (ark_mem->ark_tstopset) {
	  if ((ark_mem->ark_tn - ark_mem->ark_tstop)*ark_mem->ark_h > ZERO) 
	    ark_mem->ark_tn = ark_mem->ark_tstop;
	}

	/* Call predictor */
	ARKPredict2(ark_mem, is);  
	
	/* Set up internal data structures for evaluation of ARK residual */
	ARKSet2(ark_mem);

	/* solve implicit problem (if required) */
	if (ABS(ark_mem->ark_Ai[is][is]) > TINY) {

	  nflag = ARKNls(ark_mem, nflag);
	  kflag = ARKHandleNFlag(ark_mem, &nflag, saved_t, &ncf);

	  /* If h reduced and step needs to be retried, break loop */
	  if (kflag == PREDICT_AGAIN) break;

	  /* Return if nonlinear solve failed and recovery not possible. */
	  /***** RENAME THIS RETURN FLAG TO SIGNIFY A SUCCESSFUL/FAILED SOLVE, AND NOT A NEED TO DO AN ERROR TEST, SINCE WE'LL PERFORM MULTIPLE STAGES ****/
	  if (kflag != DO_ERROR_TEST) return(kflag);

	/* otherwise just update necessary data structures */
	} else {

	  /* set y to be RHS data computed in ARKSet */
	  N_VScale(ONE, ark_mem->ark_Fe[is], ark_mem->ark_y);

	  /* evaluate Fi at y for non-explicit problems */
	  if (!ark_mem->ark_explicit) {
	    retval = ark_mem->ark_fi(ark_mem->ark_tn, ark_mem->ark_y,
				     ark_mem->ark_Fi[is], ark_mem->ark_user_data);
	    ark_mem->ark_nfi++; 
	    if (retval < 0)  return(ARK_RHSFUNC_FAIL);
	    if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
	  }
	}

	/* successful stage solve, update remaining stage data */
	retval = ARKCompleteStage(ark_mem);
	if (retval < 0)  return(ARK_RHSFUNC_FAIL);
	if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);

      } /* loop over stage attempts */

      /* if h changed and a new prediction is needed, break out of stage loop */
      if (kflag == PREDICT_AGAIN)  break;

    } /* loop over stages */

    /* if h has changed due to convergence failure and a new 
       prediction is needed, continue to next attempt at step */
    if (kflag == PREDICT_AGAIN)  continue;
      
    /* compute time-evolved solution and error estimate */
    dsm = ARKComputeSolutions(ark_mem);

    /* Perform time accuracy error test */
    eflag = ARKDoErrorTest2(ark_mem, &nflag, saved_t, &nef, dsm);
	
    /* Go back in loop if we need to predict again (nflag=PREV_ERR_FAIL) */
    if (eflag == TRY_AGAIN)  continue;
	
    /* Return if error test failed and recovery not possible. */
    if (eflag != ARK_SUCCESS)  return(eflag);
	
    /* Error test passed (eflag=ARK_SUCCESS), break from loop */
    break;

  } /* loop over step attempts */

  /* Nonlinear system solves and error test were all successful.
     Update data, and consider change of step. */
  ARKCompleteStep(ark_mem); 
  ARKPrepareNextStep(ark_mem, dsm); 

  /* Reset growth factor for subsequent time step */
  ark_mem->ark_etamax = ETAMX2;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKAdjustParams

 This routine is called when a change in step size was decided upon,
 and it handles the required adjustments to the history array zn.
 If there is to be a change in order, we call ARKAdjustOrder and reset
 q, L = q+1, and qwait.  Then in any case, we call ARKRescale, which
 resets h and rescales the Nordsieck array.
---------------------------------------------------------------*/
static void ARKAdjustParams(ARKodeMem ark_mem)
{
  if (ark_mem->ark_qprime != ark_mem->ark_q) {
    ARKAdjustOrder(ark_mem, ark_mem->ark_qprime-ark_mem->ark_q);
    ark_mem->ark_q = ark_mem->ark_qprime;
    ark_mem->ark_L = ark_mem->ark_q+1;
    ark_mem->ark_qwait = ark_mem->ark_L;
  }
  ARKRescale(ark_mem);
}


/*---------------------------------------------------------------
 ARKAdjustParams2

 This routine is called when a change in step size was decided upon,
 and it resets h and updates the related statistics.
---------------------------------------------------------------*/
static void ARKAdjustParams2(ARKodeMem ark_mem)
{
  ark_mem->ark_h = ark_mem->ark_hscale * ark_mem->ark_eta;
  ark_mem->ark_next_h = ark_mem->ark_h;
  ark_mem->ark_hscale = ark_mem->ark_h;
}


/*---------------------------------------------------------------
 ARKAdjustOrder

 This routine is a high level routine which handles an order
 change by an amount deltaq (= +1 or -1). If a decrease in order
 is requested and q==2, then the routine returns immediately.
 Otherwise ARKAdjustBDF is called to handle the order change.
---------------------------------------------------------------*/
static void ARKAdjustOrder(ARKodeMem ark_mem, int deltaq)
{
  if ((ark_mem->ark_q==2) && (deltaq != 1)) return;
  ARKAdjustBDF(ark_mem, deltaq);
}


/*---------------------------------------------------------------
 ARKAdjustBDF

 This is a high level routine which handles adjustments to the
 history array on a change of order by deltaq.  ARKAdjustBDF 
 calls ARKIncreaseBDF if deltaq = +1 and ARKDecreaseBDF if 
 deltaq = -1 to do the actual work.
---------------------------------------------------------------*/
static void ARKAdjustBDF(ARKodeMem ark_mem, int deltaq)
{
  switch(deltaq) {
  case 1 : 
    ARKIncreaseBDF(ark_mem);
    return;
  case -1: 
    ARKDecreaseBDF(ark_mem);
    return;
  }
}


/*---------------------------------------------------------------
 ARKIncreaseBDF

 This routine adjusts the history array on an increase in the 
 order q.  A new column zn[q+1] is set equal to a multiple of 
 the saved vector (= acor) in zn[indx_acor].  Then each zn[j] is
 adjusted by a multiple of zn[q+1].  The coefficients in the 
 adjustment are the coefficients of the polynomial 
 x*x*(x+xi_1)*...*(x+xi_j), where xi_j = [t_n - t_(n-j)]/h. 
---------------------------------------------------------------*/
static void ARKIncreaseBDF(ARKodeMem ark_mem)
{
  realtype alpha0, alpha1, prod, xi, xiold, hsum, A1;
  int i, j;
  
  for (i=0; i <= ark_mem->ark_qmax; i++) ark_mem->ark_l[i] = ZERO;
  ark_mem->ark_l[2] = alpha1 = prod = xiold = ONE;
  alpha0 = -ONE;
  hsum = ark_mem->ark_hscale;
  if (ark_mem->ark_q > 1) {
    for (j=1; j < ark_mem->ark_q; j++) {
      hsum += ark_mem->ark_tau[j+1];
      xi = hsum / ark_mem->ark_hscale;
      prod *= xi;
      alpha0 -= ONE / (j+1);
      alpha1 += ONE / xi;
      for (i=j+2; i >= 2; i--) 
	ark_mem->ark_l[i] = ark_mem->ark_l[i]*xiold + ark_mem->ark_l[i-1];
      xiold = xi;
    }
  }
  A1 = (-alpha0 - alpha1) / prod;
  N_VScale(A1, ark_mem->ark_zn[ark_mem->ark_indx_acor], 
	   ark_mem->ark_zn[ark_mem->ark_L]);
  for (j=2; j <= ark_mem->ark_q; j++) {
    N_VLinearSum(ark_mem->ark_l[j], ark_mem->ark_zn[ark_mem->ark_L], 
		 ONE, ark_mem->ark_zn[j], ark_mem->ark_zn[j]);
  }  
}


/*---------------------------------------------------------------
 ARKDecreaseBDF

 This routine adjusts the history array on a decrease in the 
 order q.  Each zn[j] is adjusted by a multiple of zn[q].  The 
 coefficients in the adjustment are the coefficients of the 
 polynomial
   x*x*(x+xi_1)*...*(x+xi_j), where xi_j = [t_n - t_(n-j)]/h. 
---------------------------------------------------------------*/
static void ARKDecreaseBDF(ARKodeMem ark_mem)
{
  realtype hsum, xi;
  int i, j;
  
  for (i=0; i <= ark_mem->ark_qmax; i++) ark_mem->ark_l[i] = ZERO;
  ark_mem->ark_l[2] = ONE;
  hsum = ZERO;
  for(j=1; j <= ark_mem->ark_q-2; j++) {
    hsum += ark_mem->ark_tau[j];
    xi = hsum /ark_mem->ark_hscale;
    for (i=j+2; i >= 2; i--) 
      ark_mem->ark_l[i] = ark_mem->ark_l[i]*xi + ark_mem->ark_l[i-1];
  }
  
  for(j=2; j < ark_mem->ark_q; j++)
    N_VLinearSum(-ark_mem->ark_l[j], ark_mem->ark_zn[ark_mem->ark_q], 
		 ONE, ark_mem->ark_zn[j], ark_mem->ark_zn[j]);
}


/*---------------------------------------------------------------
 ARKRescale

 This routine rescales the Nordsieck array by multiplying the
 jth column zn[j] by eta^j, j = 1, ..., q.  Then the value of
 h is rescaled by eta, and hscale is reset to h.
---------------------------------------------------------------*/
static void ARKRescale(ARKodeMem ark_mem)
{
  int j;
  realtype factor;
  
  factor = ark_mem->ark_eta;
  for (j=1; j <= ark_mem->ark_q; j++) {
    N_VScale(factor, ark_mem->ark_zn[j], ark_mem->ark_zn[j]);
    factor *= ark_mem->ark_eta;
  }
  ark_mem->ark_h = ark_mem->ark_hscale * ark_mem->ark_eta;
  ark_mem->ark_next_h = ark_mem->ark_h;
  ark_mem->ark_hscale = ark_mem->ark_h;
}


/*---------------------------------------------------------------
 ARKPredict

 This routine computes the predicted array z_n(0), which is 
 overwritten on zn.  The prediction of zn is done by repeated 
 additions.  
---------------------------------------------------------------*/
static void ARKPredict(ARKodeMem ark_mem)
{
  int j, k;
  for (k = 1; k <= ark_mem->ark_q; k++)
    for (j = ark_mem->ark_q; j >= k; j--) 
      N_VLinearSum(ONE, ark_mem->ark_zn[j-1], ONE, 
		   ark_mem->ark_zn[j], ark_mem->ark_zn[j-1]); 
}


/*---------------------------------------------------------------
 ARKPredict2

 This routine computes the prediction for a specific internal 
 stage solution, storing the result in ark_mem->ark_zn[0].  The 
 prediction is done using the dense output routine in 
 extrapolation mode, hence stages "far" from the previous time 
 interval are predicted using lower order polynomials than the
 "nearby" stages.
---------------------------------------------------------------*/
static void ARKPredict2(ARKodeMem ark_mem, int istage)
{
  int retval, ord;
  realtype tau, tau_tol = 1.5;
  N_Vector yguess = ark_mem->ark_zn[0];

#define PRED0

  /* if this is the first step, use initial condition as guess */
  if (ark_mem->ark_nst == 0) {
    N_VScale(ONE, ark_mem->ark_ynew, yguess);
    return;
  }

#ifdef PRED0
  /***** Trivial Predictor *****/
  N_VScale(ONE, ark_mem->ark_ynew, yguess);
#endif
  
#ifdef PRED1
  /***** Dense Output Predictor 1 -- all to max order *****/
  tau = ONE + ark_mem->ark_c[istage]*ark_mem->ark_h/ark_mem->ark_hold;
  retval = ARKDenseEval(ark_mem, tau, 0, ark_mem->ark_dense_q, yguess);

  /* if predictor fails, use trivial prediction (shouldn't happen) */
  if (retval != ARK_SUCCESS) 
    N_VScale(ONE, ark_mem->ark_ynew, yguess);
#endif

#ifdef PRED2
  /***** Dense Output Predictor 2 -- decrease order w/ increasing stage *****/
  tau = ONE + ark_mem->ark_c[istage]*ark_mem->ark_h/ark_mem->ark_hold;
  ord = MAX(ark_mem->ark_dense_q - istage, 1);
  retval = ARKDenseEval(ark_mem, tau, 0, ord, yguess);

  /* if predictor fails, use trivial prediction (shouldn't happen) */
  if (retval != ARK_SUCCESS) 
    N_VScale(ONE, ark_mem->ark_ynew, yguess);
#endif

#ifdef PRED3
  /***** Dense Output Predictor for stages "close" to previous step, 
	 existing solution for subsequent stages *****/
  tau = ONE + ark_mem->ark_c[istage]*ark_mem->ark_h/ark_mem->ark_hold;
  if (tau < tau_tol) {
    ord = MAX(ark_mem->ark_dense_q - istage, 1);
    retval = ARKDenseEval(ark_mem, tau, 0, ord, yguess);
  } else {
    N_VScale(ONE, ark_mem->ark_y, yguess);
  }    

  /* if predictor fails, use trivial prediction (shouldn't happen) */
  if (retval != ARK_SUCCESS) 
    N_VScale(ONE, ark_mem->ark_ynew, yguess);
#endif

}


/*---------------------------------------------------------------
 ARKSet

 This routine is a high level routine which calls ARKSetBDF to 
 set the polynomial l, the test quantity array tq, 
 and the related variables  rl1, gamma, and gamrat.
---------------------------------------------------------------*/
static void ARKSet(ARKodeMem ark_mem)
{
  ARKSetBDF(ark_mem);
  ark_mem->ark_rl1 = ONE / ark_mem->ark_l[1];
  ark_mem->ark_gamma = ark_mem->ark_h * ark_mem->ark_rl1;
  if (ark_mem->ark_nst == 0) ark_mem->ark_gammap = ark_mem->ark_gamma;
  ark_mem->ark_gamrat = (ark_mem->ark_nst > 0) ? 
    ark_mem->ark_gamma / ark_mem->ark_gammap : ONE;  /* protect x / x != 1.0 */
}


/*---------------------------------------------------------------
 ARKSet2

 This routine sets up the RHS data for a given RK stage 
 (currently stored in ark_Fe[istage]), as well as the time-step 
 and method-related factors gamma, gammap and gamrat.
---------------------------------------------------------------*/
static void ARKSet2(ARKodeMem ark_mem)
{
  /* local data */
  realtype hA;
  int j, i = ark_mem->ark_istage;
  N_Vector RHS = ark_mem->ark_Fe[i];

  /* Initialize RHS to yn */
  N_VScale(ONE, ark_mem->ark_ynew, RHS);

  /* Iterate over each prior stage updating RHS */
  /*    Explicit pieces */
  if (!ark_mem->ark_implicit)
    for (j=0; j<i; j++) {
      hA = ark_mem->ark_h * ark_mem->ark_Ae[i][j];
      N_VLinearSum(hA, ark_mem->ark_Fe[j], ONE, RHS, RHS);
    }
  /*    Implicit pieces */
  if (!ark_mem->ark_explicit)
    for (j=0; j<i; j++) {
      hA = ark_mem->ark_h * ark_mem->ark_Ai[i][j];
      N_VLinearSum(hA, ark_mem->ark_Fi[j], ONE, RHS, RHS);
    }

  /* Update gamma */
  ark_mem->ark_gamma = ark_mem->ark_h * ark_mem->ark_Ai[i][i];
  if (ark_mem->ark_nst == 0)  
    ark_mem->ark_gammap = ark_mem->ark_gamma;
  ark_mem->ark_gamrat = (ark_mem->ark_nst > 0) ? 
    ark_mem->ark_gamma / ark_mem->ark_gammap : ONE;  /* protect x/x != 1.0 */
}


/*---------------------------------------------------------------
 ARKSetBDF

 This routine computes the coefficients l and tq.  ARKSetBDF 
 calls ARKSetTqBDF to set the test quantity array tq. 
 
 The components of the array l are the coefficients of a
 polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
                                 q-1
 Lambda(x) = (1 + x / xi*_q) * PRODUCT (1 + x / xi_i) , where
                                 i=1
  xi_i = [t_n - t_(n-i)] / h.

 The array tq is set to test quantities used in the convergence
 test, the error test, and the selection of h at a new order. 
---------------------------------------------------------------*/
static void ARKSetBDF(ARKodeMem ark_mem)
{
  realtype alpha0, alpha0_hat, xi_inv, xistar_inv, hsum;
  int i,j;
  
  ark_mem->ark_l[0] = ark_mem->ark_l[1] = xi_inv = xistar_inv = ONE;
  for (i=2; i <= ark_mem->ark_q; i++) ark_mem->ark_l[i] = ZERO;
  alpha0 = alpha0_hat = -ONE;
  hsum = ark_mem->ark_h;
  if (ark_mem->ark_q > 1) {
    for (j=2; j < ark_mem->ark_q; j++) {
      hsum += ark_mem->ark_tau[j-1];
      xi_inv = ark_mem->ark_h / hsum;
      alpha0 -= ONE / j;
      for(i=j; i >= 1; i--) 
	ark_mem->ark_l[i] += ark_mem->ark_l[i-1]*xi_inv;
      /* The l[i] are coefficients of product(1 to j) (1 + x/xi_i) */
    }
    
    /* j = q */
    alpha0 -= ONE / ark_mem->ark_q;
    xistar_inv = -ark_mem->ark_l[1] - alpha0;
    hsum += ark_mem->ark_tau[ark_mem->ark_q-1];
    xi_inv = ark_mem->ark_h / hsum;
    alpha0_hat = -ark_mem->ark_l[1] - xi_inv;
    for (i=ark_mem->ark_q; i >= 1; i--) 
      ark_mem->ark_l[i] += ark_mem->ark_l[i-1]*xistar_inv;
  }

  ARKSetTqBDF(ark_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv);
}


/*---------------------------------------------------------------
 ARKSetTqBDF

 This routine sets the test quantity array tq, that is loaded 
 with constants used in the control of estimated local errors 
 and in the nonlinear convergence test.  Specifically, while
 running at order q, the components of tq are as follows:
   tq[1] = a coefficient used to get the est. local error at order q-1
   tq[2] = a coefficient used to get the est. local error at order q
   tq[3] = a coefficient used to get the est. local error at order q+1
   tq[4] = constant used in nonlinear iteration convergence test
   tq[5] = coefficient used to get the order q+2 derivative vector used in
           the est. local error at order q+1
---------------------------------------------------------------*/
static void ARKSetTqBDF(ARKodeMem ark_mem, realtype hsum, realtype alpha0,
			realtype alpha0_hat, realtype xi_inv, realtype xistar_inv)
{
  realtype A1, A2, A3, A4, A5, A6;
  realtype C, Cpinv, Cppinv;
  
  A1 = ONE - alpha0_hat + alpha0;
  A2 = ONE + ark_mem->ark_q * A1;
  ark_mem->ark_tq[2] = ABS(A1 / (alpha0 * A2));
  ark_mem->ark_tq[5] = ABS(A2 * xistar_inv / 
			   (ark_mem->ark_l[ark_mem->ark_q] * xi_inv));
  if (ark_mem->ark_qwait == 1) {
    if (ark_mem->ark_q > 1) {
      C = xistar_inv / ark_mem->ark_l[ark_mem->ark_q];
      A3 = alpha0 + ONE / ark_mem->ark_q;
      A4 = alpha0_hat + xi_inv;
      Cpinv = (ONE - A4 + A3) / A3;
      ark_mem->ark_tq[1] = ABS(C * Cpinv);
    }
    else ark_mem->ark_tq[1] = ONE;
    hsum += ark_mem->ark_tau[ark_mem->ark_q];
    xi_inv = ark_mem->ark_h / hsum;
    A5 = alpha0 - (ONE / (ark_mem->ark_q+1));
    A6 = alpha0_hat - xi_inv;
    Cppinv = (ONE - A6 + A5) / A2;
    ark_mem->ark_tq[3] = ABS(Cppinv / 
			     (xi_inv * (ark_mem->ark_q+2) * A5));
  }
  ark_mem->ark_tq[4] = ark_mem->ark_nlscoef / ark_mem->ark_tq[2];
}


/*---------------------------------------------------------------
 ARKNls

 This routine attempts to solve the nonlinear system associated
 with a single implicit step of the linear multistep method.
 It calls ARKNlsNewton to do the work.

 Upon a successful solve, the solution is held in ark_mem->ark_y.
---------------------------------------------------------------*/
static int ARKNls(ARKodeMem ark_mem, int nflag)
{
  return(ARKNlsNewton(ark_mem, nflag));
}


/*---------------------------------------------------------------
 ARKNlsResid

 This routine evaluates the negative nonlinear residual.  It 
 assumes that any data from previous time steps/stages is 
 contained in ark_mem, and merely combines this old data with the
 current RHS vector, fy = f(t,y), to compute the nonlinear 
 residual r.

 The array ark_mem->ark_zn[1] contains both old data and the
 initial solution guess, so the current value of y only requires 
 addition of acor.

 Possible return values:  ARK_SUCCESS  or  ARK_RHSFUNC_FAIL
---------------------------------------------------------------*/
static int ARKNlsResid(ARKodeMem ark_mem, N_Vector acor, 
		       N_Vector fy, N_Vector r)
{
  /* Evaluate the negative residual of the nonlinear system */
  N_VLinearSum(ark_mem->ark_rl1, ark_mem->ark_zn[1], ONE, acor, r);
  N_VLinearSum(ark_mem->ark_gamma, fy, -ONE, r, r);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKNlsResid2

 This routine evaluates the negative nonlinear residual for the 
 additive Runge-Kutta method.  It assumes that any data from 
 previous time steps/stages is contained in ark_mem, and merely 
 combines this old data with the current RHS vector, fy = f(t,y), 
 to compute the nonlinear residual r.

 Possible return values:  ARK_SUCCESS  or  ARK_RHSFUNC_FAIL
---------------------------------------------------------------*/
static int ARKNlsResid2(ARKodeMem ark_mem, N_Vector y, 
			N_Vector fy, N_Vector r)
{
  /* At the ith stage, we compute 
       r = -zi + yn + h*sum_{j=0}^{i-1} Ae(i,j)*Fe(j) + h*sum_{j=0}^{i} Ai(i,j)*Fi(j)
       r = -zi + gamma*Fi(zi) + yn + data
     where
       zi is stored in y
       Fi(zi) is stored in fy
       (yn+data) is stored in ark_mem->ark_Fe[i],
     so we really just compute
       r = -y + gamma*fy + ark_mem->ark_Fe[i]
   */
  N_VLinearSum(-ONE, y, ONE, ark_mem->ark_Fe[ark_mem->ark_istage], r);
  N_VLinearSum(ark_mem->ark_gamma, fy, ONE, r, r);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKNlsNewton

 This routine handles the Newton iteration for the BDF method. It 
 calls lsetup if indicated, performs a modified Newton iteration 
 (using a fixed Jacobian / precond), and retries a failed attempt 
 at Newton iteration if that is indicated. 

 Upon entry, the predicted solution is held in ark_mem->ark_zn[0]; 
 this array is never changed throughout this routine.  If an 
 initial attempt at solving the nonlinear system fail (e.g. due to 
 a stale Jacobian), this allows for new attempts that first revert 
 ark_mem->ark_y back to this initial guess before trying again. 

 Upon a successful solve, the solution is held in ark_mem->ark_y.

 Possible return values:

   ARK_SUCCESS       ---> continue with error test

   ARK_RHSFUNC_FAIL  -+  
   ARK_LSETUP_FAIL    |-> halt the integration 
   ARK_LSOLVE_FAIL   -+

   CONV_FAIL        -+
   RHSFUNC_RECVR    -+-> predict again or stop if too many
---------------------------------------------------------------*/
static int ARKNlsNewton(ARKodeMem ark_mem, int nflag)
{
  N_Vector vtemp1, vtemp2, vtemp3, ycur, b;
  int convfail, retval, ier, m;
  booleantype callSetup;
  realtype del, delp, dcon;
  
  vtemp1 = ark_mem->ark_acor;  /* rename acor as vtemp1 for readability  */
  vtemp2 = ark_mem->ark_y;     /* rename y as vtemp2 for readability     */
  vtemp3 = ark_mem->ark_tempv; /* rename tempv as vtemp3 for readability */
  b      = ark_mem->ark_tempv; /* also rename tempv as b for readability */
  ycur   = ark_mem->ark_zn[0]; /* rename zn[0] as ycur for readability   */

  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ?
    ARK_NO_FAILURES : ARK_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (ark_mem->ark_setupNonNull) {      
    callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
      (ark_mem->ark_nst == 0) || 
      (ark_mem->ark_nst >= ark_mem->ark_nstlp + MSBP) || 
      (ABS(ark_mem->ark_gamrat-ONE) > DGMAX);
  } else {  
    ark_mem->ark_crate = ONE;
    callSetup = FALSE;
  }
  
  /* Looping point for attempts at solution of the nonlinear system:
       Evaluate f at predicted y, store result in ark_mem->ark_ftemp.
       Call lsetup if indicated, setting statistics and gamma factors.
       Zero out the correction array (ark_mem->ark_acor).
       Copy the predicted y (ycur) into the output (ark_mem->ark_y).
       Performs the modified Newton iteration using the existing lsetup.
       Repeat process if a recoverable failure occurred (convergence
	  failure with stale Jacobian). */
  for(;;) {

    if (!ark_mem->ark_explicit) {
      retval = ark_mem->ark_fi(ark_mem->ark_tn, ycur, 
			       ark_mem->ark_ftemp, ark_mem->ark_user_data);
      ark_mem->ark_nfi++; 
      if (retval < 0) return(ARK_RHSFUNC_FAIL);
      if (retval > 0) return(RHSFUNC_RECVR);
    }

    if (callSetup) {
      ier = ark_mem->ark_lsetup(ark_mem, convfail, ycur, 
				ark_mem->ark_ftemp, &ark_mem->ark_jcur, 
				vtemp1, vtemp2, vtemp3);
      ark_mem->ark_nsetups++;
      callSetup = FALSE;
      ark_mem->ark_gamrat = ark_mem->ark_crate = ONE; 
      ark_mem->ark_gammap = ark_mem->ark_gamma;
      ark_mem->ark_nstlp  = ark_mem->ark_nst;

      /* Return if lsetup failed */
      if (ier < 0) return(ARK_LSETUP_FAIL);
      if (ier > 0) return(CONV_FAIL);
    }

    /* Set acor to zero and load prediction into y vector */
    N_VConst(ZERO, ark_mem->ark_acor);
    N_VScale(ONE, ycur, ark_mem->ark_y);


    /***********************************
     Do the modified Newton iteration:
     - Reuses a single call to ark_mem->ark_lsetup (called by the enclosing loop). 
     - Upon entry, the predicted solution is held in ark_mem->ark_zn[0].
     - Here, ark_mem->ark_acor accumulates the difference between the predictor 
       and the final solution to the time step.
     - Upon a successful solve, the solution is held in ark_mem->ark_y. 
    */

    /* Initialize temporary variables for use in iteration */
    ark_mem->ark_mnewt = m = 0;
    del = delp = ZERO;

    /* Looping point for Newton iteration */
    for(;;) {

      /* Evaluate the nonlinear system residual, put result into b */
      retval = ARKNlsResid(ark_mem, ark_mem->ark_acor, 
			   ark_mem->ark_ftemp, b);
      if (retval != ARK_SUCCESS) {
	ier = ARK_RHSFUNC_FAIL;
	break;
      }

      /* Call the lsolve function */
      retval = ark_mem->ark_lsolve(ark_mem, b, ark_mem->ark_ewt, 
				   ark_mem->ark_y, ark_mem->ark_ftemp); 
      ark_mem->ark_nni++;
    
      if (retval < 0) {
	ier = ARK_LSOLVE_FAIL;
	break;
      }
    
      /* If lsolve had a recoverable failure and Jacobian data is
	 not current, signal to try the solution again */
      if (retval > 0) { 
	if ((!ark_mem->ark_jcur) && (ark_mem->ark_setupNonNull)) 
	  ier = TRY_AGAIN;
	else 
	  ier = CONV_FAIL;
	break;
      }

      /* Get WRMS norm of correction; add correction to acor and y */
      del = N_VWrmsNorm(b, ark_mem->ark_ewt);
      N_VLinearSum(ONE, ark_mem->ark_acor, ONE, b, ark_mem->ark_acor);
      N_VLinearSum(ONE, ycur, ONE, ark_mem->ark_acor, ark_mem->ark_y);
    
      /* Test for convergence.  If m > 0, an estimate of the convergence
	 rate constant is stored in crate, and used in the test */
      if (m > 0) {
	ark_mem->ark_crate = MAX(CRDOWN * ark_mem->ark_crate, del/delp);
      }
      ark_mem->ark_eLTE = ark_mem->ark_nlscoef / ark_mem->ark_tq[2]; /* BDF */
      /* if (ark_mem->ark_crate < ONE)                                  /\* RK *\/ */
      /* 	ark_mem->ark_eLTE = ark_mem->ark_nlscoef * (ONE - ark_mem->ark_crate); */
      /* else */
      /* 	ark_mem->ark_eLTE = ark_mem->ark_nlscoef * RCONST(0.1); */
      dcon = del * MIN(ONE, ark_mem->ark_crate) / ark_mem->ark_eLTE;
    
      if (dcon <= ONE) {
	ark_mem->ark_acnrm = (m==0) ? 
	  del : N_VWrmsNorm(ark_mem->ark_acor, ark_mem->ark_ewt);
	ark_mem->ark_jcur = FALSE;
	ier = ARK_SUCCESS;
	break;
      }

      /* update Newton iteration counter */
      ark_mem->ark_mnewt = ++m;
    
      /* Stop at maxcor iterations or if iteration seems to be diverging.
	 If still not converged and Jacobian data is not current, signal 
	 to try the solution again */
      if ((m == ark_mem->ark_maxcor) || ((m >= 2) && (del > RDIV*delp))) {
	if ((!ark_mem->ark_jcur) && (ark_mem->ark_setupNonNull)) 
	  ier = TRY_AGAIN;
	else
	  ier = CONV_FAIL;
	break;
      }
    
      /* Save norm of correction, evaluate fi, and loop again */
      delp = del;
      if (!ark_mem->ark_explicit) {
	retval = ark_mem->ark_fi(ark_mem->ark_tn, ark_mem->ark_y, 
				 ark_mem->ark_ftemp, ark_mem->ark_user_data);
	ark_mem->ark_nfi++;
	if (retval < 0) {
	  ier = ARK_RHSFUNC_FAIL;
	  break;
	}
	if (retval > 0) {
	  if ((!ark_mem->ark_jcur) && (ark_mem->ark_setupNonNull)) 
	    ier = TRY_AGAIN;
	  else
	    ier = RHSFUNC_RECVR;
	  break;
	}
      }
    } 
    /* end modified Newton iteration */
    /*********************************/
    
    /* If there is a convergence failure and the Jacobian-related 
       data appears not to be current, loop again with a call to lsetup
       in which convfail=ARK_FAIL_BAD_J.  Otherwise return. */
    if (ier != TRY_AGAIN) return(ier);
    
    callSetup = TRUE;
    convfail = ARK_FAIL_BAD_J;
  }
}


/*---------------------------------------------------------------
 ARKLs

 This routine attempts to solve the linear system associated
 with a single implicit step of the ARK method.  This should 
 only be called if the user has specified that the implicit 
 problem is linear.  In this routine, we assume that the problem 
 depends linearly on the solution, and that linear dependence 
 does not change as a function of time.  It therefore calls
 lsetup on only changes to gamma.  It then calls the 
 user-specified linear solver (with a tight tolerance) to compute 
 the time-evolved solution.  

 Upon entry, the predicted solution is held in ark_mem->ark_zn[0].

 Upon a successful solve, the solution is held in ark_mem->ark_y.

 Possible return values:

   ARK_SUCCESS       ---> continue with error test

   ARK_RHSFUNC_FAIL  -+  
   ARK_LSETUP_FAIL    |-> halt the integration 
   ARK_LSOLVE_FAIL   -+

   RHSFUNC_RECVR     --> predict again or stop if too many
---------------------------------------------------------------*/
static int ARKLs(ARKodeMem ark_mem, int nflag)
{
  N_Vector vtemp1, vtemp2, vtemp3, ycur, b;
  int convfail, retval, ier;
  booleantype callSetup;
  
  vtemp1 = ark_mem->ark_acor;  /* rename acor as vtemp1 for readability  */
  vtemp2 = ark_mem->ark_y;     /* rename y as vtemp2 for readability     */
  vtemp3 = ark_mem->ark_tempv; /* rename tempv as vtemp3 for readability */
  b      = ark_mem->ark_tempv; /* also rename tempv as b for readability */
  ycur   = ark_mem->ark_zn[0]; /* rename zn[0] as ycur for readability   */

  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ?
    ARK_NO_FAILURES : ARK_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (ark_mem->ark_setupNonNull) {      
    callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
      (ark_mem->ark_nst == 0) || 
      (ABS(ark_mem->ark_gamrat-ONE) > TINY);
  } else {  
    ark_mem->ark_crate = ONE;
    callSetup = FALSE;
  }
  
  /* Evaluate fi at predicted y, store result in ark_mem->ark_ftemp, and 
     call lsetup if indicated, setting statistics and gamma factors. */
  if (!ark_mem->ark_explicit) {
    retval = ark_mem->ark_fi(ark_mem->ark_tn, ycur, 
			     ark_mem->ark_ftemp, ark_mem->ark_user_data);
    ark_mem->ark_nfi++; 
    if (retval < 0) return(ARK_RHSFUNC_FAIL);
    if (retval > 0) return(RHSFUNC_RECVR);
  }

  if (callSetup) {
    ier = ark_mem->ark_lsetup(ark_mem, convfail, ycur, 
			      ark_mem->ark_ftemp, &ark_mem->ark_jcur, 
			      vtemp1, vtemp2, vtemp3);
    ark_mem->ark_nsetups++;
    callSetup = FALSE;
    ark_mem->ark_gamrat = ark_mem->ark_crate = ONE; 
    ark_mem->ark_gammap = ark_mem->ark_gamma;
    
    /* Return if lsetup failed */
    if (ier < 0) return(ARK_LSETUP_FAIL);
    if (ier > 0) return(CONV_FAIL);
  }

  /* Set acor to zero and load prediction into y vector */
  N_VConst(ZERO, ark_mem->ark_acor);
  N_VScale(ONE, ycur, ark_mem->ark_y);


  /* Initialize temporary variable for communication with linear solvers */
  ark_mem->ark_mnewt = 0;
  
  /* Evaluate the linear system residual, put result into b */
  retval = ARKNlsResid(ark_mem, ark_mem->ark_acor, 
		       ark_mem->ark_ftemp, b);
  if (retval != ARK_SUCCESS)  return(ARK_RHSFUNC_FAIL);
  
  /* Call the lsolve function */
  retval = ark_mem->ark_lsolve(ark_mem, b, ark_mem->ark_ewt, 
			       ark_mem->ark_y, ark_mem->ark_ftemp); 
  ark_mem->ark_nni++;
  if (retval < 0)  return(ARK_LSOLVE_FAIL);
  
  /* Get WRMS norm of correction; add correction to acor and y */
  ark_mem->ark_acnrm = N_VWrmsNorm(b, ark_mem->ark_ewt);
  N_VLinearSum(ONE, ark_mem->ark_acor, ONE, b, ark_mem->ark_acor);
  N_VLinearSum(ONE, ycur, ONE, ark_mem->ark_acor, ark_mem->ark_y);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKHandleFlag

 This routine takes action on the return value nflag = *nflagPtr
 returned by ARKNls, as follows:

 If ARKNls succeeded in solving the nonlinear system, then
 ARKHandleNFlag returns the constant DO_ERROR_TEST, which tells ARKStep
 to perform the error test.

 If the nonlinear system was not solved successfully, then ncfn and
 ncf = *ncfPtr are incremented and Nordsieck array zn is restored.

 If the solution of the nonlinear system failed due to an
 unrecoverable failure by setup, we return the value ARK_LSETUP_FAIL.
 
 If it failed due to an unrecoverable failure in solve, then we return
 the value ARK_LSOLVE_FAIL.

 If it failed due to an unrecoverable failure in rhs, then we return
 the value ARK_RHSFUNC_FAIL.

 Otherwise, a recoverable failure occurred when solving the 
 nonlinear system (ARKNls returned nflag == CONV_FAIL or RHSFUNC_RECVR). 
 In this case, if ncf is now equal to maxncf or |h| = hmin, 
 we return the value ARK_CONV_FAILURE (if nflag=CONV_FAIL) or
 ARK_REPTD_RHSFUNC_ERR (if nflag=RHSFUNC_RECVR).
 If not, we set *nflagPtr = PREV_CONV_FAIL and return the value
 PREDICT_AGAIN, telling ARKStep to reattempt the step.
---------------------------------------------------------------*/
static int ARKHandleNFlag(ARKodeMem ark_mem, int *nflagPtr, 
			  realtype saved_t, int *ncfPtr)
{
  int nflag;
  
  nflag = *nflagPtr;
  
  if (nflag == ARK_SUCCESS) return(DO_ERROR_TEST);

  /* The nonlinear soln. failed; increment ncfn and restore tn, zn */
  ark_mem->ark_ncfn++;
  ark_mem->ark_tn = saved_t;
  ARKRestore(ark_mem);
  
  /* Return if lsetup, lsolve, or rhs failed unrecoverably */
  if (nflag == ARK_LSETUP_FAIL)  return(ARK_LSETUP_FAIL);
  if (nflag == ARK_LSOLVE_FAIL)  return(ARK_LSOLVE_FAIL);
  if (nflag == ARK_RHSFUNC_FAIL) return(ARK_RHSFUNC_FAIL);
  
  /* At this point, nflag = CONV_FAIL or RHSFUNC_RECVR; increment ncf */
  (*ncfPtr)++;
  ark_mem->ark_etamax = ONE;

  /* If we had maxncf failures or |h| = hmin, 
     return ARK_CONV_FAILURE or ARK_REPTD_RHSFUNC_ERR. */
  if ((ABS(ark_mem->ark_h) <= ark_mem->ark_hmin*ONEPSM) || 
      (*ncfPtr == ark_mem->ark_maxncf)) {
    if (nflag == CONV_FAIL)     return(ARK_CONV_FAILURE);
    if (nflag == RHSFUNC_RECVR) return(ARK_REPTD_RHSFUNC_ERR);    
  }

  /* Reduce step size; return to reattempt the step */
  ark_mem->ark_eta = MAX(ETACF, ark_mem->ark_hmin / ABS(ark_mem->ark_h));
  *nflagPtr = PREV_CONV_FAIL;
  ARKRescale(ark_mem);
  ark_mem->ark_nst_con++;

  return(PREDICT_AGAIN);
}


/*---------------------------------------------------------------
 ARKRestore

 This routine undoes the prediction.  After execution of 
 ARKRestore, the Nordsieck array zn has the same values as 
 before the call to ARKPredict. 

 *** When I switch from BDF to ARK, remove ***
 *** this routine and all calls to it.     ***
---------------------------------------------------------------*/
static void ARKRestore(ARKodeMem ark_mem)
{
  int j, k;
  for (k = 1; k <= ark_mem->ark_q; k++)
    for (j = ark_mem->ark_q; j >= k; j--)
      N_VLinearSum(ONE, ark_mem->ark_zn[j-1], -ONE, 
		   ark_mem->ark_zn[j], ark_mem->ark_zn[j-1]);
}


/*---------------------------------------------------------------
 ARKComputeSolutions

 This routine calculates the final RK solution using the existing 
 data.  This solution is placed directly in ark_y.  This routine 
 also computes the error estimate ||y-ytilde||_WRMS, where ytilde 
 is the embedded solution, and the norm weights come from ark_ewt.  
 This norm value is returned.  The vector form of this estimated 
 error (y-ytilde) is stored in ark_tempv, in case the calling 
 routine wishes to examine the error locations.

 For now, we assume that the ODE is of the form 
    y' = fe(t,y) + fi(t,y)
 so that both y and ytilde can be formed from existing data.  
 However, once we extend the solver to allow for non-identity
 mass matrices, this routine will perform two additional solves
 to compute the solution and embedding.
---------------------------------------------------------------*/
static realtype ARKComputeSolutions(ARKodeMem ark_mem)
{
  /* local data */
  realtype hb;
  int j;
  N_Vector y    = ark_mem->ark_y;
  N_Vector yerr = ark_mem->ark_tempv;

  /* Initialize solution to yn, error estimate to zero */
  N_VScale(ONE, ark_mem->ark_ynew, y);
  N_VConst(ZERO, yerr);

  /* Iterate over each stage, updating solution and error estimate */
  for (j=0; j<ark_mem->ark_stages; j++) {

    /* solution */
    hb = ark_mem->ark_h * ark_mem->ark_b[j];
    N_VLinearSum(hb, ark_mem->ark_Fe[j], ONE, y, y);
    N_VLinearSum(hb, ark_mem->ark_Fi[j], ONE, y, y);

    /* error */
    hb = ark_mem->ark_h * (ark_mem->ark_b[j] - ark_mem->ark_b2[j]);
    N_VLinearSum(hb, ark_mem->ark_Fe[j], ONE, yerr, yerr);
    N_VLinearSum(hb, ark_mem->ark_Fi[j], ONE, yerr, yerr);

  }

  /* return error norm */
  return(N_VWrmsNorm(yerr, ark_mem->ark_ewt));
  
}


/*---------------------------------------------------------------
 ARKDoErrorTest

 This routine performs the local error test. 
 The weighted local error norm dsm is loaded into *dsmPtr, and 
 the test dsm ?<= 1 is made.

 If the test passes, ARKDoErrorTest returns ARK_SUCCESS. 

 If the test fails, we undo the step just taken (call ARKRestore) and 

   - if maxnef error test failures have occurred or if ABS(h) = hmin,
     we return ARK_ERR_FAILURE.

   - if more than MXNEF1 error test failures have occurred, an order
     reduction is forced. If already at order 1, restart by reloading 
     zn from scratch. If f() fails we return either ARK_RHSFUNC_FAIL
     or ARK_UNREC_RHSFUNC_ERR (no recovery is possible at this stage).

   - otherwise: update time step factor eta based on local error 
     estimate, reduce the method order if possible, and only if it is
     already at order 1 then reduce h.  Then set *nflagPtr to 
     PREV_ERR_FAIL, and return TRY_AGAIN. 
---------------------------------------------------------------*/
static booleantype ARKDoErrorTest(ARKodeMem ark_mem, int *nflagPtr,
				  realtype saved_t, int *nefPtr, realtype *dsmPtr)
{
  realtype dsm;
  int retval;
  
  dsm = ark_mem->ark_acnrm * ark_mem->ark_tq[2];

  /* If est. local error norm dsm passes test, return ARK_SUCCESS */  
  *dsmPtr = dsm; 
  if (dsm <= ONE) return(ARK_SUCCESS);
  
  /* Test failed; increment counters, set nflag, and restore tn and zn */
  (*nefPtr)++;
  ark_mem->ark_netf++;
  *nflagPtr = PREV_ERR_FAIL;
  ark_mem->ark_tn = saved_t;
  ARKRestore(ark_mem);

  /* At maxnef failures or |h| = hmin, return ARK_ERR_FAILURE */
  if ((ABS(ark_mem->ark_h) <= ark_mem->ark_hmin*ONEPSM) || 
      (*nefPtr == ark_mem->ark_maxnef)) return(ARK_ERR_FAILURE);

  /* Set etamax = 1 to prevent step size increase at end of this step */
  ark_mem->ark_etamax = ONE;

  /* Set h ratio eta from dsm, rescale, and return for retry of step */
  if (*nefPtr <= MXNEF1) {
    ark_mem->ark_eta = ONE / (RPowerR(BIAS2*dsm,ONE/ark_mem->ark_L) + ADDON);
    ark_mem->ark_eta = MAX(ETAMIN, MAX(ark_mem->ark_eta, 
				       ark_mem->ark_hmin / ABS(ark_mem->ark_h)));
    if (*nefPtr >= SMALL_NEF) ark_mem->ark_eta = MIN(ark_mem->ark_eta, ETAMXF);
    ARKRescale(ark_mem);
    return(TRY_AGAIN);
  }
  
  /* After MXNEF1 failures, force an order reduction and retry step */
  if (ark_mem->ark_q > 1) {
    ark_mem->ark_eta = MAX(ETAMIN, ark_mem->ark_hmin / ABS(ark_mem->ark_h));
    ARKAdjustOrder(ark_mem,-1);
    ark_mem->ark_L = ark_mem->ark_q;
    ark_mem->ark_q--;
    ark_mem->ark_qwait = ark_mem->ark_L;
    ARKRescale(ark_mem);
    return(TRY_AGAIN);
  }

  /* If already at order 1, restart: reload zn from scratch */
  ark_mem->ark_eta = MAX(ETAMIN, ark_mem->ark_hmin / ABS(ark_mem->ark_h));
  ark_mem->ark_h *= ark_mem->ark_eta;
  ark_mem->ark_next_h = ark_mem->ark_h;
  ark_mem->ark_hscale = ark_mem->ark_h;
  ark_mem->ark_qwait = LONG_WAIT;

  retval = ARKFullRHS(ark_mem, ark_mem->ark_tn, ark_mem->ark_zn[0], 
		      ark_mem->ark_ftemp, ark_mem->ark_tempv);
  if (retval != 0) return(ARK_RHSFUNC_FAIL);
  N_VScale(ark_mem->ark_h, ark_mem->ark_tempv, ark_mem->ark_zn[1]);

  return(TRY_AGAIN);
}


/*---------------------------------------------------------------
 ARKDoErrorTest2

 This routine performs the local error test for the ARK method. 
 The weighted local error norm dsm is passed in, and 
 the test dsm ?<= 1 is made.

 If the test passes, ARKDoErrorTest returns ARK_SUCCESS. 

 If the test fails, we revert to the last successful solution 
 time, and:
   - if maxnef error test failures have occurred or if 
     ABS(h) = hmin, we return ARK_ERR_FAILURE.
   - otherwise: update time step factor eta based on local error 
     estimate and reduce h.  Then set *nflagPtr to PREV_ERR_FAIL, 
     and return TRY_AGAIN. 
---------------------------------------------------------------*/
static booleantype ARKDoErrorTest2(ARKodeMem ark_mem, int *nflagPtr,
				   realtype saved_t, int *nefPtr, realtype dsm)
{
  /* If est. local error norm dsm passes test, return ARK_SUCCESS */  
  if (dsm <= ONE) return(ARK_SUCCESS);
  
  /* Test failed; increment counters, set nflag, and restore tn and zn */
  (*nefPtr)++;
  ark_mem->ark_netf++;
  *nflagPtr = PREV_ERR_FAIL;
  ark_mem->ark_tn = saved_t;

  /* At maxnef failures or |h| = hmin, return ARK_ERR_FAILURE */
  if ((ABS(ark_mem->ark_h) <= ark_mem->ark_hmin*ONEPSM) || 
      (*nefPtr == ark_mem->ark_maxnef)) return(ARK_ERR_FAILURE);

  /* Set etamax = 1 to prevent step size increase at end of this step */
  ark_mem->ark_etamax = ONE;

  /* Set h ratio eta from dsm, rescale, and return for retry of step */
  ark_mem->ark_eta = ONE / (RPowerR(BIAS2*dsm,ONE/ark_mem->ark_L) + ADDON);
  ark_mem->ark_eta = MAX(ETAMIN, MAX(ark_mem->ark_eta, 
				     ark_mem->ark_hmin / ABS(ark_mem->ark_h)));
  if (*nefPtr >= SMALL_NEF) ark_mem->ark_eta = MIN(ark_mem->ark_eta, ETAMXF);
  ARKRescale(ark_mem);
  return(TRY_AGAIN);
}



/*===============================================================
  Private Functions Implementation after succesful stage/step
===============================================================*/

/*---------------------------------------------------------------
 ARKCompleteStep

 This routine performs various update operations when the step 
 solution has passed the local error test. We increment the step 
 counter nst, record the values hold, told and qold, update the 
 yold and fold arrays, fill in the ynew and fnew arrays, update 
 the tau array, and apply the corrections to the zn array.  The 
 tau[i] are the last q values of h, with tau[1] the most recent.
 The counter qwait is decremented, and if 
 qwait == 1 (and q < qmax) we save acor and tq[5] for a 
 possible order increase. 
---------------------------------------------------------------*/
static void ARKCompleteStep(ARKodeMem ark_mem)
{
  int i, j;
  N_Vector tempvec;

  /* Set current time to the end of the step (in case the last 
     stage does not coincide with the step solution).  If tstop 
     is enabled, it is possible for tn + h to be past tstop by 
     roundoff, and in that case, we reset tn (after incrementing 
     by h) to tstop. */
  ark_mem->ark_tn = ark_mem->ark_told + ark_mem->ark_h;
  if (ark_mem->ark_tstopset) {
    if ((ark_mem->ark_tn - ark_mem->ark_tstop)*ark_mem->ark_h > ZERO) 
      ark_mem->ark_tn = ark_mem->ark_tstop;
  }

  /* update scalar quantities */
  ark_mem->ark_nst++;
  ark_mem->ark_hold = ark_mem->ark_h;
  ark_mem->ark_told = ark_mem->ark_tn;
  ark_mem->ark_qold = ark_mem->ark_q;
  
  /* swap yold and ynew arrays */
  tempvec = ark_mem->ark_yold;
  ark_mem->ark_yold = ark_mem->ark_ynew;
  ark_mem->ark_ynew = tempvec;

  /* swap fold and fnew arrays */
  tempvec = ark_mem->ark_fold;
  ark_mem->ark_fold = ark_mem->ark_fnew;
  ark_mem->ark_fnew = tempvec;

  /* update fnew array using explicit and implicit RHS */
  N_VConst(ZERO, ark_mem->ark_fnew);
  if (!ark_mem->ark_explicit) 
    N_VLinearSum(ONE, ark_mem->ark_Fi[ark_mem->ark_stages-1], 
		 ONE, ark_mem->ark_fnew, ark_mem->ark_fnew);
  if (!ark_mem->ark_implicit) 
    N_VLinearSum(ONE, ark_mem->ark_Fe[ark_mem->ark_stages-1], 
		 ONE, ark_mem->ark_fnew, ark_mem->ark_fnew);

  /* update ynew (cannot swap since ark_y is a pointer to user-supplied data) */
  N_VScale(ONE, ark_mem->ark_y, ark_mem->ark_ynew);

  /* update LMM data and statistics */
  for (i=ark_mem->ark_q; i >= 2; i--)  
    ark_mem->ark_tau[i] = ark_mem->ark_tau[i-1];
  if ((ark_mem->ark_q==1) && (ark_mem->ark_nst > 1)) 
    ark_mem->ark_tau[2] = ark_mem->ark_tau[1];
  ark_mem->ark_tau[1] = ark_mem->ark_h;

  for (j=0; j <= ark_mem->ark_q; j++) 
    N_VLinearSum(ark_mem->ark_l[j], ark_mem->ark_acor, ONE, 
		 ark_mem->ark_zn[j], ark_mem->ark_zn[j]);
  ark_mem->ark_qwait--;
  if ((ark_mem->ark_qwait == 1) && (ark_mem->ark_q != ark_mem->ark_qmax)) {
    N_VScale(ONE, ark_mem->ark_acor, ark_mem->ark_zn[ark_mem->ark_qmax]);
    ark_mem->ark_saved_tq5 = ark_mem->ark_tq[5];
    ark_mem->ark_indx_acor = ark_mem->ark_qmax;
  }
}


/*---------------------------------------------------------------
 ARKCompleteStage

 This routine performs various update operations when the solution
 to the nonlinear system for a given stage has completed.  
---------------------------------------------------------------*/
static int ARKCompleteStage(ARKodeMem ark_mem)
{
  int retval;

  /* fill in explicit RHS at current time and solution */
  if (!ark_mem->ark_implicit) {
    retval = ark_mem->ark_fe(ark_mem->ark_tn, ark_mem->ark_y, 
			     ark_mem->ark_tempv, ark_mem->ark_user_data);
    ark_mem->ark_nfe++;
    if (retval < 0)  return(ARK_RHSFUNC_FAIL);
    if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKprepareNextStep

 This routine handles the setting of stepsize and order for the
 next step -- hprime and qprime.  Along with hprime, it sets the
 ratio eta = hprime/h.  It also updates other state variables 
 related to a change of step size or order.
---------------------------------------------------------------*/
 static void ARKPrepareNextStep(ARKodeMem ark_mem, realtype dsm)
{
  /* If etamax = 1, defer step size or order changes */
  if (ark_mem->ark_etamax == ONE) {
    ark_mem->ark_qwait = MAX(ark_mem->ark_qwait, 2);
    ark_mem->ark_qprime = ark_mem->ark_q;
    ark_mem->ark_hprime = ark_mem->ark_h;
    ark_mem->ark_eta = ONE;
    return;
  }

  /* etaq is the ratio of new to old h at the current order */  
  ark_mem->ark_etaq = ONE /(RPowerR(BIAS2*dsm,ONE/ark_mem->ark_L) + ADDON);
  
  /* If no order change, adjust eta and acor in ARKSetEta and return */
  if (ark_mem->ark_qwait != 0) {
    ark_mem->ark_eta = ark_mem->ark_etaq;
    ark_mem->ark_qprime = ark_mem->ark_q;
    ARKSetEta(ark_mem);
    return;
  }
  
  /* If ark_mem->ark_qwait = 0, consider an order change.   etaqm1 and etaqp1 are 
     the ratios of new to old h at orders q-1 and q+1, respectively.
     ARKChooseEta selects the largest; ARKSetEta adjusts eta and acor */
  ark_mem->ark_qwait = 2;
  ark_mem->ark_etaqm1 = ARKComputeEtaqm1(ark_mem);
  ark_mem->ark_etaqp1 = ARKComputeEtaqp1(ark_mem);  
  ARKChooseEta(ark_mem); 
  ARKSetEta(ark_mem);
}


/*---------------------------------------------------------------
 ARKsetEta

 This routine adjusts the value of eta according to the various
 heuristic limits and the optional input hmax.  It also resets
 etamax to be the estimated local error vector. 
---------------------------------------------------------------*/
static void ARKSetEta(ARKodeMem ark_mem)
{

  /* If eta below the threshhold THRESH, reject a change of step size */
  if (ark_mem->ark_eta < THRESH) {
    ark_mem->ark_eta = ONE;
    ark_mem->ark_hprime = ark_mem->ark_h;
  } else {
    /* Limit eta by etamax and hmax, then set hprime */
    ark_mem->ark_eta = MIN(ark_mem->ark_eta, ark_mem->ark_etamax);
    ark_mem->ark_eta /= MAX(ONE, ABS(ark_mem->ark_h) * 
			         ark_mem->ark_hmax_inv*ark_mem->ark_eta);
    ark_mem->ark_hprime = ark_mem->ark_h * ark_mem->ark_eta;
  }
  
  /* Reset etamax for the next step size change, and scale acor */
}


/*---------------------------------------------------------------
 ARKComputeEtaqm1

 This routine computes and returns the value of etaqm1 for a
 possible decrease in order by 1.
---------------------------------------------------------------*/
static realtype ARKComputeEtaqm1(ARKodeMem ark_mem)
{
  realtype ddn;
  
  ark_mem->ark_etaqm1 = ZERO;
  if (ark_mem->ark_q > 1) {
    ddn = N_VWrmsNorm(ark_mem->ark_zn[ark_mem->ark_q], 
		      ark_mem->ark_ewt) * ark_mem->ark_tq[1];
    ark_mem->ark_etaqm1 = ONE / 
      (RPowerR(BIAS1*ddn, ONE/ark_mem->ark_q) + ADDON);
  }
  return(ark_mem->ark_etaqm1);
}


/*---------------------------------------------------------------
 ARKComputeEtaqp1

 This routine computes and returns the value of etaqp1 for a
 possible increase in order by 1.
---------------------------------------------------------------*/
static realtype ARKComputeEtaqp1(ARKodeMem ark_mem)
{
  realtype dup, cquot;
  
  ark_mem->ark_etaqp1 = ZERO;
  if (ark_mem->ark_q != ark_mem->ark_qmax) {
    if (ark_mem->ark_saved_tq5 == ZERO) return(ark_mem->ark_etaqp1);
    cquot = (ark_mem->ark_tq[5] / ark_mem->ark_saved_tq5) * 
      RPowerI(ark_mem->ark_h/ark_mem->ark_tau[2], ark_mem->ark_L);
    N_VLinearSum(-cquot, ark_mem->ark_zn[ark_mem->ark_qmax], ONE, 
		 ark_mem->ark_acor, ark_mem->ark_tempv);
    dup = N_VWrmsNorm(ark_mem->ark_tempv, ark_mem->ark_ewt) * 
      ark_mem->ark_tq[3];
    ark_mem->ark_etaqp1 = ONE / 
      (RPowerR(BIAS3*dup, ONE/(ark_mem->ark_L+1)) + ADDON);
  }
  return(ark_mem->ark_etaqp1);
}


/*---------------------------------------------------------------
 ARKChooseEta
 Given etaqm1, etaq, etaqp1 (the values of eta for qprime =
 q - 1, q, or q + 1, respectively), this routine chooses the 
 maximum eta value, sets eta to that value, and sets qprime to the
 corresponding value of q.  If there is a tie, the preference
 order is to (1) keep the same order, then (2) decrease the order,
 and finally (3) increase the order.  If the maximum eta value
 is below the threshhold THRESH, the order is kept unchanged and
 eta is set to 1.
---------------------------------------------------------------*/
static void ARKChooseEta(ARKodeMem ark_mem)
{
  realtype etam;
  
  etam = MAX(ark_mem->ark_etaqm1, MAX(ark_mem->ark_etaq, 
				      ark_mem->ark_etaqp1));
  
  if (etam < THRESH) {
    ark_mem->ark_eta = ONE;
    ark_mem->ark_qprime = ark_mem->ark_q;
    return;
  }

  if (etam == ark_mem->ark_etaq) {
    ark_mem->ark_eta = ark_mem->ark_etaq;
    ark_mem->ark_qprime = ark_mem->ark_q;

  } else if (etam == ark_mem->ark_etaqm1) {
    ark_mem->ark_eta = ark_mem->ark_etaqm1;
    ark_mem->ark_qprime = ark_mem->ark_q - 1;

  } else {
    ark_mem->ark_eta = ark_mem->ark_etaqp1;
    ark_mem->ark_qprime = ark_mem->ark_q + 1;
    
    /* Store Delta_n in zn[qmax] to be used in order increase.
       This happens at the last step of order q before an increase
       to order q+1, so it represents Delta_n in the ELTE at q+1 */
    N_VScale(ONE, ark_mem->ark_acor, ark_mem->ark_zn[ark_mem->ark_qmax]);

  }
}


/*---------------------------------------------------------------
 ARKHandleFailure

 This routine prints error messages for all cases of failure by
 ARKHin and ARKStep. It returns to ARKode the value that ARKode is 
 to return to the user.
---------------------------------------------------------------*/
static int ARKHandleFailure(ARKodeMem ark_mem, int flag)
{

  /* Depending on flag, print error message and return error flag */
  switch (flag) {
  case ARK_ERR_FAILURE: 
    ARKProcessError(ark_mem, ARK_ERR_FAILURE, "ARKODE", "ARKode", 
		    MSGARK_ERR_FAILS, ark_mem->ark_tn, ark_mem->ark_h);
    break;
  case ARK_CONV_FAILURE:
    ARKProcessError(ark_mem, ARK_CONV_FAILURE, "ARKODE", "ARKode", 
		    MSGARK_CONV_FAILS, ark_mem->ark_tn, ark_mem->ark_h);
    break;
  case ARK_LSETUP_FAIL:
    ARKProcessError(ark_mem, ARK_LSETUP_FAIL, "ARKODE", "ARKode", 
		    MSGARK_SETUP_FAILED, ark_mem->ark_tn);
    break;
  case ARK_LSOLVE_FAIL:
    ARKProcessError(ark_mem, ARK_LSOLVE_FAIL, "ARKODE", "ARKode", 
		    MSGARK_SOLVE_FAILED, ark_mem->ark_tn);
    break;
  case ARK_RHSFUNC_FAIL:
    ARKProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", "ARKode", 
		    MSGARK_RHSFUNC_FAILED, ark_mem->ark_tn);
    break;
  case ARK_UNREC_RHSFUNC_ERR:
    ARKProcessError(ark_mem, ARK_UNREC_RHSFUNC_ERR, "ARKODE", "ARKode", 
		    MSGARK_RHSFUNC_UNREC, ark_mem->ark_tn);
    break;
  case ARK_REPTD_RHSFUNC_ERR:
    ARKProcessError(ark_mem, ARK_REPTD_RHSFUNC_ERR, "ARKODE", "ARKode", 
		    MSGARK_RHSFUNC_REPTD, ark_mem->ark_tn);
    break;
  case ARK_RTFUNC_FAIL:    
    ARKProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "ARKode", 
		    MSGARK_RTFUNC_FAILED, ark_mem->ark_tn);
    break;
  case ARK_TOO_CLOSE:
    ARKProcessError(ark_mem, ARK_TOO_CLOSE, "ARKODE", "ARKode", 
		    MSGARK_TOO_CLOSE);
    break;
  default:
    return(ARK_SUCCESS);   
  }

  return(flag);
}


/*===============================================================
 Root finding   
===============================================================*/

/*---------------------------------------------------------------
 ARKRootCheck1

 This routine completes the initialization of rootfinding memory
 information, and checks whether g has a zero both at and very near
 the initial point of the IVP.

 This routine returns an int equal to:
  ARK_RTFUNC_FAIL = -12  if the g function failed, or
  ARK_SUCCESS     =   0  otherwise.
---------------------------------------------------------------*/
static int ARKRootCheck1(ARKodeMem ark_mem)
{
  int i, retval;
  realtype smallh, hratio, tplus;
  booleantype zroot;

  for (i = 0; i < ark_mem->ark_nrtfn; i++) 
    ark_mem->ark_iroots[i] = 0;
  ark_mem->ark_tlo = ark_mem->ark_tn;
  ark_mem->ark_ttol = (ABS(ark_mem->ark_tn) + 
		       ABS(ark_mem->ark_h))*ark_mem->ark_uround*HUN;

  /* Evaluate g at initial t and check for zero values. */
  retval = ark_mem->ark_gfun(ark_mem->ark_tlo, ark_mem->ark_zn[0], 
			     ark_mem->ark_glo, ark_mem->ark_user_data);
  ark_mem->ark_nge = 1;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  zroot = FALSE;
  for (i = 0; i < ark_mem->ark_nrtfn; i++) {
    if (ABS(ark_mem->ark_glo[i]) == ZERO) {
      zroot = TRUE;
      ark_mem->ark_gactive[i] = FALSE;
    }
  }
  if (!zroot) return(ARK_SUCCESS);

  /* Some g_i is zero at t0; look at g at t0+(small increment). */
  hratio = MAX(ark_mem->ark_ttol/ABS(ark_mem->ark_h), TENTH);
  smallh = hratio*ark_mem->ark_h;
  tplus = ark_mem->ark_tlo + smallh;
  N_VLinearSum(ONE, ark_mem->ark_zn[0], hratio, 
	       ark_mem->ark_zn[1], ark_mem->ark_y);
  retval = ark_mem->ark_gfun(tplus, ark_mem->ark_y, ark_mem->ark_ghi, 
			     ark_mem->ark_user_data);
  ark_mem->ark_nge++;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  /* We check now only the components of g which were exactly 0.0 at t0
   * to see if we can 'activate' them. */
  for (i = 0; i < ark_mem->ark_nrtfn; i++) {
    if (!ark_mem->ark_gactive[i] && ABS(ark_mem->ark_ghi[i]) != ZERO) {
      ark_mem->ark_gactive[i] = TRUE;
      ark_mem->ark_glo[i] = ark_mem->ark_ghi[i];
    }
  }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKRootCheck2

 This routine checks for exact zeros of g at the last root found,
 if the last return was a root.  It then checks for a close pair of
 zeros (an error condition), and for a new root at a nearby point.
 The array glo = g(tlo) at the left endpoint of the search interval
 is adjusted if necessary to assure that all g_i are nonzero
 there, before returning to do a root search in the interval.

 On entry, tlo = tretlast is the last value of tret returned by
 ARKode.  This may be the previous tn, the previous tout value, or
 the last root location.

 This routine returns an int equal to:
      ARK_RTFUNC_FAIL = -12 if the g function failed, or
      CLOSERT         =  3  if a close pair of zeros was found, or
      RTFOUND         =  1  if a new zero of g was found near tlo, or
      ARK_SUCCESS     =  0  otherwise.
---------------------------------------------------------------*/
static int ARKRootCheck2(ARKodeMem ark_mem)
{
  int i, retval;
  realtype smallh, hratio, tplus;
  booleantype zroot;

  if (ark_mem->ark_irfnd == 0) return(ARK_SUCCESS);

  (void) ARKodeGetDky(ark_mem, ark_mem->ark_tlo, 0, ark_mem->ark_y);
  retval = ark_mem->ark_gfun(ark_mem->ark_tlo, ark_mem->ark_y, 
			     ark_mem->ark_glo, ark_mem->ark_user_data);
  ark_mem->ark_nge++;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  zroot = FALSE;
  for (i = 0; i < ark_mem->ark_nrtfn; i++) 
    ark_mem->ark_iroots[i] = 0;
  for (i = 0; i < ark_mem->ark_nrtfn; i++) {
    if (!ark_mem->ark_gactive[i]) continue;
    if (ABS(ark_mem->ark_glo[i]) == ZERO) {
      zroot = TRUE;
      ark_mem->ark_iroots[i] = 1;
    }
  }
  if (!zroot) return(ARK_SUCCESS);

  /* One or more g_i has a zero at tlo.  Check g at tlo+smallh. */
  ark_mem->ark_ttol = (ABS(ark_mem->ark_tn) + 
		       ABS(ark_mem->ark_h))*ark_mem->ark_uround*HUN;
  smallh = (ark_mem->ark_h > ZERO) ? ark_mem->ark_ttol : -ark_mem->ark_ttol;
  tplus = ark_mem->ark_tlo + smallh;
  if ( (tplus - ark_mem->ark_tn)*ark_mem->ark_h >= ZERO) {
    hratio = smallh/ark_mem->ark_h;
    N_VLinearSum(ONE, ark_mem->ark_y, hratio, 
		 ark_mem->ark_zn[1], ark_mem->ark_y);
  } else {
    (void) ARKodeGetDky(ark_mem, tplus, 0, ark_mem->ark_y);
  }
  retval = ark_mem->ark_gfun(tplus, ark_mem->ark_y, ark_mem->ark_ghi, 
			     ark_mem->ark_user_data);
  ark_mem->ark_nge++;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  /* Check for close roots (error return), for a new zero at tlo+smallh,
  and for a g_i that changed from zero to nonzero. */
  zroot = FALSE;
  for (i = 0; i < ark_mem->ark_nrtfn; i++) {
    if (ABS(ark_mem->ark_ghi[i]) == ZERO) {
      if (!ark_mem->ark_gactive[i]) continue;
      if (ark_mem->ark_iroots[i] == 1) return(CLOSERT);
      zroot = TRUE;
      ark_mem->ark_iroots[i] = 1;
    } else {
      if (ark_mem->ark_iroots[i] == 1) 
	ark_mem->ark_glo[i] = ark_mem->ark_ghi[i];
    }
  }
  if (zroot) return(RTFOUND);
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKRootCheck3

 This routine interfaces to ARKRootfind to look for a root of g
 between tlo and either tn or tout, whichever comes first.
 Only roots beyond tlo in the direction of integration are sought.

 This routine returns an int equal to:
      ARK_RTFUNC_FAIL = -12 if the g function failed, or
      RTFOUND        =  1  if a root of g was found, or
      ARK_SUCCESS     =  0  otherwise.
---------------------------------------------------------------*/
static int ARKRootCheck3(ARKodeMem ark_mem)
{
  int i, retval, ier;

  /* Set thi = tn or tout, whichever comes first; set y = y(thi). */
  if (ark_mem->ark_taskc == ARK_ONE_STEP) {
    ark_mem->ark_thi = ark_mem->ark_tn;
    N_VScale(ONE, ark_mem->ark_zn[0], ark_mem->ark_y);
  }
  if (ark_mem->ark_taskc == ARK_NORMAL) {
    if ( (ark_mem->ark_toutc - ark_mem->ark_tn)*ark_mem->ark_h >= ZERO) {
      ark_mem->ark_thi = ark_mem->ark_tn; 
      N_VScale(ONE, ark_mem->ark_zn[0], ark_mem->ark_y);
    } else {
      ark_mem->ark_thi = ark_mem->ark_toutc;
      (void) ARKodeGetDky(ark_mem, ark_mem->ark_thi, 0, ark_mem->ark_y);
    }
  }

  /* Set ark_mem->ark_ghi = g(thi) and call ARKRootfind to search (tlo,thi) for roots. */
  retval = ark_mem->ark_gfun(ark_mem->ark_thi, ark_mem->ark_y, 
			     ark_mem->ark_ghi, ark_mem->ark_user_data);
  ark_mem->ark_nge++;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  ark_mem->ark_ttol = (ABS(ark_mem->ark_tn) + 
		       ABS(ark_mem->ark_h))*ark_mem->ark_uround*HUN;
  ier = ARKRootfind(ark_mem);
  if (ier == ARK_RTFUNC_FAIL) return(ARK_RTFUNC_FAIL);
  for(i=0; i<ark_mem->ark_nrtfn; i++) {
    if (!ark_mem->ark_gactive[i] && ark_mem->ark_grout[i] != ZERO) 
      ark_mem->ark_gactive[i] = TRUE;
  }
  ark_mem->ark_tlo = ark_mem->ark_trout;
  for (i = 0; i < ark_mem->ark_nrtfn; i++) 
    ark_mem->ark_glo[i] = ark_mem->ark_grout[i];

  /* If no root found, return ARK_SUCCESS. */  
  if (ier == ARK_SUCCESS) return(ARK_SUCCESS);

  /* If a root was found, interpolate to get y(trout) and return.  */
  (void) ARKodeGetDky(ark_mem, ark_mem->ark_trout, 0, ark_mem->ark_y);
  return(RTFOUND);

}


/*---------------------------------------------------------------
 ARKRootfind

 This routine solves for a root of g(t) between tlo and thi, if
 one exists.  Only roots of odd multiplicity (i.e. with a change
 of sign in one of the g_i), or exact zeros, are found.
 Here the sign of tlo - thi is arbitrary, but if multiple roots
 are found, the one closest to tlo is returned.

 The method used is the Illinois algorithm, a modified secant method.
 Reference: Kathie L. Hiebert and Lawrence F. Shampine, Implicitly
 Defined Output Points for Solutions of ODEs, Sandia National
 Laboratory Report SAND80-0180, February 1980.

 This routine uses the following parameters for communication:

 nrtfn    = number of functions g_i, or number of components of
            the vector-valued function g(t).  Input only.

 gfun     = user-defined function for g(t).  Its form is
            (void) gfun(t, y, gt, user_data)

 rootdir  = in array specifying the direction of zero-crossings.
            If rootdir[i] > 0, search for roots of g_i only if
            g_i is increasing; if rootdir[i] < 0, search for
            roots of g_i only if g_i is decreasing; otherwise
            always search for roots of g_i.

 gactive  = array specifying whether a component of g should
            or should not be monitored. gactive[i] is initially
            set to TRUE for all i=0,...,nrtfn-1, but it may be
            reset to FALSE if at the first step g[i] is 0.0
            both at the I.C. and at a small perturbation of them.
            gactive[i] is then set back on TRUE only after the 
            corresponding g function moves away from 0.0.

 nge      = cumulative counter for gfun calls.

 ttol     = a convergence tolerance for trout.  Input only.
            When a root at trout is found, it is located only to
            within a tolerance of ttol.  Typically, ttol should
            be set to a value on the order of
               100 * UROUND * max (ABS(tlo), ABS(thi))
            where UROUND is the unit roundoff of the machine.

 tlo, thi = endpoints of the interval in which roots are sought.
            On input, and must be distinct, but tlo - thi may
            be of either sign.  The direction of integration is
            assumed to be from tlo to thi.  On return, tlo and thi
            are the endpoints of the final relevant interval.

 glo, ghi = arrays of length nrtfn containing the vectors g(tlo)
            and g(thi) respectively.  Input and output.  On input,
            none of the glo[i] should be zero.

 trout    = root location, if a root was found, or thi if not.
            Output only.  If a root was found other than an exact
            zero of g, trout is the endpoint thi of the final
            interval bracketing the root, with size at most ttol.

 grout    = array of length nrtfn containing g(trout) on return.

 iroots   = int array of length nrtfn with root information.
            Output only.  If a root was found, iroots indicates
            which components g_i have a root at trout.  For
            i = 0, ..., nrtfn-1, iroots[i] = 1 if g_i has a root
            and g_i is increasing, iroots[i] = -1 if g_i has a
            root and g_i is decreasing, and iroots[i] = 0 if g_i
            has no roots or g_i varies in the direction opposite
            to that indicated by rootdir[i].

 This routine returns an int equal to:
      ARK_RTFUNC_FAIL = -12 if the g function failed, or
      RTFOUND        =  1  if a root of g was found, or
      ARK_SUCCESS     =  0  otherwise.
---------------------------------------------------------------*/
static int ARKRootfind(ARKodeMem ark_mem)
{
  realtype alpha, tmid, gfrac, maxfrac, fracint, fracsub;
  int i, retval, imax, side, sideprev;
  booleantype zroot, sgnchg;

  imax = 0;

  /* First check for change in sign in ghi or for a zero in ghi. */
  maxfrac = ZERO;
  zroot = FALSE;
  sgnchg = FALSE;
  for (i = 0;  i < ark_mem->ark_nrtfn; i++) {
    if (!ark_mem->ark_gactive[i]) continue;
    if (ABS(ark_mem->ark_ghi[i]) == ZERO) {
      if (ark_mem->ark_rootdir[i]*ark_mem->ark_glo[i] <= ZERO) {
        zroot = TRUE;
      }
    } else {
      if ( (ark_mem->ark_glo[i]*ark_mem->ark_ghi[i] < ZERO) && 
	   (ark_mem->ark_rootdir[i]*ark_mem->ark_glo[i] <= ZERO) ) {
        gfrac = ABS(ark_mem->ark_ghi[i]/(ark_mem->ark_ghi[i] - ark_mem->ark_glo[i]));
        if (gfrac > maxfrac) {
          sgnchg = TRUE;
          maxfrac = gfrac;
          imax = i;
        }
      }
    }
  }

  /* If no sign change was found, reset trout and grout.  Then return
     ARK_SUCCESS if no zero was found, or set iroots and return RTFOUND.  */ 
  if (!sgnchg) {
    ark_mem->ark_trout = ark_mem->ark_thi;
    for (i = 0; i < ark_mem->ark_nrtfn; i++) 
      ark_mem->ark_grout[i] = ark_mem->ark_ghi[i];
    if (!zroot) return(ARK_SUCCESS);
    for (i = 0; i < ark_mem->ark_nrtfn; i++) {
      ark_mem->ark_iroots[i] = 0;
      if (!ark_mem->ark_gactive[i]) continue;
      if (ABS(ark_mem->ark_ghi[i]) == ZERO) 
	ark_mem->ark_iroots[i] = ark_mem->ark_glo[i] > 0 ? -1:1;
    }
    return(RTFOUND);
  }

  /* Initialize alpha to avoid compiler warning */
  alpha = ONE;

  /* A sign change was found.  Loop to locate nearest root. */
  side = 0;  sideprev = -1;
  for(;;) {                                    /* Looping point */

    /* Set weight alpha.
       On the first two passes, set alpha = 1.  Thereafter, reset alpha
       according to the side (low vs high) of the subinterval in which
       the sign change was found in the previous two passes.
       If the sides were opposite, set alpha = 1.
       If the sides were the same, then double alpha (if high side),
       or halve alpha (if low side).
       The next guess tmid is the secant method value if alpha = 1, but
       is closer to tlo if alpha < 1, and closer to thi if alpha > 1.    */
    if (sideprev == side) {
      alpha = (side == 2) ? alpha*TWO : alpha*HALF;
    } else {
      alpha = ONE;
    }

    /* Set next root approximation tmid and get g(tmid).
       If tmid is too close to tlo or thi, adjust it inward,
       by a fractional distance that is between 0.1 and 0.5.  */
    tmid = ark_mem->ark_thi - (ark_mem->ark_thi - ark_mem->ark_tlo) *
      ark_mem->ark_ghi[imax]/(ark_mem->ark_ghi[imax] - alpha*ark_mem->ark_glo[imax]);
    if (ABS(tmid - ark_mem->ark_tlo) < HALF*ark_mem->ark_ttol) {
      fracint = ABS(ark_mem->ark_thi - ark_mem->ark_tlo)/ark_mem->ark_ttol;
      fracsub = (fracint > FIVE) ? TENTH : HALF/fracint;
      tmid = ark_mem->ark_tlo + fracsub*(ark_mem->ark_thi - ark_mem->ark_tlo);
    }
    if (ABS(ark_mem->ark_thi - tmid) < HALF*ark_mem->ark_ttol) {
      fracint = ABS(ark_mem->ark_thi - ark_mem->ark_tlo)/ark_mem->ark_ttol;
      fracsub = (fracint > FIVE) ? TENTH : HALF/fracint;
      tmid = ark_mem->ark_thi - fracsub*(ark_mem->ark_thi - ark_mem->ark_tlo);
    }

    (void) ARKodeGetDky(ark_mem, tmid, 0, ark_mem->ark_y);
    retval = ark_mem->ark_gfun(tmid, ark_mem->ark_y, ark_mem->ark_grout, 
			       ark_mem->ark_user_data);
    ark_mem->ark_nge++;
    if (retval != 0) return(ARK_RTFUNC_FAIL);

    /* Check to see in which subinterval g changes sign, and reset imax.
       Set side = 1 if sign change is on low side, or 2 if on high side.  */  
    maxfrac = ZERO;
    zroot = FALSE;
    sgnchg = FALSE;
    sideprev = side;
    for (i = 0;  i < ark_mem->ark_nrtfn; i++) {
      if (!ark_mem->ark_gactive[i]) continue;
      if (ABS(ark_mem->ark_grout[i]) == ZERO) {
        if (ark_mem->ark_rootdir[i]*ark_mem->ark_glo[i] <= ZERO) {
          zroot = TRUE;
        }
      } else {
        if ( (ark_mem->ark_glo[i]*ark_mem->ark_grout[i] < ZERO) && 
	     (ark_mem->ark_rootdir[i]*ark_mem->ark_glo[i] <= ZERO) ) {
          gfrac = ABS(ark_mem->ark_grout[i]/(ark_mem->ark_grout[i] - ark_mem->ark_glo[i]));
          if (gfrac > maxfrac) {
            sgnchg = TRUE;
            maxfrac = gfrac;
            imax = i;
          }
        }
      }
    }
    if (sgnchg) {
      /* Sign change found in (tlo,tmid); replace thi with tmid. */
      ark_mem->ark_thi = tmid;
      for (i = 0; i < ark_mem->ark_nrtfn; i++) 
	ark_mem->ark_ghi[i] = ark_mem->ark_grout[i];
      side = 1;
      /* Stop at root thi if converged; otherwise loop. */
      if (ABS(ark_mem->ark_thi - ark_mem->ark_tlo) <= ark_mem->ark_ttol) break;
      continue;  /* Return to looping point. */
    }

    if (zroot) {
      /* No sign change in (tlo,tmid), but g = 0 at tmid; return root tmid. */
      ark_mem->ark_thi = tmid;
      for (i = 0; i < ark_mem->ark_nrtfn; i++) 
	ark_mem->ark_ghi[i] = ark_mem->ark_grout[i];
      break;
    }

    /* No sign change in (tlo,tmid), and no zero at tmid.
       Sign change must be in (tmid,thi).  Replace tlo with tmid. */
    ark_mem->ark_tlo = tmid;
    for (i = 0; i < ark_mem->ark_nrtfn; i++) 
      ark_mem->ark_glo[i] = ark_mem->ark_grout[i];
    side = 2;
    /* Stop at root thi if converged; otherwise loop back. */
    if (ABS(ark_mem->ark_thi - ark_mem->ark_tlo) <= ark_mem->ark_ttol) 
      break;

  } /* End of root-search loop */

  /* Reset trout and grout, set iroots, and return RTFOUND. */
  ark_mem->ark_trout = ark_mem->ark_thi;
  for (i = 0; i < ark_mem->ark_nrtfn; i++) {
    ark_mem->ark_grout[i] = ark_mem->ark_ghi[i];
    ark_mem->ark_iroots[i] = 0;
    if (!ark_mem->ark_gactive[i]) continue;
    if ( (ABS(ark_mem->ark_ghi[i]) == ZERO) && 
	 (ark_mem->ark_rootdir[i]*ark_mem->ark_glo[i] <= ZERO) ) 
      ark_mem->ark_iroots[i] = ark_mem->ark_glo[i] > 0 ? -1:1;
    if ( (ark_mem->ark_glo[i]*ark_mem->ark_ghi[i] < ZERO) && 
	 (ark_mem->ark_rootdir[i]*ark_mem->ark_glo[i] <= ZERO) ) 
      ark_mem->ark_iroots[i] = ark_mem->ark_glo[i] > 0 ? -1:1;
  }
  return(RTFOUND);
}


/*===============================================================
  Internal EWT function
===============================================================*/

/*---------------------------------------------------------------
 ARKEwtSet

 This routine is responsible for setting the error weight vector ewt,
 according to tol_type, as follows:

 (1) ewt[i] = 1 / (reltol * ABS(ycur[i]) + *abstol), i=0,...,neq-1
     if tol_type = ARK_SS
 (2) ewt[i] = 1 / (reltol * ABS(ycur[i]) + abstol[i]), i=0,...,neq-1
     if tol_type = ARK_SV

 ARKEwtSet returns 0 if ewt is successfully set as above to a
 positive vector and -1 otherwise. In the latter case, ewt is
 considered undefined.

 All the real work is done in the routines ARKEwtSetSS, ARKEwtSetSV.
---------------------------------------------------------------*/
int ARKEwtSet(N_Vector ycur, N_Vector weight, void *data)
{
  ARKodeMem ark_mem;
  int flag = 0;

  /* data points to ark_mem here */
  ark_mem = (ARKodeMem) data;

  switch(ark_mem->ark_itol) {
  case ARK_SS: 
    flag = ARKEwtSetSS(ark_mem, ycur, weight);
    break;
  case ARK_SV: 
    flag = ARKEwtSetSV(ark_mem, ycur, weight);
    break;
  }
  
  return(flag);
}


/*---------------------------------------------------------------
 ARKEwtSetSS

 This routine sets ewt as decribed above in the case tol_type = ARK_SS.
 It tests for non-positive components before inverting. ARKEwtSetSS
 returns 0 if ewt is successfully set to a positive vector
 and -1 otherwise. In the latter case, ewt is considered undefined.
---------------------------------------------------------------*/
static int ARKEwtSetSS(ARKodeMem ark_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, ark_mem->ark_tempv);
  N_VScale(ark_mem->ark_reltol, ark_mem->ark_tempv, ark_mem->ark_tempv);
  N_VAddConst(ark_mem->ark_tempv, ark_mem->ark_Sabstol, ark_mem->ark_tempv);
  if (N_VMin(ark_mem->ark_tempv) <= ZERO) return(-1);
  N_VInv(ark_mem->ark_tempv, weight);
  return(0);
}


/*---------------------------------------------------------------
 ARKEwtSetSV

 This routine sets ewt as decribed above in the case tol_type = ARK_SV.
 It tests for non-positive components before inverting. ARKEwtSetSV
 returns 0 if ewt is successfully set to a positive vector
 and -1 otherwise. In the latter case, ewt is considered undefined.
---------------------------------------------------------------*/
static int ARKEwtSetSV(ARKodeMem ark_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, ark_mem->ark_tempv);
  N_VLinearSum(ark_mem->ark_reltol, ark_mem->ark_tempv, ONE, 
	       ark_mem->ark_Vabstol, ark_mem->ark_tempv);
  if (N_VMin(ark_mem->ark_tempv) <= ZERO) return(-1);
  N_VInv(ark_mem->ark_tempv, weight);
  return(0);
}


/*===============================================================
  ARKODE Error Handling function   
===============================================================*/

/*---------------------------------------------------------------
 ARKProcessError is a high level error handling function
 - if ark_mem==NULL it prints the error message to stderr
 - otherwise, it sets-up and calls the error handling function 
   pointed to by ark_ehfun 
---------------------------------------------------------------*/
void ARKProcessError(ARKodeMem ark_mem, 
                    int error_code, const char *module, const char *fname, 
                    const char *msgfmt, ...)
{
  va_list ap;
  char msg[256];

  /* Initialize the argument pointer variable 
     (msgfmt is the last required argument to ARKProcessError) */
  va_start(ap, msgfmt);

  if (ark_mem == NULL) {    /* We write to stderr */

#ifndef NO_FPRINTF_OUTPUT
    fprintf(stderr, "\n[%s ERROR]  %s\n  ", module, fname);
    fprintf(stderr, msgfmt);
    fprintf(stderr, "\n\n");
#endif

  } else {                 /* We can call ehfun */

    /* Compose the message */
    vsprintf(msg, msgfmt, ap);

    /* Call ehfun */
    ark_mem->ark_ehfun(error_code, module, fname, msg, 
		       ark_mem->ark_eh_data);

  }

  /* Finalize argument processing */
  va_end(ap);

  return;
}

/*---------------------------------------------------------------
 ARKErrHandler is the default error handling function.
   It sends the error message to the stream pointed to by ark_errfp 
---------------------------------------------------------------*/
void ARKErrHandler(int error_code, const char *module,
		   const char *function, char *msg, void *data)
{
  ARKodeMem ark_mem;
  char err_type[10];

  /* data points to ark_mem here */
  ark_mem = (ARKodeMem) data;

  if (error_code == ARK_WARNING)
    sprintf(err_type,"WARNING");
  else
    sprintf(err_type,"ERROR");

#ifndef NO_FPRINTF_OUTPUT
  if (ark_mem->ark_errfp!=NULL) {
    fprintf(ark_mem->ark_errfp,"\n[%s %s]  %s\n",module,err_type,function);
    fprintf(ark_mem->ark_errfp,"  %s\n\n",msg);
  }
#endif

  return;
}

/*---------------------------------------------------------------
 ARKAdapt is the time step adaptivity wrapper function.  This 
 should be called after ARKCompleteStep, as it depends on the 
 values updated by that routine.
---------------------------------------------------------------*/
static int ARKAdapt(ARKodeMem ark_mem, realtype *hnew)
{
  int ier;
  realtype h_acc, h_cfl, cflfac, safety, growth, lbound, ubound;
  cflfac = ark_mem->ark_hadapt_params[0];
  safety = ark_mem->ark_hadapt_params[1];
  growth = ark_mem->ark_hadapt_params[2];
  lbound = ark_mem->ark_hadapt_params[3];
  ubound = ark_mem->ark_hadapt_params[4];

  /* Call algorithm-specific error adaptivity method */
  switch (ark_mem->ark_hadapt_imethod) {
  case(0):    /* PID controller */
    ier = ARKAdaptPID(ark_mem, &h_acc);
    break;
  case(1):    /* PI controller */
    ier = ARKAdaptPI(ark_mem, &h_acc);
    break;
  case(2):    /* I controller */
    ier = ARKAdaptI(ark_mem, &h_acc);
    break;
  case(3):    /* explicit Gustafsson controller */
    ier = ARKAdaptExpGus(ark_mem, &h_acc);
    break;
  case(4):    /* implicit Gustafsson controller */
    ier = ARKAdaptImpGus(ark_mem, &h_acc);
    break;
  case(5):    /* imex Gustafsson controller */
    ier = ARKAdaptImExGus(ark_mem, &h_acc);
    break;
  case(-1):   /* user-supplied controller */
    ier = ark_mem->ark_hadapt(ark_mem->ark_zn[0], 
			      ark_mem->ark_tn, ark_mem->ark_h, 
			      ark_mem->ark_hadapt_ehist[0],
			      ark_mem->ark_hadapt_ehist[1],
			      ark_mem->ark_hadapt_ehist[2],
			      ark_mem->ark_q, ark_mem->ark_p, 
			      &h_acc, ark_mem->ark_hadapt_data);
    break;
  default:
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKAdapt", 
		    "Illegal adapt_imethod.");
    return (ARK_ILL_INPUT);
  }
  if (ier != ARK_SUCCESS) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKAdapt", 
		    "Error in accuracy-based adaptivity function.");
    return (ARK_ILL_INPUT);
  }


  /* Call explicit stability function */
  ier = ark_mem->ark_expstab(ark_mem->ark_zn[0], ark_mem->ark_tn, 
			     &h_cfl, ark_mem->ark_estab_data);
  if (ier != ARK_SUCCESS) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKAdapt", 
		    "Error in explicit stability function.");
    return (ARK_ILL_INPUT);
  }

  /* enforce safety factors, growth bound */
  h_acc *= safety;
  h_cfl *= cflfac;
  h_acc = MIN(h_acc, growth*ark_mem->ark_h);

  /* increment the relevant step counter, set desired step */
  if (h_acc < h_cfl)
    ark_mem->ark_nst_acc++;
  else
    ark_mem->ark_nst_exp++;
  h_acc = MIN(h_acc, h_cfl);

  /* enforce adaptivity bounds */
  if ((h_acc > lbound*ark_mem->ark_h) && (h_acc < ubound*ark_mem->ark_h))
    h_acc = ark_mem->ark_h;
  h_acc = MAX(h_acc, ark_mem->ark_hmin);

  /* set the output and return */
  *hnew = h_acc;

  return(ier);
}


/*---------------------------------------------------------------
 ARKAdaptPID implements a PID time step control algorithm.
---------------------------------------------------------------*/
static int ARKAdaptPID(ARKodeMem ark_mem, realtype *hnew)
{
  realtype k1, k2, k3, e1, e2, e3, hcur, h_acc;

  /* set usable time-step adaptivity parameters */
  k1 = -ark_mem->ark_hadapt_params[5] / ark_mem->ark_p;
  k2 =  ark_mem->ark_hadapt_params[6] / ark_mem->ark_p;
  k3 = -ark_mem->ark_hadapt_params[7] / ark_mem->ark_p;
  e1 = MAX(ark_mem->ark_hadapt_ehist[0], TINY);
  e2 = MAX(ark_mem->ark_hadapt_ehist[1], TINY);
  e3 = MAX(ark_mem->ark_hadapt_ehist[2], TINY);
  hcur = ark_mem->ark_h;
  
  /* compute estimated optimal time step size, set into output */
  h_acc = hcur * RPowerR(e1,k1) * RPowerR(e2,k2) * RPowerR(e3,k3);
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKAdaptPI implements a PI time step control algorithm.
---------------------------------------------------------------*/
static int ARKAdaptPI(ARKodeMem ark_mem, realtype *hnew)
{
  realtype k1, k2, e1, e2, hcur, h_acc;

  /* set usable time-step adaptivity parameters */
  k1 = -ark_mem->ark_hadapt_params[5] / ark_mem->ark_p;
  k2 =  ark_mem->ark_hadapt_params[6] / ark_mem->ark_p;
  e1 = MAX(ark_mem->ark_hadapt_ehist[0], TINY);
  e2 = MAX(ark_mem->ark_hadapt_ehist[1], TINY);
  hcur = ark_mem->ark_h;
  
  /* compute estimated optimal time step size, set into output */
  h_acc = hcur * RPowerR(e1,k1) * RPowerR(e2,k2);
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKAdaptI implements an I time step control algorithm.
---------------------------------------------------------------*/
static int ARKAdaptI(ARKodeMem ark_mem, realtype *hnew)
{
  realtype k1, e1, hcur, h_acc;

  /* set usable time-step adaptivity parameters */
  k1 = -ark_mem->ark_hadapt_params[5] / ark_mem->ark_p;
  e1 = MAX(ark_mem->ark_hadapt_ehist[0], TINY);
  hcur = ark_mem->ark_h;
  
  /* compute estimated optimal time step size, set into output */
  h_acc = hcur * RPowerR(e1,k1);
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKAdaptExpGus implements the explicit Gustafsson time step
 control algorithm.
---------------------------------------------------------------*/
static int ARKAdaptExpGus(ARKodeMem ark_mem, realtype *hnew)
{
  realtype k1, k2, e1, e2, hcur, h_acc;

  /* modified method for first step */
  if (ark_mem->ark_nst < 2) {

    k1 = -ONE / ark_mem->ark_p;
    hcur = ark_mem->ark_h;
    e1 = MAX(ark_mem->ark_hadapt_ehist[0], TINY);
    h_acc = hcur * RPowerR(e1,k1);

  /* general estimate */
  } else {

    k1 = -ark_mem->ark_hadapt_params[5] / ark_mem->ark_p;
    k2 = -ark_mem->ark_hadapt_params[6] / ark_mem->ark_p;
    e1 = MAX(ark_mem->ark_hadapt_ehist[0], TINY);
    e2 = e1 / MAX(ark_mem->ark_hadapt_ehist[1], TINY);
    hcur = ark_mem->ark_h;
    h_acc = hcur * RPowerR(e1,k1) * RPowerR(e2,k2);

  }
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKAdaptImpGus implements the implicit Gustafsson time step 
 control algorithm.
---------------------------------------------------------------*/
static int ARKAdaptImpGus(ARKodeMem ark_mem, realtype *hnew)
{
  realtype k1, k2, e1, e2, hcur, hrat, h_acc;

  /* modified method for first step */
  if (ark_mem->ark_nst < 2) {

    k1 = -ONE / ark_mem->ark_p;
    hcur = ark_mem->ark_h;
    e1 = MAX(ark_mem->ark_hadapt_ehist[0], TINY);
    h_acc = hcur * RPowerR(e1,k1);

  /* general estimate */
  } else {

    k1 = -ark_mem->ark_hadapt_params[5] / ark_mem->ark_p;
    k2 = -ark_mem->ark_hadapt_params[6] / ark_mem->ark_p;
    e1 = MAX(ark_mem->ark_hadapt_ehist[0], TINY);
    e2 = e1 / MAX(ark_mem->ark_hadapt_ehist[1], TINY);
    hcur = ark_mem->ark_h;
    hrat = hcur / ark_mem->ark_hold;
    h_acc = hcur * hrat * RPowerR(e1,k1) * RPowerR(e2,k2);

  }
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKAdaptImExGus implements a combination implicit/explicit 
 Gustafsson time step control algorithm.
---------------------------------------------------------------*/
static int ARKAdaptImExGus(ARKodeMem ark_mem, realtype *hnew)
{
  realtype k1, k2, k3, e1, e2, hcur, hrat, h_acc;

  /* modified method for first step */
  if (ark_mem->ark_nst < 2) {

    k1 = -ONE / ark_mem->ark_p;
    hcur = ark_mem->ark_h;
    e1 = MAX(ark_mem->ark_hadapt_ehist[0], TINY);
    h_acc = hcur * RPowerR(e1,k1);

  /* general estimate */
  } else {

    k1 = -ark_mem->ark_hadapt_params[5] / ark_mem->ark_p;
    k2 = -ark_mem->ark_hadapt_params[6] / ark_mem->ark_p;
    k3 = -ark_mem->ark_hadapt_params[7] / ark_mem->ark_p;
    e1 = MAX(ark_mem->ark_hadapt_ehist[0], TINY);
    e2 = e1 / MAX(ark_mem->ark_hadapt_ehist[1], TINY);
    hcur = ark_mem->ark_h;
    hrat = hcur / ark_mem->ark_hold;
    /* implicit estimate */
    h_acc = hcur * hrat * RPowerR(e1,k3) * RPowerR(e2,k3);
    /* explicit estimate */
    h_acc = MIN(h_acc, hcur * RPowerR(e1,k1) * RPowerR(e2,k2));

  }
  *hnew = h_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKExpStab is the default explicit stability estimation function
---------------------------------------------------------------*/
int ARKExpStab(N_Vector y, realtype t, realtype *hstab, void *data)
{
  /* FILL THIS IN!!! */
  *hstab = RCONST(1.0e200);

  return(ARK_SUCCESS);
}

/*---------------------------------------------------------------
 ARKSetButcherTables

 This routine determines the ERK/DIRK/ARK method to use, based 
 on the desired accuracy and information on whether the problem
 is explicit, implicit or imex.

 VERIFIED TO CORRECTLY FILL IN TABLES!
---------------------------------------------------------------*/
static int ARKSetButcherTables(ARKodeMem ark_mem)
{

  /* if tables have already been specified, just return */
  int i,j;
  booleantype A_set = FALSE;
  for (i=0; i<ARK_S_MAX; i++)
    for (j=0; j<ARK_S_MAX; j++) {
      if (ABS(ark_mem->ark_Ae[i][j]) > TINY)  A_set = TRUE;
      if (ABS(ark_mem->ark_Ai[i][j]) > TINY)  A_set = TRUE;
    }
  if (A_set)  return (ARK_SUCCESS);


  /* TEMPORARY */
  int qsave = ark_mem->ark_q;
  ark_mem->ark_q = 4;
  /*************/


  /**** explicit methods ****/
  if (ark_mem->ark_explicit) {

    switch (ark_mem->ark_q) {
    case(2):    /* Heun-Euler: q=2, p=1, s=2 */
      ARKodeLoadButcherTable(0, &ark_mem->ark_stagesE, 
			     &ark_mem->ark_qE, 
			     &ark_mem->ark_pE, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_bE, 
			     ark_mem->ark_cE, 
			     ark_mem->ark_b2E);
      break;

    case(3):    /* ERK-3-2: q=3, p=2, s=3 */
      ARKodeLoadButcherTable(1, &ark_mem->ark_stagesE, 
			     &ark_mem->ark_qE, 
			     &ark_mem->ark_pE, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_bE, 
			     ark_mem->ark_cE, 
			     ark_mem->ark_b2E);
      break;

    case(4):    /* Merson: q=4, p=3, s=5 */
      ARKodeLoadButcherTable(4, &ark_mem->ark_stagesE, 
			     &ark_mem->ark_qE, 
			     &ark_mem->ark_pE, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_bE, 
			     ark_mem->ark_cE, 
			     ark_mem->ark_b2E);
      break;

    case(5):    /* Dormand-Prince: q=5, p=4, s=7 */
      ARKodeLoadButcherTable(10, &ark_mem->ark_stagesE, 
			     &ark_mem->ark_qE, 
			     &ark_mem->ark_pE, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_bE, 
			     ark_mem->ark_cE, 
			     ark_mem->ark_b2E);
      break;

    case(6):    /* Verner: q=6, p=5, s=8 */
      ARKodeLoadButcherTable(12, &ark_mem->ark_stagesE, 
			     &ark_mem->ark_qE, 
			     &ark_mem->ark_pE, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_bE, 
			     ark_mem->ark_cE, 
			     ark_mem->ark_b2E);
      break;

    default:    /* no available method, set default */
      ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "ARKSetButcherTables", 
		      "No explicit method at requested order, using q=6.");
      ARKodeLoadButcherTable(12, &ark_mem->ark_stagesE, 
			     &ark_mem->ark_qE, 
			     &ark_mem->ark_pE, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_bE, 
			     ark_mem->ark_cE, 
			     ark_mem->ark_b2E);
      break;
    }

    /**** implicit methods ****/
  } else if (ark_mem->ark_implicit) {

    switch (ark_mem->ark_q) {

    case(2):
    case(3):    /* TRBDF2 ESDIRK: q=3, p=2, s=3 */
      ARKodeLoadButcherTable(14, &ark_mem->ark_stages, 
			     &ark_mem->ark_q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    case(4):    /* SDIRK-5-4: q=4, p=3, s=5, A/L stable */
      ARKodeLoadButcherTable(20, &ark_mem->ark_stages, 
			     &ark_mem->ark_q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    case(5):    /* Kvaerno(7,4,5): q=5, p=4, s=7, A/L stable */
      ARKodeLoadButcherTable(24, &ark_mem->ark_stages, 
			     &ark_mem->ark_q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    default:    /* no available method, set default */
      ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "ARKSetButcherTables", 
		      "No implicit method at requested order, using q=5.");
      ARKodeLoadButcherTable(24, &ark_mem->ark_stages, 
			     &ark_mem->ark_q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;
    }
 
    /**** imex methods ****/
  } else {

    switch (ark_mem->ark_q) {

    case(2):
    case(3):    /* ARK3(2)4L[2]SA: q=3, p=2, s=4, DIRK is A/L stable */
      ARKodeLoadButcherTable(3, &ark_mem->ark_stagesE, 
			     &ark_mem->ark_qE, 
			     &ark_mem->ark_pE, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_bE, 
			     ark_mem->ark_cE, 
			     ark_mem->ark_b2E);
      ARKodeLoadButcherTable(16, &ark_mem->ark_stages, 
			     &ark_mem->ark_q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    case(4):    /* ARK4(3)6L[2]SA: q=4, p=3, s=6, DIRK is A/L stable */
      ARKodeLoadButcherTable(6, &ark_mem->ark_stagesE, 
			     &ark_mem->ark_qE, 
			     &ark_mem->ark_pE, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_bE, 
			     ark_mem->ark_cE, 
			     ark_mem->ark_b2E);
      ARKodeLoadButcherTable(22, &ark_mem->ark_stages, 
			     &ark_mem->ark_q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    case(5):    /* ARK5(4)8L[2]SA: q=5, p=4, s=8, DIRK is A stable */
      ARKodeLoadButcherTable(11, &ark_mem->ark_stagesE, 
			     &ark_mem->ark_qE, 
			     &ark_mem->ark_pE, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_bE, 
			     ark_mem->ark_cE, 
			     ark_mem->ark_b2E);
      ARKodeLoadButcherTable(26, &ark_mem->ark_stages, 
			     &ark_mem->ark_q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;

    default:    /* no available method, set default */
      ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "ARKSetButcherTables", 
		      "No imex method at requested order, using q=5.");
      ARKodeLoadButcherTable(11, &ark_mem->ark_stagesE, 
			     &ark_mem->ark_qE, 
			     &ark_mem->ark_pE, 
			     ark_mem->ark_Ae, 
			     ark_mem->ark_bE, 
			     ark_mem->ark_cE, 
			     ark_mem->ark_b2E);
      ARKodeLoadButcherTable(26, &ark_mem->ark_stages, 
			     &ark_mem->ark_q, 
			     &ark_mem->ark_p, 
			     ark_mem->ark_Ai, 
			     ark_mem->ark_b, 
			     ark_mem->ark_c, 
			     ark_mem->ark_b2);
      break;
    }

  }

  /* TEMPORARY */
  ark_mem->ark_q = qsave;
  /*************/

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKCheckButcherTables

 This routine runs through the explicit and/or implicit Butcher 
 tables to ensure that they meet all necessary requirements, 
 including:
     strictly lower-triangular (ERK)
     lower-triangular with some nonzeros on diagonal (IRK)
     method order q > 0 (all)
     embedding order q > 0 (all)
     stages > 0 (all)
     matching stages (ARK)
     matching p (ARK)
     matching q (ARK)
     matching c (ARK)
     matching b (ARK)
     matching b2 (ARK)

Returns ARK_SUCCESS if tables pass, ARK_ILL_INPUT otherwise.
---------------------------------------------------------------*/
static int ARKCheckButcherTables(ARKodeMem ark_mem)
{
  int i, j;
  booleantype okay;
  realtype tol = RCONST(1.0e-12);

  /* check that ERK table is strictly lower triangular */
  if (!ark_mem->ark_implicit) {
    okay = TRUE;
    for (i=0; i<ark_mem->ark_stagesE; i++)
      for (j=i; j<ark_mem->ark_stagesE; j++)
	if (ABS(ark_mem->ark_Ae[i][j]) > tol)  
	  okay = FALSE;
    if (!okay) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKCheckButcherTables",
		      "Ae Butcher table is implicit!");
      return(ARK_ILL_INPUT);
    }
  }

  /* check that IRK table is implicit */
  if (!ark_mem->ark_explicit) {
    okay = FALSE;
    for (i=0; i<ark_mem->ark_stages; i++)
      if (ABS(ark_mem->ark_Ai[i][i]) > tol)  
	okay = TRUE;
    if (!okay) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKCheckButcherTables",
		      "Ai Butcher table is explicit!");
      return(ARK_ILL_INPUT);
    }
  }

  /* check that IRK table is lower triangular */
  if (!ark_mem->ark_explicit) {
    okay = TRUE;
    for (i=0; i<ark_mem->ark_stages; i++)
      for (j=i+1; j<ark_mem->ark_stages; j++)
	if (ABS(ark_mem->ark_Ai[i][j]) > tol)  
	  okay = FALSE;
    if (!okay) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKCheckButcherTables",
		      "Ai Butcher table has entries above diagonal!");
      return(ARK_ILL_INPUT);
    }
  }

  /* check that method order q > 0 */
  if (!ark_mem->ark_implicit) 
    if (ark_mem->ark_qE < 1) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKCheckButcherTables",
		      "ERK method order < 1!");
      return(ARK_ILL_INPUT);
    }
  if (!ark_mem->ark_explicit) 
    if (ark_mem->ark_q < 1) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKCheckButcherTables",
		      "IRK method order < 1!");
      return(ARK_ILL_INPUT);
    }

  /* check that embedding order p > 0 */
  if (!ark_mem->ark_implicit) 
    if (ark_mem->ark_pE < 1) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKCheckButcherTables",
		      "ERK embedding order < 1!");
      return(ARK_ILL_INPUT);
    }
  if (!ark_mem->ark_explicit) 
    if (ark_mem->ark_p < 1) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKCheckButcherTables",
		      "IRK embedding order < 1!");
      return(ARK_ILL_INPUT);
    }

  /* check that stages > 0 */
  if (!ark_mem->ark_implicit) 
    if (ark_mem->ark_stagesE < 1) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKCheckButcherTables",
		      "ERK stages < 1!");
      return(ARK_ILL_INPUT);
    }
  if (!ark_mem->ark_explicit) 
    if (ark_mem->ark_stages < 1) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKCheckButcherTables",
		      "IRK stages < 1!");
      return(ARK_ILL_INPUT);
    }

  /* if purely explicit, update ark_stages */
  if (ark_mem->ark_explicit)
    ark_mem->ark_stages = ark_mem->ark_stagesE;

  /* if purely explicit or purely implicit, return */
  if (ark_mem->ark_implicit || ark_mem->ark_explicit)
    return(ARK_SUCCESS);


  /* check that shared ERK/IRK table data match */
  okay = TRUE;

  /* stages */
  if (ark_mem->ark_stages != ark_mem->ark_stagesE) 
    okay = FALSE;

  /* q */
  if (ark_mem->ark_q != ark_mem->ark_qE) 
    okay = FALSE;

  /* p */
  if (ark_mem->ark_p != ark_mem->ark_pE) 
    okay = FALSE;

  /* c, b, b2 */
  for (i=0; i<ark_mem->ark_stages; i++) {
    if (ABS(ark_mem->ark_c[i]  - ark_mem->ark_cE[i])  > tol)  
      okay = FALSE;
    if (ABS(ark_mem->ark_b[i]  - ark_mem->ark_bE[i])  > tol)  
      okay = FALSE;
    if (ABS(ark_mem->ark_b2[i] - ark_mem->ark_b2E[i]) > tol)  
      okay = FALSE;
  }

  if (!okay) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKCheckButcherTables",
		    "incompatible ERK/IRK Butcher tables!");
    return(ARK_ILL_INPUT);
  }
  
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKDenseEval

 This routine evaluates the dense output formula for the method,
 using a cubic Hermite interpolating formula with the data 
 {yn,yp,fn,fp}, where yn,yp are the stored solution values and 
 fn,fp are the sotred RHS function values at the endpoints of the
 last successful time step [tn,tp].  If greater polynomial order 
 than 3 is requested (and the method order allows it, then we 
 can bootstrap up to a 5th-order accurate interpolant.  For lower 
 order interpolants than cubic, we use {yn,yp,fp} for 
 quadratic, {yn,yp} for linear, and {0.5*(yn+yp)} for constant.
 
 Derivatives have reduced accuracy than the dense output function
 itself, losing one order per derivative.  We will provide 
 derivatives up to d = min(5,q).

 The input 'tau' specifies the time at which to return derivative
 information, the formula is 
    t = tn + tau*h,
 where h = tp-tn, i.e. values 0<tau<1 provide interpolation, 
 other values result in extrapolation.
---------------------------------------------------------------*/
static int ARKDenseEval(ARKodeMem ark_mem, realtype tau, 
			int d, int order, N_Vector yout)
{
  /* local variables */
  int q, retval;
  realtype h, tval, a0, a1, a2, a3, a4, a5, tau2, tau3, tau4, tau5;
  N_Vector fa, fb, ftmp;

  h = ark_mem->ark_hold;
  tau2 = tau*tau;
  tau3 = tau*tau2;
  tau4 = tau*tau3;
  tau5 = tau*tau4;

  /* determine polynomial order q */
  q = MIN(order, ark_mem->ark_dense_q);   /* respect Set routine  */
  q = MIN(q, ark_mem->ark_q);             /* respect method order */
  q = MAX(q, 0);                          /* respect lower bound  */

  /* check that d is possible */
  if ((d > MIN(5,q)) || (d < 0)) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKDenseEval", "Requested illegal derivative.");
    return (ARK_ILL_INPUT);
  }

  /* build polynomial based on order */
  switch (q) {

  case(0):    /* constant interpolant, yout = 0.5*(yn+yp) */
    N_VLinearSum(HALF, ark_mem->ark_yold, HALF, ark_mem->ark_ynew, yout);
    break;

  case(1):    /* linear interpolant */
    if (d == 0) {
      a0 = ONE-tau;
      a1 = tau;
    } else {  /* d=1 */
      a0 =  ONE/h;
      a1 = -ONE/h;
    }
    N_VLinearSum(a0, ark_mem->ark_yold, a1, ark_mem->ark_ynew, yout);
    break;

  case(2):    /* quadratic interpolant */
    if (d == 0) {
      a0 = ONE - TWO*tau + tau2;
      a1 = TWO*tau - tau2;
      a2 = tau2 - tau;
    } else if (d == 1) {
      a0 = (-TWO + TWO*tau)/h;
      a1 = (TWO - TWO*tau)/h;
      a2 = (TWO*tau - ONE)/h;
    } else {  /* d == 2 */
      a0 = TWO/h/h;
      a1 = -TWO/h/h;
      a2 = TWO/h/h;
    }
    N_VLinearSum(a0, ark_mem->ark_yold, a1, ark_mem->ark_ynew, yout);
    N_VLinearSum(a2, ark_mem->ark_fnew, ONE, yout, yout);
    break;

  case(3):    /* cubic interpolant */
    if (d == 0) {
      a0 = ONE - THREE*tau2 + TWO*tau3;
      a1 = THREE*tau2 - TWO*tau3;
      a2 = h*(tau3 - TWO*tau2 + tau);
      a3 = h*(tau3 - tau2);
    } else if (d == 1) {
      a0 = (-SIX*tau + SIX*tau2)/h;
      a1 = (SIX*tau - SIX*tau2)/h;
      a2 = THREE*tau2 - FOUR*tau + ONE;
      a3 = THREE*tau2 - TWO*tau;
    } else if (d == 2) {
      a0 = (-SIX + TWELVE*tau)/h/h;
      a1 = (SIX - TWELVE*tau)/h/h;
      a2 = (SIX*tau - FOUR)/h;
      a3 = (SIX*tau - TWO)/h;
    } else {  /* d == 3 */
      a0 = TWELVE/h/h/h;
      a1 = -TWELVE/h/h/h;
      a2 = SIX/h/h;
      a3 = SIX/h/h;
    }
    N_VLinearSum(a0, ark_mem->ark_yold, a1, ark_mem->ark_ynew, yout);
    N_VLinearSum(a2, ark_mem->ark_fold, ONE, yout, yout);
    N_VLinearSum(a3, ark_mem->ark_fnew, ONE, yout, yout);
    break;

  case(4):    /* quartic interpolant */
    /* first, evaluate cubic interpolant at tau=1/3 */
    tval = ONE/THREE;
    retval = ARKDenseEval(ark_mem, tval, 0, 3, yout);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* second, evaluate RHS at tau=1/3, storing temporary values 
       in ftemp and the result in fa */
    fa = ark_mem->ark_fa;
    ftmp = ark_mem->ark_ftemp;
    tval = ark_mem->ark_told + h/THREE;
    retval = ARKFullRHS(ark_mem, tval, yout, ftmp, fa);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* evaluate desired function */
    if (d == 0) {
      a0 = ONE + SIX*tau2 - RCONST(16.0)*tau3 + RCONST(9.0)*tau4;
      a1 = -SIX*tau2 + RCONST(16.0)*tau3 - RCONST(9.0)*tau4;
      a2 = h*(tau - TWO*tau2 + tau3);
      a3 = h*FOURTH*(FIVE*tau2 - RCONST(14.0)*tau3 + RCONST(9.0)*tau4);
      a4 = h*RCONST(27.0)*FOURTH*(tau2 - TWO*tau3 + tau4);
    } else if (d == 1) {
      a0 = (TWELVE*tau - RCONST(48.0)*tau2 + RCONST(36.0)*tau3)/h;
      a1 = (-TWELVE*tau + RCONST(48.0)*tau2 - RCONST(36.0)*tau3)/h;
      a2 = (ONE - FOUR*tau + THREE*tau2);
      a3 = HALF*(FIVE*tau - RCONST(21.0)*tau2 + RCONST(18.0)*tau3);
      a4 = RCONST(27.0)*HALF*(tau - THREE*tau2 + TWO*tau3);
    } else if (d == 2) {
      a0 = (TWELVE - RCONST(96.0)*tau + RCONST(108.0)*tau2)/h/h;
      a1 = (-TWELVE + RCONST(96.0)*tau - RCONST(108.0)*tau2)/h/h;
      a2 = (-FOUR + SIX*tau)/h;
      a3 = (FIVE*HALF - RCONST(21.0)*tau + RCONST(27.0)*tau2)/h;
      a4 = (RCONST(27.0)*HALF - RCONST(81.0)*tau + RCONST(81.0)*tau2)/h;
    } else if (d == 3) {
      a0 = (-RCONST(96.0) + RCONST(216.0)*tau)/h/h/h;
      a1 = (RCONST(96.0) - RCONST(216.0)*tau)/h/h/h;
      a2 = SIX/h/h;
      a3 = (-RCONST(21.0) + RCONST(54.0)*tau)/h/h;
      a4 = (-RCONST(81.0) + RCONST(162.0)*tau)/h/h;
    } else {  /* d == 4 */
      a0 = RCONST(216.0)/h/h/h/h;
      a1 = -RCONST(216.0)/h/h/h/h;
      a2 = ZERO;
      a3 = RCONST(54.0)/h/h/h;
      a4 = RCONST(162.0)/h/h/h;
    }
    N_VLinearSum(a0, ark_mem->ark_yold, a1, ark_mem->ark_ynew, yout);
    N_VLinearSum(a2, ark_mem->ark_fold, ONE, yout, yout);
    N_VLinearSum(a3, ark_mem->ark_fnew, ONE, yout, yout);
    N_VLinearSum(a4, fa, ONE, yout, yout);
    break;

  case(5):    /* quintic interpolant */
    /* first, evaluate quartic interpolant at tau=2/3 */
    tval = TWO/THREE;
    retval = ARKDenseEval(ark_mem, tval, 0, 4, yout);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* second, evaluate RHS at tau=2/3, storing temporary values 
       in ftemp and the result in fb */
    fb = ark_mem->ark_fb;
    ftmp = ark_mem->ark_ftemp;
    tval = ark_mem->ark_told + h*TWO/THREE;
    retval = ARKFullRHS(ark_mem, tval, yout, ftmp, fb);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* third, evaluate quartic interpolant at tau=1/3 */
    tval = ONE/THREE;
    retval = ARKDenseEval(ark_mem, tval, 0, 4, yout);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* fourth, evaluate RHS at tau=1/3, storing temporary values 
       in ftemp and the result in fb */
    fa = ark_mem->ark_fa;
    ftmp = ark_mem->ark_ftemp;
    tval = ark_mem->ark_told + h/THREE;
    retval = ARKFullRHS(ark_mem, tval, yout, ftmp, fa);
    if (retval != 0)  return(ARK_RHSFUNC_FAIL);

    /* evaluate desired function */
    if (d == 0) {
      a0 = ONE - RCONST(30.0)*tau2 + RCONST(110.0)*tau3 - 
	RCONST(135.0)*tau4 + RCONST(54.0)*tau5;
      a1 = RCONST(30.0)*tau2 - RCONST(110.0)*tau3 + 
	RCONST(135.0)*tau4 - RCONST(54.0)*tau5;
      a2 = h*FOURTH*(FOUR*tau - RCONST(26.0)*tau2 + RCONST(67.0)*tau3 - 
		     RCONST(72.0)*tau4 + RCONST(27.0)*tau5);
      a3 = h*FOURTH*(-RCONST(13.0)*tau2 + RCONST(49.0)*tau3 - 
		     RCONST(63.0)*tau4 + RCONST(27.0)*tau5);
      a4 = h*FOURTH*RCONST(27.0)*(-tau2 + FIVE*tau3 - 
				  SEVEN*tau4 + THREE*tau5);
      a5 = h*FOURTH*RCONST(27.0)*(-TWO*tau2 + SEVEN*tau3 - 
				  RCONST(8.0)*tau4 + THREE*tau5);
    } else if (d == 1) {
      a0 = (-RCONST(60.0)*tau + RCONST(330.0)*tau2 - 
	    RCONST(540.0)*tau3 + RCONST(270.0)*tau4)/h;
      a1 = (RCONST(60.0)*tau - RCONST(330.0)*tau2 + 
	    RCONST(540.0)*tau3 - RCONST(270.0)*tau4)/h;
      a2 = FOURTH*(FOUR - RCONST(52.0)*tau + RCONST(201.0)*tau2 - 
		   RCONST(288.0)*tau3 + RCONST(135.0)*tau4);
      a3 = FOURTH*(-RCONST(26.0)*tau + RCONST(147.0)*tau2 - 
		   RCONST(252.0)*tau3 + RCONST(135.0)*tau4);
      a4 = FOURTH*RCONST(27.0)*(-TWO*tau + RCONST(15.0)*tau2 - 
		      RCONST(28.0)*tau3 + RCONST(15.0)*tau4);
      a5 = FOURTH*RCONST(27.0)*(-FOUR*tau + RCONST(21.0)*tau2 - 
		      RCONST(32.0)*tau3 + RCONST(15.0)*tau4);
    } else if (d == 2) {
      a0 = (-RCONST(60.0) + RCONST(660.0)*tau - 
	    RCONST(1620.0)*tau2 + RCONST(1080.0)*tau3)/h/h;
      a1 = (RCONST(60.0) - RCONST(660.0)*tau + 
	    RCONST(1620.0)*tau2 - RCONST(1080.0)*tau3)/h/h;
      a2 = FOURTH*(-RCONST(52+402.0)*tau - RCONST(864.0)*tau2 + 
		   RCONST(540.0)*tau3)/h;
      a3 = FOURTH*(-RCONST(26.0) + RCONST(294.0)*tau - 
		   RCONST(756.0)*tau2 + RCONST(540.0)*tau3)/h;
      a4 = FOURTH*RCONST(27.0)*(-TWO + RCONST(30.0)*tau - 
				RCONST(84.0)*tau2 + RCONST(60.0)*tau3)/h;
      a5 = FOURTH*RCONST(27.0)*(-FOUR + RCONST(42.0)*tau - 
				RCONST(96.0)*tau2 + RCONST(60.0)*tau3)/h;
    } else if (d == 3) {
      a0 = (RCONST(660.0) - RCONST(3240.0)*tau + RCONST(3240.0)*tau2)/h/h/h;
      a1 = (-RCONST(660.0) + RCONST(3240.0)*tau - RCONST(3240.0)*tau2)/h/h/h;
      a2 = FOURTH*(RCONST(402.0) - RCONST(1728.0)*tau + 
		   RCONST(1620.0)*tau2)/h/h;
      a3 = FOURTH*(RCONST(294.0) - RCONST(1512.0)*tau + 
		   RCONST(1620.0)*tau2)/h/h;
      a4 = FOURTH*RCONST(27.0)*(RCONST(30.0) - RCONST(168.0)*tau + 
				RCONST(180.0)*tau2)/h/h;
      a5 = FOURTH*RCONST(27.0)*(RCONST(42.0) - RCONST(192.0)*tau + 
				RCONST(180.0)*tau2)/h/h;
    } else if (d == 4) {
      a0 = (-RCONST(3240.0) + RCONST(6480.0)*tau)/h/h/h/h;
      a1 = (RCONST(3240.0) - RCONST(6480.0)*tau)/h/h/h/h;
      a2 = (-RCONST(432.0) + RCONST(810.0)*tau)/h/h/h;
      a3 = (-RCONST(378.0) + RCONST(810.0)*tau)/h/h/h;
      a4 = RCONST(27.0)*(-RCONST(42.0) + RCONST(90.0)*tau)/h/h/h;
      a5 = RCONST(27.0)*(-RCONST(48.0) + RCONST(90.0)*tau)/h/h/h;
    } else {  /* d == 5 */
      a0 = RCONST(6480.0)/h/h/h/h/h;
      a1 = -RCONST(6480.0)/h/h/h/h/h;
      a2 = RCONST(810.0)/h/h/h/h;
      a3 = RCONST(810.0)/h/h/h/h;
      a4 = RCONST(2430.0)/h/h/h/h;
      a5 = RCONST(2430.0)/h/h/h/h;
    }
    N_VLinearSum(a0, ark_mem->ark_yold, a1, ark_mem->ark_ynew, yout);
    N_VLinearSum(a2, ark_mem->ark_fold, ONE, yout, yout);
    N_VLinearSum(a3, ark_mem->ark_fnew, ONE, yout, yout);
    N_VLinearSum(a4, fa, ONE, yout, yout);
    N_VLinearSum(a5, fb, ONE, yout, yout);
    break;

  default:
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKDenseEval", 
		    "Illegal polynomial order.");
    return (ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*===============================================================
   EOF
===============================================================*/
