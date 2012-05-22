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
             ARKODE Private Constants                             
===============================================================*/
#define ZERO   RCONST(0.0)      /* real 0.0     */
#define TINY   RCONST(1.0e-10)  /* small number */
#define TENTH  RCONST(0.1)      /* real 0.1     */
#define POINT2 RCONST(0.2)      /* real 0.2     */
#define FOURTH RCONST(0.25)     /* real 0.25    */
#define HALF   RCONST(0.5)      /* real 0.5     */
#define ONE    RCONST(1.0)      /* real 1.0     */
#define TWO    RCONST(2.0)      /* real 2.0     */
#define THREE  RCONST(3.0)      /* real 3.0     */
#define FOUR   RCONST(4.0)      /* real 4.0     */
#define FIVE   RCONST(5.0)      /* real 5.0     */
#define TWELVE RCONST(12.0)     /* real 12.0    */
#define HUN    RCONST(100.0)    /* real 100.0   */


/*===============================================================
             ARKODE Routine-Specific Constants                   
===============================================================*/

/*---------------------------------------------------------------
 Control constants for lower-level functions used by ARKStep:
-----------------------------------------------------------------
 ARKHin return values:
    ARK_SUCCESS
    ARK_RHSFUNC_FAIL
    ARK_TOO_CLOSE

 ARKStep control constants:
    DO_ERROR_TEST
    PREDICT_AGAIN

 ARKStep return values: 
    ARK_SUCCESS
    ARK_LSETUP_FAIL
    ARK_LSOLVE_FAIL
    ARK_RHSFUNC_FAIL
    ARK_RTFUNC_FAIL
    ARK_CONV_FAILURE
    ARK_ERR_FAILURE
    ARK_FIRST_RHSFUNC_ERR

 ARKNls input nflag values:
    FIRST_CALL
    PREV_CONV_FAIL
    PREV_ERR_FAIL
    
 ARKNls return values: 
    ARK_SUCCESS
    ARK_LSETUP_FAIL 
    ARK_LSOLVE_FAIL 
    ARK_RHSFUNC_FAIL
    CONV_FAIL
    RHSFUNC_RECVR
 
 ARKNewtonIteration return values:
    ARK_SUCCESS
    ARK_LSOLVE_FAIL
    ARK_RHSFUNC_FAIL
    CONV_FAIL
    RHSFUNC_RECVR
    TRY_AGAIN
---------------------------------------------------------------*/
#define DO_ERROR_TEST    +2
#define PREDICT_AGAIN    +3

#define CONV_FAIL        +4 
#define TRY_AGAIN        +5

#define FIRST_CALL       +6
#define PREV_CONV_FAIL   +7
#define PREV_ERR_FAIL    +8

#define RHSFUNC_RECVR    +9


/*---------------------------------------------------------------
 Control constants for lower-level rootfinding functions
-----------------------------------------------------------------
 ARKRootCheck1 return values:
    ARK_SUCCESS
    ARK_RTFUNC_FAIL

 ARKRootCheck2 return values:
    ARK_SUCCESS
    ARK_RTFUNC_FAIL
    CLOSERT
    RTFOUND

 ARKRootCheck3 return values:
    ARK_SUCCESS
    ARK_RTFUNC_FAIL
    RTFOUND

 ARKRootfind return values:
    ARK_SUCCESS
    ARK_RTFUNC_FAIL
    RTFOUND
---------------------------------------------------------------*/

#define RTFOUND          +1
#define CLOSERT          +3

/*---------------------------------------------------------------
 Control constants for tolerances
---------------------------------------------------------------*/
#define ARK_NN  0
#define ARK_SS  1
#define ARK_SV  2
#define ARK_WF  3

/*---------------------------------------------------------------
 Algorithmic constants
-----------------------------------------------------------------
 ARKodeGetDky and ARKStep
    FUZZ_FACTOR

 ARKHin
    HLB_FACTOR
    HUB_FACTOR
    H_BIAS,
    MAX_ITERS

 ARKodeCreate 
   CORTES

 ARKStep
    THRESH
    ETAMX1
    ETAMX2
    ETAMX3
    ETAMXF
    ETAMIN
    ETACF
    ADDON
    BIAS1
    BIAS2
    BIAS3
    ONEPSM

    SMALL_NST   nst > SMALL_NST => use ETAMX3 
    MXNCF       max no. of convergence failures during one step 
                try
    MXNEF       max no. of error test failures during one step 
                try
    MXNEF1      max no. of error test failures before forcing a 
                reduction of order
    SMALL_NEF   if an error failure occurs and 
                SMALL_NEF <= nef <= MXNEF1, then reset 
                eta = MIN(eta, ETAMXF)
    LONG_WAIT   number of steps to wait before considering an 
                order change when q==1 and MXNEF1 error test 
                failures have occurred

 ARKNls
    NLS_MAXCOR  maximum no. of corrector iterations for the 
                nonlinear solver
    CRDOWN      constant used in the estimation of the 
                convergence rate (crate) of the iterates for 
                the nonlinear equation
    DGMAX       if |gamma/gammap-1| > DGMAX then call lsetup
    RDIV        declare divergence if ratio del/delp > RDIV
    MSBP        max no. of steps between lsetup calls
---------------------------------------------------------------*/
#define FUZZ_FACTOR RCONST(100.0)

#define HLB_FACTOR  RCONST(100.0)
#define HUB_FACTOR  RCONST(0.1)
#define H_BIAS      HALF
#define MAX_ITERS   4

#define CORTES      RCONST(0.1)

#define THRESH      RCONST(1.5)
#define ETAMX1      RCONST(10000.0) 
#define ETAMX2      RCONST(10.0)
#define ETAMX3      RCONST(10.0)
#define ETAMXF      RCONST(0.2)
#define ETAMIN      RCONST(0.1)
#define ETACF       RCONST(0.25)
#define ADDON       RCONST(0.000001)
#define BIAS1       RCONST(6.0)
#define BIAS2       RCONST(6.0)
#define BIAS3       RCONST(10.0)
#define ONEPSM      RCONST(1.000001)

#define SMALL_NST   10
#define MXNCF       10
#define MXNEF       7
#define MXNEF1      3
#define SMALL_NEF   2
#define LONG_WAIT   10

#define NLS_MAXCOR  3
#define CRDOWN      RCONST(0.3)
#define DGMAX       RCONST(0.3)

#define RDIV        TWO
#define MSBP        20


/*===============================================================
             Private Helper Functions Prototypes
===============================================================*/
static booleantype ARKCheckNvector(N_Vector tmpl);

static int ARKInitialSetup(ARKodeMem ark_mem);

static booleantype ARKAllocVectors(ARKodeMem ark_mem, 
				   N_Vector tmpl);
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

static void ARKAdjustParams(ARKodeMem ark_mem);
static void ARKAdjustOrder(ARKodeMem ark_mem, int deltaq);
static void ARKAdjustBDF(ARKodeMem ark_mem, int deltaq);
static void ARKIncreaseBDF(ARKodeMem ark_mem);
static void ARKDecreaseBDF(ARKodeMem ark_mem);

static void ARKRescale(ARKodeMem ark_mem);

static void ARKPredict(ARKodeMem ark_mem);

static void ARKSet(ARKodeMem ark_mem);
static void ARKSetBDF(ARKodeMem ark_mem);
static void ARKSetTqBDF(ARKodeMem ark_mem, realtype hsum, 
			realtype alpha0, realtype alpha0_hat, 
			realtype xi_inv, realtype xistar_inv);

static int ARKNls(ARKodeMem ark_mem, int nflag);
static int ARKNlsNewton(ARKodeMem ark_mem, int nflag);
static int ARKNewtonIteration(ARKodeMem ark_mem);

static int ARKHandleNFlag(ARKodeMem ark_mem, int *nflagPtr, 
			  realtype saved_t, int *ncfPtr);

static void ARKRestore(ARKodeMem ark_mem, realtype saved_t);

static int ARKDoErrorTest(ARKodeMem ark_mem, int *nflagPtr,
                         realtype saved_t, int *nefPtr, 
			  realtype *dsmPtr);

static void ARKCompleteStep(ARKodeMem ark_mem);

static void ARKPrepareNextStep(ARKodeMem ark_mem, realtype dsm);
static void ARKSetEta(ARKodeMem ark_mem);
static realtype ARKComputeEtaqm1(ARKodeMem ark_mem);
static realtype ARKComputeEtaqp1(ARKodeMem ark_mem);
static void ARKChooseEta(ARKodeMem ark_mem);

static int  ARKHandleFailure(ARKodeMem ark_mem,int flag);

static int ARKRootCheck1(ARKodeMem ark_mem);
static int ARKRootCheck2(ARKodeMem ark_mem);
static int ARKRootCheck3(ARKodeMem ark_mem);
static int ARKRootfind(ARKodeMem ark_mem);


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
  int maxord;
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

  maxord = Q_MAX;

  /* Set uround */
  ark_mem->ark_uround = UNIT_ROUNDOFF;

  /* Set default values for integrator optional inputs */
  ark_mem->ark_fe          = NULL;
  ark_mem->ark_fi          = NULL;
  ark_mem->ark_expstab     = ARKExpStab;
  ark_mem->ark_estab_data  = ark_mem;
  ark_mem->ark_hadapt      = ARKAdapt;
  ark_mem->ark_hadapt_data = ark_mem;
  ark_mem->ark_user_data   = NULL;
  ark_mem->ark_itol        = ARK_NN;
  ark_mem->ark_user_efun   = FALSE;
  ark_mem->ark_efun        = NULL;
  ark_mem->ark_e_data      = NULL;
  ark_mem->ark_ehfun       = ARKErrHandler;
  ark_mem->ark_eh_data     = ark_mem;
  ark_mem->ark_errfp       = stderr;
  ark_mem->ark_qmax        = maxord;
  ark_mem->ark_mxstep      = MXSTEP_DEFAULT;
  ark_mem->ark_mxhnil      = MXHNIL_DEFAULT;
  ark_mem->ark_hin         = ZERO;
  ark_mem->ark_hmin        = HMIN_DEFAULT;
  ark_mem->ark_hmax_inv    = HMAX_INV_DEFAULT;
  ark_mem->ark_tstopset    = FALSE;
  ark_mem->ark_maxcor      = NLS_MAXCOR;
  ark_mem->ark_maxnef      = MXNEF;
  ark_mem->ark_maxncf      = MXNCF;
  ark_mem->ark_nlscoef     = CORTES;

  /* Initialize root finding variables */
  ark_mem->ark_glo        = NULL;
  ark_mem->ark_ghi        = NULL;
  ark_mem->ark_grout      = NULL;
  ark_mem->ark_iroots     = NULL;
  ark_mem->ark_rootdir    = NULL;
  ark_mem->ark_gfun       = NULL;
  ark_mem->ark_nrtfn      = 0;
  ark_mem->ark_gactive    = NULL;
  ark_mem->ark_mxgnull    = 1;

  /* Set the saved value qmax_alloc */
  ark_mem->ark_qmax_alloc = maxord;
  
  /* Initialize lrw and liw */
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
  if (fe == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeInit", MSGARK_NULL_F);
    return(ARK_ILL_INPUT);
  }
  if (fi == NULL) {
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

  /* Allocate the vectors (using y0 as a template) */
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

  ark_mem->ark_qu    = 0;
  ark_mem->ark_hu    = ZERO;
  ark_mem->ark_tolsf = ONE;

  /* Set the linear solver addresses to NULL.
     (We check != NULL later, in ARKode.) */
  ark_mem->ark_linit  = NULL;
  ark_mem->ark_lsetup = NULL;
  ark_mem->ark_lsolve = NULL;
  ark_mem->ark_lfree  = NULL;
  ark_mem->ark_lmem   = NULL;

  /* Initialize zn[0] in the history array */
  N_VScale(ONE, y0, ark_mem->ark_zn[0]);

  /* Initialize all the counters */
  ark_mem->ark_nst     = 0;
  ark_mem->ark_nfe     = 0;
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

  ark_mem->ark_qu    = 0;
  ark_mem->ark_hu    = ZERO;
  ark_mem->ark_tolsf = ONE;

  /* Initialize zn[0] in the history array */
  N_VScale(ONE, y0, ark_mem->ark_zn[0]);
 
  /* Initialize all the counters */
  ark_mem->ark_nst     = 0;
  ark_mem->ark_nfe     = 0;
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
  ark_mem->ark_reltol = reltol;
  ark_mem->ark_Sabstol = abstol;
  ark_mem->ark_itol = ARK_SS;
  ark_mem->ark_user_efun = FALSE;
  ark_mem->ark_efun = ARKEwtSet;
  ark_mem->ark_e_data = NULL; /* set to arkode_mem in InitialSetup */

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

  ark_mem->ark_reltol = reltol;
  N_VScale(ONE, abstol, ark_mem->ark_Vabstol);
  ark_mem->ark_itol = ARK_SV;
  ark_mem->ark_user_efun = FALSE;
  ark_mem->ark_efun = ARKEwtSet;
  ark_mem->ark_e_data = NULL; /* set to arkode_mem in InitialSetup */

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

  ark_mem->ark_itol = ARK_WF;
  ark_mem->ark_user_efun = TRUE;
  ark_mem->ark_efun = efun;
  ark_mem->ark_e_data = NULL; /* set to user_data in InitialSetup */
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
 initialized problem, it computes a tentative initial step size h.

 ARKode supports three modes, specified by itask: ARK_NORMAL, 
 ARK_ONE_STEP and ARK_FIXED_STEP.  In the ARK_NORMAL mode, the 
 solver steps until it reaches or passes tout and then 
 interpolates to obtain y(tout).  In the ARK_ONE_STEP mode, it 
 takes one internal step and returns.  In the ARK_FIXED_STEP mode, 
 it takes fixed step sizes (given by hmin) and turns of all time 
 adaptivity.
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
    
    /* Call f at (t0,y0), set zn[1] = y'(t0), 
       set initial h (from H0 or ARKHin), and scale zn[1] by h.
       Also check for zeros of root function g at and near t0.    */
    retval = ark_mem->ark_fi(ark_mem->ark_tn, ark_mem->ark_zn[0], 
			     ark_mem->ark_zn[1], ark_mem->ark_user_data); 
    ark_mem->ark_nfe++;
    if (retval < 0) {
      ARKProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", 
		      "ARKode", MSGARK_RHSFUNC_FAILED, ark_mem->ark_tn);
      return(ARK_RHSFUNC_FAIL);
    }
    if (retval > 0) {
      ARKProcessError(ark_mem, ARK_FIRST_RHSFUNC_ERR, "ARKODE", 
		      "ARKode", MSGARK_RHSFUNC_FIRST);
      return(ARK_FIRST_RHSFUNC_ERR);
    }

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
    (ABS(ark_mem->ark_tn) + ABS(ark_mem->ark_hu));
  if (ark_mem->ark_hu < ZERO) tfuzz = -tfuzz;
  tp = ark_mem->ark_tn - ark_mem->ark_hu - tfuzz;
  tn1 = ark_mem->ark_tn + tfuzz;
  if ((t-tp)*(t-tn1) > ZERO) {
    ARKProcessError(ark_mem, ARK_BAD_T, "ARKODE", "ARKodeGetDky", 
		    MSGARK_BAD_T, t, ark_mem->ark_tn-ark_mem->ark_hu, 
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
 ARKAllocVectors, and the memory lmem for the linear solver 
 (deallocated by a call to lfree).
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

 This routine allocates the ARKODE vectors ewt, acor, tempv, ftemp, and
 zn[0], ..., zn[maxord].
 If all memory allocations are successful, ARKAllocVectors returns TRUE. 
 Otherwise all allocated memory is freed and ARKAllocVectors returns FALSE.
 This routine also sets the optional outputs lrw and liw, which are
 (respectively) the lengths of the real and integer work spaces
 allocated here.
---------------------------------------------------------------*/
static booleantype ARKAllocVectors(ARKodeMem ark_mem, N_Vector tmpl)
{
  int i, j;

  /* Allocate ewt, acor, tempv, ftemp */
  ark_mem->ark_ewt = N_VClone(tmpl);
  if (ark_mem->ark_ewt == NULL) return(FALSE);
  ark_mem->ark_acor = N_VClone(tmpl);
  if (ark_mem->ark_acor == NULL) {
    N_VDestroy(ark_mem->ark_ewt);
    return(FALSE);
  }
  ark_mem->ark_tempv = N_VClone(tmpl);
  if (ark_mem->ark_tempv == NULL) {
    N_VDestroy(ark_mem->ark_ewt);
    N_VDestroy(ark_mem->ark_acor);
    return(FALSE);
  }
  ark_mem->ark_ftemp = N_VClone(tmpl);
  if (ark_mem->ark_ftemp == NULL) {
    N_VDestroy(ark_mem->ark_tempv);
    N_VDestroy(ark_mem->ark_ewt);
    N_VDestroy(ark_mem->ark_acor);
    return(FALSE);
  }

  /* Allocate zn[0] ... zn[qmax] */
  for (j=0; j <= ark_mem->ark_qmax; j++) {
    ark_mem->ark_zn[j] = N_VClone(tmpl);
    if (ark_mem->ark_zn[j] == NULL) {
      N_VDestroy(ark_mem->ark_ewt);
      N_VDestroy(ark_mem->ark_acor);
      N_VDestroy(ark_mem->ark_tempv);
      N_VDestroy(ark_mem->ark_ftemp);
      for (i=0; i < j; i++) N_VDestroy(ark_mem->ark_zn[i]);
      return(FALSE);
    }
  }

  /* Update solver workspace lengths  */
  ark_mem->ark_lrw += (ark_mem->ark_qmax + 5)*ark_mem->ark_lrw1;
  ark_mem->ark_liw += (ark_mem->ark_qmax + 5)*ark_mem->ark_liw1;

  /* Store the value of qmax used here */
  ark_mem->ark_qmax_alloc = ark_mem->ark_qmax;

  return(TRUE);
}


/*---------------------------------------------------------------
 ARKFreeVectors

 This routine frees the ARKODE vectors allocated in ARKAllocVectors.
---------------------------------------------------------------*/
static void ARKFreeVectors(ARKodeMem ark_mem)
{
  int j, maxord;
  
  maxord = ark_mem->ark_qmax_alloc;

  N_VDestroy(ark_mem->ark_ewt);
  N_VDestroy(ark_mem->ark_acor);
  N_VDestroy(ark_mem->ark_tempv);
  N_VDestroy(ark_mem->ark_ftemp);
  for(j=0; j <= maxord; j++) N_VDestroy(ark_mem->ark_zn[j]);

  ark_mem->ark_lrw -= (maxord + 5)*ark_mem->ark_lrw1;
  ark_mem->ark_liw -= (maxord + 5)*ark_mem->ark_liw1;

  if (ark_mem->ark_VabstolMallocDone) {
    N_VDestroy(ark_mem->ark_Vabstol);
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
  int ier;

  /* Did the user specify tolerances? */
  if (ark_mem->ark_itol == ARK_NN) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKInitialSetup", MSGARK_NO_TOLS);
    return(ARK_ILL_INPUT);
  }

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
 where ydd = estimated second derivative of y.

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

  N_VAbs(ark_mem->ark_zn[0], temp2);
  ark_mem->ark_efun(ark_mem->ark_zn[0], temp1, ark_mem->ark_e_data);
  N_VInv(temp1, temp1);
  N_VLinearSum(HUB_FACTOR, temp2, ONE, temp1, temp1);

  N_VAbs(ark_mem->ark_zn[1], temp2);

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
 ARKYddNorm

 This routine computes an estimate of the second derivative of y
 using a difference quotient, and returns its WRMS norm.
---------------------------------------------------------------*/
static int ARKYddNorm(ARKodeMem ark_mem, realtype hg, realtype *yddnrm)
{
  int retval;

  N_VLinearSum(hg, ark_mem->ark_zn[1], ONE, 
	       ark_mem->ark_zn[0], ark_mem->ark_y);
  retval = ark_mem->ark_fi(ark_mem->ark_tn+hg, ark_mem->ark_y, 
			   ark_mem->ark_tempv, ark_mem->ark_user_data);
  ark_mem->ark_nfe++;
  if (retval < 0) return(ARK_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  N_VLinearSum(ONE, ark_mem->ark_tempv, -ONE, 
	       ark_mem->ark_zn[1], ark_mem->ark_tempv);
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

    ARKPredict(ark_mem);  
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

  ark_mem->ark_etamax = (ark_mem->ark_nst <= SMALL_NST) ? ETAMX2 : ETAMX3;

  /*  Finally, we rescale the acor array to be the 
      estimated local error vector. */
  N_VScale(ark_mem->ark_tq[2], ark_mem->ark_acor, ark_mem->ark_acor);
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

 This routine advances tn by the tentative step size h, and computes
 the predicted array z_n(0), which is overwritten on zn.  The
 prediction of zn is done by repeated additions.
 If tstop is enabled, it is possible for tn + h to be past tstop by roundoff,
 and in that case, we reset tn (after incrementing by h) to tstop.
---------------------------------------------------------------*/
static void ARKPredict(ARKodeMem ark_mem)
{
  int j, k;
  
  ark_mem->ark_tn += ark_mem->ark_h;
  if (ark_mem->ark_tstopset) {
    if ((ark_mem->ark_tn - ark_mem->ark_tstop)*ark_mem->ark_h > ZERO) 
      ark_mem->ark_tn = ark_mem->ark_tstop;
  }
  for (k = 1; k <= ark_mem->ark_q; k++)
    for (j = ark_mem->ark_q; j >= k; j--) 
      N_VLinearSum(ONE, ark_mem->ark_zn[j-1], ONE, 
		   ark_mem->ark_zn[j], ark_mem->ark_zn[j-1]); 
}


/*---------------------------------------------------------------
 ARKSet

 This routine is a high level routine which calls ARKSetBDF to 
 set the polynomial l, the test quantity array tq, 
 and the related variables  rl1, gamma, and gamrat.

 The array tq is loaded with constants used in the control of estimated
 local errors and in the nonlinear convergence test.  Specifically, while
 running at order q, the components of tq are as follows:
   tq[1] = a coefficient used to get the est. local error at order q-1
   tq[2] = a coefficient used to get the est. local error at order q
   tq[3] = a coefficient used to get the est. local error at order q+1
   tq[4] = constant used in nonlinear iteration convergence test
   tq[5] = coefficient used to get the order q+2 derivative vector used in
           the est. local error at order q+1
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

 This routine sets the test quantity array tq.
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
---------------------------------------------------------------*/
static int ARKNls(ARKodeMem ark_mem, int nflag)
{
  return(ARKNlsNewton(ark_mem, nflag));
}


/*---------------------------------------------------------------
 ARKNlsNewton

 This routine handles the Newton iteration. It calls lsetup if 
 indicated, calls ARKNewtonIteration to perform the iteration, and 
 retries a failed attempt at Newton iteration if that is indicated.

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
  N_Vector vtemp1, vtemp2, vtemp3;
  int convfail, retval, ier;
  booleantype callSetup;
  
  vtemp1 = ark_mem->ark_acor;  /* rename acor as vtemp1 for readability  */
  vtemp2 = ark_mem->ark_y;     /* rename y as vtemp2 for readability     */
  vtemp3 = ark_mem->ark_tempv; /* rename tempv as vtemp3 for readability */
  
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
  
  /* Looping point for the solution of the nonlinear system.
     Evaluate f at the predicted y, call lsetup if indicated, and
     call ARKNewtonIteration for the Newton iteration itself.      */
  for(;;) {

    retval = ark_mem->ark_fi(ark_mem->ark_tn, ark_mem->ark_zn[0], 
			     ark_mem->ark_ftemp, ark_mem->ark_user_data);
    ark_mem->ark_nfe++; 
    if (retval < 0) return(ARK_RHSFUNC_FAIL);
    if (retval > 0) return(RHSFUNC_RECVR);

    if (callSetup) {
      ier = ark_mem->ark_lsetup(ark_mem, convfail, ark_mem->ark_zn[0], 
				ark_mem->ark_ftemp, &ark_mem->ark_jcur, 
				vtemp1, vtemp2, vtemp3);
      ark_mem->ark_nsetups++;
      callSetup = FALSE;
      ark_mem->ark_gamrat = ark_mem->ark_crate = ONE; 
      ark_mem->ark_gammap = ark_mem->ark_gamma;
      ark_mem->ark_nstlp = ark_mem->ark_nst;

      /* Return if lsetup failed */
      if (ier < 0) return(ARK_LSETUP_FAIL);
      if (ier > 0) return(CONV_FAIL);
    }

    /* Set acor to zero and load prediction into y vector */
    N_VConst(ZERO, ark_mem->ark_acor);
    N_VScale(ONE, ark_mem->ark_zn[0], ark_mem->ark_y);

    /* Do the Newton iteration */
    ier = ARKNewtonIteration(ark_mem);

    /* If there is a convergence failure and the Jacobian-related 
       data appears not to be current, loop again with a call to lsetup
       in which convfail=ARK_FAIL_BAD_J.  Otherwise return. */
    if (ier != TRY_AGAIN) return(ier);
    
    callSetup = TRUE;
    convfail = ARK_FAIL_BAD_J;
  }
}


/*---------------------------------------------------------------
 ARKNewtonIteration

 This routine performs the Newton iteration. If the iteration succeeds,
 it returns the value ARK_SUCCESS. If not, it may signal the ARKNlsNewton 
 routine to call lsetup again and reattempt the iteration, by
 returning the value TRY_AGAIN. (In this case, ARKNlsNewton must set 
 convfail to ARK_FAIL_BAD_J before calling setup again). 
 Otherwise, this routine returns one of the appropriate values 
 ARK_LSOLVE_FAIL, ARK_RHSFUNC_FAIL, CONV_FAIL, or RHSFUNC_RECVR back 
 to ARKNlsNewton.
---------------------------------------------------------------*/
static int ARKNewtonIteration(ARKodeMem ark_mem)
{
  int m, retval;
  realtype del, delp, dcon;
  N_Vector b;

  ark_mem->ark_mnewt = m = 0;

  /* Initialize delp to avoid compiler warning message */
  del = delp = ZERO;

  /* Looping point for Newton iteration */
  for(;;) {

    /* Evaluate the residual of the nonlinear system*/
    N_VLinearSum(ark_mem->ark_rl1, ark_mem->ark_zn[1], ONE, 
		 ark_mem->ark_acor, ark_mem->ark_tempv);
    N_VLinearSum(ark_mem->ark_gamma, ark_mem->ark_ftemp, -ONE, 
		 ark_mem->ark_tempv, ark_mem->ark_tempv);

    /* Call the lsolve function */
    b = ark_mem->ark_tempv;
    retval = ark_mem->ark_lsolve(ark_mem, b, ark_mem->ark_ewt, 
				 ark_mem->ark_y, ark_mem->ark_ftemp); 
    ark_mem->ark_nni++;
    
    if (retval < 0) return(ARK_LSOLVE_FAIL);
    
    /* If lsolve had a recoverable failure and Jacobian data is
       not current, signal to try the solution again */
    if (retval > 0) { 
      if ((!ark_mem->ark_jcur) && (ark_mem->ark_setupNonNull)) 
	return(TRY_AGAIN);
      else
	return(CONV_FAIL);
    }

    /* Get WRMS norm of correction; add correction to acor and y */
    del = N_VWrmsNorm(b, ark_mem->ark_ewt);
    N_VLinearSum(ONE, ark_mem->ark_acor, ONE, b, ark_mem->ark_acor);
    N_VLinearSum(ONE, ark_mem->ark_zn[0], ONE, 
		 ark_mem->ark_acor, ark_mem->ark_y);
    
    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test */
    if (m > 0) {
      ark_mem->ark_crate = MAX(CRDOWN * ark_mem->ark_crate, del/delp);
    }
    dcon = del * MIN(ONE, ark_mem->ark_crate) / ark_mem->ark_tq[4];
    
    if (dcon <= ONE) {
      ark_mem->ark_acnrm = (m==0) ? 
	del : N_VWrmsNorm(ark_mem->ark_acor, ark_mem->ark_ewt);
      ark_mem->ark_jcur = FALSE;
      return(ARK_SUCCESS); /* Nonlinear system was solved successfully */
    }
    
    ark_mem->ark_mnewt = ++m;
    
    /* Stop at maxcor iterations or if iter. seems to be diverging.
       If still not converged and Jacobian data is not current, 
       signal to try the solution again */
    if ((m == ark_mem->ark_maxcor) || ((m >= 2) && (del > RDIV*delp))) {
      if ((!ark_mem->ark_jcur) && (ark_mem->ark_setupNonNull)) 
	return(TRY_AGAIN);
      else
	return(CONV_FAIL);
    }
    
    /* Save norm of correction, evaluate f, and loop again */
    delp = del;
    retval = ark_mem->ark_fi(ark_mem->ark_tn, ark_mem->ark_y, 
			     ark_mem->ark_ftemp, ark_mem->ark_user_data);
    ark_mem->ark_nfe++;
    if (retval < 0) return(ARK_RHSFUNC_FAIL);
    if (retval > 0) {
      if ((!ark_mem->ark_jcur) && (ark_mem->ark_setupNonNull)) 
	return(TRY_AGAIN);
      else
	return(RHSFUNC_RECVR);
    }

  } /* end loop */
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

  /* The nonlinear soln. failed; increment ncfn and restore zn */
  ark_mem->ark_ncfn++;
  ARKRestore(ark_mem, saved_t);
  
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

  return(PREDICT_AGAIN);
}


/*---------------------------------------------------------------
 ARKRestore

 This routine restores the value of tn to saved_t and undoes the
 prediction.  After execution of ARKRestore, the Nordsieck array zn has
 the same values as before the call to ARKPredict. 
---------------------------------------------------------------*/
static void ARKRestore(ARKodeMem ark_mem, realtype saved_t)
{
  int j, k;
  
  ark_mem->ark_tn = saved_t;
  for (k = 1; k <= ark_mem->ark_q; k++)
    for (j = ark_mem->ark_q; j >= k; j--)
      N_VLinearSum(ONE, ark_mem->ark_zn[j-1], -ONE, 
		   ark_mem->ark_zn[j], ark_mem->ark_zn[j-1]);
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

   - otherwise, set *nflagPtr to PREV_ERR_FAIL, and return TRY_AGAIN. 
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
  
  /* Test failed; increment counters, set nflag, and restore zn array */
  (*nefPtr)++;
  ark_mem->ark_netf++;
  *nflagPtr = PREV_ERR_FAIL;
  ARKRestore(ark_mem, saved_t);

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

  retval = ark_mem->ark_fi(ark_mem->ark_tn, ark_mem->ark_zn[0], 
			   ark_mem->ark_tempv, ark_mem->ark_user_data);
  ark_mem->ark_nfe++;
  if (retval < 0)  return(ARK_RHSFUNC_FAIL);
  if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);

  N_VScale(ark_mem->ark_h, ark_mem->ark_tempv, ark_mem->ark_zn[1]);

  return(TRY_AGAIN);
}


/*===============================================================
  Private Functions Implementation after succesful step
===============================================================*/

/*---------------------------------------------------------------
 ARKCompleteStep

 This routine performs various update operations when the solution
 to the nonlinear system has passed the local error test. 
 We increment the step counter nst, record the values hu and qu,
 update the tau array, and apply the corrections to the zn array.
 The tau[i] are the last q values of h, with tau[1] the most recent.
 The counter qwait is decremented, and if qwait == 1 (and q < qmax)
 we save acor and tq[5] for a possible order increase. 
---------------------------------------------------------------*/
static void ARKCompleteStep(ARKodeMem ark_mem)
{
  int i, j;
  
  ark_mem->ark_nst++;
  ark_mem->ark_hu = ark_mem->ark_h;
  ark_mem->ark_qu = ark_mem->ark_q;

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
 ARKAdapt is the default time step adaptivity function.
---------------------------------------------------------------*/
int ARKAdapt(N_Vector y, realtype t, realtype h, 
	     realtype e1, realtype e2, 
	     realtype e3, int q, int p, 
	     realtype *hnew, void *data)
{
  /* FILL THIS IN!!! */
  *hnew = h;

  return(0);
}


/*---------------------------------------------------------------
 ARKExpStab is the default explicit stability estimation function
---------------------------------------------------------------*/
int ARKExpStab(N_Vector y, realtype t, realtype *hstab, void *data)
{
  /* FILL THIS IN!!! */
  *hstab = RCONST(1.0e200);

  return(0);
}

/*===============================================================
   EOF
===============================================================*/

