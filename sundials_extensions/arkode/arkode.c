/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * This is the implementation file for the main ARKODE integrator.
 * It is independent of the ARKODE linear solver in use.
 * -----------------------------------------------------------------
 ***** UNTOUCHED *****
 * -----------------------------------------------------------------
 */

/*=================================================================*/
/*             Import Header Files                                 */
/*=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

/*=================================================================*/
/*             Macros                                              */
/*=================================================================*/

/* Macro: loop */
#define loop for(;;)

/*=================================================================*/
/*             ARKODE Private Constants                             */
/*=================================================================*/

#define ZERO   RCONST(0.0)     /* real 0.0     */
#define TINY   RCONST(1.0e-10) /* small number */
#define TENTH  RCONST(0.1)     /* real 0.1     */
#define POINT2 RCONST(0.2)     /* real 0.2     */
#define FOURTH RCONST(0.25)    /* real 0.25    */
#define HALF   RCONST(0.5)     /* real 0.5     */
#define ONE    RCONST(1.0)     /* real 1.0     */
#define TWO    RCONST(2.0)     /* real 2.0     */
#define THREE  RCONST(3.0)     /* real 3.0     */
#define FOUR   RCONST(4.0)     /* real 4.0     */
#define FIVE   RCONST(5.0)     /* real 5.0     */
#define TWELVE RCONST(12.0)    /* real 12.0    */
#define HUN    RCONST(100.0)   /* real 100.0   */

/*=================================================================*/
/*             ARKODE Routine-Specific Constants                   */
/*=================================================================*/

/* 
 * Control constants for lower-level functions used by ARKStep 
 * ----------------------------------------------------------
 *
 * ARKHin return values:
 *    ARK_SUCCESS
 *    ARK_RHSFUNC_FAIL
 *    ARK_TOO_CLOSE
 *
 * ARKStep control constants:
 *    DO_ERROR_TEST
 *    PREDICT_AGAIN
 *
 * ARKStep return values: 
 *    ARK_SUCCESS,
 *    ARK_LSETUP_FAIL,  ARK_LSOLVE_FAIL, 
 *    ARK_RHSFUNC_FAIL, ARK_RTFUNC_FAIL
 *    ARK_CONV_FAILURE, ARK_ERR_FAILURE,
 *    ARK_FIRST_RHSFUNC_ERR
 *
 * ARKNls input nflag values:
 *    FIRST_CALL
 *    PREV_CONV_FAIL
 *    PREV_ERR_FAIL
 *    
 * ARKNls return values: 
 *    ARK_SUCCESS,
 *    ARK_LSETUP_FAIL, ARK_LSOLVE_FAIL, ARK_RHSFUNC_FAIL,
 *    CONV_FAIL, RHSFUNC_RECVR
 * 
 * ARKNewtonIteration return values:
 *    ARK_SUCCESS, 
 *    ARK_LSOLVE_FAIL, ARK_RHSFUNC_FAIL
 *    CONV_FAIL, RHSFUNC_RECVR,
 *    TRY_AGAIN
 * 
 */

#define DO_ERROR_TEST    +2
#define PREDICT_AGAIN    +3

#define CONV_FAIL        +4 
#define TRY_AGAIN        +5

#define FIRST_CALL       +6
#define PREV_CONV_FAIL   +7
#define PREV_ERR_FAIL    +8

#define RHSFUNC_RECVR    +9

/*
 * Control constants for lower-level rootfinding functions
 * -------------------------------------------------------
 *
 * ARKRcheck1 return values:
 *    ARK_SUCCESS,
 *    ARK_RTFUNC_FAIL,
 * ARKRcheck2 return values:
 *    ARK_SUCCESS
 *    ARK_RTFUNC_FAIL,
 *    CLOSERT
 *    RTFOUND
 * ARKRcheck3 return values:
 *    ARK_SUCCESS
 *    ARK_RTFUNC_FAIL,
 *    RTFOUND
 * ARKRootfind return values:
 *    ARK_SUCCESS
 *    ARK_RTFUNC_FAIL,
 *    RTFOUND
 */

#define RTFOUND          +1
#define CLOSERT          +3

/*
 * Control constants for tolerances
 * --------------------------------
 */

#define ARK_NN  0
#define ARK_SS  1
#define ARK_SV  2
#define ARK_WF  3

/*
 * Algorithmic constants
 * ---------------------
 *
 * ARKodeGetDky and ARKStep
 *
 *    FUZZ_FACTOR
 *
 * ARKHin
 *
 *    HLB_FACTOR
 *    HUB_FACTOR
 *    H_BIAS
 *    MAX_ITERS
 *
 * ARKodeCreate 
 *
 *   CORTES
 *
 * ARKStep
 *
 *    THRESH
 *    ETAMX1
 *    ETAMX2
 *    ETAMX3
 *    ETAMXF
 *    ETAMIN
 *    ETACF
 *    ADDON
 *    BIAS1
 *    BIAS2
 *    BIAS3
 *    ONEPSM
 *
 *    SMALL_NST   nst > SMALL_NST => use ETAMX3 
 *    MXNCF       max no. of convergence failures during one step try
 *    MXNEF       max no. of error test failures during one step try
 *    MXNEF1      max no. of error test failures before forcing a reduction of order
 *    SMALL_NEF   if an error failure occurs and SMALL_NEF <= nef <= MXNEF1, then
 *                reset eta =  MIN(eta, ETAMXF)
 *    LONG_WAIT   number of steps to wait before considering an order change when
 *                q==1 and MXNEF1 error test failures have occurred
 *
 * ARKNls
 *    
 *    NLS_MAXCOR  maximum no. of corrector iterations for the nonlinear solver
 *    CRDOWN      constant used in the estimation of the convergence rate (crate)
 *                of the iterates for the nonlinear equation
 *    DGMAX       iter == ARK_NEWTON, |gamma/gammap-1| > DGMAX => call lsetup
 *    RDIV        declare divergence if ratio del/delp > RDIV
 *    MSBP        max no. of steps between lsetup calls
 *    
 */


#define FUZZ_FACTOR RCONST(100.0)

#define HLB_FACTOR RCONST(100.0)
#define HUB_FACTOR RCONST(0.1)
#define H_BIAS     HALF
#define MAX_ITERS  4

#define CORTES RCONST(0.1)

#define THRESH RCONST(1.5)
#define ETAMX1 RCONST(10000.0) 
#define ETAMX2 RCONST(10.0)
#define ETAMX3 RCONST(10.0)
#define ETAMXF RCONST(0.2)
#define ETAMIN RCONST(0.1)
#define ETACF  RCONST(0.25)
#define ADDON  RCONST(0.000001)
#define BIAS1  RCONST(6.0)
#define BIAS2  RCONST(6.0)
#define BIAS3  RCONST(10.0)
#define ONEPSM RCONST(1.000001)

#define SMALL_NST    10
#define MXNCF        10
#define MXNEF         7
#define MXNEF1        3
#define SMALL_NEF     2
#define LONG_WAIT    10

#define NLS_MAXCOR 3
#define CRDOWN RCONST(0.3)
#define DGMAX  RCONST(0.3)

#define RDIV      TWO
#define MSBP       20

/*=================================================================*/
/*             Private Helper Functions Prototypes                 */
/*=================================================================*/

static booleantype ARKCheckNvector(N_Vector tmpl);

static int ARKInitialSetup(ARKodeMem ark_mem);

static booleantype ARKAlloarkectors(ARKodeMem ark_mem, N_Vector tmpl);
static void ARKFreeVectors(ARKodeMem ark_mem);

static int ARKEwtSetSS(ARKodeMem ark_mem, N_Vector ycur, N_Vector weight);
static int ARKEwtSetSV(ARKodeMem ark_mem, N_Vector ycur, N_Vector weight);

static int ARKHin(ARKodeMem ark_mem, realtype tout);
static realtype ARKUpperBoundH0(ARKodeMem ark_mem, realtype tdist);
static int ARKYddNorm(ARKodeMem ark_mem, realtype hg, realtype *yddnrm);

static int ARKStep(ARKodeMem ark_mem);

static int ARKsldet(ARKodeMem ark_mem);

static void ARKAdjustParams(ARKodeMem ark_mem);
static void ARKAdjustOrder(ARKodeMem ark_mem, int deltaq);
static void ARKAdjustAdams(ARKodeMem ark_mem, int deltaq);
static void ARKAdjustBDF(ARKodeMem ark_mem, int deltaq);
static void ARKIncreaseBDF(ARKodeMem ark_mem);
static void ARKDecreaseBDF(ARKodeMem ark_mem);

static void ARKRescale(ARKodeMem ark_mem);

static void ARKPredict(ARKodeMem ark_mem);

static void ARKSet(ARKodeMem ark_mem);
static void ARKSetAdams(ARKodeMem ark_mem);
static realtype ARKAdamsStart(ARKodeMem ark_mem, realtype m[]);
static void ARKAdamsFinish(ARKodeMem ark_mem, realtype m[], realtype M[], realtype hsum);
static realtype ARKAltSum(int iend, realtype a[], int k);
static void ARKSetBDF(ARKodeMem ark_mem);
static void ARKSetTqBDF(ARKodeMem ark_mem, realtype hsum, realtype alpha0,
                       realtype alpha0_hat, realtype xi_inv, realtype xistar_inv);

static int ARKNls(ARKodeMem ark_mem, int nflag);
static int ARKNlsFunctional(ARKodeMem ark_mem);
static int ARKNlsNewton(ARKodeMem ark_mem, int nflag);
static int ARKNewtonIteration(ARKodeMem ark_mem);

static int ARKHandleNFlag(ARKodeMem ark_mem, int *nflagPtr, realtype saved_t,
                         int *ncfPtr);

static void ARKRestore(ARKodeMem ark_mem, realtype saved_t);

static int ARKDoErrorTest(ARKodeMem ark_mem, int *nflagPtr,
                         realtype saved_t, int *nefPtr, realtype *dsmPtr);

static void ARKCompleteStep(ARKodeMem ark_mem);

static void ARKPrepareNextStep(ARKodeMem ark_mem, realtype dsm);
static void ARKSetEta(ARKodeMem ark_mem);
static realtype ARKComputeEtaqm1(ARKodeMem ark_mem);
static realtype ARKComputeEtaqp1(ARKodeMem ark_mem);
static void ARKChooseEta(ARKodeMem ark_mem);
static void ARKBDFStab(ARKodeMem ark_mem);

static int  ARKHandleFailure(ARKodeMem ark_mem,int flag);

static int ARKRcheck1(ARKodeMem ark_mem);
static int ARKRcheck2(ARKodeMem ark_mem);
static int ARKRcheck3(ARKodeMem ark_mem);
static int ARKRootfind(ARKodeMem ark_mem);

/* 
 * =================================================================
 * EXPORTED FUNCTIONS IMPLEMENTATION
 * =================================================================
 */

/* 
 * ARKodeCreate
 *
 * ARKodeCreate creates an internal memory block for a problem to 
 * be solved by ARKODE.
 * If successful, ARKodeCreate returns a pointer to the problem memory. 
 * This pointer should be passed to ARKodeInit.  
 * If an initialization error occurs, ARKodeCreate prints an error 
 * message to standard err and returns NULL. 
 */

void *ARKodeCreate(int lmm, int iter)
{
  int maxord;
  ARKodeMem ark_mem;

  /* Test inputs */

  if ((lmm != ARK_ADAMS) && (lmm != ARK_BDF)) {
    ARKProcessError(NULL, 0, "ARKODE", "ARKodeCreate", MSGARK_BAD_LMM);
    return(NULL);
  }
  
  if ((iter != ARK_FUNCTIONAL) && (iter != ARK_NEWTON)) {
    ARKProcessError(NULL, 0, "ARKODE", "ARKodeCreate", MSGARK_BAD_ITER);
    return(NULL);
  }

  ark_mem = NULL;
  ark_mem = (ARKodeMem) malloc(sizeof(struct ARKodeMemRec));
  if (ark_mem == NULL) {
    ARKProcessError(NULL, 0, "ARKODE", "ARKodeCreate", MSGARK_ARKMEM_FAIL);
    return(NULL);
  }

  /* Zero out ark_mem */
  memset(ark_mem, 0, sizeof(struct ARKodeMemRec));

  maxord = (lmm == ARK_ADAMS) ? ADAMS_Q_MAX : BDF_Q_MAX;

  /* copy input parameters into ark_mem */
  ark_mem->ark_lmm  = lmm;
  ark_mem->ark_iter = iter;

  /* Set uround */
  ark_mem->ark_uround = UNIT_ROUNDOFF;

  /* Set default values for integrator optional inputs */
  ark_mem->ark_f          = NULL;
  ark_mem->ark_user_data  = NULL;
  ark_mem->ark_itol       = ARK_NN;
  ark_mem->ark_user_efun  = FALSE;
  ark_mem->ark_efun       = NULL;
  ark_mem->ark_e_data     = NULL;
  ark_mem->ark_ehfun      = ARKErrHandler;
  ark_mem->ark_eh_data    = ark_mem;
  ark_mem->ark_errfp      = stderr;
  ark_mem->ark_qmax       = maxord;
  ark_mem->ark_mxstep     = MXSTEP_DEFAULT;
  ark_mem->ark_mxhnil     = MXHNIL_DEFAULT;
  ark_mem->ark_sldeton    = FALSE;
  ark_mem->ark_hin        = ZERO;
  ark_mem->ark_hmin       = HMIN_DEFAULT;
  ark_mem->ark_hmax_inv   = HMAX_INV_DEFAULT;
  ark_mem->ark_tstopset   = FALSE;
  ark_mem->ark_maxcor     = NLS_MAXCOR;
  ark_mem->ark_maxnef     = MXNEF;
  ark_mem->ark_maxncf     = MXNCF;
  ark_mem->ark_nlscoef    = CORTES;

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

/*-----------------------------------------------------------------*/

#define iter (ark_mem->ark_iter)  
#define lmm  (ark_mem->ark_lmm) 
#define lrw  (ark_mem->ark_lrw)
#define liw  (ark_mem->ark_liw)

/*-----------------------------------------------------------------*/

/*
 * ARKodeInit
 * 
 * ARKodeInit allocates and initializes memory for a problem. All 
 * problem inputs are checked for errors. If any error occurs during 
 * initialization, it is reported to the file whose file pointer is 
 * errfp and an error flag is returned. Otherwise, it returns ARK_SUCCESS
 */

int ARKodeInit(void *arkode_mem, ARKRhsFn fe, ARKRhsFn fi, 
	       ARKExpStabFn EStab, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  booleantype nvectorOK, allocOK;
  long int lrw1, liw1;
  int i,k;

  /* Check arkode_mem */

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKodeInit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check for legal input parameters */

  if (y0==NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKodeInit", MSGARK_NULL_Y0);
    return(ARK_ILL_INPUT);
  }

  if (fe == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKodeInit", MSGARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  if (fi == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKodeInit", MSGARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  if (EStab == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKodeInit", MSGARK_NULL_F);
    return(ARK_ILL_INPUT);
  }

  /* Test if all required vector operations are implemented */

  nvectorOK = ARKCheckNvector(y0);
  if(!nvectorOK) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKodeInit", MSGARK_BAD_NVECTOR);
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

  allocOK = ARKAlloarkectors(ark_mem, y0);
  if (!allocOK) {
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "ARKodeInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  /* All error checking is complete at this point */

  /* Copy the input parameters into ARKODE state */

  ark_mem->ark_f  = fi;
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
     (We check != NULL later, in ARKode, if using ARK_NEWTON.) */

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
  ark_mem->ark_nscon   = 0;
  ark_mem->ark_nge     = 0;

  ark_mem->ark_irfnd   = 0;

  /* Initialize other integrator optional outputs */

  ark_mem->ark_h0u      = ZERO;
  ark_mem->ark_next_h   = ZERO;
  ark_mem->ark_next_q   = 0;

  /* Initialize Stablilty Limit Detection data */
  /* NOTE: We do this even if stab lim det was not
     turned on yet. This way, the user can turn it
     on at any time */

  ark_mem->ark_nor = 0;
  for (i = 1; i <= 5; i++)
    for (k = 1; k <= 3; k++) 
      ark_mem->ark_ssdat[i-1][k-1] = ZERO;

  /* Problem has been successfully initialized */

  ark_mem->ark_MallocDone = TRUE;

  return(ARK_SUCCESS);
}

/*-----------------------------------------------------------------*/

#define lrw1 (ark_mem->ark_lrw1)
#define liw1 (ark_mem->ark_liw1)

/*-----------------------------------------------------------------*/

/*
 * ARKodeReInit
 *
 * ARKodeReInit re-initializes ARKODE's memory for a problem, assuming
 * it has already been allocated in a prior ARKodeInit call.
 * All problem specification inputs are checked for errors.
 * If any error occurs during initialization, it is reported to the
 * file whose file pointer is errfp.
 * The return value is ARK_SUCCESS = 0 if no errors occurred, or
 * a negative value otherwise.
 */

int ARKodeReInit(void *arkode_mem, realtype t0, N_Vector y0)
{
  ARKodeMem ark_mem;
  int i,k;
 
  /* Check arkode_mem */

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKodeReInit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check if arkode_mem was allocated */

  if (ark_mem->ark_MallocDone == FALSE) {
    ARKProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", "ARKodeReInit", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check for legal input parameters */

  if (y0 == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKodeReInit", MSGARK_NULL_Y0);
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
  ark_mem->ark_nscon   = 0;
  ark_mem->ark_nge     = 0;

  ark_mem->ark_irfnd   = 0;

  /* Initialize other integrator optional outputs */

  ark_mem->ark_h0u      = ZERO;
  ark_mem->ark_next_h   = ZERO;
  ark_mem->ark_next_q   = 0;

  /* Initialize Stablilty Limit Detection data */

  ark_mem->ark_nor = 0;
  for (i = 1; i <= 5; i++)
    for (k = 1; k <= 3; k++) 
      ark_mem->ark_ssdat[i-1][k-1] = ZERO;
  
  /* Problem has been successfully re-initialized */

  return(ARK_SUCCESS);
}

/*-----------------------------------------------------------------*/

/*
 * ARKodeSStolerances
 * ARKodeSVtolerances
 * ARKodeWFtolerances
 *
 * These functions specify the integration tolerances. One of them
 * MUST be called before the first call to ARKode.
 *
 * ARKodeSStolerances specifies scalar relative and absolute tolerances.
 * ARKodeSVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance (a potentially different absolute tolerance 
 *   for each vector component).
 * ARKodeWFtolerances specifies a user-provides function (of type ARKEwtFn)
 *   which will be called to set the error weight vector.
 */

int ARKodeSStolerances(void *arkode_mem, realtype reltol, realtype abstol)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKodeSStolerances", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_MallocDone == FALSE) {
    ARKProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", "ARKodeSStolerances", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check inputs */

  if (reltol < ZERO) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKodeSStolerances", MSGARK_BAD_RELTOL);
    return(ARK_ILL_INPUT);
  }

  if (abstol < ZERO) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKodeSStolerances", MSGARK_BAD_ABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Copy tolerances into memory */
  
  ark_mem->ark_reltol = reltol;
  ark_mem->ark_Sabstol = abstol;

  ark_mem->ark_itol = ARK_SS;

  ark_mem->ark_user_efun = FALSE;
  ark_mem->ark_efun = ARKEwtSet;
  ark_mem->ark_e_data = NULL; /* will be set to arkode_mem in InitialSetup */

  return(ARK_SUCCESS);
}


int ARKodeSVtolerances(void *arkode_mem, realtype reltol, N_Vector abstol)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKodeSVtolerances", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_MallocDone == FALSE) {
    ARKProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", "ARKodeSVtolerances", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  /* Check inputs */

  if (reltol < ZERO) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKodeSVtolerances", MSGARK_BAD_RELTOL);
    return(ARK_ILL_INPUT);
  }

  if (N_VMin(abstol) < ZERO) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKodeSVtolerances", MSGARK_BAD_ABSTOL);
    return(ARK_ILL_INPUT);
  }

  /* Copy tolerances into memory */
  
  if ( !(ark_mem->ark_VabstolMallocDone) ) {
    ark_mem->ark_Vabstol = N_VClone(ark_mem->ark_ewt);
    lrw += lrw1;
    liw += liw1;
    ark_mem->ark_VabstolMallocDone = TRUE;
  }

  ark_mem->ark_reltol = reltol;
  N_VScale(ONE, abstol, ark_mem->ark_Vabstol);

  ark_mem->ark_itol = ARK_SV;

  ark_mem->ark_user_efun = FALSE;
  ark_mem->ark_efun = ARKEwtSet;
  ark_mem->ark_e_data = NULL; /* will be set to arkode_mem in InitialSetup */

  return(ARK_SUCCESS);
}


int ARKodeWFtolerances(void *arkode_mem, ARKEwtFn efun)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKodeWFtolerances", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (ark_mem->ark_MallocDone == FALSE) {
    ARKProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", "ARKodeWFtolerances", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }

  ark_mem->ark_itol = ARK_WF;

  ark_mem->ark_user_efun = TRUE;
  ark_mem->ark_efun = efun;
  ark_mem->ark_e_data = NULL; /* will be set to user_data in InitialSetup */

  return(ARK_SUCCESS);
}

/*-----------------------------------------------------------------*/

#define gfun    (ark_mem->ark_gfun)
#define glo     (ark_mem->ark_glo)
#define ghi     (ark_mem->ark_ghi)
#define grout   (ark_mem->ark_grout)
#define iroots  (ark_mem->ark_iroots)
#define rootdir (ark_mem->ark_rootdir)
#define gactive (ark_mem->ark_gactive)

/*-----------------------------------------------------------------*/

/*
 * ARKodeRootInit
 *
 * ARKodeRootInit initializes a rootfinding problem to be solved
 * during the integration of the ODE system.  It loads the root
 * function pointer and the number of root functions, and allocates
 * workspace memory.  The return value is ARK_SUCCESS = 0 if no errors
 * occurred, or a negative value otherwise.
 */

int ARKodeRootInit(void *arkode_mem, int nrtfn, ARKRootFn g)
{
  ARKodeMem ark_mem;
  int i, nrt;

  /* Check arkode_mem pointer */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKodeRootInit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  nrt = (nrtfn < 0) ? 0 : nrtfn;

  /* If rerunning ARKodeRootInit() with a different number of root
     functions (changing number of gfun components), then free
     currently held memory resources */
  if ((nrt != ark_mem->ark_nrtfn) && (ark_mem->ark_nrtfn > 0)) {
    free(glo); glo = NULL;
    free(ghi); ghi = NULL;
    free(grout); grout = NULL;
    free(iroots); iroots = NULL;
    free(rootdir); rootdir = NULL;
    free(gactive); gactive = NULL;

    lrw -= 3 * (ark_mem->ark_nrtfn);
    liw -= 3 * (ark_mem->ark_nrtfn);
  }

  /* If ARKodeRootInit() was called with nrtfn == 0, then set ark_nrtfn to
     zero and ark_gfun to NULL before returning */
  if (nrt == 0) {
    ark_mem->ark_nrtfn = nrt;
    gfun = NULL;
    return(ARK_SUCCESS);
  }

  /* If rerunning ARKodeRootInit() with the same number of root functions
     (not changing number of gfun components), then check if the root
     function argument has changed */
  /* If g != NULL then return as currently reserved memory resources
     will suffice */
  if (nrt == ark_mem->ark_nrtfn) {
    if (g != gfun) {
      if (g == NULL) {
        free(glo); glo = NULL;
        free(ghi); ghi = NULL;
        free(grout); grout = NULL;
        free(iroots); iroots = NULL;
        free(rootdir); rootdir = NULL;
        free(gactive); gactive = NULL;

        lrw -= 3*nrt;
        liw -= 3*nrt;

        ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKodeRootInit", MSGARK_NULL_G);
        return(ARK_ILL_INPUT);
      }
      else {
        gfun = g;
        return(ARK_SUCCESS);
      }
    }
    else return(ARK_SUCCESS);
  }

  /* Set variable values in ARKode memory block */
  ark_mem->ark_nrtfn = nrt;
  if (g == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKodeRootInit", MSGARK_NULL_G);
    return(ARK_ILL_INPUT);
  }
  else gfun = g;

  /* Allocate necessary memory and return */
  glo = NULL;
  glo = (realtype *) malloc(nrt*sizeof(realtype));
  if (glo == NULL) {
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  ghi = NULL;
  ghi = (realtype *) malloc(nrt*sizeof(realtype));
  if (ghi == NULL) {
    free(glo); glo = NULL;
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  grout = NULL;
  grout = (realtype *) malloc(nrt*sizeof(realtype));
  if (grout == NULL) {
    free(glo); glo = NULL;
    free(ghi); ghi = NULL;
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  iroots = NULL;
  iroots = (int *) malloc(nrt*sizeof(int));
  if (iroots == NULL) {
    free(glo); glo = NULL; 
    free(ghi); ghi = NULL;
    free(grout); grout = NULL;
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  rootdir = NULL;
  rootdir = (int *) malloc(nrt*sizeof(int));
  if (rootdir == NULL) {
    free(glo); glo = NULL; 
    free(ghi); ghi = NULL;
    free(grout); grout = NULL;
    free(iroots); iroots = NULL;
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODE", "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  gactive = NULL;
  gactive = (booleantype *) malloc(nrt*sizeof(booleantype));
  if (gactive == NULL) {
    free(glo); glo = NULL; 
    free(ghi); ghi = NULL;
    free(grout); grout = NULL;
    free(iroots); iroots = NULL;
    free(rootdir); rootdir = NULL;
    ARKProcessError(ark_mem, ARK_MEM_FAIL, "ARKODES", "ARKodeRootInit", MSGARK_MEM_FAIL);
    return(ARK_MEM_FAIL);
  }

  /* Set default values for rootdir (both directions) */
  for(i=0; i<nrt; i++) rootdir[i] = 0;

  /* Set default values for gactive (all active) */
  for(i=0; i<nrt; i++) gactive[i] = TRUE;

  lrw += 3*nrt;
  liw += 3*nrt;

  return(ARK_SUCCESS);
}

/* 
 * =================================================================
 * Readibility Constants
 * =================================================================
 */

#define f              (ark_mem->ark_f)      
#define user_data      (ark_mem->ark_user_data) 
#define efun           (ark_mem->ark_efun)
#define e_data         (ark_mem->ark_e_data) 
#define qmax           (ark_mem->ark_qmax)
#define mxstep         (ark_mem->ark_mxstep)
#define mxhnil         (ark_mem->ark_mxhnil)
#define sldeton        (ark_mem->ark_sldeton)
#define hin            (ark_mem->ark_hin)
#define hmin           (ark_mem->ark_hmin)
#define hmax_inv       (ark_mem->ark_hmax_inv)
#define tstop          (ark_mem->ark_tstop)
#define tstopset       (ark_mem->ark_tstopset)
#define maxnef         (ark_mem->ark_maxnef)
#define maxncf         (ark_mem->ark_maxncf)
#define maxcor         (ark_mem->ark_maxcor)
#define nlscoef        (ark_mem->ark_nlscoef)
#define itol           (ark_mem->ark_itol)         
#define reltol         (ark_mem->ark_reltol)       
#define Sabstol        (ark_mem->ark_Sabstol)
#define Vabstol        (ark_mem->ark_Vabstol)

#define uround         (ark_mem->ark_uround)  
#define zn             (ark_mem->ark_zn) 
#define ewt            (ark_mem->ark_ewt)  
#define y              (ark_mem->ark_y)
#define acor           (ark_mem->ark_acor)
#define tempv          (ark_mem->ark_tempv)
#define ftemp          (ark_mem->ark_ftemp) 
#define q              (ark_mem->ark_q)
#define qprime         (ark_mem->ark_qprime)
#define next_q         (ark_mem->ark_next_q)
#define qwait          (ark_mem->ark_qwait)
#define L              (ark_mem->ark_L)
#define h              (ark_mem->ark_h)
#define hprime         (ark_mem->ark_hprime)
#define next_h         (ark_mem->ark_next_h)
#define eta            (ark_mem->ark_eta) 
#define etaqm1         (ark_mem->ark_etaqm1) 
#define etaq           (ark_mem->ark_etaq) 
#define etaqp1         (ark_mem->ark_etaqp1) 
#define nscon          (ark_mem->ark_nscon)
#define hscale         (ark_mem->ark_hscale)
#define tn             (ark_mem->ark_tn)
#define tau            (ark_mem->ark_tau)
#define tq             (ark_mem->ark_tq)
#define l              (ark_mem->ark_l)
#define rl1            (ark_mem->ark_rl1)
#define gamma          (ark_mem->ark_gamma) 
#define gammap         (ark_mem->ark_gammap) 
#define gamrat         (ark_mem->ark_gamrat)
#define crate          (ark_mem->ark_crate)
#define acnrm          (ark_mem->ark_acnrm)
#define mnewt          (ark_mem->ark_mnewt)
#define etamax         (ark_mem->ark_etamax)
#define nst            (ark_mem->ark_nst)
#define nfe            (ark_mem->ark_nfe)
#define ncfn           (ark_mem->ark_ncfn)
#define netf           (ark_mem->ark_netf)
#define nni            (ark_mem->ark_nni)
#define nsetups        (ark_mem->ark_nsetups)
#define nhnil          (ark_mem->ark_nhnil)
#define linit          (ark_mem->ark_linit)
#define lsetup         (ark_mem->ark_lsetup)
#define lsolve         (ark_mem->ark_lsolve) 
#define lfree          (ark_mem->ark_lfree) 
#define lmem           (ark_mem->ark_lmem) 
#define qu             (ark_mem->ark_qu)          
#define nstlp          (ark_mem->ark_nstlp)  
#define h0u            (ark_mem->ark_h0u)
#define hu             (ark_mem->ark_hu)         
#define saved_tq5      (ark_mem->ark_saved_tq5)  
#define indx_acor      (ark_mem->ark_indx_acor)
#define jcur           (ark_mem->ark_jcur)         
#define tolsf          (ark_mem->ark_tolsf)      
#define setupNonNull   (ark_mem->ark_setupNonNull) 
#define nor            (ark_mem->ark_nor)
#define ssdat          (ark_mem->ark_ssdat)

#define nrtfn          (ark_mem->ark_nrtfn)
#define tlo            (ark_mem->ark_tlo)
#define thi            (ark_mem->ark_thi)
#define tretlast       (ark_mem->ark_tretlast)
#define toutc          (ark_mem->ark_toutc)
#define trout          (ark_mem->ark_trout)
#define ttol           (ark_mem->ark_ttol)
#define taskc          (ark_mem->ark_taskc)
#define irfnd          (ark_mem->ark_irfnd)
#define nge            (ark_mem->ark_nge)


/*-----------------------------------------------------------------*/

/*
 * ARKode
 *
 * This routine is the main driver of the ARKODE package. 
 *
 * It integrates over a time interval defined by the user, by calling
 * ARKStep to do internal time steps.
 *
 * The first time that ARKode is called for a successfully initialized
 * problem, it computes a tentative initial step size h.
 *
 * ARKode supports two modes, specified by itask: ARK_NORMAL, ARK_ONE_STEP.
 * In the ARK_NORMAL mode, the solver steps until it reaches or passes tout
 * and then interpolates to obtain y(tout).
 * In the ARK_ONE_STEP mode, it takes one internal step and returns.
 */

int ARKode(void *arkode_mem, realtype tout, N_Vector yout, 
          realtype *tret, int itask)
{
  ARKodeMem ark_mem;
  long int nstloc;
  int retval, hflag, kflag, istate, ir, ier, irfndp;
  int ewtsetOK;
  realtype troundoff, tout_hin, rh, nrm;
  booleantype inactive_roots;

  /*
   * -------------------------------------
   * 1. Check and process inputs
   * -------------------------------------
   */

  /* Check if arkode_mem exists */
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKode", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Check if arkode_mem was allocated */
  if (ark_mem->ark_MallocDone == FALSE) {
    ARKProcessError(ark_mem, ARK_NO_MALLOC, "ARKODE", "ARKode", MSGARK_NO_MALLOC);
    return(ARK_NO_MALLOC);
  }
  
  /* Check for yout != NULL */
  if ((y = yout) == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", MSGARK_YOUT_NULL);
    return(ARK_ILL_INPUT);
  }

  /* Check for tret != NULL */
  if (tret == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", MSGARK_TRET_NULL);
    return(ARK_ILL_INPUT);
  }

  /* Check for valid itask */
  if ( (itask != ARK_NORMAL) && (itask != ARK_ONE_STEP) ) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", MSGARK_BAD_ITASK);
    return(ARK_ILL_INPUT);
  }

  if (itask == ARK_NORMAL) toutc = tout;
  taskc = itask;

  /*
   * ----------------------------------------
   * 2. Initializations performed only at
   *    the first step (nst=0):
   *    - initial setup
   *    - initialize Nordsieck history array
   *    - compute initial step size
   *    - check for approach to tstop
   *    - check for approach to a root
   * ----------------------------------------
   */

  if (nst == 0) {

    ier = ARKInitialSetup(ark_mem);
    if (ier!= ARK_SUCCESS) return(ier);
    
    /* Call f at (t0,y0), set zn[1] = y'(t0), 
       set initial h (from H0 or ARKHin), and scale zn[1] by h.
       Also check for zeros of root function g at and near t0.    */
    
    retval = f(tn, zn[0], zn[1], user_data); 
    nfe++;
    if (retval < 0) {
      ARKProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", "ARKode", MSGARK_RHSFUNC_FAILED, tn);
      return(ARK_RHSFUNC_FAIL);
    }
    if (retval > 0) {
      ARKProcessError(ark_mem, ARK_FIRST_RHSFUNC_ERR, "ARKODE", "ARKode", MSGARK_RHSFUNC_FIRST);
      return(ARK_FIRST_RHSFUNC_ERR);
    }

    /* Set initial h (from H0 or ARKHin). */

    h = hin;
    if ( (h != ZERO) && ((tout-tn)*h < ZERO) ) {
      ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", MSGARK_BAD_H0);
      return(ARK_ILL_INPUT);
    }
    if (h == ZERO) {
      tout_hin = tout;
      if ( tstopset && (tout-tn)*(tout-tstop) > 0 ) tout_hin = tstop; 
      hflag = ARKHin(ark_mem, tout_hin);
      if (hflag != ARK_SUCCESS) {
        istate = ARKHandleFailure(ark_mem, hflag);
        return(istate);
      }
    }
    rh = ABS(h)*hmax_inv;
    if (rh > ONE) h /= rh;
    if (ABS(h) < hmin) h *= hmin/ABS(h);

    /* Check for approach to tstop */

    if (tstopset) {
      if ( (tstop - tn)*h < ZERO ) {
        ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", MSGARK_BAD_TSTOP, tstop, tn);
        return(ARK_ILL_INPUT);
      }
      if ( (tn + h - tstop)*h > ZERO ) 
        h = (tstop - tn)*(ONE-FOUR*uround);
    }

    /* Scale zn[1] by h.*/

    hscale = h; 
    h0u    = h;
    hprime = h;

    N_VScale(h, zn[1], zn[1]);

    /* Check for zeros of root function g at and near t0. */

    if (nrtfn > 0) {

      retval = ARKRcheck1(ark_mem);

      if (retval == ARK_RTFUNC_FAIL) {
        ARKProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "ARKRcheck1", MSGARK_RTFUNC_FAILED, tn);
        return(ARK_RTFUNC_FAIL);
      }

    }

  } /* end of first call block */

  /*
   * ------------------------------------------------------
   * 3. At following steps, perform stop tests:
   *    - check for root in last step
   *    - check if we passed tstop
   *    - check if we passed tout (NORMAL mode)
   *    - check if current tn was returned (ONE_STEP mode)
   *    - check if we are close to tstop
   *      (adjust step size if needed)
   * -------------------------------------------------------
   */

  if (nst > 0) {

    /* Estimate an infinitesimal time interval to be used as
       a roundoff for time quantities (based on current time 
       and step size) */
    troundoff = FUZZ_FACTOR*uround*(ABS(tn) + ABS(h));

    /* First, check for a root in the last step taken, other than the
       last root found, if any.  If itask = ARK_ONE_STEP and y(tn) was not
       returned because of an intervening root, return y(tn) now.     */
    if (nrtfn > 0) {

      irfndp = irfnd;
      
      retval = ARKRcheck2(ark_mem);

      if (retval == CLOSERT) {
        ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKRcheck2", MSGARK_CLOSE_ROOTS, tlo);
        return(ARK_ILL_INPUT);
      } else if (retval == ARK_RTFUNC_FAIL) {
        ARKProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "ARKRcheck2", MSGARK_RTFUNC_FAILED, tlo);
        return(ARK_RTFUNC_FAIL);
      } else if (retval == RTFOUND) {
        tretlast = *tret = tlo;
        return(ARK_ROOT_RETURN);
      }

      /* If tn is distinct from tretlast (within roundoff),
         check remaining interval for roots */
      if ( ABS(tn - tretlast) > troundoff ) {

        retval = ARKRcheck3(ark_mem);

        if (retval == ARK_SUCCESS) {     /* no root found */
          irfnd = 0;
          if ((irfndp == 1) && (itask == ARK_ONE_STEP)) {
            tretlast = *tret = tn;
            N_VScale(ONE, zn[0], yout);
            return(ARK_SUCCESS);
          }
        } else if (retval == RTFOUND) {  /* a new root was found */
          irfnd = 1;
          tretlast = *tret = tlo;
          return(ARK_ROOT_RETURN);
        } else if (retval == ARK_RTFUNC_FAIL) {  /* g failed */
          ARKProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "ARKRcheck3", MSGARK_RTFUNC_FAILED, tlo);
          return(ARK_RTFUNC_FAIL);
        }

      }

    } /* end of root stop check */

    /* In ARK_NORMAL mode, test if tout was reached */
    if ( (itask == ARK_NORMAL) && ((tn-tout)*h >= ZERO) ) {
      tretlast = *tret = tout;
      ier =  ARKodeGetDky(ark_mem, tout, 0, yout);
      if (ier != ARK_SUCCESS) {
        ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", MSGARK_BAD_TOUT, tout);
        return(ARK_ILL_INPUT);
      }
      return(ARK_SUCCESS);
    }

    /* In ARK_ONE_STEP mode, test if tn was returned */
    if ( itask == ARK_ONE_STEP && ABS(tn - tretlast) > troundoff ) {
      tretlast = *tret = tn;
      N_VScale(ONE, zn[0], yout);
      return(ARK_SUCCESS);
    }

    /* Test for tn at tstop or near tstop */
    if ( tstopset ) {

      if ( ABS(tn - tstop) <= troundoff) {
        ier =  ARKodeGetDky(ark_mem, tstop, 0, yout);
        if (ier != ARK_SUCCESS) {
          ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", MSGARK_BAD_TSTOP, tstop, tn);
          return(ARK_ILL_INPUT);
        }
        tretlast = *tret = tstop;
        tstopset = FALSE;
        return(ARK_TSTOP_RETURN);
      }
      
      /* If next step would overtake tstop, adjust stepsize */
      if ( (tn + hprime - tstop)*h > ZERO ) {
        hprime = (tstop - tn)*(ONE-FOUR*uround);
        eta = hprime/h;
      }

    }
    
  } /* end stopping tests block */  

  /*
   * --------------------------------------------------
   * 4. Looping point for internal steps
   *
   *    4.1. check for errors (too many steps, too much
   *         accuracy requested, step size too small)
   *    4.2. take a new step (call ARKStep)
   *    4.3. stop on error 
   *    4.4. perform stop tests:
   *         - check for root in last step
   *         - check if tout was passed
   *         - check if close to tstop
   *         - check if in ONE_STEP mode (must return)
   * --------------------------------------------------
   */

  nstloc = 0;
  loop {
   
    next_h = h;
    next_q = q;
    
    /* Reset and check ewt */
    if (nst > 0) {

      ewtsetOK = efun(zn[0], ewt, e_data);

      if (ewtsetOK != 0) {

        if (itol == ARK_WF) 
          ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", MSGARK_EWT_NOW_FAIL, tn);
        else 
          ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKode", MSGARK_EWT_NOW_BAD, tn);
	
        istate = ARK_ILL_INPUT;
        tretlast = *tret = tn;
        N_VScale(ONE, zn[0], yout);
        break;

      }
    }
    
    /* Check for too many steps */
    if ( (mxstep>0) && (nstloc >= mxstep) ) {
      ARKProcessError(ark_mem, ARK_TOO_MUCH_WORK, "ARKODE", "ARKode", MSGARK_MAX_STEPS, tn);
      istate = ARK_TOO_MUCH_WORK;
      tretlast = *tret = tn;
      N_VScale(ONE, zn[0], yout);
      break;
    }

    /* Check for too much accuracy requested */
    nrm = N_VWrmsNorm(zn[0], ewt);
    tolsf = uround * nrm;
    if (tolsf > ONE) {
      ARKProcessError(ark_mem, ARK_TOO_MUCH_ACC, "ARKODE", "ARKode", MSGARK_TOO_MUCH_ACC, tn);
      istate = ARK_TOO_MUCH_ACC;
      tretlast = *tret = tn;
      N_VScale(ONE, zn[0], yout);
      tolsf *= TWO;
      break;
    } else {
      tolsf = ONE;
    }

    /* Check for h below roundoff level in tn */
    if (tn + h == tn) {
      nhnil++;
      if (nhnil <= mxhnil) 
        ARKProcessError(ark_mem, ARK_WARNING, "ARKODE", "ARKode", MSGARK_HNIL, tn, h);
      if (nhnil == mxhnil) 
        ARKProcessError(ark_mem, ARK_WARNING, "ARKODE", "ARKode", MSGARK_HNIL_DONE);
    }

    /* Call ARKStep to take a step */
    kflag = ARKStep(ark_mem);

    /* Process failed step cases, and exit loop */
    if (kflag != ARK_SUCCESS) {
      istate = ARKHandleFailure(ark_mem, kflag);
      tretlast = *tret = tn;
      N_VScale(ONE, zn[0], yout);
      break;
    }
    
    nstloc++;

    /* Check for root in last step taken. */
    if (nrtfn > 0) {

      retval = ARKRcheck3(ark_mem);

      if (retval == RTFOUND) {  /* A new root was found */
        irfnd = 1;
        istate = ARK_ROOT_RETURN;
        tretlast = *tret = tlo;
        break;
      } else if (retval == ARK_RTFUNC_FAIL) { /* g failed */
        ARKProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "ARKRcheck3", MSGARK_RTFUNC_FAILED, tlo);
        istate = ARK_RTFUNC_FAIL;
        break;
      }

      /* If we are at the end of the first step and we still have
       * some event functions that are inactive, issue a warning
       * as this may indicate a user error in the implementation
       * of the root function. */

      if (nst==1) {
        inactive_roots = FALSE;
        for (ir=0; ir<nrtfn; ir++) { 
          if (!gactive[ir]) {
            inactive_roots = TRUE;
            break;
          }
        }
        if ((ark_mem->ark_mxgnull > 0) && inactive_roots) {
          ARKProcessError(ark_mem, ARK_WARNING, "ARKODES", "ARKode", MSGARK_INACTIVE_ROOTS);
        }
      }

    }

    /* In NORMAL mode, check if tout reached */
    if ( (itask == ARK_NORMAL) &&  (tn-tout)*h >= ZERO ) {
      istate = ARK_SUCCESS;
      tretlast = *tret = tout;
      (void) ARKodeGetDky(ark_mem, tout, 0, yout);
      next_q = qprime;
      next_h = hprime;
      break;
    }

    /* Check if tn is at tstop or near tstop */
    if ( tstopset ) {

      troundoff = FUZZ_FACTOR*uround*(ABS(tn) + ABS(h));
      if ( ABS(tn - tstop) <= troundoff) {
        (void) ARKodeGetDky(ark_mem, tstop, 0, yout);
        tretlast = *tret = tstop;
        tstopset = FALSE;
        istate = ARK_TSTOP_RETURN;
        break;
      }

      if ( (tn + hprime - tstop)*h > ZERO ) {
        hprime = (tstop - tn)*(ONE-FOUR*uround);
        eta = hprime/h;
      }

    }

    /* In ONE_STEP mode, copy y and exit loop */
    if (itask == ARK_ONE_STEP) {
      istate = ARK_SUCCESS;
      tretlast = *tret = tn;
      N_VScale(ONE, zn[0], yout);
      next_q = qprime;
      next_h = hprime;
      break;
    }

  } /* end looping for internal steps */

  return(istate);
}

/*-----------------------------------------------------------------*/

/*
 * ARKodeGetDky
 *
 * This routine computes the k-th derivative of the interpolating
 * polynomial at the time t and stores the result in the vector dky.
 * The formula is:
 *         q 
 *  dky = SUM c(j,k) * (t - tn)^(j-k) * h^(-j) * zn[j] , 
 *        j=k 
 * where c(j,k) = j*(j-1)*...*(j-k+1), q is the current order, and
 * zn[j] is the j-th column of the Nordsieck history array.
 *
 * This function is called by ARKode with k = 0 and t = tout, but
 * may also be called directly by the user.
 */

int ARKodeGetDky(void *arkode_mem, realtype t, int k, N_Vector dky)
{
  realtype s, c, r;
  realtype tfuzz, tp, tn1;
  int i, j;
  ARKodeMem ark_mem;
  
  /* Check all inputs for legality */
 
  if (arkode_mem == NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", "ARKodeGetDky", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  if (dky == NULL) {
    ARKProcessError(ark_mem, ARK_BAD_DKY, "ARKODE", "ARKodeGetDky", MSGARK_NULL_DKY);
    return(ARK_BAD_DKY);
  }

  if ((k < 0) || (k > q)) {
    ARKProcessError(ark_mem, ARK_BAD_K, "ARKODE", "ARKodeGetDky", MSGARK_BAD_K);
    return(ARK_BAD_K);
  }
  
  /* Allow for some slack */
  tfuzz = FUZZ_FACTOR * uround * (ABS(tn) + ABS(hu));
  if (hu < ZERO) tfuzz = -tfuzz;
  tp = tn - hu - tfuzz;
  tn1 = tn + tfuzz;
  if ((t-tp)*(t-tn1) > ZERO) {
    ARKProcessError(ark_mem, ARK_BAD_T, "ARKODE", "ARKodeGetDky", MSGARK_BAD_T, t, tn-hu, tn);
    return(ARK_BAD_T);
  }

  /* Sum the differentiated interpolating polynomial */

  s = (t - tn) / h;
  for (j=q; j >= k; j--) {
    c = ONE;
    for (i=j; i >= j-k+1; i--) c *= i;
    if (j == q) {
      N_VScale(c, zn[q], dky);
    } else {
      N_VLinearSum(c, zn[j], s, dky, dky);
    }
  }
  if (k == 0) return(ARK_SUCCESS);
  r = RPowerI(h,-k);
  N_VScale(r, dky, dky);
  return(ARK_SUCCESS);
}

/*
 * ARKodeFree
 *
 * This routine frees the problem memory allocated by ARKodeInit.
 * Such memory includes all the vectors allocated by ARKAlloarkectors,
 * and the memory lmem for the linear solver (deallocated by a call
 * to lfree).
 */

void ARKodeFree(void **arkode_mem)
{
  ARKodeMem ark_mem;

  if (*arkode_mem == NULL) return;

  ark_mem = (ARKodeMem) (*arkode_mem);
  
  ARKFreeVectors(ark_mem);

  if (iter == ARK_NEWTON && lfree != NULL) lfree(ark_mem);

  if (nrtfn > 0) {
    free(glo); glo = NULL;
    free(ghi); ghi = NULL;
    free(grout); grout = NULL;
    free(iroots); iroots = NULL;
    free(rootdir); rootdir = NULL;
    free(gactive); gactive = NULL;
  }

  free(*arkode_mem);
  *arkode_mem = NULL;
}

/* 
 * =================================================================
 *  Private Functions Implementation
 * =================================================================
 */

/*
 * ARKCheckNvector
 * This routine checks if all required vector operations are present.
 * If any of them is missing it returns FALSE.
 */

static booleantype ARKCheckNvector(N_Vector tmpl)
{
  if((tmpl->ops->nvclone     == NULL) ||
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

/*
 * ARKAlloarkectors
 *
 * This routine allocates the ARKODE vectors ewt, acor, tempv, ftemp, and
 * zn[0], ..., zn[maxord].
 * If all memory allocations are successful, ARKAlloarkectors returns TRUE. 
 * Otherwise all allocated memory is freed and ARKAlloarkectors returns FALSE.
 * This routine also sets the optional outputs lrw and liw, which are
 * (respectively) the lengths of the real and integer work spaces
 * allocated here.
 */

static booleantype ARKAlloarkectors(ARKodeMem ark_mem, N_Vector tmpl)
{
  int i, j;

  /* Allocate ewt, acor, tempv, ftemp */
  
  ewt = N_VClone(tmpl);
  if (ewt == NULL) return(FALSE);

  acor = N_VClone(tmpl);
  if (acor == NULL) {
    N_VDestroy(ewt);
    return(FALSE);
  }

  tempv = N_VClone(tmpl);
  if (tempv == NULL) {
    N_VDestroy(ewt);
    N_VDestroy(acor);
    return(FALSE);
  }

  ftemp = N_VClone(tmpl);
  if (ftemp == NULL) {
    N_VDestroy(tempv);
    N_VDestroy(ewt);
    N_VDestroy(acor);
    return(FALSE);
  }

  /* Allocate zn[0] ... zn[qmax] */

  for (j=0; j <= qmax; j++) {
    zn[j] = N_VClone(tmpl);
    if (zn[j] == NULL) {
      N_VDestroy(ewt);
      N_VDestroy(acor);
      N_VDestroy(tempv);
      N_VDestroy(ftemp);
      for (i=0; i < j; i++) N_VDestroy(zn[i]);
      return(FALSE);
    }
  }

  /* Update solver workspace lengths  */
  lrw += (qmax + 5)*lrw1;
  liw += (qmax + 5)*liw1;

  /* Store the value of qmax used here */
  ark_mem->ark_qmax_alloc = qmax;

  return(TRUE);
}

/*  
 * ARKFreeVectors
 *
 * This routine frees the ARKODE vectors allocated in ARKAlloarkectors.
 */

static void ARKFreeVectors(ARKodeMem ark_mem)
{
  int j, maxord;
  
  maxord = ark_mem->ark_qmax_alloc;

  N_VDestroy(ewt);
  N_VDestroy(acor);
  N_VDestroy(tempv);
  N_VDestroy(ftemp);
  for(j=0; j <= maxord; j++) N_VDestroy(zn[j]);

  lrw -= (maxord + 5)*lrw1;
  liw -= (maxord + 5)*liw1;

  if (ark_mem->ark_VabstolMallocDone) {
    N_VDestroy(Vabstol);
    lrw -= lrw1;
    liw -= liw1;
  }
}

/*  
 * ARKInitialSetup
 *
 * This routine performs input consistency checks at the first step.
 * If needed, it also checks the linear solver module and calls the
 * linear solver initialization routine.
 */

static int ARKInitialSetup(ARKodeMem ark_mem)
{
  int ier;

  /* Did the user specify tolerances? */
  if (itol == ARK_NN) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKInitialSetup", MSGARK_NO_TOLS);
    return(ARK_ILL_INPUT);
  }

  /* Set data for efun */
  if (ark_mem->ark_user_efun) e_data = user_data;
  else                      e_data = ark_mem;

  /* Load initial error weights */
  ier = efun(zn[0], ewt, e_data);
  if (ier != 0) {
    if (itol == ARK_WF) 
      ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKInitialSetup", MSGARK_EWT_FAIL);
    else 
      ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKInitialSetup", MSGARK_BAD_EWT);
    return(ARK_ILL_INPUT);
  }
  
  /* Check if lsolve function exists (if needed) and call linit function (if it exists) */
  if (iter == ARK_NEWTON) {
    if (lsolve == NULL) {
      ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", "ARKInitialSetup", MSGARK_LSOLVE_NULL);
      return(ARK_ILL_INPUT);
    }
    if (linit != NULL) {
      ier = linit(ark_mem);
      if (ier != 0) {
        ARKProcessError(ark_mem, ARK_LINIT_FAIL, "ARKODE", "ARKInitialSetup", MSGARK_LINIT_FAIL);
        return(ARK_LINIT_FAIL);
      }
    }
  }

  return(ARK_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * PRIVATE FUNCTIONS FOR ARKODE
 * -----------------------------------------------------------------
 */

/*
 * ARKHin
 *
 * This routine computes a tentative initial step size h0. 
 * If tout is too close to tn (= t0), then ARKHin returns ARK_TOO_CLOSE
 * and h remains uninitialized. Note that here tout is either the value
 * passed to ARKode at the first call or the value of tstop (if tstop is 
 * enabled and it is closer to t0=tn than tout).
 * If the RHS function fails unrecoverably, ARKHin returns ARK_RHSFUNC_FAIL.
 * If the RHS function fails recoverably too many times and recovery is
 * not possible, ARKHin returns ARK_REPTD_RHSFUNC_ERR.
 * Otherwise, ARKHin sets h to the chosen value h0 and returns ARK_SUCCESS.
 *
 * The algorithm used seeks to find h0 as a solution of
 *       (WRMS norm of (h0^2 ydd / 2)) = 1, 
 * where ydd = estimated second derivative of y.
 *
 * We start with an initial estimate equal to the geometric mean of the
 * lower and upper bounds on the step size.
 *
 * Loop up to MAX_ITERS times to find h0.
 * Stop if new and previous values differ by a factor < 2.
 * Stop if hnew/hg > 2 after one iteration, as this probably means
 * that the ydd value is bad because of cancellation error.        
 *  
 * For each new proposed hg, we allow MAX_ITERS attempts to
 * resolve a possible recoverable failure from f() by reducing
 * the proposed stepsize by a factor of 0.2. If a legal stepsize
 * still cannot be found, fall back on a previous value if possible,
 * or else return ARK_REPTD_RHSFUNC_ERR.
 *
 * Finally, we apply a bias (0.5) and verify that h0 is within bounds.
 */

static int ARKHin(ARKodeMem ark_mem, realtype tout)
{
  int retval, sign, count1, count2;
  realtype tdiff, tdist, tround, hlb, hub;
  realtype hg, hgs, hs, hnew, hrat, h0, yddnrm;
  booleantype hgOK, hnewOK;

  /* If tout is too close to tn, give up */
  
  if ((tdiff = tout-tn) == ZERO) return(ARK_TOO_CLOSE);
  
  sign = (tdiff > ZERO) ? 1 : -1;
  tdist = ABS(tdiff);
  tround = uround * MAX(ABS(tn), ABS(tout));

  if (tdist < TWO*tround) return(ARK_TOO_CLOSE);
  
  /* 
     Set lower and upper bounds on h0, and take geometric mean 
     as first trial value.
     Exit with this value if the bounds cross each other.
  */

  hlb = HLB_FACTOR * tround;
  hub = ARKUpperBoundH0(ark_mem, tdist);

  hg  = RSqrt(hlb*hub);

  if (hub < hlb) {
    if (sign == -1) h = -hg;
    else            h =  hg;
    return(ARK_SUCCESS);
  }
  
  /* Outer loop */

  hnewOK = FALSE;
  hs = hg;         /* safeguard against 'uninitialized variable' warning */

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

    /* If the stopping criteria was met, or if this is the last pass, stop */
    if ( (hnewOK) || (count1 == MAX_ITERS))  {hnew = hg; break;}

    /* Propose new step size */
    hnew = (yddnrm*hub*hub > TWO) ? RSqrt(TWO/yddnrm) : RSqrt(hg*hub);
    hrat = hnew/hg;
    
    /* Accept hnew if it does not differ from hg by more than a factor of 2 */
    if ((hrat > HALF) && (hrat < TWO)) {
      hnewOK = TRUE;
    }

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
  h = h0;

  return(ARK_SUCCESS);
}

/*
 * ARKUpperBoundH0
 *
 * This routine sets an upper bound on abs(h0) based on
 * tdist = tn - t0 and the values of y[i]/y'[i].
 */

static realtype ARKUpperBoundH0(ARKodeMem ark_mem, realtype tdist)
{
  realtype hub_inv, hub;
  N_Vector temp1, temp2;

  /* 
   * Bound based on |y0|/|y0'| -- allow at most an increase of
   * HUB_FACTOR in y0 (based on a forward Euler step). The weight 
   * factor is used as a safeguard against zero components in y0. 
   */

  temp1 = tempv;
  temp2 = acor;

  N_VAbs(zn[0], temp2);
  efun(zn[0], temp1, e_data);
  N_VInv(temp1, temp1);
  N_VLinearSum(HUB_FACTOR, temp2, ONE, temp1, temp1);

  N_VAbs(zn[1], temp2);

  N_VDiv(temp2, temp1, temp1);
  hub_inv = N_VMaxNorm(temp1);

  /*
   * bound based on tdist -- allow at most a step of magnitude
   * HUB_FACTOR * tdist
   */

  hub = HUB_FACTOR*tdist;

  /* Use the smaler of the two */

  if (hub*hub_inv > ONE) hub = ONE/hub_inv;

  return(hub);
}

/*
 * ARKYddNorm
 *
 * This routine computes an estimate of the second derivative of y
 * using a difference quotient, and returns its WRMS norm.
 */

static int ARKYddNorm(ARKodeMem ark_mem, realtype hg, realtype *yddnrm)
{
  int retval;

  N_VLinearSum(hg, zn[1], ONE, zn[0], y);
  retval = f(tn+hg, y, tempv, user_data);
  nfe++;
  if (retval < 0) return(ARK_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  N_VLinearSum(ONE, tempv, -ONE, zn[1], tempv);
  N_VScale(ONE/hg, tempv, tempv);

  *yddnrm = N_VWrmsNorm(tempv, ewt);

  return(ARK_SUCCESS);
}

/* 
 * ARKStep
 *
 * This routine performs one internal arkode step, from tn to tn + h.
 * It calls other routines to do all the work.
 *
 * The main operations done here are as follows:
 * - preliminary adjustments if a new step size was chosen;
 * - prediction of the Nordsieck history array zn at tn + h;
 * - setting of multistep method coefficients and test quantities;
 * - solution of the nonlinear system;
 * - testing the local error;
 * - updating zn and other state data if successful;
 * - resetting stepsize and order for the next step.
 * - if SLDET is on, check for stability, reduce order if necessary.
 * On a failure in the nonlinear system solution or error test, the
 * step may be reattempted, depending on the nature of the failure.
 */

static int ARKStep(ARKodeMem ark_mem)
{
  realtype saved_t, dsm;
  int ncf, nef;
  int nflag, kflag, eflag;
  
  saved_t = tn;
  ncf = nef = 0;
  nflag = FIRST_CALL;

  if ((nst > 0) && (hprime != h)) ARKAdjustParams(ark_mem);
  
  /* Looping point for attempts to take a step */
  loop {  

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

  /* If Stablilty Limit Detection is turned on, call stability limit
     detection routine for possible order reduction. */

  if (sldeton) ARKBDFStab(ark_mem);

  etamax = (nst <= SMALL_NST) ? ETAMX2 : ETAMX3;

  /*  Finally, we rescale the acor array to be the 
      estimated local error vector. */

  N_VScale(tq[2], acor, acor);
  return(ARK_SUCCESS);
      
}

/*
 * ARKAdjustParams
 *
 * This routine is called when a change in step size was decided upon,
 * and it handles the required adjustments to the history array zn.
 * If there is to be a change in order, we call ARKAdjustOrder and reset
 * q, L = q+1, and qwait.  Then in any case, we call ARKRescale, which
 * resets h and rescales the Nordsieck array.
 */

static void ARKAdjustParams(ARKodeMem ark_mem)
{
  if (qprime != q) {
    ARKAdjustOrder(ark_mem, qprime-q);
    q = qprime;
    L = q+1;
    qwait = L;
  }
  ARKRescale(ark_mem);
}

/*
 * ARKAdjustOrder
 *
 * This routine is a high level routine which handles an order
 * change by an amount deltaq (= +1 or -1). If a decrease in order
 * is requested and q==2, then the routine returns immediately.
 * Otherwise ARKAdjustAdams or ARKAdjustBDF is called to handle the
 * order change (depending on the value of lmm).
 */

static void ARKAdjustOrder(ARKodeMem ark_mem, int deltaq)
{
  if ((q==2) && (deltaq != 1)) return;
  
  switch(lmm){
  case ARK_ADAMS: 
    ARKAdjustAdams(ark_mem, deltaq);
    break;
  case ARK_BDF:   
    ARKAdjustBDF(ark_mem, deltaq);
    break;
  }
}

/*
 * ARKAdjustAdams
 *
 * This routine adjusts the history array on a change of order q by
 * deltaq, in the case that lmm == ARK_ADAMS.
 */

static void ARKAdjustAdams(ARKodeMem ark_mem, int deltaq)
{
  int i, j;
  realtype xi, hsum;

  /* On an order increase, set new column of zn to zero and return */
  
  if (deltaq==1) {
    N_VConst(ZERO, zn[L]);
    return;
  }

  /*
   * On an order decrease, each zn[j] is adjusted by a multiple of zn[q].
   * The coeffs. in the adjustment are the coeffs. of the polynomial:
   *        x
   * q * INT { u * ( u + xi_1 ) * ... * ( u + xi_{q-2} ) } du 
   *        0
   * where xi_j = [t_n - t_(n-j)]/h => xi_0 = 0
   */

  for (i=0; i <= qmax; i++) l[i] = ZERO;
  l[1] = ONE;
  hsum = ZERO;
  for (j=1; j <= q-2; j++) {
    hsum += tau[j];
    xi = hsum / hscale;
    for (i=j+1; i >= 1; i--) l[i] = l[i]*xi + l[i-1];
  }
  
  for (j=1; j <= q-2; j++) l[j+1] = q * (l[j] / (j+1));
  
  for (j=2; j < q; j++)
    N_VLinearSum(-l[j], zn[q], ONE, zn[j], zn[j]);
}

/*
 * ARKAdjustBDF
 *
 * This is a high level routine which handles adjustments to the
 * history array on a change of order by deltaq in the case that 
 * lmm == ARK_BDF.  ARKAdjustBDF calls ARKIncreaseBDF if deltaq = +1 and 
 * ARKDecreaseBDF if deltaq = -1 to do the actual work.
 */

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

/*
 * ARKIncreaseBDF
 *
 * This routine adjusts the history array on an increase in the 
 * order q in the case that lmm == ARK_BDF.  
 * A new column zn[q+1] is set equal to a multiple of the saved 
 * vector (= acor) in zn[indx_acor].  Then each zn[j] is adjusted by
 * a multiple of zn[q+1].  The coefficients in the adjustment are the 
 * coefficients of the polynomial x*x*(x+xi_1)*...*(x+xi_j),
 * where xi_j = [t_n - t_(n-j)]/h.
 */

static void ARKIncreaseBDF(ARKodeMem ark_mem)
{
  realtype alpha0, alpha1, prod, xi, xiold, hsum, A1;
  int i, j;
  
  for (i=0; i <= qmax; i++) l[i] = ZERO;
  l[2] = alpha1 = prod = xiold = ONE;
  alpha0 = -ONE;
  hsum = hscale;
  if (q > 1) {
    for (j=1; j < q; j++) {
      hsum += tau[j+1];
      xi = hsum / hscale;
      prod *= xi;
      alpha0 -= ONE / (j+1);
      alpha1 += ONE / xi;
      for (i=j+2; i >= 2; i--) l[i] = l[i]*xiold + l[i-1];
      xiold = xi;
    }
  }
  A1 = (-alpha0 - alpha1) / prod;
  N_VScale(A1, zn[indx_acor], zn[L]);
  for (j=2; j <= q; j++) {
    N_VLinearSum(l[j], zn[L], ONE, zn[j], zn[j]);
  }  
}

/*
 * ARKDecreaseBDF
 *
 * This routine adjusts the history array on a decrease in the 
 * order q in the case that lmm == ARK_BDF.  
 * Each zn[j] is adjusted by a multiple of zn[q].  The coefficients
 * in the adjustment are the coefficients of the polynomial
 *   x*x*(x+xi_1)*...*(x+xi_j), where xi_j = [t_n - t_(n-j)]/h.
 */

static void ARKDecreaseBDF(ARKodeMem ark_mem)
{
  realtype hsum, xi;
  int i, j;
  
  for (i=0; i <= qmax; i++) l[i] = ZERO;
  l[2] = ONE;
  hsum = ZERO;
  for(j=1; j <= q-2; j++) {
    hsum += tau[j];
    xi = hsum /hscale;
    for (i=j+2; i >= 2; i--) l[i] = l[i]*xi + l[i-1];
  }
  
  for(j=2; j < q; j++)
    N_VLinearSum(-l[j], zn[q], ONE, zn[j], zn[j]);
}

/*
 * ARKRescale
 *
 * This routine rescales the Nordsieck array by multiplying the
 * jth column zn[j] by eta^j, j = 1, ..., q.  Then the value of
 * h is rescaled by eta, and hscale is reset to h.
 */

static void ARKRescale(ARKodeMem ark_mem)
{
  int j;
  realtype factor;
  
  factor = eta;
  for (j=1; j <= q; j++) {
    N_VScale(factor, zn[j], zn[j]);
    factor *= eta;
  }
  h = hscale * eta;
  next_h = h;
  hscale = h;
  nscon = 0;
}

/*
 * ARKPredict
 *
 * This routine advances tn by the tentative step size h, and computes
 * the predicted array z_n(0), which is overwritten on zn.  The
 * prediction of zn is done by repeated additions.
 * If tstop is enabled, it is possible for tn + h to be past tstop by roundoff,
 * and in that case, we reset tn (after incrementing by h) to tstop.
 */

static void ARKPredict(ARKodeMem ark_mem)
{
  int j, k;
  
  tn += h;
  if (tstopset) {
    if ((tn - tstop)*h > ZERO) tn = tstop;
  }
  for (k = 1; k <= q; k++)
    for (j = q; j >= k; j--) 
      N_VLinearSum(ONE, zn[j-1], ONE, zn[j], zn[j-1]); 
}

/*
 * ARKSet
 *
 * This routine is a high level routine which calls ARKSetAdams or
 * ARKSetBDF to set the polynomial l, the test quantity array tq, 
 * and the related variables  rl1, gamma, and gamrat.
 *
 * The array tq is loaded with constants used in the control of estimated
 * local errors and in the nonlinear convergence test.  Specifically, while
 * running at order q, the components of tq are as follows:
 *   tq[1] = a coefficient used to get the est. local error at order q-1
 *   tq[2] = a coefficient used to get the est. local error at order q
 *   tq[3] = a coefficient used to get the est. local error at order q+1
 *   tq[4] = constant used in nonlinear iteration convergence test
 *   tq[5] = coefficient used to get the order q+2 derivative vector used in
 *           the est. local error at order q+1
 */

static void ARKSet(ARKodeMem ark_mem)
{
  switch(lmm) {
  case ARK_ADAMS: 
    ARKSetAdams(ark_mem);
    break;
  case ARK_BDF:
    ARKSetBDF(ark_mem);
    break;
  }
  rl1 = ONE / l[1];
  gamma = h * rl1;
  if (nst == 0) gammap = gamma;
  gamrat = (nst > 0) ? gamma / gammap : ONE;  /* protect x / x != 1.0 */
}

/*
 * ARKSetAdams
 *
 * This routine handles the computation of l and tq for the
 * case lmm == ARK_ADAMS.
 *
 * The components of the array l are the coefficients of a
 * polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
 *                          q-1
 * (d/dx) Lambda(x) = c * PRODUCT (1 + x / xi_i) , where
 *                          i=1
 *  Lambda(-1) = 0, Lambda(0) = 1, and c is a normalization factor.
 * Here xi_i = [t_n - t_(n-i)] / h.
 *
 * The array tq is set to test quantities used in the convergence
 * test, the error test, and the selection of h at a new order.
 */

static void ARKSetAdams(ARKodeMem ark_mem)
{
  realtype m[L_MAX], M[3], hsum;
  
  if (q == 1) {
    l[0] = l[1] = tq[1] = tq[5] = ONE;
    tq[2] = HALF;
    tq[3] = ONE/TWELVE;
    tq[4] = nlscoef / tq[2];       /* = 0.1 / tq[2] */
    return;
  }
  
  hsum = ARKAdamsStart(ark_mem, m);
  
  M[0] = ARKAltSum(q-1, m, 1);
  M[1] = ARKAltSum(q-1, m, 2);
  
  ARKAdamsFinish(ark_mem, m, M, hsum);
}

/*
 * ARKAdamsStart
 *
 * This routine generates in m[] the coefficients of the product
 * polynomial needed for the Adams l and tq coefficients for q > 1.
 */

static realtype ARKAdamsStart(ARKodeMem ark_mem, realtype m[])
{
  realtype hsum, xi_inv, sum;
  int i, j;
  
  hsum = h;
  m[0] = ONE;
  for (i=1; i <= q; i++) m[i] = ZERO;
  for (j=1; j < q; j++) {
    if ((j==q-1) && (qwait == 1)) {
      sum = ARKAltSum(q-2, m, 2);
      tq[1] = q * sum / m[q-2];
    }
    xi_inv = h / hsum;
    for (i=j; i >= 1; i--) m[i] += m[i-1] * xi_inv;
    hsum += tau[j];
    /* The m[i] are coefficients of product(1 to j) (1 + x/xi_i) */
  }
  return(hsum);
}

/*
 * ARKAdamsFinish
 *
 * This routine completes the calculation of the Adams l and tq.
 */

static void ARKAdamsFinish(ARKodeMem ark_mem, realtype m[], realtype M[], realtype hsum)
{
  int i;
  realtype M0_inv, xi, xi_inv;
  
  M0_inv = ONE / M[0];
  
  l[0] = ONE;
  for (i=1; i <= q; i++) l[i] = M0_inv * (m[i-1] / i);
  xi = hsum / h;
  xi_inv = ONE / xi;
  
  tq[2] = M[1] * M0_inv / xi;
  tq[5] = xi / l[q];

  if (qwait == 1) {
    for (i=q; i >= 1; i--) m[i] += m[i-1] * xi_inv;
    M[2] = ARKAltSum(q, m, 2);
    tq[3] = M[2] * M0_inv / L;
  }

  tq[4] = nlscoef / tq[2];
}

/*  
 * ARKAltSum
 *
 * ARKAltSum returns the value of the alternating sum
 *   sum (i= 0 ... iend) [ (-1)^i * (a[i] / (i + k)) ].
 * If iend < 0 then ARKAltSum returns 0.
 * This operation is needed to compute the integral, from -1 to 0,
 * of a polynomial x^(k-1) M(x) given the coefficients of M(x).
 */

static realtype ARKAltSum(int iend, realtype a[], int k)
{
  int i, sign;
  realtype sum;
  
  if (iend < 0) return(ZERO);
  
  sum = ZERO;
  sign = 1;
  for (i=0; i <= iend; i++) {
    sum += sign * (a[i] / (i+k));
    sign = -sign;
  }
  return(sum);
}

/*
 * ARKSetBDF
 *
 * This routine computes the coefficients l and tq in the case
 * lmm == ARK_BDF.  ARKSetBDF calls ARKSetTqBDF to set the test
 * quantity array tq. 
 * 
 * The components of the array l are the coefficients of a
 * polynomial Lambda(x) = l_0 + l_1 x + ... + l_q x^q, given by
 *                                 q-1
 * Lambda(x) = (1 + x / xi*_q) * PRODUCT (1 + x / xi_i) , where
 *                                 i=1
 *  xi_i = [t_n - t_(n-i)] / h.
 *
 * The array tq is set to test quantities used in the convergence
 * test, the error test, and the selection of h at a new order.
 */

static void ARKSetBDF(ARKodeMem ark_mem)
{
  realtype alpha0, alpha0_hat, xi_inv, xistar_inv, hsum;
  int i,j;
  
  l[0] = l[1] = xi_inv = xistar_inv = ONE;
  for (i=2; i <= q; i++) l[i] = ZERO;
  alpha0 = alpha0_hat = -ONE;
  hsum = h;
  if (q > 1) {
    for (j=2; j < q; j++) {
      hsum += tau[j-1];
      xi_inv = h / hsum;
      alpha0 -= ONE / j;
      for(i=j; i >= 1; i--) l[i] += l[i-1]*xi_inv;
      /* The l[i] are coefficients of product(1 to j) (1 + x/xi_i) */
    }
    
    /* j = q */
    alpha0 -= ONE / q;
    xistar_inv = -l[1] - alpha0;
    hsum += tau[q-1];
    xi_inv = h / hsum;
    alpha0_hat = -l[1] - xi_inv;
    for (i=q; i >= 1; i--) l[i] += l[i-1]*xistar_inv;
  }

  ARKSetTqBDF(ark_mem, hsum, alpha0, alpha0_hat, xi_inv, xistar_inv);
}

/*
 * ARKSetTqBDF
 *
 * This routine sets the test quantity array tq in the case
 * lmm == ARK_BDF.
 */

static void ARKSetTqBDF(ARKodeMem ark_mem, realtype hsum, realtype alpha0,
                       realtype alpha0_hat, realtype xi_inv, realtype xistar_inv)
{
  realtype A1, A2, A3, A4, A5, A6;
  realtype C, Cpinv, Cppinv;
  
  A1 = ONE - alpha0_hat + alpha0;
  A2 = ONE + q * A1;
  tq[2] = ABS(A1 / (alpha0 * A2));
  tq[5] = ABS(A2 * xistar_inv / (l[q] * xi_inv));
  if (qwait == 1) {
    if (q > 1) {
      C = xistar_inv / l[q];
      A3 = alpha0 + ONE / q;
      A4 = alpha0_hat + xi_inv;
      Cpinv = (ONE - A4 + A3) / A3;
      tq[1] = ABS(C * Cpinv);
    }
    else tq[1] = ONE;
    hsum += tau[q];
    xi_inv = h / hsum;
    A5 = alpha0 - (ONE / (q+1));
    A6 = alpha0_hat - xi_inv;
    Cppinv = (ONE - A6 + A5) / A2;
    tq[3] = ABS(Cppinv / (xi_inv * (q+2) * A5));
  }
  tq[4] = nlscoef / tq[2];
}

/*
 * ARKNls
 *
 * This routine attempts to solve the nonlinear system associated
 * with a single implicit step of the linear multistep method.
 * Depending on iter, it calls ARKNlsFunctional or ARKNlsNewton
 * to do the work.
 */

static int ARKNls(ARKodeMem ark_mem, int nflag)
{
  int flag = ARK_SUCCESS;

  switch(iter) {
  case ARK_FUNCTIONAL: 
    flag = ARKNlsFunctional(ark_mem);
    break;
  case ARK_NEWTON:
    flag = ARKNlsNewton(ark_mem, nflag);
    break;
  }

  return(flag);
}

/*
 * ARKNlsFunctional
 *
 * This routine attempts to solve the nonlinear system using 
 * functional iteration (no matrices involved).
 *
 * Possible return values are:
 *
 *   ARK_SUCCESS      --->  continue with error test
 *
 *   ARK_RHSFUNC_FAIL --->  halt the integration
 *
 *   CONV_FAIL       -+
 *   RHSFUNC_RECVR   -+->  predict again or stop if too many
 *
 */

static int ARKNlsFunctional(ARKodeMem ark_mem)
{
  int retval, m;
  realtype del, delp, dcon;

  /* Initialize counter and evaluate f at predicted y */
  
  crate = ONE;
  m = 0;

  retval = f(tn, zn[0], tempv, user_data);
  nfe++;
  if (retval < 0) return(ARK_RHSFUNC_FAIL);
  if (retval > 0) return(RHSFUNC_RECVR);

  N_VConst(ZERO, acor);

  /* Initialize delp to avoid compiler warning message */
  del = delp = ZERO;

  /* Loop until convergence; accumulate corrections in acor */

  loop {

    nni++;

    /* Correct y directly from the last f value */
    N_VLinearSum(h, tempv, -ONE, zn[1], tempv);
    N_VScale(rl1, tempv, tempv);
    N_VLinearSum(ONE, zn[0], ONE, tempv, y);
    /* Get WRMS norm of current correction to use in convergence test */
    N_VLinearSum(ONE, tempv, -ONE, acor, acor);
    del = N_VWrmsNorm(acor, ewt);
    N_VScale(ONE, tempv, acor);
    
    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.        */
    if (m > 0) crate = MAX(CRDOWN * crate, del / delp);
    dcon = del * MIN(ONE, crate) / tq[4];
    if (dcon <= ONE) {
      acnrm = (m == 0) ? del : N_VWrmsNorm(acor, ewt);
      return(ARK_SUCCESS);  /* Convergence achieved */
    }

    /* Stop at maxcor iterations or if iter. seems to be diverging */
    m++;
    if ((m==maxcor) || ((m >= 2) && (del > RDIV * delp))) return(CONV_FAIL);

    /* Save norm of correction, evaluate f, and loop again */
    delp = del;

    retval = f(tn, y, tempv, user_data);
    nfe++;
    if (retval < 0) return(ARK_RHSFUNC_FAIL);
    if (retval > 0) return(RHSFUNC_RECVR);

  }
}

/*
 * ARKNlsNewton
 *
 * This routine handles the Newton iteration. It calls lsetup if 
 * indicated, calls ARKNewtonIteration to perform the iteration, and 
 * retries a failed attempt at Newton iteration if that is indicated.
 *
 * Possible return values:
 *
 *   ARK_SUCCESS       ---> continue with error test
 *
 *   ARK_RHSFUNC_FAIL  -+  
 *   ARK_LSETUP_FAIL    |-> halt the integration 
 *   ARK_LSOLVE_FAIL   -+
 *
 *   CONV_FAIL        -+
 *   RHSFUNC_RECVR    -+-> predict again or stop if too many
 *
 */

static int ARKNlsNewton(ARKodeMem ark_mem, int nflag)
{
  N_Vector vtemp1, vtemp2, vtemp3;
  int convfail, retval, ier;
  booleantype callSetup;
  
  vtemp1 = acor;  /* rename acor as vtemp1 for readability  */
  vtemp2 = y;     /* rename y as vtemp2 for readability     */
  vtemp3 = tempv; /* rename tempv as vtemp3 for readability */
  
  /* Set flag convfail, input to lsetup for its evaluation decision */
  convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL)) ?
    ARK_NO_FAILURES : ARK_FAIL_OTHER;

  /* Decide whether or not to call setup routine (if one exists) */
  if (setupNonNull) {      
    callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL) ||
      (nst == 0) || (nst >= nstlp + MSBP) || (ABS(gamrat-ONE) > DGMAX);
  } else {  
    crate = ONE;
    callSetup = FALSE;
  }
  
  /* Looping point for the solution of the nonlinear system.
     Evaluate f at the predicted y, call lsetup if indicated, and
     call ARKNewtonIteration for the Newton iteration itself.      */
  
  loop {

    retval = f(tn, zn[0], ftemp, user_data);
    nfe++; 
    if (retval < 0) return(ARK_RHSFUNC_FAIL);
    if (retval > 0) return(RHSFUNC_RECVR);

    if (callSetup) {
      ier = lsetup(ark_mem, convfail, zn[0], ftemp, &jcur, 
                   vtemp1, vtemp2, vtemp3);
      nsetups++;
      callSetup = FALSE;
      gamrat = crate = ONE; 
      gammap = gamma;
      nstlp = nst;
      /* Return if lsetup failed */
      if (ier < 0) return(ARK_LSETUP_FAIL);
      if (ier > 0) return(CONV_FAIL);
    }

    /* Set acor to zero and load prediction into y vector */
    N_VConst(ZERO, acor);
    N_VScale(ONE, zn[0], y);

    /* Do the Newton iteration */
    ier = ARKNewtonIteration(ark_mem);

    /* If there is a convergence failure and the Jacobian-related 
       data appears not to be current, loop again with a call to lsetup
       in which convfail=ARK_FAIL_BAD_J.  Otherwise return.                 */
    if (ier != TRY_AGAIN) return(ier);
    
    callSetup = TRUE;
    convfail = ARK_FAIL_BAD_J;
  }
}

/*
 * ARKNewtonIteration
 *
 * This routine performs the Newton iteration. If the iteration succeeds,
 * it returns the value ARK_SUCCESS. If not, it may signal the ARKNlsNewton 
 * routine to call lsetup again and reattempt the iteration, by
 * returning the value TRY_AGAIN. (In this case, ARKNlsNewton must set 
 * convfail to ARK_FAIL_BAD_J before calling setup again). 
 * Otherwise, this routine returns one of the appropriate values 
 * ARK_LSOLVE_FAIL, ARK_RHSFUNC_FAIL, CONV_FAIL, or RHSFUNC_RECVR back 
 * to ARKNlsNewton.
 */

static int ARKNewtonIteration(ARKodeMem ark_mem)
{
  int m, retval;
  realtype del, delp, dcon;
  N_Vector b;

  mnewt = m = 0;

  /* Initialize delp to avoid compiler warning message */
  del = delp = ZERO;

  /* Looping point for Newton iteration */
  loop {

    /* Evaluate the residual of the nonlinear system*/
    N_VLinearSum(rl1, zn[1], ONE, acor, tempv);
    N_VLinearSum(gamma, ftemp, -ONE, tempv, tempv);

    /* Call the lsolve function */
    b = tempv;
    retval = lsolve(ark_mem, b, ewt, y, ftemp); 
    nni++;
    
    if (retval < 0) return(ARK_LSOLVE_FAIL);
    
    /* If lsolve had a recoverable failure and Jacobian data is
       not current, signal to try the solution again            */
    if (retval > 0) { 
      if ((!jcur) && (setupNonNull)) return(TRY_AGAIN);
      else                           return(CONV_FAIL);
    }

    /* Get WRMS norm of correction; add correction to acor and y */
    del = N_VWrmsNorm(b, ewt);
    N_VLinearSum(ONE, acor, ONE, b, acor);
    N_VLinearSum(ONE, zn[0], ONE, acor, y);
    
    /* Test for convergence.  If m > 0, an estimate of the convergence
       rate constant is stored in crate, and used in the test.        */
    if (m > 0) {
      crate = MAX(CRDOWN * crate, del/delp);
    }
    dcon = del * MIN(ONE, crate) / tq[4];
    
    if (dcon <= ONE) {
      acnrm = (m==0) ? del : N_VWrmsNorm(acor, ewt);
      jcur = FALSE;
      return(ARK_SUCCESS); /* Nonlinear system was solved successfully */
    }
    
    mnewt = ++m;
    
    /* Stop at maxcor iterations or if iter. seems to be diverging.
       If still not converged and Jacobian data is not current, 
       signal to try the solution again                            */
    if ((m == maxcor) || ((m >= 2) && (del > RDIV*delp))) {
      if ((!jcur) && (setupNonNull)) return(TRY_AGAIN);
      else                           return(CONV_FAIL);
    }
    
    /* Save norm of correction, evaluate f, and loop again */
    delp = del;
    retval = f(tn, y, ftemp, user_data);
    nfe++;
    if (retval < 0) return(ARK_RHSFUNC_FAIL);
    if (retval > 0) {
      if ((!jcur) && (setupNonNull)) return(TRY_AGAIN);
      else                           return(RHSFUNC_RECVR);
    }

  } /* end loop */
}

/*
 * ARKHandleFlag
 *
 * This routine takes action on the return value nflag = *nflagPtr
 * returned by ARKNls, as follows:
 *
 * If ARKNls succeeded in solving the nonlinear system, then
 * ARKHandleNFlag returns the constant DO_ERROR_TEST, which tells ARKStep
 * to perform the error test.
 *
 * If the nonlinear system was not solved successfully, then ncfn and
 * ncf = *ncfPtr are incremented and Nordsieck array zn is restored.
 *
 * If the solution of the nonlinear system failed due to an
 * unrecoverable failure by setup, we return the value ARK_LSETUP_FAIL.
 * 
 * If it failed due to an unrecoverable failure in solve, then we return
 * the value ARK_LSOLVE_FAIL.
 *
 * If it failed due to an unrecoverable failure in rhs, then we return
 * the value ARK_RHSFUNC_FAIL.
 *
 * Otherwise, a recoverable failure occurred when solving the 
 * nonlinear system (ARKNls returned nflag == CONV_FAIL or RHSFUNC_RECVR). 
 * In this case, if ncf is now equal to maxncf or |h| = hmin, 
 * we return the value ARK_CONV_FAILURE (if nflag=CONV_FAIL) or
 * ARK_REPTD_RHSFUNC_ERR (if nflag=RHSFUNC_RECVR).
 * If not, we set *nflagPtr = PREV_CONV_FAIL and return the value
 * PREDICT_AGAIN, telling ARKStep to reattempt the step.
 *
 */

static int ARKHandleNFlag(ARKodeMem ark_mem, int *nflagPtr, realtype saved_t,
                         int *ncfPtr)
{
  int nflag;
  
  nflag = *nflagPtr;
  
  if (nflag == ARK_SUCCESS) return(DO_ERROR_TEST);

  /* The nonlinear soln. failed; increment ncfn and restore zn */
  ncfn++;
  ARKRestore(ark_mem, saved_t);
  
  /* Return if lsetup, lsolve, or rhs failed unrecoverably */
  if (nflag == ARK_LSETUP_FAIL)  return(ARK_LSETUP_FAIL);
  if (nflag == ARK_LSOLVE_FAIL)  return(ARK_LSOLVE_FAIL);
  if (nflag == ARK_RHSFUNC_FAIL) return(ARK_RHSFUNC_FAIL);
  
  /* At this point, nflag = CONV_FAIL or RHSFUNC_RECVR; increment ncf */
  
  (*ncfPtr)++;
  etamax = ONE;

  /* If we had maxncf failures or |h| = hmin, 
     return ARK_CONV_FAILURE or ARK_REPTD_RHSFUNC_ERR. */

  if ((ABS(h) <= hmin*ONEPSM) || (*ncfPtr == maxncf)) {
    if (nflag == CONV_FAIL)     return(ARK_CONV_FAILURE);
    if (nflag == RHSFUNC_RECVR) return(ARK_REPTD_RHSFUNC_ERR);    
  }

  /* Reduce step size; return to reattempt the step */

  eta = MAX(ETACF, hmin / ABS(h));
  *nflagPtr = PREV_CONV_FAIL;
  ARKRescale(ark_mem);

  return(PREDICT_AGAIN);
}

/*
 * ARKRestore
 *
 * This routine restores the value of tn to saved_t and undoes the
 * prediction.  After execution of ARKRestore, the Nordsieck array zn has
 * the same values as before the call to ARKPredict.
 */

static void ARKRestore(ARKodeMem ark_mem, realtype saved_t)
{
  int j, k;
  
  tn = saved_t;
  for (k = 1; k <= q; k++)
    for (j = q; j >= k; j--)
      N_VLinearSum(ONE, zn[j-1], -ONE, zn[j], zn[j-1]);
}

/*
 * ARKDoErrorTest
 *
 * This routine performs the local error test. 
 * The weighted local error norm dsm is loaded into *dsmPtr, and 
 * the test dsm ?<= 1 is made.
 *
 * If the test passes, ARKDoErrorTest returns ARK_SUCCESS. 
 *
 * If the test fails, we undo the step just taken (call ARKRestore) and 
 *
 *   - if maxnef error test failures have occurred or if ABS(h) = hmin,
 *     we return ARK_ERR_FAILURE.
 *
 *   - if more than MXNEF1 error test failures have occurred, an order
 *     reduction is forced. If already at order 1, restart by reloading 
 *     zn from scratch. If f() fails we return either ARK_RHSFUNC_FAIL
 *     or ARK_UNREC_RHSFUNC_ERR (no recovery is possible at this stage).
 *
 *   - otherwise, set *nflagPtr to PREV_ERR_FAIL, and return TRY_AGAIN. 
 *
 */

static booleantype ARKDoErrorTest(ARKodeMem ark_mem, int *nflagPtr,
                                realtype saved_t, int *nefPtr, realtype *dsmPtr)
{
  realtype dsm;
  int retval;
  
  dsm = acnrm * tq[2];

  /* If est. local error norm dsm passes test, return ARK_SUCCESS */  
  *dsmPtr = dsm; 
  if (dsm <= ONE) return(ARK_SUCCESS);
  
  /* Test failed; increment counters, set nflag, and restore zn array */
  (*nefPtr)++;
  netf++;
  *nflagPtr = PREV_ERR_FAIL;
  ARKRestore(ark_mem, saved_t);

  /* At maxnef failures or |h| = hmin, return ARK_ERR_FAILURE */
  if ((ABS(h) <= hmin*ONEPSM) || (*nefPtr == maxnef)) return(ARK_ERR_FAILURE);

  /* Set etamax = 1 to prevent step size increase at end of this step */
  etamax = ONE;

  /* Set h ratio eta from dsm, rescale, and return for retry of step */
  if (*nefPtr <= MXNEF1) {
    eta = ONE / (RPowerR(BIAS2*dsm,ONE/L) + ADDON);
    eta = MAX(ETAMIN, MAX(eta, hmin / ABS(h)));
    if (*nefPtr >= SMALL_NEF) eta = MIN(eta, ETAMXF);
    ARKRescale(ark_mem);
    return(TRY_AGAIN);
  }
  
  /* After MXNEF1 failures, force an order reduction and retry step */
  if (q > 1) {
    eta = MAX(ETAMIN, hmin / ABS(h));
    ARKAdjustOrder(ark_mem,-1);
    L = q;
    q--;
    qwait = L;
    ARKRescale(ark_mem);
    return(TRY_AGAIN);
  }

  /* If already at order 1, restart: reload zn from scratch */

  eta = MAX(ETAMIN, hmin / ABS(h));
  h *= eta;
  next_h = h;
  hscale = h;
  qwait = LONG_WAIT;
  nscon = 0;

  retval = f(tn, zn[0], tempv, user_data);
  nfe++;
  if (retval < 0)  return(ARK_RHSFUNC_FAIL);
  if (retval > 0)  return(ARK_UNREC_RHSFUNC_ERR);

  N_VScale(h, tempv, zn[1]);

  return(TRY_AGAIN);
}

/* 
 * =================================================================
 *  Private Functions Implementation after succesful step
 * =================================================================
 */

/*
 * ARKCompleteStep
 *
 * This routine performs various update operations when the solution
 * to the nonlinear system has passed the local error test. 
 * We increment the step counter nst, record the values hu and qu,
 * update the tau array, and apply the corrections to the zn array.
 * The tau[i] are the last q values of h, with tau[1] the most recent.
 * The counter qwait is decremented, and if qwait == 1 (and q < qmax)
 * we save acor and tq[5] for a possible order increase.
 */

static void ARKCompleteStep(ARKodeMem ark_mem)
{
  int i, j;
  
  nst++;
  nscon++;
  hu = h;
  qu = q;

  for (i=q; i >= 2; i--)  tau[i] = tau[i-1];
  if ((q==1) && (nst > 1)) tau[2] = tau[1];
  tau[1] = h;

  for (j=0; j <= q; j++) 
    N_VLinearSum(l[j], acor, ONE, zn[j], zn[j]);
  qwait--;
  if ((qwait == 1) && (q != qmax)) {
    N_VScale(ONE, acor, zn[qmax]);
    saved_tq5 = tq[5];
    indx_acor = qmax;
  }
}

/*
 * ARKprepareNextStep
 *
 * This routine handles the setting of stepsize and order for the
 * next step -- hprime and qprime.  Along with hprime, it sets the
 * ratio eta = hprime/h.  It also updates other state variables 
 * related to a change of step size or order. 
 */

 static void ARKPrepareNextStep(ARKodeMem ark_mem, realtype dsm)
{
  /* If etamax = 1, defer step size or order changes */
  if (etamax == ONE) {
    qwait = MAX(qwait, 2);
    qprime = q;
    hprime = h;
    eta = ONE;
    return;
  }

  /* etaq is the ratio of new to old h at the current order */  
  etaq = ONE /(RPowerR(BIAS2*dsm,ONE/L) + ADDON);
  
  /* If no order change, adjust eta and acor in ARKSetEta and return */
  if (qwait != 0) {
    eta = etaq;
    qprime = q;
    ARKSetEta(ark_mem);
    return;
  }
  
  /* If qwait = 0, consider an order change.   etaqm1 and etaqp1 are 
     the ratios of new to old h at orders q-1 and q+1, respectively.
     ARKChooseEta selects the largest; ARKSetEta adjusts eta and acor */
  qwait = 2;
  etaqm1 = ARKComputeEtaqm1(ark_mem);
  etaqp1 = ARKComputeEtaqp1(ark_mem);  
  ARKChooseEta(ark_mem); 
  ARKSetEta(ark_mem);
}

/*
 * ARKsetEta
 *
 * This routine adjusts the value of eta according to the various
 * heuristic limits and the optional input hmax.  It also resets
 * etamax to be the estimated local error vector.
 */

static void ARKSetEta(ARKodeMem ark_mem)
{

  /* If eta below the threshhold THRESH, reject a change of step size */
  if (eta < THRESH) {
    eta = ONE;
    hprime = h;
  } else {
    /* Limit eta by etamax and hmax, then set hprime */
    eta = MIN(eta, etamax);
    eta /= MAX(ONE, ABS(h)*hmax_inv*eta);
    hprime = h * eta;
    if (qprime < q) nscon = 0;
  }
  
  /* Reset etamax for the next step size change, and scale acor */
}

/*
 * ARKComputeEtaqm1
 *
 * This routine computes and returns the value of etaqm1 for a
 * possible decrease in order by 1.
 */

static realtype ARKComputeEtaqm1(ARKodeMem ark_mem)
{
  realtype ddn;
  
  etaqm1 = ZERO;
  if (q > 1) {
    ddn = N_VWrmsNorm(zn[q], ewt) * tq[1];
    etaqm1 = ONE/(RPowerR(BIAS1*ddn, ONE/q) + ADDON);
  }
  return(etaqm1);
}

/*
 * ARKComputeEtaqp1
 *
 * This routine computes and returns the value of etaqp1 for a
 * possible increase in order by 1.
 */

static realtype ARKComputeEtaqp1(ARKodeMem ark_mem)
{
  realtype dup, cquot;
  
  etaqp1 = ZERO;
  if (q != qmax) {
    if (saved_tq5 == ZERO) return(etaqp1);
    cquot = (tq[5] / saved_tq5) * RPowerI(h/tau[2], L);
    N_VLinearSum(-cquot, zn[qmax], ONE, acor, tempv);
    dup = N_VWrmsNorm(tempv, ewt) * tq[3];
    etaqp1 = ONE / (RPowerR(BIAS3*dup, ONE/(L+1)) + ADDON);
  }
  return(etaqp1);
}

/*
 * ARKChooseEta
 * Given etaqm1, etaq, etaqp1 (the values of eta for qprime =
 * q - 1, q, or q + 1, respectively), this routine chooses the 
 * maximum eta value, sets eta to that value, and sets qprime to the
 * corresponding value of q.  If there is a tie, the preference
 * order is to (1) keep the same order, then (2) decrease the order,
 * and finally (3) increase the order.  If the maximum eta value
 * is below the threshhold THRESH, the order is kept unchanged and
 * eta is set to 1.
 */

static void ARKChooseEta(ARKodeMem ark_mem)
{
  realtype etam;
  
  etam = MAX(etaqm1, MAX(etaq, etaqp1));
  
  if (etam < THRESH) {
    eta = ONE;
    qprime = q;
    return;
  }

  if (etam == etaq) {

    eta = etaq;
    qprime = q;

  } else if (etam == etaqm1) {

    eta = etaqm1;
    qprime = q - 1;

  } else {

    eta = etaqp1;
    qprime = q + 1;

    if (lmm == ARK_BDF) {

      /* 
       * Store Delta_n in zn[qmax] to be used in order increase 
       *
       * This happens at the last step of order q before an increase
       * to order q+1, so it represents Delta_n in the ELTE at q+1
       */

      N_VScale(ONE, acor, zn[qmax]);

    }

  }

}

/*
 * ARKHandleFailure
 *
 * This routine prints error messages for all cases of failure by
 * ARKHin and ARKStep. It returns to ARKode the value that ARKode is 
 * to return to the user.
 */

static int ARKHandleFailure(ARKodeMem ark_mem, int flag)
{

  /* Set vector of  absolute weighted local errors */
  /*
  N_VProd(acor, ewt, tempv);
  N_VAbs(tempv, tempv);
  */

  /* Depending on flag, print error message and return error flag */
  switch (flag) {
  case ARK_ERR_FAILURE: 
    ARKProcessError(ark_mem, ARK_ERR_FAILURE, "ARKODE", "ARKode", MSGARK_ERR_FAILS, tn, h);
    break;
  case ARK_CONV_FAILURE:
    ARKProcessError(ark_mem, ARK_CONV_FAILURE, "ARKODE", "ARKode", MSGARK_CONV_FAILS, tn, h);
    break;
  case ARK_LSETUP_FAIL:
    ARKProcessError(ark_mem, ARK_LSETUP_FAIL, "ARKODE", "ARKode", MSGARK_SETUP_FAILED, tn);
    break;
  case ARK_LSOLVE_FAIL:
    ARKProcessError(ark_mem, ARK_LSOLVE_FAIL, "ARKODE", "ARKode", MSGARK_SOLVE_FAILED, tn);
    break;
  case ARK_RHSFUNC_FAIL:
    ARKProcessError(ark_mem, ARK_RHSFUNC_FAIL, "ARKODE", "ARKode", MSGARK_RHSFUNC_FAILED, tn);
    break;
  case ARK_UNREC_RHSFUNC_ERR:
    ARKProcessError(ark_mem, ARK_UNREC_RHSFUNC_ERR, "ARKODE", "ARKode", MSGARK_RHSFUNC_UNREC, tn);
    break;
  case ARK_REPTD_RHSFUNC_ERR:
    ARKProcessError(ark_mem, ARK_REPTD_RHSFUNC_ERR, "ARKODE", "ARKode", MSGARK_RHSFUNC_REPTD, tn);
    break;
  case ARK_RTFUNC_FAIL:    
    ARKProcessError(ark_mem, ARK_RTFUNC_FAIL, "ARKODE", "ARKode", MSGARK_RTFUNC_FAILED, tn);
    break;
  case ARK_TOO_CLOSE:
    ARKProcessError(ark_mem, ARK_TOO_CLOSE, "ARKODE", "ARKode", MSGARK_TOO_CLOSE);
    break;
  default:
    return(ARK_SUCCESS);   
  }

  return(flag);

}

/* 
 * =================================================================
 * BDF Stability Limit Detection                       
 * =================================================================
 */

/*
 * ARKBDFStab
 *
 * This routine handles the BDF Stability Limit Detection Algorithm
 * STALD.  It is called if lmm = ARK_BDF and the SLDET option is on.
 * If the order is 3 or more, the required norm data is saved.
 * If a decision to reduce order has not already been made, and
 * enough data has been saved, ARKsldet is called.  If it signals
 * a stability limit violation, the order is reduced, and the step
 * size is reset accordingly.
 */

void ARKBDFStab(ARKodeMem ark_mem)
{
  int i,k, ldflag, factorial;
  realtype sq, sqm1, sqm2;
      
  /* If order is 3 or greater, then save scaled derivative data,
     push old data down in i, then add current values to top.    */

  if (q >= 3) {
    for (k = 1; k <= 3; k++)
      { for (i = 5; i >= 2; i--) ssdat[i][k] = ssdat[i-1][k]; }
    factorial = 1;
    for (i = 1; i <= q-1; i++) factorial *= i;
    sq = factorial*q*(q+1)*acnrm/MAX(tq[5],TINY);
    sqm1 = factorial*q*N_VWrmsNorm(zn[q], ewt);
    sqm2 = factorial*N_VWrmsNorm(zn[q-1], ewt);
    ssdat[1][1] = sqm2*sqm2;
    ssdat[1][2] = sqm1*sqm1;
    ssdat[1][3] = sq*sq;
  }  

  if (qprime >= q) {

    /* If order is 3 or greater, and enough ssdat has been saved,
       nscon >= q+5, then call stability limit detection routine.  */

    if ( (q >= 3) && (nscon >= q+5) ) {
      ldflag = ARKsldet(ark_mem);
      if (ldflag > 3) {
        /* A stability limit violation is indicated by
           a return flag of 4, 5, or 6.
           Reduce new order.                     */
        qprime = q-1;
        eta = etaqm1; 
        eta = MIN(eta,etamax);
        eta = eta/MAX(ONE,ABS(h)*hmax_inv*eta);
        hprime = h*eta;
        nor = nor + 1;
      }
    }
  }
  else {
    /* Otherwise, let order increase happen, and 
       reset stability limit counter, nscon.     */
    nscon = 0;
  }
}

/*
 * ARKsldet
 *
 * This routine detects stability limitation using stored scaled 
 * derivatives data. ARKsldet returns the magnitude of the
 * dominate characteristic root, rr. The presents of a stability
 * limit is indicated by rr > "something a little less then 1.0",  
 * and a positive kflag. This routine should only be called if
 * order is greater than or equal to 3, and data has been collected
 * for 5 time steps. 
 * 
 * Returned values:
 *    kflag = 1 -> Found stable characteristic root, normal matrix case
 *    kflag = 2 -> Found stable characteristic root, quartic solution
 *    kflag = 3 -> Found stable characteristic root, quartic solution,
 *                 with Newton correction
 *    kflag = 4 -> Found stability violation, normal matrix case
 *    kflag = 5 -> Found stability violation, quartic solution
 *    kflag = 6 -> Found stability violation, quartic solution,
 *                 with Newton correction
 *
 *    kflag < 0 -> No stability limitation, 
 *                 or could not compute limitation.
 *
 *    kflag = -1 -> Min/max ratio of ssdat too small.
 *    kflag = -2 -> For normal matrix case, vmax > vrrt2*vrrt2
 *    kflag = -3 -> For normal matrix case, The three ratios
 *                  are inconsistent.
 *    kflag = -4 -> Small coefficient prevents elimination of quartics.  
 *    kflag = -5 -> R value from quartics not consistent.
 *    kflag = -6 -> No corrected root passes test on qk values
 *    kflag = -7 -> Trouble solving for sigsq.
 *    kflag = -8 -> Trouble solving for B, or R via B.
 *    kflag = -9 -> R via sigsq[k] disagrees with R from data.
 */

static int ARKsldet(ARKodeMem ark_mem)
{
  int i, k, j, it, kmin, kflag = 0;
  realtype rat[5][4], rav[4], qkr[4], sigsq[4], smax[4], ssmax[4];
  realtype drr[4], rrc[4],sqmx[4], qjk[4][4], vrat[5], qc[6][4], qco[6][4];
  realtype rr, rrcut, vrrtol, vrrt2, sqtol, rrtol;
  realtype smink, smaxk, sumrat, sumrsq, vmin, vmax, drrmax, adrr;
  realtype tem, sqmax, saqk, qp, s, sqmaxk, saqj, sqmin;
  realtype rsa, rsb, rsc, rsd, rd1a, rd1b, rd1c;
  realtype rd2a, rd2b, rd3a, cest1, corr1; 
  realtype ratp, ratm, qfac1, qfac2, bb, rrb;

  /* The following are cutoffs and tolerances used by this routine */

  rrcut  = RCONST(0.98);
  vrrtol = RCONST(1.0e-4);
  vrrt2  = RCONST(5.0e-4);
  sqtol  = RCONST(1.0e-3);
  rrtol  = RCONST(1.0e-2);
  
  rr = ZERO;
  
  /*  Index k corresponds to the degree of the interpolating polynomial. */
  /*      k = 1 -> q-1          */
  /*      k = 2 -> q            */
  /*      k = 3 -> q+1          */
  
  /*  Index i is a backward-in-time index, i = 1 -> current time, */
  /*      i = 2 -> previous step, etc    */
  
  /* get maxima, minima, and variances, and form quartic coefficients  */
  
  for (k=1; k<=3; k++) {
    smink = ssdat[1][k];
    smaxk = ZERO;
    
    for (i=1; i<=5; i++) {
      smink = MIN(smink,ssdat[i][k]);
      smaxk = MAX(smaxk,ssdat[i][k]);
    }
    
    if (smink < TINY*smaxk) {
      kflag = -1;  
      return(kflag);
    }
    smax[k] = smaxk;
    ssmax[k] = smaxk*smaxk;
    
    sumrat = ZERO;
    sumrsq = ZERO;
    for (i=1; i<=4; i++) {
      rat[i][k] = ssdat[i][k]/ssdat[i+1][k];
      sumrat = sumrat + rat[i][k];
      sumrsq = sumrsq + rat[i][k]*rat[i][k];
    } 
    rav[k] = FOURTH*sumrat;
    vrat[k] = ABS(FOURTH*sumrsq - rav[k]*rav[k]);
    
    qc[5][k] = ssdat[1][k]*ssdat[3][k] - ssdat[2][k]*ssdat[2][k];
    qc[4][k] = ssdat[2][k]*ssdat[3][k] - ssdat[1][k]*ssdat[4][k];
    qc[3][k] = ZERO;
    qc[2][k] = ssdat[2][k]*ssdat[5][k] - ssdat[3][k]*ssdat[4][k];
    qc[1][k] = ssdat[4][k]*ssdat[4][k] - ssdat[3][k]*ssdat[5][k];
    
    for (i=1; i<=5; i++) {
      qco[i][k] = qc[i][k];
    }
  }                            /* End of k loop */
  
  /* Isolate normal or nearly-normal matrix case. Three quartic will
     have common or nearly-common roots in this case. 
     Return a kflag = 1 if this procedure works. If three root 
     differ more than vrrt2, return error kflag = -3.    */
  
  vmin = MIN(vrat[1],MIN(vrat[2],vrat[3]));
  vmax = MAX(vrat[1],MAX(vrat[2],vrat[3]));
  
  if(vmin < vrrtol*vrrtol) {
    if (vmax > vrrt2*vrrt2) {
      kflag = -2;  
      return(kflag);
    } else {
      rr = (rav[1] + rav[2] + rav[3])/THREE;
      
      drrmax = ZERO;
      for(k = 1;k<=3;k++) {
        adrr = ABS(rav[k] - rr);
        drrmax = MAX(drrmax, adrr);
      }
      if (drrmax > vrrt2) {
        kflag = -3;    
      }
      
      kflag = 1;

      /*  can compute charactistic root, drop to next section   */
      
    }
  } else {

    /* use the quartics to get rr. */
    
    if (ABS(qco[1][1]) < TINY*ssmax[1]) {
      kflag = -4;    
      return(kflag);
    }
    
    tem = qco[1][2]/qco[1][1];
    for(i=2; i<=5; i++) {
      qco[i][2] = qco[i][2] - tem*qco[i][1];
    }

    qco[1][2] = ZERO;
    tem = qco[1][3]/qco[1][1];
    for(i=2; i<=5; i++) {
      qco[i][3] = qco[i][3] - tem*qco[i][1];
    }
    qco[1][3] = ZERO;
    
    if (ABS(qco[2][2]) < TINY*ssmax[2]) {
      kflag = -4;    
      return(kflag);
    }
    
    tem = qco[2][3]/qco[2][2];
    for(i=3; i<=5; i++) {
      qco[i][3] = qco[i][3] - tem*qco[i][2];
    }
    
    if (ABS(qco[4][3]) < TINY*ssmax[3]) {
      kflag = -4;    
      return(kflag);
    }
    
    rr = -qco[5][3]/qco[4][3];
    
    if (rr < TINY || rr > HUN) {
      kflag = -5;   
      return(kflag);
    }
    
    for(k=1; k<=3; k++) {
      qkr[k] = qc[5][k] + rr*(qc[4][k] + rr*rr*(qc[2][k] + rr*qc[1][k]));
    }  
    
    sqmax = ZERO;
    for(k=1; k<=3; k++) {
      saqk = ABS(qkr[k])/ssmax[k];
      if (saqk > sqmax) sqmax = saqk;
    } 
    
    if (sqmax < sqtol) {
      kflag = 2;
      
      /*  can compute charactistic root, drop to "given rr,etc"   */
      
    } else {

      /* do Newton corrections to improve rr.  */
      
      for(it=1; it<=3; it++) {
        for(k=1; k<=3; k++) {
          qp = qc[4][k] + rr*rr*(THREE*qc[2][k] + rr*FOUR*qc[1][k]);
          drr[k] = ZERO;
          if (ABS(qp) > TINY*ssmax[k]) drr[k] = -qkr[k]/qp;
          rrc[k] = rr + drr[k];
        } 
        
        for(k=1; k<=3; k++) {
          s = rrc[k];
          sqmaxk = ZERO;
          for(j=1; j<=3; j++) {
            qjk[j][k] = qc[5][j] + s*(qc[4][j] + 
                                      s*s*(qc[2][j] + s*qc[1][j]));
            saqj = ABS(qjk[j][k])/ssmax[j];
            if (saqj > sqmaxk) sqmaxk = saqj;
          } 
          sqmx[k] = sqmaxk;
        } 

        sqmin = sqmx[1]; kmin = 1;
        for(k=2; k<=3; k++) {
          if (sqmx[k] < sqmin) {
            kmin = k;
            sqmin = sqmx[k];
          }
        } 
        rr = rrc[kmin];
        
        if (sqmin < sqtol) {
          kflag = 3;
          /*  can compute charactistic root   */
          /*  break out of Newton correction loop and drop to "given rr,etc" */ 
          break;
        } else {
          for(j=1; j<=3; j++) {
            qkr[j] = qjk[j][kmin];
          }
        }     
      }          /*  end of Newton correction loop  */ 
      
      if (sqmin > sqtol) {
        kflag = -6;
        return(kflag);
      }
    }     /*  end of if (sqmax < sqtol) else   */
  }      /*  end of if(vmin < vrrtol*vrrtol) else, quartics to get rr. */
  
  /* given rr, find sigsq[k] and verify rr.  */
  /* All positive kflag drop to this section  */
  
  for(k=1; k<=3; k++) {
    rsa = ssdat[1][k];
    rsb = ssdat[2][k]*rr;
    rsc = ssdat[3][k]*rr*rr;
    rsd = ssdat[4][k]*rr*rr*rr;
    rd1a = rsa - rsb;
    rd1b = rsb - rsc;
    rd1c = rsc - rsd;
    rd2a = rd1a - rd1b;
    rd2b = rd1b - rd1c;
    rd3a = rd2a - rd2b;
    
    if (ABS(rd1b) < TINY*smax[k]) {
      kflag = -7;
      return(kflag);
    }
    
    cest1 = -rd3a/rd1b;
    if (cest1 < TINY || cest1 > FOUR) {
      kflag = -7;
      return(kflag);
    }
    corr1 = (rd2b/cest1)/(rr*rr);
    sigsq[k] = ssdat[3][k] + corr1;
  }
  
  if (sigsq[2] < TINY) {
    kflag = -8;
    return(kflag);
  }
  
  ratp = sigsq[3]/sigsq[2];
  ratm = sigsq[1]/sigsq[2];
  qfac1 = FOURTH*(q*q - ONE);
  qfac2 = TWO/(q - ONE);
  bb = ratp*ratm - ONE - qfac1*ratp;
  tem = ONE - qfac2*bb;
  
  if (ABS(tem) < TINY) {
    kflag = -8;
    return(kflag);
  }
  
  rrb = ONE/tem;
  
  if (ABS(rrb - rr) > rrtol) {
    kflag = -9;
    return(kflag);
  }
  
  /* Check to see if rr is above cutoff rrcut  */
  if (rr > rrcut) {
    if (kflag == 1) kflag = 4;
    if (kflag == 2) kflag = 5;
    if (kflag == 3) kflag = 6;
  }
  
  /* All positive kflag returned at this point  */
  
  return(kflag);
  
}

/* 
 * =================================================================
 * Root finding   
 * =================================================================
 */

/*-----------------------------------------------------------------*/

/* 
 * ARKRcheck1
 *
 * This routine completes the initialization of rootfinding memory
 * information, and checks whether g has a zero both at and very near
 * the initial point of the IVP.
 *
 * This routine returns an int equal to:
 *  ARK_RTFUNC_FAIL = -12  if the g function failed, or
 *  ARK_SUCCESS     =   0  otherwise.
 */

static int ARKRcheck1(ARKodeMem ark_mem)
{
  int i, retval;
  realtype smallh, hratio, tplus;
  booleantype zroot;

  for (i = 0; i < nrtfn; i++) iroots[i] = 0;
  tlo = tn;
  ttol = (ABS(tn) + ABS(h))*uround*HUN;

  /* Evaluate g at initial t and check for zero values. */
  retval = gfun(tlo, zn[0], glo, user_data);
  nge = 1;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) {
    if (ABS(glo[i]) == ZERO) {
      zroot = TRUE;
      gactive[i] = FALSE;
    }
  }
  if (!zroot) return(ARK_SUCCESS);

  /* Some g_i is zero at t0; look at g at t0+(small increment). */
  hratio = MAX(ttol/ABS(h), TENTH);
  smallh = hratio*h;
  tplus = tlo + smallh;
  N_VLinearSum(ONE, zn[0], hratio, zn[1], y);
  retval = gfun(tplus, y, ghi, user_data);
  nge++;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  /* We check now only the components of g which were exactly 0.0 at t0
   * to see if we can 'activate' them. */
  for (i = 0; i < nrtfn; i++) {
    if (!gactive[i] && ABS(ghi[i]) != ZERO) {
      gactive[i] = TRUE;
      glo[i] = ghi[i];
    }
  }
  return(ARK_SUCCESS);
}

/*
 * ARKRcheck2
 *
 * This routine checks for exact zeros of g at the last root found,
 * if the last return was a root.  It then checks for a close pair of
 * zeros (an error condition), and for a new root at a nearby point.
 * The array glo = g(tlo) at the left endpoint of the search interval
 * is adjusted if necessary to assure that all g_i are nonzero
 * there, before returning to do a root search in the interval.
 *
 * On entry, tlo = tretlast is the last value of tret returned by
 * ARKode.  This may be the previous tn, the previous tout value, or
 * the last root location.
 *
 * This routine returns an int equal to:
 *      ARK_RTFUNC_FAIL = -12 if the g function failed, or
 *      CLOSERT        =  3  if a close pair of zeros was found, or
 *      RTFOUND        =  1  if a new zero of g was found near tlo, or
 *      ARK_SUCCESS     =  0  otherwise.
 */

static int ARKRcheck2(ARKodeMem ark_mem)
{
  int i, retval;
  realtype smallh, hratio, tplus;
  booleantype zroot;

  if (irfnd == 0) return(ARK_SUCCESS);

  (void) ARKodeGetDky(ark_mem, tlo, 0, y);
  retval = gfun(tlo, y, glo, user_data);
  nge++;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) iroots[i] = 0;
  for (i = 0; i < nrtfn; i++) {
    if (!gactive[i]) continue;
    if (ABS(glo[i]) == ZERO) {
      zroot = TRUE;
      iroots[i] = 1;
    }
  }
  if (!zroot) return(ARK_SUCCESS);

  /* One or more g_i has a zero at tlo.  Check g at tlo+smallh. */
  ttol = (ABS(tn) + ABS(h))*uround*HUN;
  smallh = (h > ZERO) ? ttol : -ttol;
  tplus = tlo + smallh;
  if ( (tplus - tn)*h >= ZERO) {
    hratio = smallh/h;
    N_VLinearSum(ONE, y, hratio, zn[1], y);
  } else {
    (void) ARKodeGetDky(ark_mem, tplus, 0, y);
  }
  retval = gfun(tplus, y, ghi, user_data);
  nge++;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  /* Check for close roots (error return), for a new zero at tlo+smallh,
  and for a g_i that changed from zero to nonzero. */
  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) {
    if (ABS(ghi[i]) == ZERO) {
      if (!gactive[i]) continue;
      if (iroots[i] == 1) return(CLOSERT);
      zroot = TRUE;
      iroots[i] = 1;
    } else {
      if (iroots[i] == 1) glo[i] = ghi[i];
    }
  }
  if (zroot) return(RTFOUND);
  return(ARK_SUCCESS);
}

/*
 * ARKRcheck3
 *
 * This routine interfaces to ARKRootfind to look for a root of g
 * between tlo and either tn or tout, whichever comes first.
 * Only roots beyond tlo in the direction of integration are sought.
 *
 * This routine returns an int equal to:
 *      ARK_RTFUNC_FAIL = -12 if the g function failed, or
 *      RTFOUND        =  1  if a root of g was found, or
 *      ARK_SUCCESS     =  0  otherwise.
 */

static int ARKRcheck3(ARKodeMem ark_mem)
{
  int i, retval, ier;

  /* Set thi = tn or tout, whichever comes first; set y = y(thi). */
  if (taskc == ARK_ONE_STEP) {
    thi = tn;
    N_VScale(ONE, zn[0], y);
  }
  if (taskc == ARK_NORMAL) {
    if ( (toutc - tn)*h >= ZERO) {
      thi = tn; 
      N_VScale(ONE, zn[0], y);
    } else {
      thi = toutc;
      (void) ARKodeGetDky(ark_mem, thi, 0, y);
    }
  }

  /* Set ghi = g(thi) and call ARKRootfind to search (tlo,thi) for roots. */
  retval = gfun(thi, y, ghi, user_data);
  nge++;
  if (retval != 0) return(ARK_RTFUNC_FAIL);

  ttol = (ABS(tn) + ABS(h))*uround*HUN;
  ier = ARKRootfind(ark_mem);
  if (ier == ARK_RTFUNC_FAIL) return(ARK_RTFUNC_FAIL);
  for(i=0; i<nrtfn; i++) {
    if(!gactive[i] && grout[i] != ZERO) gactive[i] = TRUE;
  }
  tlo = trout;
  for (i = 0; i < nrtfn; i++) glo[i] = grout[i];

  /* If no root found, return ARK_SUCCESS. */  
  if (ier == ARK_SUCCESS) return(ARK_SUCCESS);

  /* If a root was found, interpolate to get y(trout) and return.  */
  (void) ARKodeGetDky(ark_mem, trout, 0, y);
  return(RTFOUND);

}

/*
 * ARKRootfind
 *
 * This routine solves for a root of g(t) between tlo and thi, if
 * one exists.  Only roots of odd multiplicity (i.e. with a change
 * of sign in one of the g_i), or exact zeros, are found.
 * Here the sign of tlo - thi is arbitrary, but if multiple roots
 * are found, the one closest to tlo is returned.
 *
 * The method used is the Illinois algorithm, a modified secant method.
 * Reference: Kathie L. Hiebert and Lawrence F. Shampine, Implicitly
 * Defined Output Points for Solutions of ODEs, Sandia National
 * Laboratory Report SAND80-0180, February 1980.
 *
 * This routine uses the following parameters for communication:
 *
 * nrtfn    = number of functions g_i, or number of components of
 *            the vector-valued function g(t).  Input only.
 *
 * gfun     = user-defined function for g(t).  Its form is
 *            (void) gfun(t, y, gt, user_data)
 *
 * rootdir  = in array specifying the direction of zero-crossings.
 *            If rootdir[i] > 0, search for roots of g_i only if
 *            g_i is increasing; if rootdir[i] < 0, search for
 *            roots of g_i only if g_i is decreasing; otherwise
 *            always search for roots of g_i.
 *
 * gactive  = array specifying whether a component of g should
 *            or should not be monitored. gactive[i] is initially
 *            set to TRUE for all i=0,...,nrtfn-1, but it may be
 *            reset to FALSE if at the first step g[i] is 0.0
 *            both at the I.C. and at a small perturbation of them.
 *            gactive[i] is then set back on TRUE only after the 
 *            corresponding g function moves away from 0.0.
 *
 * nge      = cumulative counter for gfun calls.
 *
 * ttol     = a convergence tolerance for trout.  Input only.
 *            When a root at trout is found, it is located only to
 *            within a tolerance of ttol.  Typically, ttol should
 *            be set to a value on the order of
 *               100 * UROUND * max (ABS(tlo), ABS(thi))
 *            where UROUND is the unit roundoff of the machine.
 *
 * tlo, thi = endpoints of the interval in which roots are sought.
 *            On input, and must be distinct, but tlo - thi may
 *            be of either sign.  The direction of integration is
 *            assumed to be from tlo to thi.  On return, tlo and thi
 *            are the endpoints of the final relevant interval.
 *
 * glo, ghi = arrays of length nrtfn containing the vectors g(tlo)
 *            and g(thi) respectively.  Input and output.  On input,
 *            none of the glo[i] should be zero.
 *
 * trout    = root location, if a root was found, or thi if not.
 *            Output only.  If a root was found other than an exact
 *            zero of g, trout is the endpoint thi of the final
 *            interval bracketing the root, with size at most ttol.
 *
 * grout    = array of length nrtfn containing g(trout) on return.
 *
 * iroots   = int array of length nrtfn with root information.
 *            Output only.  If a root was found, iroots indicates
 *            which components g_i have a root at trout.  For
 *            i = 0, ..., nrtfn-1, iroots[i] = 1 if g_i has a root
 *            and g_i is increasing, iroots[i] = -1 if g_i has a
 *            root and g_i is decreasing, and iroots[i] = 0 if g_i
 *            has no roots or g_i varies in the direction opposite
 *            to that indicated by rootdir[i].
 *
 * This routine returns an int equal to:
 *      ARK_RTFUNC_FAIL = -12 if the g function failed, or
 *      RTFOUND        =  1  if a root of g was found, or
 *      ARK_SUCCESS     =  0  otherwise.
 */

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
  for (i = 0;  i < nrtfn; i++) {
    if(!gactive[i]) continue;
    if (ABS(ghi[i]) == ZERO) {
      if(rootdir[i]*glo[i] <= ZERO) {
        zroot = TRUE;
      }
    } else {
      if ( (glo[i]*ghi[i] < ZERO) && (rootdir[i]*glo[i] <= ZERO) ) {
        gfrac = ABS(ghi[i]/(ghi[i] - glo[i]));
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
    trout = thi;
    for (i = 0; i < nrtfn; i++) grout[i] = ghi[i];
    if (!zroot) return(ARK_SUCCESS);
    for (i = 0; i < nrtfn; i++) {
      iroots[i] = 0;
      if(!gactive[i]) continue;
      if (ABS(ghi[i]) == ZERO) iroots[i] = glo[i] > 0 ? -1:1;
    }
    return(RTFOUND);
  }

  /* Initialize alpha to avoid compiler warning */
  alpha = ONE;

  /* A sign change was found.  Loop to locate nearest root. */

  side = 0;  sideprev = -1;
  loop {                                    /* Looping point */

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
    tmid = thi - (thi - tlo)*ghi[imax]/(ghi[imax] - alpha*glo[imax]);
    if (ABS(tmid - tlo) < HALF*ttol) {
      fracint = ABS(thi - tlo)/ttol;
      fracsub = (fracint > FIVE) ? TENTH : HALF/fracint;
      tmid = tlo + fracsub*(thi - tlo);
    }
    if (ABS(thi - tmid) < HALF*ttol) {
      fracint = ABS(thi - tlo)/ttol;
      fracsub = (fracint > FIVE) ? TENTH : HALF/fracint;
      tmid = thi - fracsub*(thi - tlo);
    }

    (void) ARKodeGetDky(ark_mem, tmid, 0, y);
    retval = gfun(tmid, y, grout, user_data);
    nge++;
    if (retval != 0) return(ARK_RTFUNC_FAIL);

    /* Check to see in which subinterval g changes sign, and reset imax.
       Set side = 1 if sign change is on low side, or 2 if on high side.  */  
    maxfrac = ZERO;
    zroot = FALSE;
    sgnchg = FALSE;
    sideprev = side;
    for (i = 0;  i < nrtfn; i++) {
      if(!gactive[i]) continue;
      if (ABS(grout[i]) == ZERO) {
        if(rootdir[i]*glo[i] <= ZERO) {
          zroot = TRUE;
        }
      } else {
        if ( (glo[i]*grout[i] < ZERO) && (rootdir[i]*glo[i] <= ZERO) ) {
          gfrac = ABS(grout[i]/(grout[i] - glo[i]));
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
      thi = tmid;
      for (i = 0; i < nrtfn; i++) ghi[i] = grout[i];
      side = 1;
      /* Stop at root thi if converged; otherwise loop. */
      if (ABS(thi - tlo) <= ttol) break;
      continue;  /* Return to looping point. */
    }

    if (zroot) {
      /* No sign change in (tlo,tmid), but g = 0 at tmid; return root tmid. */
      thi = tmid;
      for (i = 0; i < nrtfn; i++) ghi[i] = grout[i];
      break;
    }

    /* No sign change in (tlo,tmid), and no zero at tmid.
       Sign change must be in (tmid,thi).  Replace tlo with tmid. */
    tlo = tmid;
    for (i = 0; i < nrtfn; i++) glo[i] = grout[i];
    side = 2;
    /* Stop at root thi if converged; otherwise loop back. */
    if (ABS(thi - tlo) <= ttol) break;

  } /* End of root-search loop */

  /* Reset trout and grout, set iroots, and return RTFOUND. */
  trout = thi;
  for (i = 0; i < nrtfn; i++) {
    grout[i] = ghi[i];
    iroots[i] = 0;
    if(!gactive[i]) continue;
    if ( (ABS(ghi[i]) == ZERO) && (rootdir[i]*glo[i] <= ZERO) ) 
      iroots[i] = glo[i] > 0 ? -1:1;
    if ( (glo[i]*ghi[i] < ZERO) && (rootdir[i]*glo[i] <= ZERO) ) 
      iroots[i] = glo[i] > 0 ? -1:1;
  }
  return(RTFOUND);
}

/* 
 * =================================================================
 * Internal EWT function
 * =================================================================
 */

/*
 * ARKEwtSet
 *
 * This routine is responsible for setting the error weight vector ewt,
 * according to tol_type, as follows:
 *
 * (1) ewt[i] = 1 / (reltol * ABS(ycur[i]) + *abstol), i=0,...,neq-1
 *     if tol_type = ARK_SS
 * (2) ewt[i] = 1 / (reltol * ABS(ycur[i]) + abstol[i]), i=0,...,neq-1
 *     if tol_type = ARK_SV
 *
 * ARKEwtSet returns 0 if ewt is successfully set as above to a
 * positive vector and -1 otherwise. In the latter case, ewt is
 * considered undefined.
 *
 * All the real work is done in the routines ARKEwtSetSS, ARKEwtSetSV.
 */

int ARKEwtSet(N_Vector ycur, N_Vector weight, void *data)
{
  ARKodeMem ark_mem;
  int flag = 0;

  /* data points to ark_mem here */

  ark_mem = (ARKodeMem) data;

  switch(itol) {
  case ARK_SS: 
    flag = ARKEwtSetSS(ark_mem, ycur, weight);
    break;
  case ARK_SV: 
    flag = ARKEwtSetSV(ark_mem, ycur, weight);
    break;
  }
  
  return(flag);
}

/*
 * ARKEwtSetSS
 *
 * This routine sets ewt as decribed above in the case tol_type = ARK_SS.
 * It tests for non-positive components before inverting. ARKEwtSetSS
 * returns 0 if ewt is successfully set to a positive vector
 * and -1 otherwise. In the latter case, ewt is considered undefined.
 */

static int ARKEwtSetSS(ARKodeMem ark_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, tempv);
  N_VScale(reltol, tempv, tempv);
  N_VAddConst(tempv, Sabstol, tempv);
  if (N_VMin(tempv) <= ZERO) return(-1);
  N_VInv(tempv, weight);
  return(0);
}

/*
 * ARKEwtSetSV
 *
 * This routine sets ewt as decribed above in the case tol_type = ARK_SV.
 * It tests for non-positive components before inverting. ARKEwtSetSV
 * returns 0 if ewt is successfully set to a positive vector
 * and -1 otherwise. In the latter case, ewt is considered undefined.
 */

static int ARKEwtSetSV(ARKodeMem ark_mem, N_Vector ycur, N_Vector weight)
{
  N_VAbs(ycur, tempv);
  N_VLinearSum(reltol, tempv, ONE, Vabstol, tempv);
  if (N_VMin(tempv) <= ZERO) return(-1);
  N_VInv(tempv, weight);
  return(0);
}

/* 
 * =================================================================
 * ARKODE Error Handling function   
 * =================================================================
 */

/* 
 * ARKProcessError is a high level error handling function
 * - if ark_mem==NULL it prints the error message to stderr
 * - otherwise, it sets-up and calls the error hadling function 
 *   pointed to by ark_ehfun
 */

#define ehfun    (ark_mem->ark_ehfun)
#define eh_data  (ark_mem->ark_eh_data)

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

    ehfun(error_code, module, fname, msg, eh_data);

  }

  /* Finalize argument processing */
  
  va_end(ap);

  return;

}

/* ARKErrHandler is the default error handling function.
   It sends the error message to the stream pointed to by ark_errfp */

#define errfp    (ark_mem->ark_errfp)

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
  if (errfp!=NULL) {
    fprintf(errfp,"\n[%s %s]  %s\n",module,err_type,function);
    fprintf(errfp,"  %s\n\n",msg);
  }
#endif

  return;
}
