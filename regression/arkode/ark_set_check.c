/*---------------------------------------------------------------
 $Revision: $
 $Date: $
-----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Set routine checking routine:
 
 Runs through all possible "set" routines in ARKODE to see 
 whether the relevant parameters get set appropriately within 
 ark_mem.
---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <math.h>
#include <arkode/arkode.h>
#include <arkode/arkode_impl.h>
#include <nvector/nvector_serial.h>
#include <arkode/arkode_dense.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

/* User-supplied Functions Called by the Solver */
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static void efun(int error_code, const char *module, const char *function, 
		 char *msg, void *user_data);
static int hfun(N_Vector y, realtype t, realtype h1, realtype h2, 
		realtype h3, realtype e1, realtype e2, realtype e3, 
		int q, int p, realtype *hnew, void *user_data);
static int sfun(N_Vector y, realtype t, realtype *hstab, void *user_data);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);


/* Main Program */
int main() {

  /* Set up basic ARKode solver structure */
  int flag;
  realtype T0 = RCONST(0.0);
  void *arkode_mem = NULL;
  long int NEQ = 1;
  N_Vector y = NULL;
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;
  N_VConst(0.0, y);
  arkode_mem = ARKodeCreate();
  if (check_flag((void *) arkode_mem, "ARKodeCreate", 0)) return 1;
  flag = ARKodeInit(arkode_mem, fe, fi, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return 1;

  /* Set ARKodeMem pointer to solver structure (for querying) */
  ARKodeMem ark_mem;
  ark_mem = (ARKodeMem) arkode_mem;

  /* Set default ARKode parameters */
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;
  

  /***** set routine checks *****/
  /* error handler */
  flag = ARKodeSetErrHandlerFn(arkode_mem, efun, NULL);
  if (check_flag(&flag, "ARKodeSetErrHandlerFn", 1)) return 1;
  if (ark_mem->ark_ehfun != efun) {
    printf("Error in ARKodeSetErrHandlerFn: did not set function\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* error file */
  FILE *EFID = NULL;
  flag = ARKodeSetErrFile(arkode_mem, EFID);
  if (check_flag(&flag, "ARKodeSetErrFile", 1)) return 1;
  if (ark_mem->ark_errfp != NULL) {
    printf("Error in ARKodeSetErrFile: did not set file pointer\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* user data */
  double data;
  flag = ARKodeSetUserData(arkode_mem, &data);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;
  if (ark_mem->ark_user_data != &data) {
    printf("Error in ARKodeSetUserData: did not set memory reference\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* solver diagnostics */
  FILE *DFID = stderr;
  flag = ARKodeSetDiagnostics(arkode_mem, DFID);
  if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return 1;
  if (ark_mem->ark_diagfp != DFID) {
    printf("Error in ARKodeSetDiagnostics: did not set file pointer\n");
    return 1;
  }
  if (ark_mem->ark_report != TRUE) {
    printf("Error in ARKodeSetDiagnostics: did not enable reporting\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* solver order */
  flag = ARKodeSetOrder(arkode_mem, 1);
  if (check_flag(&flag, "ARKodeSetOrder", 1)) return 1;
  if (ark_mem->ark_q != 1) {
    printf("Error in ARKodeSetOrder: did not set desired order\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* dense output order */
  flag = ARKodeSetDenseOrder(arkode_mem, 0);
  if (check_flag(&flag, "ARKodeSetDenseOrder", 1)) return 1;
  if (ark_mem->ark_dense_q != 0) {
    printf("Error in ARKodeSetDenseOrder: did not set desired order\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* linear problem */
  flag = ARKodeSetLinear(arkode_mem);
  if (check_flag(&flag, "ARKodeSetLinear", 1)) return 1;
  if (ark_mem->ark_linear != TRUE) {
    printf("Error in ARKodeSetLinear: did not set solver linearity\n");
    return 1;
  }

  /* nonlinear problem */
  flag = ARKodeSetNonlinear(arkode_mem);
  if (check_flag(&flag, "ARKodeSetNonlinear", 1)) return 1;
  if (ark_mem->ark_linear != FALSE) {
    printf("Error in ARKodeSetNonlinear: did not set solver nonlinearity\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* explicit problem */
  flag = ARKodeSetExplicit(arkode_mem);
  if (check_flag(&flag, "ARKodeSetExplicit", 1)) return 1;
  if (ark_mem->ark_explicit != TRUE) {
    printf("Error in ARKodeSetExplicit: did not set ark_explicit\n");
    return 1;
  }
  if (ark_mem->ark_implicit != FALSE) {
    printf("Error in ARKodeSetExplicit: did not set ark_implicit\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* implicit problem */
  flag = ARKodeSetImplicit(arkode_mem);
  if (check_flag(&flag, "ARKodeSetImplicit", 1)) return 1;
  if (ark_mem->ark_explicit != FALSE) {
    printf("Error in ARKodeSetImplicit: did not set ark_explicit\n");
    return 1;
  }
  if (ark_mem->ark_implicit != TRUE) {
    printf("Error in ARKodeSetImplicit: did not set ark_implicit\n");
    return 1;
  }

  /* imex problem */
  flag = ARKodeSetImEx(arkode_mem);
  if (check_flag(&flag, "ARKodeSetImEx", 1)) return 1;
  if (ark_mem->ark_explicit != FALSE) {
    printf("Error in ARKodeSetImEx: did not set ark_explicit\n");
    return 1;
  }
  if (ark_mem->ark_implicit != FALSE) {
    printf("Error in ARKodeSetImEx: did not set ark_implicit\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* ERK table */
  int s = 2;
  int q = 7;
  int p = 8;
  realtype c[]  = {0.1, 1.0};
  realtype A[]  = {0.2, 0.3, 0.4, 0.5};
  realtype b[]  = {0.6, 0.7};
  realtype be[] = {0.8, 0.9};
  flag = ARKodeSetERKTable(arkode_mem, s, q, p, c, A, b, be);
  if (check_flag(&flag, "ARKodeSetERKTable", 1)) return 1;
  if (ark_mem->ark_stages != s) {
    printf("Error in ARKodeSetERKTable: did not set ark_stages\n");
    return 1;
  }
  if (ark_mem->ark_q != q) {
    printf("Error in ARKodeSetERKTable: did not set ark_q\n");
    return 1;
  }
  if (ark_mem->ark_p != p) {
    printf("Error in ARKodeSetERKTable: did not set ark_p\n");
    return 1;
  }
  if (ark_mem->ark_c[0] != c[0]) {
    printf("Error in ARKodeSetERKTable: did not set ark_c\n");
    return 1;
  }
  if (ark_mem->ark_b[0] != b[0]) {
    printf("Error in ARKodeSetERKTable: did not set ark_b\n");
    return 1;
  }
  if (ark_mem->ark_b2[0] != be[0]) {
    printf("Error in ARKodeSetERKTable: did not set ark_b2\n");
    return 1;
  }
  if (ark_mem->ark_Ae[0] != A[0]) {
    printf("Error in ARKodeSetERKTable: did not set ark_Ae\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* IRK table */
  realtype A_[]  = {0.4, 0.5, 0.6, 0.7};
  flag = ARKodeSetIRKTable(arkode_mem, s, q, p, c, A_, b, be);
  if (check_flag(&flag, "ARKodeSetIRKTable", 1)) return 1;
  if (ark_mem->ark_stages != s) {
    printf("Error in ARKodeSetIRKTable: did not set ark_stages\n");
    return 1;
  }
  if (ark_mem->ark_q != q) {
    printf("Error in ARKodeSetIRKTable: did not set ark_q\n");
    return 1;
  }
  if (ark_mem->ark_p != p) {
    printf("Error in ARKodeSetIRKTable: did not set ark_p\n");
    return 1;
  }
  if (ark_mem->ark_c[0] != c[0]) {
    printf("Error in ARKodeSetIRKTable: did not set ark_c\n");
    return 1;
  }
  if (ark_mem->ark_b[0] != b[0]) {
    printf("Error in ARKodeSetIRKTable: did not set ark_b\n");
    return 1;
  }
  if (ark_mem->ark_b2[0] != be[0]) {
    printf("Error in ARKodeSetIRKTable: did not set ark_b2\n");
    return 1;
  }
  if (ark_mem->ark_Ai[0] != A_[0]) {
    printf("Error in ARKodeSetIRKTable: did not set ark_Ai\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* ARK tables */
  flag = ARKodeSetARKTables(arkode_mem, s, q, p, c, A_, A, b, be);
  if (check_flag(&flag, "ARKodeSetARKTable", 1)) return 1;
  if (ark_mem->ark_stages != s) {
    printf("Error in ARKodeSetARKTable: did not set ark_stages\n");
    return 1;
  }
  if (ark_mem->ark_q != q) {
    printf("Error in ARKodeSetARKTable: did not set ark_q\n");
    return 1;
  }
  if (ark_mem->ark_p != p) {
    printf("Error in ARKodeSetARKTable: did not set ark_p\n");
    return 1;
  }
  if (ark_mem->ark_c[0] != c[0]) {
    printf("Error in ARKodeSetARKTable: did not set ark_c\n");
    return 1;
  }
  if (ark_mem->ark_b[0] != b[0]) {
    printf("Error in ARKodeSetARKTable: did not set ark_b\n");
    return 1;
  }
  if (ark_mem->ark_b2[0] != be[0]) {
    printf("Error in ARKodeSetARKTable: did not set ark_b2\n");
    return 1;
  }
  if (ark_mem->ark_Ae[0] != A[0]) {
    printf("Error in ARKodeSetARKTable: did not set ark_Ae\n");
    return 1;
  }
  if (ark_mem->ark_Ai[0] != A_[0]) {
    printf("Error in ARKodeSetARKTable: did not set ark_Ai\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* ERK table num */
  flag = ARKodeSetERKTableNum(arkode_mem, 4);
  if (check_flag(&flag, "ARKodeSetERKTableNum", 1)) return 1;
  if (ark_mem->ark_stages == 0) {
    printf("Error in ARKodeSetERKTableNum: did not set ark_stages\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* IRK table num */
  flag = ARKodeSetIRKTableNum(arkode_mem, 24);
  if (check_flag(&flag, "ARKodeSetIRKTableNum", 1)) return 1;
  if (ark_mem->ark_stages == 0) {
    printf("Error in ARKodeSetIRKTableNum: did not set ark_stages\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* ARK table num */
  flag = ARKodeSetARKTableNum(arkode_mem, 22, 6);
  if (check_flag(&flag, "ARKodeSetARKTableNum", 1)) return 1;
  if (ark_mem->ark_stages == 0) {
    printf("Error in ARKodeSetARKTableNum: did not set ark_stages\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* max num steps */
  flag = ARKodeSetMaxNumSteps(arkode_mem, 95487);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;
  if (ark_mem->ark_mxstep != 95487) {
    printf("Error in ARKodeSetMaxNumSteps: did not set value\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* max hnil warns */
  flag = ARKodeSetMaxHnilWarns(arkode_mem, 24);
  if (check_flag(&flag, "ARKodeSetMaxHnilWarns", 1)) return 1;
  if (ark_mem->ark_mxhnil != 24) {
    printf("Error in ARKodeSetMaxHnilWarns: did not set value\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* init step */
  realtype hin=12.34;
  flag = ARKodeSetInitStep(arkode_mem, hin);
  if (check_flag(&flag, "ARKodeSetInitStep", 1)) return 1;
  if (ark_mem->ark_hin != hin) {
    printf("Error in ARKodeSetInitStep: did not set value\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* min step */
  realtype hmin=0.00123;
  flag = ARKodeSetMinStep(arkode_mem, hmin);
  if (check_flag(&flag, "ARKodeSetMinStep", 1)) return 1;
  if (ark_mem->ark_hmin != hmin) {
    printf("Error in ARKodeSetMinStep: did not set value\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* max step */
  realtype hmax=2.0;
  realtype hmaxinv=1.0/hmax;
  flag = ARKodeSetMaxStep(arkode_mem, hmax);
  if (check_flag(&flag, "ARKodeSetMaxStep", 1)) return 1;
  if (ark_mem->ark_hmax_inv != hmaxinv) {
    printf("Error in ARKodeSetMaxStep: did not set value\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* tstop */
  realtype tstop=8.0;
  flag = ARKodeSetStopTime(arkode_mem, tstop);
  if (check_flag(&flag, "ARKodeSetStopTime", 1)) return 1;
  if (ark_mem->ark_tstop != tstop) {
    printf("Error in ARKodeSetStopTime: did not set tstop\n");
    return 1;
  }
  if (ark_mem->ark_tstopset != TRUE) {
    printf("Error in ARKodeSetStopTime: did not set flag\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* adaptivity parameters */
  flag = ARKodeSetCFLFraction(arkode_mem, 0.2);
  if (check_flag(&flag, "ARKodeSetCFLFraction", 1)) return 1;
  if (ark_mem->ark_hadapt_cfl != 0.2) {
    printf("Error in ARKodeSetCFLFraction: did not set cfl fraction\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  flag = ARKodeSetSafetyFactor(arkode_mem, 0.25);
  if (check_flag(&flag, "ARKodeSetSafetyFactor", 1)) return 1;
  if (ark_mem->ark_hadapt_safety != 0.25) {
    printf("Error in ARKodeSetSafetyFactor: did not set correctly\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  flag = ARKodeSetErrorBias(arkode_mem, 12.5);
  if (check_flag(&flag, "ARKodeSetErrorBias", 1)) return 1;
  if (ark_mem->ark_hadapt_bias != 12.5) {
    printf("Error in ARKodeSetErrorBias: did not set correctly\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  flag = ARKodeSetMaxGrowth(arkode_mem, 20.5);
  if (check_flag(&flag, "ARKodeSetMaxGrowth", 1)) return 1;
  if (ark_mem->ark_hadapt_growth != 20.5) {
    printf("Error in ARKodeSetMaxGrowth: did not set correctly\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  flag = ARKodeSetFixedStepBounds(arkode_mem, 0.5, 2.5);
  if (check_flag(&flag, "ARKodeSetFixedStepBounds", 1)) return 1;
  if (ark_mem->ark_hadapt_lbound != 0.5) {
    printf("Error in ARKodeSetFixedStepBounds: did not set lb correctly\n");
    return 1;
  }
  if (ark_mem->ark_hadapt_ubound != 2.5) {
    printf("Error in ARKodeSetFixedStepBounds: did not set ub correctly\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* adaptivity method itself */
  flag = ARKodeSetAdaptivityMethod(arkode_mem, 2, 1, 1, NULL);
  if (check_flag(&flag, "ARKodeSetAdaptivityMethod", 1)) return 1;
  if (ark_mem->ark_hadapt_imethod != 2) {
    printf("Error in ARKodeSetAdaptivityMethod: did not set method\n");
    return 1;
  }
  if (!ark_mem->ark_hadapt_pq) {
    printf("Error in ARKodeSetAdaptivityMethod: did not set pq flag\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* adaptivity function */
  flag = ARKodeSetAdaptivityFn(arkode_mem, hfun, NULL);
  if (check_flag(&flag, "ARKodeSetAdaptivityFn", 1)) return 1;
  if (ark_mem->ark_hadapt_imethod != -1) {
    printf("Error in ARKodeSetAdaptivityFn: did not set method number\n");
    return 1;
  }
  if (ark_mem->ark_hadapt != hfun) {
    printf("Error in ARKodeSetAdaptivityFn: did not set method function\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* adaptivity constants */
  flag = ARKodeSetMaxFirstGrowth(arkode_mem, 15.0);
  if (check_flag(&flag, "ARKodeSetMaxFirstGrowth", 1)) return 1;
  if (ark_mem->ark_etamx1 != 15.0) {
    printf("Error in ARKodeSetMaxFirstGrowth: did not set etamx1\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  flag = ARKodeSetMaxEFailGrowth(arkode_mem, 0.2);
  if (check_flag(&flag, "ARKodeSetMaxEFailGrowth", 1)) return 1;
  if (ark_mem->ark_etamxf != 0.2) {
    printf("Error in ARKodeSetMaxEFailGrowth: did not set etamxf\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  flag = ARKodeSetMaxCFailGrowth(arkode_mem, 0.3);
  if (check_flag(&flag, "ARKodeSetMaxCFailGrowth", 1)) return 1;
  if (ark_mem->ark_etacf != 0.3) {
    printf("Error in ARKodeSetMaxCFailGrowth: did not set etacf\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  flag = ARKodeSetSmallNumEFails(arkode_mem, 4);
  if (check_flag(&flag, "ARKodeSetSmallNumEFails", 1)) return 1;
  if (ark_mem->ark_small_nef != 4) {
    printf("Error in ARKodeSetSmallNumEFails: did not set small_nef\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* Newton constants */
  flag = ARKodeSetNewtonCRDown(arkode_mem, 0.1);
  if (check_flag(&flag, "ARKodeSetNewtonCRDown", 1)) return 1;
  if (ark_mem->ark_crdown != 0.1) {
    printf("Error in ARKodeSetNewtonCRDown: did not set crdown\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  flag = ARKodeSetNewtonRDiv(arkode_mem, 0.2);
  if (check_flag(&flag, "ARKodeSetNewtonRDiv", 1)) return 1;
  if (ark_mem->ark_rdiv != 0.2) {
    printf("Error in ARKodeSetNewtonRDiv: did not set rdiv\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* LSetup constants */
  flag = ARKodeSetDeltaGammaMax(arkode_mem, 0.1);
  if (check_flag(&flag, "ARKodeSetDeltaGammaMax", 1)) return 1;
  if (ark_mem->ark_dgmax != 0.1) {
    printf("Error in ARKodeSetDeltaGammaMax: did not set dgmax\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  flag = ARKodeSetMaxStepsBetweenLSet(arkode_mem, 13);
  if (check_flag(&flag, "ARKodeSetMaxStepsBetweenLSet", 1)) return 1;
  if (ark_mem->ark_msbp != 13) {
    printf("Error in ARKodeSetMaxStepsBetweenLSet: did not set msbp\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* predictor method */
  flag = ARKodeSetPredictorMethod(arkode_mem, 3);
  if (check_flag(&flag, "ARKodeSetPredictorMethod", 1)) return 1;
  if (ark_mem->ark_predictor != 3) {
    printf("Error in ARKodeSetPredictorMethod: did not set value\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* stability function */
  flag = ARKodeSetStabilityFn(arkode_mem, sfun, NULL);
  if (check_flag(&flag, "ARKodeSetStabilityFn", 1)) return 1;
  if (ark_mem->ark_expstab != sfun) {
    printf("Error in ARKodeSetAdaptivityFn: did not set function\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* max err test fails */
  flag = ARKodeSetMaxErrTestFails(arkode_mem, 17);
  if (check_flag(&flag, "ARKodeSetMaxErrTestFails", 1)) return 1;
  if (ark_mem->ark_maxnef != 17) {
    printf("Error in ARKodeSetMaxErrTestFails: did not set value\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* max conv fails */
  flag = ARKodeSetMaxConvFails(arkode_mem, 12);
  if (check_flag(&flag, "ARKodeSetMaxConvFails", 1)) return 1;
  if (ark_mem->ark_maxncf != 12) {
    printf("Error in ARKodeSetMaxConvFails: did not set value\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* max nonlin iters fails */
  flag = ARKodeSetMaxNonlinIters(arkode_mem, 7);
  if (check_flag(&flag, "ARKodeSetMaxNonlinIters", 1)) return 1;
  if (ark_mem->ark_maxcor != 7) {
    printf("Error in ARKodeSetMaxNonlinIters: did not set value\n");
    return 1;
  }
  flag = ARKodeSetDefaults(arkode_mem);
  if (check_flag(&flag, "ARKodeSetDefaults", 1)) return 1;

  /* nonlin conv coef */
  realtype nlscoef = 0.01;
  flag = ARKodeSetNonlinConvCoef(arkode_mem, nlscoef);
  if (check_flag(&flag, "ARKodeSetNonlinConvCoef", 1)) return 1;
  if (ark_mem->ark_nlscoef != nlscoef) {
    printf("Error in ARKodeSetNonlinConvCoef: did not set value\n");
    return 1;
  }

  /* if we made it here, all tests passed */
  printf("All tests passed\n");


  /* Free y vector */
  N_VDestroy_Serial(y);

  /* Free integrator memory */
  ARKodeFree(&arkode_mem);

  return 0;
}


/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* fe routine template */
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  N_VConst(0.0, ydot);
  return 0;
}

/* fi routine template */
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  N_VConst(0.0, ydot);
  return 0;
}

/* error handler routine template */
static void efun(int error_code, const char *module, 
		 const char *function, char *msg, void *user_data) {
  return;
}

/* adaptivity routine template */
static int hfun(N_Vector y, realtype t, realtype h1, realtype h2, 
		realtype h3, realtype e1, realtype e2, realtype e3, 
		int q, int p, realtype *hnew, void *user_data) {
  *hnew = 1.0;
  return 0;
}

/* stability routine template */
static int sfun(N_Vector y, realtype t, realtype *hstab, void *user_data) {
  *hstab = 1.0;
  return 0;
}





/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer  
*/
static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return 1; }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return 1; }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return 1; }

  return 0;
}


/*---- end of file ----*/
