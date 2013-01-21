/* -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following test simulates a brusselator problem from chemical 
 * kinetics.  This is an ODE system with 3 components, Y = [u,v,w], 
 * satisfying the equations,
 *    du/dt = a - (w+1)*u + v*u^2
 *    dv/dt = w*u - v*u^2
 *    dw/dt = (b-w)/ep - w*u
 * for t in the interval [0.0, 10.0], with initial conditions 
 * Y0 = [u0,v0,w0]. 
 * 
 * This testing routine is designed to see whether ARKode can 
 * handle two scenarios:
 *   (a) a user-supplied Butcher table
 *   (b) a change in method mid-integration
 * 
 * This program solves the problem first using a DIRK method, 
 * followed with an ARK method.  The DIRK method uses a custom 
 * Butcher table, while the ARK method is the default method of 
 * order 4.
 *
 * 100 outputs are printed at equal intervals, and run statistics 
 * are printed at the end.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>

/* Header files with a description of contents used */
#include <arkode/arkode.h>           /* prototypes for ARKODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <arkode/arkode_dense.h>     /* prototype for ARKDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */



/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int JacI(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);



/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);
  realtype Tf = RCONST(10.0);
  realtype dTout = RCONST(0.5);
  int Nt = ceil(Tf/dTout);
  int Nt1 = Nt/2;
  int Nt2 = Nt-Nt1;
  realtype a, b, ep, u0, v0, w0;
  long int NEQ = 3;

  /* declare solver parameters */
  int flag = 0;

  /* general problem variables */
  N_Vector y = NULL;
  N_Vector ytrue = NULL;
  void *arkode_mem = NULL;
  void *arktrue_mem = NULL;

  /* read problem parameter and tolerances from input file:
     test   - test problem choice
     reltol - desired relative tolerance
     abstol - desired absolute tolerance */
  int test;
  double reltol_, abstol_;
  FILE *FID;
  FID=fopen("input_brusselator.txt","r");
  fscanf(FID,"  test = %i\n", &test);
  fscanf(FID,"  reltol = %lf\n", &reltol_);
  fscanf(FID,"  abstol = %lf\n", &abstol_);
  fclose(FID);

  /* convert the inputs to 'realtype' format */
  realtype reltol = reltol_;
  realtype abstol = abstol_;
  realtype reltol2 = reltol_*1.0e-3;
  realtype abstol2 = abstol_*1.0e-3;

  /* open solver diagnostics output file for writing */
  FILE *DFID;
  DFID=fopen("diags_ark_brusselator_custom.txt","w");
  
  /* set up the test problem according to the desired input */
  if (test == 1) {
    u0 = RCONST(3.9);
    v0 = RCONST(1.1);
    w0 = RCONST(2.8);
    a  = RCONST(1.2);
    b  = RCONST(2.5);
    ep = RCONST(1.0e-5);
  } else if (test == 3) {
    u0 = RCONST(3.0);
    v0 = RCONST(3.0);
    w0 = RCONST(3.5);
    a  = RCONST(0.5);
    b  = RCONST(3.0);
    ep = RCONST(5.0e-4);
  } else {
    u0 = RCONST(1.2);
    v0 = RCONST(3.1);
    w0 = RCONST(3.0);
    a  = RCONST(1.0);
    b  = RCONST(3.5);
    ep = RCONST(5.0e-6);
  }

  /* set user data to contain problem-defining parameters */
  realtype rdata[3] = {a, b, ep};

  /* Initial problem output */
  printf("\nBrusselator ODE test problem:\n");
  printf("    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n",u0,v0,w0);
  printf("    problem parameters:  a = %g,  b = %g,  ep = %g\n",a,b,ep);
  printf("    reltol = %.1e,  abstol = %.1e\n\n",reltol,abstol);


  /* Create serial vector of length NEQ for initial condition */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  ytrue = N_VNew_Serial(NEQ);
  if (check_flag((void *)ytrue, "N_VNew_Serial", 0)) return(1);

  /* Set initial conditions into y */
  NV_Ith_S(y,0) = u0;
  NV_Ith_S(y,1) = v0;
  NV_Ith_S(y,2) = w0;
  NV_Ith_S(ytrue,0) = u0;
  NV_Ith_S(ytrue,1) = v0;
  NV_Ith_S(ytrue,2) = w0;

  /* Call ARKodeCreate to create the solver memory */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return(1);
  arktrue_mem = ARKodeCreate();
  if (check_flag((void *)arktrue_mem, "ARKodeCreate", 0)) return(1);
  
  /* Call ARKodeInit to initialize the integrator memory and specify the
     user's right hand side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y.  Start out using a DIRK solver. */
  printf("Starting in purely implicit mode\n");
  flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return(1);

  /* Compute reference solution with default implicit method */
  flag = ARKodeInit(arktrue_mem, NULL, f, T0, ytrue);
  if (check_flag(&flag, "ARKodeInit", 1)) return(1);

  /* Call ARKodeSetUserData to pass rdata to user functions */
  flag = ARKodeSetUserData(arkode_mem, (void *) rdata);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return(1);
  flag = ARKodeSetUserData(arktrue_mem, (void *) rdata);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return(1);

  /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
  flag = ARKodeSetDiagnostics(arkode_mem, DFID);
  if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return(1);

  /* Set customized DIRK method to start out with */
  int s = 4;
  int q = 3;
  int p = 2;
  realtype c[] = {RCONST(0.0), RCONST(0.871733043), 
		  RCONST(1.0), RCONST(1.0)};
  realtype bmeth[] = {RCONST(0.308809969973036), RCONST(1.490563388254106), 
		      RCONST(-1.235239879727145), RCONST(0.4358665215)};
  realtype bembed[] = {RCONST(0.490563388419108), RCONST(0.073570090080892), 
		       RCONST(0.4358665215), RCONST(0.0)};
  realtype A[] = {RCONST(0.0), RCONST(0.0), RCONST(0.0), RCONST(0.0), 
		  RCONST(0.4358665215), RCONST(0.4358665215), RCONST(0.0), 
		  RCONST(0.0), RCONST(0.490563388419108), 
		  RCONST(0.073570090080892), RCONST(0.4358665215), 
		  RCONST(0.0), RCONST(0.308809969973036), 
		  RCONST(1.490563388254106), RCONST(-1.235239879727145), 
		  RCONST(0.4358665215)};
  flag = ARKodeSetIRKTable(arkode_mem, s, q, p, c, A, bmeth, bembed);
  if (check_flag(&flag, "ARKodeSetIRKTable", 1)) return(1);

  /* Set additional solver parameters for reference solution */
  flag = ARKodeSetMaxNumSteps(arktrue_mem, 100000);

  /* Call ARKodeSetMaxNumSteps to increase default (for testing) */
  flag = ARKodeSetMaxNumSteps(arkode_mem, 100000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return(1);

  /* Call ARKodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return(1);
  flag = ARKodeSStolerances(arktrue_mem, reltol2, abstol2);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return(1);

  /* Call ARKDense to specify the ARKDENSE dense linear solver */
  flag = ARKDense(arkode_mem, NEQ);
  if (check_flag(&flag, "ARKDense", 1)) return(1);
  flag = ARKDense(arktrue_mem, NEQ);
  if (check_flag(&flag, "ARKDense", 1)) return(1);

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);
  if (check_flag(&flag, "ARKDlsSetDenseJacFn", 1)) return(1);
  flag = ARKDlsSetDenseJacFn(arktrue_mem, Jac);
  if (check_flag(&flag, "ARKDlsSetDenseJacFn", 1)) return(1);

  /* Write all solver parameters to stdout */
  printf("\n");
  flag = ARKodeWriteParameters(arkode_mem, stdout);
  if (check_flag(&flag, "ARKodeWriteParameters", 1)) return(1);

  /* In loop, call ARKode, print results, and test for error.
     Break out of loop when the half of the final output time has been reached */
  realtype t = T0;
  realtype t2 = T0;
  realtype tout = T0+dTout;
  realtype u, v, w, uerr, verr, werr, errI=0.0, err2=0.0;
  printf("        t           u           v           w        uerr          verr          werr\n");
  printf("   ---------------------------------------------------------------------------------------\n");
  int iout;
  for (iout=0; iout<Nt1; iout++) {
    flag = ARKodeSetStopTime(arktrue_mem, tout);
    flag = ARKode(arktrue_mem, tout, ytrue, &t2, ARK_NORMAL);
    flag = ARKodeSetStopTime(arkode_mem, tout);
    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_flag(&flag, "ARKode", 1)) break;
    if (flag >= 0) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }

    u = NV_Ith_S(y,0);
    v = NV_Ith_S(y,1);
    w = NV_Ith_S(y,2);
    uerr = fabs(NV_Ith_S(ytrue,0) - u);
    verr = fabs(NV_Ith_S(ytrue,1) - v);
    werr = fabs(NV_Ith_S(ytrue,2) - w);
    errI = (errI > verr) ? errI : verr;
    errI = (errI > uerr) ? errI : uerr;
    errI = (errI > werr) ? errI : werr;
    err2 += uerr*uerr + verr*verr + werr*werr;
    printf("  %10.6f  %10.6f  %10.6f  %10.6f  %12.5e  %12.5e  %12.5e\n", 
	   t, u, v, w, uerr, verr, werr);
  }
  printf("   ---------------------------------------------------------------------------------------\n");


  /* Print some statistics */
  long int nst, nst_a, nst_c, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKodeGetNumSteps", 1);
  flag = ARKodeGetNumAccSteps(arkode_mem, &nst_a);
  check_flag(&flag, "ARKodeGetNumAccSteps", 1);
  flag = ARKodeGetNumConvSteps(arkode_mem, &nst_c);
  check_flag(&flag, "ARKodeGetNumConvSteps", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
  flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1);
  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKodeGetNumErrTestFails", 1);
  flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1);
  flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1);
  flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKDlsGetNumJacEvals", 1);
  flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKDlsGetNumRhsEvals", 1);

  printf("\nFirst Phase Solver Statistics:\n");
  printf("   Internal solver steps = %li (acc = %li,  conv = %li)\n", 
	 nst, nst_a, nst_c);
  printf("   RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Linear solver setups = %li\n", nsetups);
  printf("   RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Number of Jacobian evaluations = %li\n", nje);
  printf("   Number of Newton iterations = %li\n", nni);
  printf("   Number of linear solver convergence failures = %li\n", ncfn);
  printf("   Number of error test failures = %li\n\n", netf);



  /* Call ARKodeReInit to re-initialize with an ARK solver. */
  printf("Re-Starting in ImEx mode\n");
  flag = ARKodeReInit(arkode_mem, fe, fi, tout, y);
  if (check_flag(&flag, "ARKodeReInit", 1)) return(1);

  /* Set desired method order */
  printf("  Setting order = 4\n");
  flag = ARKodeSetOrder(arkode_mem, 4);
  if (check_flag(&flag, "ARKodeSetOrder", 1)) return(1);

  /* Reset the Jacobian routine */
  flag = ARKDlsSetDenseJacFn(arkode_mem, JacI);
  if (check_flag(&flag, "ARKDlsSetDenseJacFn", 1)) return(1);

  /* Write all solver parameters to stdout */
  printf("\n");
  flag = ARKodeWriteParameters(arkode_mem, stdout);
  if (check_flag(&flag, "ARKodeWriteParameters", 1)) return(1);

  /* In loop, call ARKode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  printf("        t           u           v           w        uerr          verr          werr\n");
  printf("   ---------------------------------------------------------------------------------------\n");
  tout += dTout;
  for (iout=0; iout<Nt2; iout++) {
    flag = ARKodeSetStopTime(arktrue_mem, tout);
    flag = ARKode(arktrue_mem, tout, ytrue, &t2, ARK_NORMAL);
    flag = ARKodeSetStopTime(arkode_mem, tout);
    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_flag(&flag, "ARKode", 1)) break;
    if (flag >= 0) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }

    u = NV_Ith_S(y,0);
    v = NV_Ith_S(y,1);
    w = NV_Ith_S(y,2);
    uerr = fabs(NV_Ith_S(ytrue,0) - u);
    verr = fabs(NV_Ith_S(ytrue,1) - v);
    werr = fabs(NV_Ith_S(ytrue,2) - w);
    errI = (errI > verr) ? errI : verr;
    errI = (errI > uerr) ? errI : uerr;
    errI = (errI > werr) ? errI : werr;
    err2 += uerr*uerr + verr*verr + werr*werr;
    printf("  %10.6f  %10.6f  %10.6f  %10.6f  %12.5e  %12.5e  %12.5e\n", 
	   t, u, v, w, uerr, verr, werr);
  }
  err2 = sqrt(err2 / 3.0 / Nt);
  printf("   ---------------------------------------------------------------------------------------\n");

  /* Print some final statistics */
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKodeGetNumSteps", 1);
  flag = ARKodeGetNumAccSteps(arkode_mem, &nst_a);
  check_flag(&flag, "ARKodeGetNumAccSteps", 1);
  flag = ARKodeGetNumConvSteps(arkode_mem, &nst_c);
  check_flag(&flag, "ARKodeGetNumConvSteps", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
  flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1);
  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKodeGetNumErrTestFails", 1);
  flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1);
  flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1);
  flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKDlsGetNumJacEvals", 1);
  flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKDlsGetNumRhsEvals", 1);

  printf("\nSecond Phase Solver Statistics:\n");
  printf("   Internal solver steps = %li (acc = %li,  conv = %li)\n", 
	 nst, nst_a, nst_c);
  printf("   RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Linear solver setups = %li\n", nsetups);
  printf("   RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Number of Jacobian evaluations = %li\n", nje);
  printf("   Number of Newton iterations = %li\n", nni);
  printf("   Number of linear solver convergence failures = %li\n", ncfn);
  printf("   Number of error test failures = %li\n\n", netf);

  printf("\nOverall accuracy results:\n");
  printf("   Error: max = %g, rms = %g\n", errI, err2);
  printf("   Oversolve = %g\n\n", reltol/err2);

  /* Free y vector */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(ytrue);

  /* Free integrator memory */
  ARKodeFree(&arkode_mem);
  ARKodeFree(&arktrue_mem);

  /* close solver diagnostics output file */
  fclose(DFID);

  return(0);
}


/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;
  realtype a  = rdata[0];
  realtype b  = rdata[1];
  realtype ep = rdata[2];
  realtype u = NV_Ith_S(y,0);
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);

  /* du/dt = a - (w+1)*u + v*u^2 */
  NV_Ith_S(ydot,0) = a - (w+1.0)*u + v*u*u;

  /* dv/dt = w*u - v*u^2 */
  NV_Ith_S(ydot,1) = w*u - v*u*u;

  /* dw/dt = (b-w)/ep - w*u */
  NV_Ith_S(ydot,2) = (b-w)/ep - w*u;

  return(0);
}

/* fe routine to compute the explicit portion of f(t,y). */
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;
  realtype a  = rdata[0];
  realtype u = NV_Ith_S(y,0);
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);

  /* du/dt = a - (w+1)*u + v*u^2 */
  NV_Ith_S(ydot,0) = a - (w+1.0)*u + v*u*u;

  /* dv/dt = w*u - v*u^2 */
  NV_Ith_S(ydot,1) = w*u - v*u*u;

  /* dw/dt = -w*u */
  NV_Ith_S(ydot,2) = -w*u;

  return(0);
}

/* fi routine to compute the implicit portion of f(t,y). */
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;
  realtype b  = rdata[1];
  realtype ep = rdata[2];
  realtype w = NV_Ith_S(y,2);

  /* du/dt = a - (w+1)*u + v*u^2 */
  NV_Ith_S(ydot,0) = 0.0;

  /* dv/dt = w*u - v*u^2 */
  NV_Ith_S(ydot,1) = 0.0;

  /* dw/dt = (b-w)/ep - w*u */
  NV_Ith_S(ydot,2) = (b-w)/ep;

  return(0);
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *rdata = (realtype *) user_data;
  realtype ep = rdata[2];
  realtype u = NV_Ith_S(y,0);
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);

  /* du/dt = a - (w+1)*u + v*u^2 */
  DENSE_ELEM(J,0,0) = -(w+1.0) + 2.0*u*v;
  DENSE_ELEM(J,0,1) = u*u;
  DENSE_ELEM(J,0,2) = -u;

  /* dv/dt = w*u - v*u^2 */
  DENSE_ELEM(J,1,0) = w - 2.0*u*v;
  DENSE_ELEM(J,1,1) = -u*u;
  DENSE_ELEM(J,1,2) = u;

  /* dw/dt = (b-w)/ep - w*u */
  DENSE_ELEM(J,2,0) = -w;
  DENSE_ELEM(J,2,1) = 0.0;
  DENSE_ELEM(J,2,2) = -1.0/ep - u;

  return(0);
}

/* Jacobian routine to compute J(t,y) = dfi/dy. */
static int JacI(long int N, realtype t,
		N_Vector y, N_Vector fy, DlsMat J, void *user_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *rdata = (realtype *) user_data;
  realtype ep = rdata[2];
  SetToZero(J);
  DENSE_ELEM(J,2,2) = -1.0/ep;
  return(0);
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
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}


/*---- end of file ----*/
