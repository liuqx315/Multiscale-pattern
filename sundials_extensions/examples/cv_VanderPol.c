/* -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following test simulates the Van der Pol oscillator.  This 
 * is an ODE system with 2 components, Y = [u,v], satisfying the 
 * equations,
 *    du/dt = v
 *    dv/dt = (v - v*u^2 - ep2*u)/ep.
 * for t in [0,Tf] and initial conditions Y0 = [u0,v0].
 * 
 * We have 3 different testing scenarios:
 *
 * Test 1:  Tf=1.5, u0=-2, v0=-2.3553013976081, ep=1.0e-3, ep2=1.0
 *    Here, v quickly increases to around 0.5, then both change very 
 *    slowly until around t=0.6, when v increases exponentially to a 
 *    value of near 12, before it suddenly drops down to around -1.0.  
 *    When v drops, u jumps to around 2.0.  Both proceed slowly 
 *    afterwards for the remainder of the time interval.
 *
 * Test 2:  Tf=12,  u0=2,  v0=0,  ep=ep2=0.2
 *    Here, both u and v proceed smoothly until around t=5, when v 
 *    dips to -7.0 and then returns to around 0.5 within an interval
 *    of around 1 unit time.  During that change, u drops from 
 *    around 1.0 to around -2.0.  This process then repeats but in 
 *    the opposite direction at around t=11.
 *
 * Test 3:  Tf=1000,  u0=2,  v0=0,  ep=ep2=1.0e-3
 *    Here, both u and v proceed smoothly until around t=800, when 
 *    within less than 10 time units u suddenly plunges from around 
 *    1.1 to -2.0.  Meanwhile, v remains relatively constant 
 *    throughout the entire simulation.
 *
 * These tests are selected within the input file (test = {1,2,3}), 
 * with the default set to test 1 in case the input is invalid.
 * Also in the input file, we allow specification of the desired 
 * relative and absolute tolerances.
 * 
 * This program solves the problem with the BDF method, using a
 * Newton iteration with the CVDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 *
 * 100 outputs are printed at equal intervals, and run statistics 
 * are printed at the end.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>

/* Header files with a description of contents used */
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */



/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);



/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0=0.0;
  realtype Tf, dTout, ep, ep2, u0, v0;
  int NEQ = 2;
  int Nt = 100;

  /* general problem variables */
  int flag, flag2;
  N_Vector y = NULL;
  N_Vector ytrue = NULL;
  void *cvode_mem = NULL;
  void *cvtrue_mem = NULL;

  /* read problem parameter and tolerances from input file:
     test   - test problem choice
     reltol - desired relative tolerance
     abstol - desired absolute tolerance */
  int test;
  double reltol_, abstol_;
  FILE *FID;
  FID=fopen("input_VanderPol.txt","r");
  fscanf(FID,"  test = %i\n", &test);
  fscanf(FID,"  reltol = %lf\n", &reltol_);
  fscanf(FID,"  abstol = %lf\n", &abstol_);
  fclose(FID);

  /* convert the inputs to 'realtype' format */
  realtype reltol = reltol_;
  realtype abstol = abstol_;
  realtype reltol2 = reltol_*1.0e-3;
  realtype abstol2 = abstol_*1.0e-3;

  /* set up the test problem according to the desired input */
  if (test == 2) {
    Tf = RCONST(12.0);
    u0 = RCONST(2.0);
    v0 = RCONST(0.0);
    ep = RCONST(0.2);
    ep2 = RCONST(0.2);
  } else if (test == 3) {
    Tf = RCONST(1000.0);
    u0 = RCONST(2.0);
    v0 = RCONST(0.0);
    ep = RCONST(1.0e-3);
    ep2 = RCONST(1.0e-3);
  } else {
    Tf = RCONST(1.5);
    u0 = RCONST(-2.0);
    v0 = RCONST(-2.3553013976081);
    ep = RCONST(1.0e-3);
    ep2 = RCONST(1.0);
  }
  dTout = (Tf-T0)/Nt;

  /* set user data to contain stiffness parameters */
  realtype rdata[2] = {ep, ep2};

  /* Initial problem output */
  printf("\nVan der Pol ODE test problem:\n");
  printf("    time interval:  t in [%g,%g]\n",T0,Tf);
  printf("    output interval = %g\n",dTout);
  printf("    initial conditions:  u0 = %g,  v0 = %g\n",u0,v0);
  printf("    stiffness parameters:  ep = %g,  ep2 = %g\n",ep,ep2);
  printf("    tolerances:  reltol = %.1e,  abstol = %.1e\n\n",reltol,abstol);


  /* Create serial vector of length NEQ for initial condition */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  ytrue = N_VNew_Serial(NEQ);
  if (check_flag((void *)ytrue, "N_VNew_Serial", 0)) return(1);

  /* Set initial conditions into y, ytrue */
  NV_Ith_S(y,0) = u0;
  NV_Ith_S(y,1) = v0;
  NV_Ith_S(ytrue,0) = u0;
  NV_Ith_S(ytrue,1) = v0;

  /* Call CVodeCreate to create the solver memory and specify the 
     Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
  cvtrue_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvtrue_mem, "CVodeCreate", 0)) return(1);
  
  /* Call CVodeInit to initialize the integrator memory and specify the
     user's right hand side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y */
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);
  flag = CVodeInit(cvtrue_mem, f, T0, ytrue);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  /* Call CVodeSetUserData to pass rdata to user functions */
  flag = CVodeSetUserData(cvode_mem, (void *) rdata);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);
  flag = CVodeSetUserData(cvtrue_mem, (void *) rdata);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* Set additional solver parameters for reference solution */
  flag = CVodeSetMaxNumSteps(cvtrue_mem, 100000);

  /* Call CVodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);
  flag = CVodeSStolerances(cvtrue_mem, reltol2, abstol2);
  if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)) return(1);
  flag = CVDense(cvtrue_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)) return(1);

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);
  flag = CVDlsSetDenseJacFn(cvtrue_mem, Jac);
  if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t  = T0;
  realtype t2 = T0;
  realtype tout = dTout;
  realtype u, v, uerr, verr, errI=0.0, err2=0.0;
  printf("        t           u           v        uerr          verr\n");
  printf("   --------------------------------------------------------------\n");
  int iout;
  for (iout=0; iout<Nt; iout++) {
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    u = NV_Ith_S(y,0);
    v = NV_Ith_S(y,1);
    flag2 = CVode(cvtrue_mem, tout, ytrue, &t2, CV_NORMAL);
    uerr = fabs(NV_Ith_S(ytrue,0) - u);
    verr = fabs(NV_Ith_S(ytrue,1) - v);
    errI = (errI > verr) ? errI : verr;
    errI = (errI > uerr) ? errI : uerr;
    err2 += uerr*uerr + verr*verr;
    printf("  %10.6f  %10.6f  %10.6f  %12.5e  %12.5e\n", 
	   t, u, v, uerr, verr);

    if (check_flag(&flag, "CVode", 1)) break;
    if (flag == CV_SUCCESS) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
  }
  err2 = sqrt(err2 / 2.0 / Nt);
  printf("   --------------------------------------------------------------\n");

  /* Print some final statistics */
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf;
  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Total internal solver steps = %li\n", nst);
  printf("   Total RHS evals = %li\n", nfe);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of linear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Error: max = %g, rms = %g\n", errI, err2);
  printf("   Oversolve = %g\n\n", reltol/err2);

  /* Free y vector */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(ytrue);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);
  CVodeFree(&cvtrue_mem);

  return(0);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* f routine to compute the ODE RHS function f(t,y). */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;
  realtype ep  = rdata[0];
  realtype ep2 = rdata[1];
  realtype u = NV_Ith_S(y,0);
  realtype v = NV_Ith_S(y,1);

  /* du/dt = v */
  NV_Ith_S(ydot,0) = v;

  /* dv/dt = (v - v*u^2 - ep2*u)/ep */
  NV_Ith_S(ydot,1) = (v - v*u*u - ep2*u)/ep;

  return(0);
}

/* Jacobian routine to compute J(t,y) = df/dy. */

static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *rdata = (realtype *) user_data;
  realtype ep  = rdata[0];
  realtype ep2 = rdata[1];
  realtype u = NV_Ith_S(y,0);
  realtype v = NV_Ith_S(y,1);

  /* du/dt = v */
  DENSE_ELEM(J,0,0) = 0.0;
  DENSE_ELEM(J,0,1) = 1.0;

  /* dv/dt = (v - v*u^2 - ep2*u)/ep */
  DENSE_ELEM(J,1,0) = -(2.0*v*u + ep2)/ep;
  DENSE_ELEM(J,1,1) = (1.0 - u*u)/ep;

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
