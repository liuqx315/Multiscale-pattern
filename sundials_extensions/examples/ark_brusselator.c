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
 * We have 3 different testing scenarios:
 *
 * Test 1:  u0=3.9,  v0=1.1,  w0=2.8,  a=1.2,  b=2.5,  ep=1.0e-5
 *    Here, all three components exhibit a rapid transient change 
 *    during the first 0.2 time units, followed by a slow and 
 *    smooth evolution.
 *
 * Test 2:  u0=1.2,  v0=3.1,  w0=3,  a=1,  b=3.5,  ep=5.0e-6
 *    Here, w experiences a fast initial transient, jumping 0.5 
 *    within a few steps.  All values proceed smoothly until 
 *    around t=6.5, when both u and v undergo a sharp transition, 
 *    with u increaseing from around 0.5 to 5 and v decreasing 
 *    from around 6 to 1 in less than 0.5 time units.  After this
 *    transition, both u and v continue to evolve somewhat 
 *    rapidly for another 1.4 time units, and finish off smoothly.
 *
 * Test 3:  u0=3,  v0=3,  w0=3.5,  a=0.5,  b=3,  ep=5.0e-4
 *    Here, all components undergo very rapid initial transients 
 *    during the first 0.3 time units, and all then proceed very 
 *    smoothly for the remainder of the simulation.
 *
 * These tests are selected within the input file (test = {1,2,3}), 
 * with the default set to test 2 in case the input is invalid.
 * Also in the input file, we allow specification of the desired 
 * relative and absolute tolerances.
 * 
 * This program solves the problem with the DIRK method, using a
 * Newton iteration with the ARKDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
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
static int Jac(long int N, realtype t,
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
  realtype dTout = RCONST(1.0);
  int Nt = ceil(Tf/dTout);
  realtype a, b, ep, u0, v0, w0;
  long int NEQ = 3;

  /* general problem variables */
  int flag;
  N_Vector y = NULL;
  void *arkode_mem = NULL;

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

  /* open solver diagnostics output file for writing */
  FILE *DFID;
  DFID=fopen("diags_ark_brusselator.txt","w");
  
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

  /* Set initial conditions into y */
  NV_Ith_S(y,0) = u0;
  NV_Ith_S(y,1) = v0;
  NV_Ith_S(y,2) = w0;

  /* Call ARKodeCreate to create the solver memory */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return(1);
  
  /* Call ARKodeInit to initialize the integrator memory and specify the
     user's right hand side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y */
  flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return(1);

  /* Call ARKodeSetUserData to pass rdata to user functions */
  flag = ARKodeSetUserData(arkode_mem, (void *) rdata);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return(1);

  /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
  flag = ARKodeSetDiagnostics(arkode_mem, DFID);
  if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return(1);

  /* Call ARKodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return(1);

  /* Call ARKDense to specify the ARKDENSE dense linear solver */
  flag = ARKDense(arkode_mem, NEQ);
  if (check_flag(&flag, "ARKDense", 1)) return(1);

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);
  if (check_flag(&flag, "ARKDlsSetDenseJacFn", 1)) return(1);

  /* In loop, call ARKode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t = T0;
  realtype tout = dTout;
  realtype u, v, w;
  printf("        t           u           v           w\n");
  printf("   -------------------------------------------\n");
  int iout;
  for (iout=0; iout<Nt; iout++) {
    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
    u = NV_Ith_S(y,0);
    v = NV_Ith_S(y,1);
    w = NV_Ith_S(y,2);
    printf("  %10.6f  %10.6f  %10.6f  %10.6f\n", t, u, v, w);

    if (check_flag(&flag, "ARKode", 1)) break;
    if (flag == ARK_SUCCESS) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
  }
  printf("   -------------------------------------------\n");

  /* Print some final statistics */
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

  printf("\nFinal Solver Statistics:\n");
  printf("   Total internal solver steps = %li (acc = %li,  conv = %li)\n", 
	 nst, nst_a, nst_c);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of linear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n\n", netf);

  /* Free y vector */
  N_VDestroy_Serial(y);

  /* Free integrator memory */
  ARKodeFree(&arkode_mem);

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
