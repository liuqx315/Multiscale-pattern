/*---------------------------------------------------------------
 $Revision: $
 $Date: $
-----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Example problem:
 
 The following is a simple example problem with analytical 
 solution,
     dy/dt = (t^2+1)*exp(-y)
 for t in the interval [0.0, -3.0] (goes backwards in time), 
 with initial condition: y=3. This has analytical solution 
      y(t) = log(exp(3) + t + (1/3)*t^3)
 
 This program solves the problem with the DIRK method,
 Newton iteration with the ARKDENSE dense linear solver, and a
 user-supplied Jacobian routine.
 Output is printed every 1.0 units of time (10 total).
 Run statistics (optional outputs) are printed at the end.
---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <math.h>
#include <arkode/arkode.h>
#include <nvector/nvector_serial.h>
#include <arkode/arkode_dense.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);

/* Parameter input helper function */
int init_from_file(void *ark_mem, char *fname, ARKRhsFn f, 
		   ARKRhsFn fe, ARKRhsFn fi, realtype T0, 
		   N_Vector y0, int *ImEx, int *dorder, 
		   int *fxpt, realtype *RTol, realtype *ATol);


/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);
  realtype Tf = RCONST(-3.0);
  realtype dTout = RCONST(-0.1);
  long int NEQ = 1;

  /* declare solver parameters */
  int flag, dense_order, imex, fixedpt;

  /* general problem variables */
  int idense;
  N_Vector y = NULL;
  void *arkode_mem = NULL;

  /* open solver diagnostics output file for writing */
  FILE *DFID;
  DFID=fopen("diags_ark_analytic_nonlin_backwards.txt","w");
  
  /* Initial problem output */
  printf("\nAnalytical backwards ODE test problem:\n");

  /* Create serial vector of length NEQ for initial condition */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;

  /* Initialize y to 0 */
  NV_Ith_S(y,0) = RCONST(3.0);

  /* Call ARKodeCreate to create the solver memory */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return 1;
  
  /* Call init_from_file helper routine to read and set solver parameters */
  realtype rtol, atol;
  flag = init_from_file(arkode_mem, "solve_params.txt", f, NULL, NULL, T0,
			y, &imex, &dense_order, &fixedpt, &rtol, &atol);
  if (check_flag(&flag, "init_from_file", 1)) return 1;
  if (rtol <= 0.0)  rtol = 1.e-6;
  if (atol <= 0.0)  atol = 1.e-10;
  realtype reltol = rtol;
  realtype abstol = atol;

  /* If (dense_order == -1), tell integrator to use tstop */
  if (dense_order == -1) {
    idense = 0;
  } else {    /* otherwise tell integrator to use dense output */
    idense = 1;
  }

  /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
  flag = ARKodeSetDiagnostics(arkode_mem, DFID);
  if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return 1;

  /* Call ARKodeSetMaxNumSteps to increase default (for testing) */
  flag = ARKodeSetMaxNumSteps(arkode_mem, 10000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;

  /* Call ARKodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Call ARKDense to specify the ARKDENSE dense linear solver */
  flag = ARKDense(arkode_mem, NEQ);
  if (check_flag(&flag, "ARKDense", 1)) return 1;

  /* Set the Jacobian routine to Jac (user-supplied) */
  if (imex == 0) {
    flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);
    if (check_flag(&flag, "ARKDlsSetDenseJacFn", 1)) return 1;
  }

  /* Write all solver parameters to stdout */
  printf("\n");
  flag = ARKodeWriteParameters(arkode_mem, stdout);
  if (check_flag(&flag, "ARKodeWriteParameters", 1)) return 1;

  /* In loop, call ARKode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t = T0;
  realtype tout = T0+dTout;
  realtype u, uerr, err2=0.0, errI=0.0;
  int Nt=0;
  printf("        t           u         error\n");
  printf("   ------------------------------------\n");
  while (t - Tf > 1.0e-8) {
    if (!idense) 
      flag = ARKodeSetStopTime(arkode_mem, tout);
    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_flag(&flag, "ARKode", 1)) break;
    if (flag >= 0) {
      tout += dTout;
      tout = (tout < Tf) ? Tf : tout;
    }
    u = NV_Ith_S(y,0);
    uerr = fabs(u - log(exp(3.0)+t+t*t*t/3.0));
    errI = (errI > uerr) ? errI : uerr;
    err2 += uerr*uerr;
    Nt++;
    printf("  %10.6f  %10.6f  %12.5e\n", t, u, uerr);
  }
  err2 = sqrt(err2 / Nt);
  printf("   ------------------------------------\n");

  /* Print some final statistics */
  long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKodeGetNumSteps", 1);
  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ARKodeGetNumStepAttempts", 1);
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
  printf("   Internal solver steps = %li (attempted = %li)\n", 
	 nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of nonlinear iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Error: max = %g, rms = %g\n", errI, err2);
  printf("   Oversolve = %g\n\n", reltol/(err2+1.e-10*reltol));

  /* Free y vector */
  N_VDestroy_Serial(y);

  /* Free integrator memory */
  ARKodeFree(&arkode_mem);

  /* close solver diagnostics output file */
  fclose(DFID);

  return 0;
}


/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  NV_Ith_S(ydot,0) = (t*t+1.0)*exp(-NV_Ith_S(y,0));
  return 0;
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  DENSE_ELEM(J,0,0) = -(t*t+1.0)*exp(-NV_Ith_S(y,0));
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
