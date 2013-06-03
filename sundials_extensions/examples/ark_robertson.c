/*---------------------------------------------------------------
 $Revision: $
 $Date: $
-----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Example problem:
 
 The following test simulates the Robertson problem, 
 corresponding to the kinetics of an autocatalytic reaction.  
 This is an ODE system with 3 components, Y = [u,v,w], satisfying
 the equations,
    du/dt = -0.04*u + 1e4*v*w
    dv/dt = 0.04*u - 1e4*v*w - 3e7*v^2
    dw/dt = 3e7*v^2
 for t in the interval [0.0, 1e11], with initial conditions 
 Y0 = [1,0,0]. 
 
 In the input file, input_robertson.txt, we allow specification 
 of the desired relative and absolute tolerances.
 
 This program solves the problem with one of the solvers, ERK, 
 DIRK or ARK.  For DIRK and ARK, implicit subsystems are solved 
 using a Newton iteration with the ARKDENSE dense linear solver, 
 and a user-supplied Jacobian routine.

 100 outputs are printed at equal intervals, and run statistics 
 are printed at the end.
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


/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);
  realtype Tf = RCONST(1.e11);
  realtype dTout = (Tf-T0)/100;
  int Nt = ceil(Tf/dTout);
  realtype u0, v0, w0, h0;
  long int NEQ = 3;

  /* declare solver parameters */
  int flag;

  /* general problem variables */
  int idense;
  N_Vector y = NULL;
  void *arkode_mem = NULL;

  /* open solver diagnostics output file for writing */
  FILE *DFID;
  DFID=fopen("diags_ark_robertson.txt","w");
  
  /* set up the initial conditions */
  u0 = RCONST(1.0);
  v0 = RCONST(0.0);
  w0 = RCONST(0.0);

  /* Initial problem output */
  printf("\nRobertson ODE test problem:\n");
  printf("    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n",u0,v0,w0);

  /* Create serial vectors of length NEQ for initial condition */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;

  /* Set initial conditions into y */
  NV_Ith_S(y,0) = u0;
  NV_Ith_S(y,1) = v0;
  NV_Ith_S(y,2) = w0;

  /* Call ARKodeCreate to create the solver memory */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return 1;
  
  /* Call ARKodeInit to initialize the integrator memory and specify the
     user's right hand side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y */
  flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return 1;

  /* Set tolerances */
  realtype reltol = 1.e-4;
  realtype abstol = 1.e-8;
  h0 = 1.e-4 * reltol;

  /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
  flag = ARKodeSetDiagnostics(arkode_mem, DFID);
  if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return 1;

  /* Set custom initial step */
  flag = ARKodeSetInitStep(arkode_mem, h0);
  if (check_flag(&flag, "ARKodeSetInitStep", 1)) return 1;

  /* Increase maximum number of error test failures */
  flag = ARKodeSetMaxErrTestFails(arkode_mem, 20);
  if (check_flag(&flag, "ARKodeSetMaxErrTestFails", 1)) return 1;

  /* Call ARKodeSetmaxNonlinIters to increase default for this problem*/
  flag = ARKodeSetMaxNonlinIters(arkode_mem, 8);
  if (check_flag(&flag, "ARKodeSetMaxNonlinIters", 1)) return 1;

  /* Call ARKodeSetNonlinConvCoef */
  flag = ARKodeSetNonlinConvCoef(arkode_mem, 1.e-7);
  if (check_flag(&flag, "ARKodeSetNonlinConvCoef", 1)) return 1;

  /* Call ARKodeSetMaxNumSteps to increase default (for testing) */
  flag = ARKodeSetMaxNumSteps(arkode_mem, 100000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;

  /* Call ARKodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Call ARKDense to specify the ARKDENSE dense linear solver */
  flag = ARKDense(arkode_mem, NEQ);
  if (check_flag(&flag, "ARKDense", 1)) return 1;

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);
  if (check_flag(&flag, "ARKDlsSetDenseJacFn", 1)) return 1;

  /* In loop, call ARKode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t = T0;
  realtype tout = T0+dTout;
  realtype u, v, w;
  printf("        t           u           v           w\n");
  printf("   --------------------------------------------------\n");
  printf("  %10.3e  %12.5e  %12.5e  %12.5e\n", 
	 t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
  int iout;
  for (iout=0; iout<Nt; iout++) {

    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
    printf("  %10.3e  %12.5e  %12.5e  %12.5e\n", 
	   t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));

    if (check_flag(&flag, "ARKode", 1)) break;
    if (flag >= 0) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }

  }
  printf("   --------------------------------------------------\n");

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
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);

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
  realtype u = NV_Ith_S(y,0);
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);

  /* du/dt = -0.04*u + 1.e4*v*w */
  NV_Ith_S(ydot,0) = -0.04*u + 1.e4*v*w;

  /* dv/dt = 0.04*u - 1.e4*v*w - 3.e7*v*v */
  NV_Ith_S(ydot,1) = 0.04*u - 1.e4*v*w - 3.e7*v*v;

  /* dw/dt = 3.e7*v*v */
  NV_Ith_S(ydot,2) = 3.e7*v*v;

  return 0;
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);
  SetToZero(J);

  /* du/dt = -0.04*u + 1.e4*v*w */
  DENSE_ELEM(J,0,0) = -0.04;
  DENSE_ELEM(J,0,1) = 1.e4*w;
  DENSE_ELEM(J,0,2) = 1.e4*v;

  /* dv/dt = 0.04*u - 1.e4*v*w - 3.e7*v*v */
  DENSE_ELEM(J,1,0) = 0.04;
  DENSE_ELEM(J,1,1) = -1.e4*w - 6.e7*v;
  DENSE_ELEM(J,1,2) = -1.e4*v;

  /* dw/dt = 3.e7*v*v */
  DENSE_ELEM(J,2,1) = 6.e7*v;

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
