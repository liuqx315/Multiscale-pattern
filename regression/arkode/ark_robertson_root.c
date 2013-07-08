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

 While integrating the system, we use the rootfinding feature 
 to find the times at which either u=1e-4 or w=1e-2.
 
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
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int JacI(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int g(realtype t, N_Vector y, 
	     realtype *gout, void *user_data);


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
  realtype T1 = RCONST(0.4);
  realtype TMult = RCONST(10.0);
  int Nt = 12;
  realtype u0, v0, w0;
  long int NEQ = 3;
  int rootsfound[2];
  long int nst, nst_a, nfe, nfi, nsetups;
  long int nje, nfeLS, nni, ncfn, netf, nge;
  
  /* declare solver parameters */
  int flag, rtflag, dense_order, imex, fixedpt;

  /* general problem variables */
  int idense;
  N_Vector y = NULL;
  N_Vector atols = NULL;
  void *arkode_mem = NULL;

  /* open solver diagnostics output file for writing */
  FILE *DFID;
  DFID=fopen("diags_ark_robertson_root.txt","w");
  
  /* set up the initial conditions */
  u0 = RCONST(1.0);
  v0 = RCONST(0.0);
  w0 = RCONST(0.0);

  /* Initial problem output */
  printf("\nRobertson ODE test problem (with rootfinding):\n");
  printf("    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n",u0,v0,w0);

  /* Create serial vectors of length NEQ for initial condition */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;
  atols = N_VNew_Serial(NEQ);
  if (check_flag((void *) atols, "N_VNew_Serial", 0)) return 1;

  /* Set initial conditions into y */
  NV_Ith_S(y,0) = u0;
  NV_Ith_S(y,1) = v0;
  NV_Ith_S(y,2) = w0;

  /* Call ARKodeCreate to create the solver memory */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return 1;
  
  /* Call init_from_file helper routine to read and set solver parameters */
  realtype rtol, atol;
  flag = init_from_file(arkode_mem, "solve_params.txt", f, fe, fi, T0,
			y, &imex, &dense_order, &fixedpt, &rtol, &atol);
  if (check_flag(&flag, "init_from_file", 1)) return 1;

  /* Update tolerances */
  realtype reltol = RCONST(1.0e-4);
  /* NV_Ith_S(atols,0) = RCONST(1.0e-8); */
  /* NV_Ith_S(atols,1) = RCONST(1.0e-14); */
  /* NV_Ith_S(atols,2) = RCONST(1.0e-6); */
  NV_Ith_S(atols,0) = RCONST(1.0e-8);
  NV_Ith_S(atols,1) = RCONST(1.0e-8);
  NV_Ith_S(atols,2) = RCONST(1.0e-8);

  /* If (dense_order == -1), tell integrator to use tstop */
  if (dense_order == -1) {
    idense = 0;
  } else {    /* otherwise tell integrator to use dense output */
    idense = 1;
  }

  /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
  flag = ARKodeSetDiagnostics(arkode_mem, DFID);
  if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return 1;

  /* Increase maximum number of error test failures */
  flag = ARKodeSetMaxErrTestFails(arkode_mem, 20);
  if (check_flag(&flag, "ARKodeSetMaxErrTestFails", 1)) return 1;

  /* Call ARKodeSetMaxNumSteps to increase default (for testing) */
  flag = ARKodeSetMaxNumSteps(arkode_mem, 100000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;

  /* Call ARKodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = ARKodeSVtolerances(arkode_mem, reltol, atols);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Call ARKodeRootInit to specify the root function with 2 equations */
  flag = ARKodeRootInit(arkode_mem, 2, g);
  if (check_flag(&flag, "ARKodeRootInit", 1)) return 1;

  /* Call ARKDense to specify the ARKDENSE dense linear solver */
  flag = ARKDense(arkode_mem, NEQ);
  if (check_flag(&flag, "ARKDense", 1)) return 1;

  /* Set the Jacobian routine to Jac (user-supplied) */
  switch (imex) {
  case 0:         /* purely implicit */
    flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);   break;
  case 1:         /* purely explicit */
    break;
  default:        /* imex */
    flag = ARKDlsSetDenseJacFn(arkode_mem, JacI);  break;
  }
  if (check_flag(&flag, "ARKDlsSetDenseJacFn", 1)) return 1;

  /* Write all solver parameters to stdout */
  printf("\n");
  flag = ARKodeWriteParameters(arkode_mem, stdout);
  if (check_flag(&flag, "ARKodeWriteParameters", 1)) return 1;


  /* Print out root-finding expectations */
  printf("\n Roots should be found at the following times and values:\n");
  printf("   t=2.64019e-1, u=9.89965e-1, v=3.47049e-5,  w=1.00000e-2  [ 0 1]\n");
  printf("   t=2.07956e+7, u=1.00000e-4, v=3.96207e-10, w=9.99900e-1  [-1 0]\n\n");


  /* In loop, call ARKode, print results, and test for error.
     Break out of loop when Nt preset output times have been reached */
  realtype t = T0;
  printf("        t             u             v             w\n");
  printf("   -----------------------------------------------------\n");
  printf("  %12.5e  %12.5e  %12.5e  %12.5e\n", 
	 t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
  realtype tout = T1;
  int iout=0;
  while(1) {

    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
    printf("  %12.5e  %12.5e  %12.5e  %12.5e\n",  t, 
	   NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));

    if (flag == ARK_ROOT_RETURN) {
      rtflag = ARKodeGetRootInfo(arkode_mem, rootsfound);
      if (check_flag(&flag, "ARKodeGetRootInfo", 1)) return 1;
      printf("      rootsfound[] = %3d %3d\n", 
	     rootsfound[0], rootsfound[1]);
    }

    if (check_flag(&flag, "ARKode", 1)) break;
    if (flag == ARK_SUCCESS) {
      iout++;
      tout *= TMult;
    }

    if (iout == Nt) break;
  }
  printf("   -----------------------------------------------------\n");

  /* Print some final statistics */
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

  flag = ARKodeGetNumGEvals(arkode_mem, &nge);
  check_flag(&flag, "ARKodeGetNumGEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", 
	 nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of nonlinear iterations = %li\n", nni);
  printf("   Total root-function g evals = %li\n", nge);
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

/* fe routine to compute the explicit portion of f(t,y). */
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype u = NV_Ith_S(y,0);

  /* du/dt = -0.04*u */
  NV_Ith_S(ydot,0) = -0.04*u;

  /* dv/dt = 0.04*u */
  NV_Ith_S(ydot,1) = 0.04*u;

  /* dw/dt = 0.0 */
  NV_Ith_S(ydot,2) = 0.0;

  return 0;
}

/* fi routine to compute the implicit portion of f(t,y). */
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);

  /* du/dt = 1.e4*v*w */
  NV_Ith_S(ydot,0) = 1.e4*v*w;

  /* dv/dt = - 1.e4*v*w - 3.e7*v*v */
  NV_Ith_S(ydot,1) = 1.e4*v*w - 3.e7*v*v;

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

/* Jacobian routine to compute J(t,y) = dfi/dy. */
static int JacI(long int N, realtype t,
		N_Vector y, N_Vector fy, DlsMat J, void *user_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);
  SetToZero(J);

  /* du/dt = 1.e4*v*w */
  DENSE_ELEM(J,0,1) = 1.e4*w;
  DENSE_ELEM(J,0,2) = 1.e4*v;

  /* dv/dt = -1.e4*v*w - 3.e7*v*v */
  DENSE_ELEM(J,1,1) = -1.e4*w - 6.e7*v;
  DENSE_ELEM(J,1,2) = -1.e4*v;

  /* dw/dt = 3.e7*v*v */
  DENSE_ELEM(J,2,1) = 6.e7*v;

  return 0;
}

/* g routine to compute the root-finding function g(t,y). */
static int g(realtype t, N_Vector y, realtype *gout, void *user_data)
{
  realtype u = NV_Ith_S(y,0);
  realtype w = NV_Ith_S(y,2);

  /* check for u == 1e-4 */
  gout[0] = u - RCONST(0.0001);

  /* check for w == 1e-2 */
  gout[1] = w - RCONST(0.01);

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
