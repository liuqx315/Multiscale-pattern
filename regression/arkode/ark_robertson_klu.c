/*---------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2013, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 ----------------------------------------------------------------
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
 
 This program solves the problem with one of the solvers DIRK or 
 ARK.  For DIRK and ARK, implicit subsystems are solved using a 
 Newton iteration with the ARKKLU sparse linear solver.

 100 outputs are printed at equal intervals, and run statistics 
 are printed at the end.
---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <math.h>
#include <arkode/arkode.h>
#include <nvector/nvector_serial.h>
#include <arkode/arkode_klu.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_types.h>

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(realtype t, N_Vector y, N_Vector fy, SlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int JacI(realtype t, N_Vector y, N_Vector fy, SlsMat J, void *user_data,
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
  realtype Tf = RCONST(1.e11);
  realtype dTout = (Tf-T0)/100;
  int Nt = ceil(Tf/dTout);
  realtype u0, v0, w0, h0;
  long int nnz, NEQ = 3;

  /* declare solver parameters */
  int flag, dense_order, imex, fixedpt;

  /* general problem variables */
  int idense;
  N_Vector y = NULL;
  N_Vector ytrue = NULL;
  void *arkode_mem = NULL;
  void *arktrue_mem = NULL;

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
  ytrue = N_VNew_Serial(NEQ);
  if (check_flag((void *) ytrue, "N_VNew_Serial", 0)) return 1;

  /* Set initial conditions into y */
  NV_Ith_S(y,0) = u0;
  NV_Ith_S(y,1) = v0;
  NV_Ith_S(y,2) = w0;
  NV_Ith_S(ytrue,0) = u0;
  NV_Ith_S(ytrue,1) = v0;
  NV_Ith_S(ytrue,2) = w0;

  /* Call ARKodeCreate to create the solver memory */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return 1;
  arktrue_mem = ARKodeCreate();
  if (check_flag((void *)arktrue_mem, "ARKodeCreate", 0)) return 1;
  
  /* Call init_from_file helper routine to read and set solver parameters */
  realtype rtol, atol;
  flag = init_from_file(arkode_mem, "solve_params.txt", f, fe, fi, T0, 
			y, &imex, &dense_order, &fixedpt, &rtol, &atol);
  if (check_flag(&flag, "init_from_file", 1)) return 1;
  rtol = 1.e-4;      /* Update tolerances */
  atol = 1.e-8;
  realtype reltol = rtol;
  realtype abstol = atol;
  realtype reltol2 = reltol*1.0e-3;
  realtype abstol2 = abstol*1.0e-3;
  h0 = 1.e-4 * reltol;

  /* If (dense_order == -1), tell integrator to use tstop */
  if (dense_order == -1) {
    idense = 0;
  } else {    /* otherwise tell integrator to use dense output */
    idense = 1;
  }

  /* Reference solution uses default implicit method */
  flag = ARKodeInit(arktrue_mem, NULL, f, T0, ytrue);
  if (check_flag(&flag, "ARKodeInit", 1)) return 1;

  /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
  flag = ARKodeSetDiagnostics(arkode_mem, DFID);
  if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return 1;

  /* Set custom initial step */
  flag = ARKodeSetInitStep(arkode_mem, h0);
  if (check_flag(&flag, "ARKodeSetInitStep", 1)) return 1;

  /* Increase maximum number of error test failures */
  flag = ARKodeSetMaxErrTestFails(arkode_mem, 20);
  if (check_flag(&flag, "ARKodeSetMaxErrTestFails", 1)) return 1;
  flag = ARKodeSetMaxErrTestFails(arktrue_mem, 30);
  if (check_flag(&flag, "ARKodeSetMaxErrTestFails", 1)) return 1;

  /* Call ARKodeSetMaxNumSteps to increase default (for testing) */
  flag = ARKodeSetMaxNumSteps(arkode_mem, 100000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;
  flag = ARKodeSetMaxNumSteps(arktrue_mem, 100000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;

  /* Tighten inner solver tolerances for this problem */
  flag = ARKodeSetNonlinConvCoef(arkode_mem, 1.e-3);
  if (check_flag(&flag, "ARKodeSetNonlinConvCoef", 1)) return 1;
  flag = ARKodeSetNonlinConvCoef(arktrue_mem, 1.e-3);
  if (check_flag(&flag, "ARKodeSetNonlinConvCoef", 1)) return 1;
  flag = ARKodeSetMaxNonlinIters(arkode_mem, 8);
  if (check_flag(&flag, "ARKodeSetMaxNonlinIters", 1)) return 1;
  flag = ARKodeSetMaxNonlinIters(arktrue_mem, 8);
  if (check_flag(&flag, "ARKodeSetMaxNonlinIters", 1)) return 1;

  /* Call ARKodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;
  flag = ARKodeSStolerances(arktrue_mem, reltol2, abstol2);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Call ARKKLU to specify the ARKKLU sparse linear solver */
  nnz = 7;
  flag = ARKKLU(arkode_mem, NEQ, nnz);
  if (check_flag(&flag, "ARKKLU", 1)) return 1;
  flag = ARKKLU(arktrue_mem, NEQ, nnz);
  if (check_flag(&flag, "ARKKLU", 1)) return 1;

  /* Set the Jacobian routine to Jac (user-supplied) */
  switch (imex) {
  case 0:         /* purely implicit */
    flag = ARKSlsSetSparseJacFn(arkode_mem, Jac);   break;
  case 1:         /* purely explicit */
    break;
  default:        /* imex */
    flag = ARKSlsSetSparseJacFn(arkode_mem, JacI);  break;
  }
  if (check_flag(&flag, "ARKSlsSetSparseJacFn", 1)) return 1;
  flag = ARKSlsSetSparseJacFn(arktrue_mem, Jac);
  if (check_flag(&flag, "ARKDlsSetSparseJacFn", 1)) return 1;

  /* Write all solver parameters to stdout */
  printf("\n");
  flag = ARKodeWriteParameters(arkode_mem, stdout);
  if (check_flag(&flag, "ARKodeWriteParameters", 1)) return 1;

  /* In loop, call ARKode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t = T0;
  realtype t2 = T0;
  realtype tout = T0+dTout;
  realtype u, v, w, uerr, verr, werr, errI=0.0, err2=0.0;
  printf("        t           u           v           w        uerr          verr          werr\n");
  printf("   ----------------------------------------------------------------------------------------------\n");
  printf("  %10.3e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n", 
	 t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2), 0.0, 0.0, 0.0);
  int iout;
  for (iout=0; iout<Nt; iout++) {
    if (!idense)
      flag = ARKodeSetStopTime(arktrue_mem, tout);
    flag = ARKode(arktrue_mem, tout, ytrue, &t2, ARK_NORMAL);
    if (!idense)
      flag = ARKodeSetStopTime(arkode_mem, tout);
    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (flag >= 0) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {
      fprintf(stderr,"Solver failure, stopping integration\n");
      return 1;
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
    printf("  %10.3e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n", 
	   t, u, v, w, uerr, verr, werr);
  }
  err2 = sqrt(err2 / 3.0 / Nt);
  printf("   ----------------------------------------------------------------------------------------------\n");

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
  flag = ARKSlsGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKSlsGetNumJacEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", 
	 nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of nonlinear iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Error: max = %g, rms = %g\n", errI, err2);
  printf("   Oversolve = %g\n\n", reltol/(err2+1.e-10*reltol));

  /* Free y vector */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(ytrue);

  /* Free integrator memory */
  ARKodeFree(&arkode_mem);
  ARKodeFree(&arktrue_mem);

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
  NV_Ith_S(ydot,0) = -RCONST(0.04)*u + RCONST(1.e4)*v*w;

  /* dv/dt = 0.04*u - 1.e4*v*w - 3.e7*v*v */
  NV_Ith_S(ydot,1) = RCONST(0.04)*u - RCONST(1.e4)*v*w - RCONST(3.e7)*v*v;

  /* dw/dt = 3.e7*v*v */
  NV_Ith_S(ydot,2) = RCONST(3.e7)*v*v;

  return 0;
}

/* fe routine to compute the explicit portion of f(t,y). */
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype u = NV_Ith_S(y,0);

  /* du/dt = -0.04*u */
  NV_Ith_S(ydot,0) = -RCONST(0.04)*u;

  /* dv/dt = 0.04*u */
  NV_Ith_S(ydot,1) = RCONST(0.04)*u;

  /* dw/dt = 0.0 */
  NV_Ith_S(ydot,2) = RCONST(0.0);

  return 0;
}

/* fi routine to compute the implicit portion of f(t,y). */
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);

  /* du/dt = 1.e4*v*w */
  NV_Ith_S(ydot,0) = RCONST(1.e4)*v*w;

  /* dv/dt = - 1.e4*v*w - 3.e7*v*v */
  NV_Ith_S(ydot,1) = RCONST(1.e4)*v*w - RCONST(3.e7)*v*v;

  /* dw/dt = 3.e7*v*v */
  NV_Ith_S(ydot,2) = RCONST(3.e7)*v*v;

  return 0;
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(realtype t, N_Vector y, N_Vector fy, SlsMat J, 
	       void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int nz=0;
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);
  SlsSetToZero(J);

  /* first column:
     d(du/dt)/du = -0.04
     d(dv/dt)/du = 0.04
     d(dw/dt)/du = 0 */
  J->colptrs[nz] = nz;
  J->data[nz] = RCONST(-0.04);
  J->rowvals[nz++] = 0;
  J->data[nz] = RCONST(0.04);
  J->rowvals[nz++] = 1;

  /* second column:
     d(du/dt)/dv = 1.e4*w
     d(dv/dt)/dv = -1.e4*w - 6.e7*v
     d(dw/dt)/dv = 6.e7*v */
  J->colptrs[1] = nz;
  J->data[nz] = RCONST(1.0e4)*w;
  J->rowvals[nz++] = 0;
  J->data[nz] = RCONST(-1.0e4)*w - RCONST(6.0e7)*v;
  J->rowvals[nz++] = 1;
  J->data[nz] = RCONST(6.0e7)*v;
  J->rowvals[nz++] = 2;

  /* third column:
     d(du/dt)/dw = 1.e4*v
     d(dv/dt)/dw = -1.e4*v
     d(dw/dt)/dw = 0 */
  J->colptrs[2] = nz;
  J->data[nz] = RCONST(1.0e4)*v;
  J->rowvals[nz++] = 0;
  J->data[nz] = RCONST(-1.0e4)*v;
  J->rowvals[nz++] = 1;
  
  /* end of data */
  J->colptrs[3] = nz;

  return 0;
}

/* Jacobian routine to compute J(t,y) = dfi/dy. */
static int JacI(realtype t, N_Vector y, N_Vector fy, SlsMat J, 
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int nz=0;
  realtype v = NV_Ith_S(y,1);
  realtype w = NV_Ith_S(y,2);
  SlsSetToZero(J);

  /* first column:
     d(du/dt)/du = 0
     d(dv/dt)/du = 0
     d(dw/dt)/du = 0 */
  J->colptrs[nz] = nz;

  /* second column:
     d(du/dt)/dv = 1.e4*w
     d(dv/dt)/dv = -1.e4*w - 6.e7*v
     d(dw/dt)/dv = 6.e7*v */
  J->colptrs[1] = nz;
  J->data[nz] = RCONST(1.0e4)*w;
  J->rowvals[nz++] = 0;
  J->data[nz] = RCONST(-1.0e4)*w - RCONST(6.0e7)*v;
  J->rowvals[nz++] = 1;
  J->data[nz] = RCONST(6.0e7)*v;
  J->rowvals[nz++] = 2;

  /* third column:
     d(du/dt)/dw = 1.e4*v
     d(dv/dt)/dw = -1.e4*v
     d(dw/dt)/dw = 0 */
  J->colptrs[2] = nz;
  J->data[nz] = RCONST(1.0e4)*v;
  J->rowvals[nz++] = 0;
  J->data[nz] = RCONST(-1.0e4)*v;
  J->rowvals[nz++] = 1;
  
  /* end of data */
  J->colptrs[3] = nz;

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
