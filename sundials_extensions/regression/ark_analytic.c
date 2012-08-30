/* -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following is a simple example problem with analytical 
 * solution,
 *    dy/dt = lamda*y + 1/(1+t^2) - lamda*atan(t)
 * for t in the interval [0.0, 10.0], with initial condition: y=0. 
 * 
 * The stiffness of the problem is directly proportional to the 
 * value of "lamda", which is specified through an input file, along 
 * with the desired relative and absolute tolerances.  The value of
 * lamda should be negative to result in a well-posed ODE; for values
 * with magnitude larger than 100 the problem becomes quite stiff.
 *
 * In the example input file, we choose lamda = -100.
 * 
 * This program solves the problem with the DIRK method,
 * Newton iteration with the ARKDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * Output is printed every 1.0 units of time (10 total).
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>

/* Header files with a description of contents used */

#include <arkode/arkode.h>             /* prototypes for ARKODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <arkode/arkode_dense.h>       /* prototype for ARKDense */
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
  long int NEQ = 1;

  /* declare solver parameters */
  int flag, order, dense_order, btable, adapt_method, small_nef, 
    msbp, maxcor, predictor;
  flag = order = adapt_method = small_nef = msbp = maxcor = predictor = 0;
  dense_order = btable = -1;
  double cflfac, safety, bias, growth, hfixed_lb, hfixed_ub, k1, 
    k2, k3, etamx1, etamxf, etacf, crdown, rdiv, dgmax, nlscoef;
  cflfac = safety = bias = growth = hfixed_lb = hfixed_ub = k1 = k2 = k3
    = etamx1 = etamxf = etacf = crdown = rdiv = dgmax = nlscoef = 0.0;

  /* general problem variables */
  int idense;
  N_Vector y = NULL;
  void *arkode_mem = NULL;

  /* read problem parameter and tolerances from input file:
     lamda  - problem stiffness parameter
     reltol - desired relative tolerance
     abstol - desired absolute tolerance */
  double reltol_, abstol_, lamda_;
  FILE *FID;
  FID=fopen("input_analytic.txt","r");
  fscanf(FID,"  lamda = %lf\n",  &lamda_);
  fscanf(FID,"  reltol = %lf\n", &reltol_);
  fscanf(FID,"  abstol = %lf\n", &abstol_);
  fclose(FID);

  /* convert the inputs to 'realtype' format */
  realtype reltol = reltol_;
  realtype abstol = abstol_;
  realtype lamda  = lamda_;

  /* read solver parameters from file */
  FID=fopen("solve_params.txt","r");
  fscanf(FID,"order = %i\n",  &order);
  fscanf(FID,"dense_order = %i\n", &dense_order);
  fscanf(FID,"btable = %i\n",  &btable);
  fscanf(FID,"adapt_method = %i\n", &adapt_method);
  fscanf(FID,"cflfac = %lf\n", &cflfac);
  fscanf(FID,"safety = %lf\n", &safety);
  fscanf(FID,"bias = %lf\n", &bias);
  fscanf(FID,"growth = %lf\n", &growth);
  fscanf(FID,"hfixed_lb = %lf\n", &hfixed_lb);
  fscanf(FID,"hfixed_ub = %lf\n", &hfixed_ub);
  fscanf(FID,"k1 = %lf\n", &k1);
  fscanf(FID,"k2 = %lf\n", &k2);
  fscanf(FID,"k3 = %lf\n", &k3);
  fscanf(FID,"etamx1 = %lf\n", &etamx1);
  fscanf(FID,"etamxf = %lf\n", &etamxf);
  fscanf(FID,"etacf = %lf\n", &etacf);
  fscanf(FID,"small_nef = %i\n", &small_nef);
  fscanf(FID,"crdown = %lf\n", &crdown);
  fscanf(FID,"rdiv = %lf\n", &rdiv);
  fscanf(FID,"dgmax = %lf\n", &dgmax);
  fscanf(FID,"predictor = %i\n", &predictor);
  fscanf(FID,"msbp = %i\n", &msbp);
  fscanf(FID,"maxcor = %i\n", &maxcor);
  fscanf(FID,"nlscoef = %lf\n", &nlscoef);
  fclose(FID);

  realtype adapt_params[] = {cflfac, safety, bias, growth, 
			     hfixed_lb, hfixed_ub, k1, k2, k3};

  /* open solver diagnostics output file for writing */
  FILE *DFID;
  DFID=fopen("diags_ark_analytic.txt","w");
  
  /* Initial problem output */
  printf("\nAnalytical ODE test problem:\n");
  printf("    lamda = %g\n",    lamda);
  printf("   reltol = %.1e\n",  reltol);
  printf("   abstol = %.1e\n\n",abstol);


  /* Create serial vector of length NEQ for initial condition */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  /* Initialize y to 0 */
  NV_Ith_S(y,0) = 0.0;

  /* Call ARKodeCreate to create the solver memory and specify the 
     Backward Differentiation Formula and the use of a Newton iteration */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return(1);
  
  /* Call ARKodeInit to initialize the integrator memory and specify the
     user's right hand side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y */
  flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return(1);

  /* Call ARKodeSetUserData to pass lamda to user functions */
  flag = ARKodeSetUserData(arkode_mem, (void *) &lamda);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return(1);

  /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
  flag = ARKodeSetDiagnostics(arkode_mem, DFID);
  if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return(1);

  /* Call ARKodeSet routines to insert solver parameters */
  if (order != 0) {     /* order overrides btable */
    printf("  Setting order = %i\n",order);
    flag = ARKodeSetOrder(arkode_mem, order);
    if (flag != 0) {
      fprintf(stderr,"Error in ARKodeSetOrder = %i\n",flag);
      return(1);
    }
  } else if (btable != -1) {
    printf("  Setting IRK Table number = %i\n",btable);
    flag = ARKodeSetIRKTableNum(arkode_mem, btable);
    if (flag != 0) {
      fprintf(stderr,"Error in ARKodeSetIRKTableNum = %i\n",flag);
      return(1);
    }
  }
  printf("  Setting dense order = %i\n",dense_order);
  flag = ARKodeSetDenseOrder(arkode_mem, dense_order);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetDenseOrder = %i\n",flag);
    return(1);
  }
  printf("  Setting adaptivity method = %i\n",adapt_method);
  printf("  Setting adaptivity params = %g %g %g %g %g %g %g %g %g\n",
	 adapt_params[0], adapt_params[1], adapt_params[2], 
	 adapt_params[3], adapt_params[4], adapt_params[5], 
	 adapt_params[6], adapt_params[7], adapt_params[8]);
  flag = ARKodeSetAdaptivityMethod(arkode_mem, adapt_method, adapt_params);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetAdaptMethod = %i\n",flag);
    return(1);
  }
  printf("  Setting adaptivity constants = %g %g %g %i\n",
	 etamx1, etamxf, etacf, small_nef);
  flag = ARKodeSetAdaptivityConstants(arkode_mem, etamx1, etamxf, etacf, small_nef);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetAdaptConstants = %i\n",flag);
    return(1);
  }
  printf("  Setting Newton constants = %g %g\n", crdown, rdiv);
  flag = ARKodeSetNewtonConstants(arkode_mem, crdown, rdiv);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetNewtonConstants = %i\n",flag);
    return(1);
  }
  printf("  Setting LSetup constants = %g %i\n", dgmax, msbp);
  flag = ARKodeSetLSetupConstants(arkode_mem, dgmax, msbp);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetLSetupConstants = %i\n",flag);
    return(1);
  }
  printf("  Setting predictor method = %i\n", predictor);
  flag = ARKodeSetPredictorMethod(arkode_mem, predictor);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetPredictorMethod = %i\n",flag);
    return(1);
  }
  printf("  Setting max Newton iters = %i\n", maxcor);
  flag = ARKodeSetMaxNonlinIters(arkode_mem, maxcor);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetMaxNonlinIters = %i\n",flag);
    return(1);
  }
  printf("  Setting nonlinear solver coefficient = %g\n", nlscoef);
  flag = ARKodeSetNonlinConvCoef(arkode_mem, nlscoef);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetMaxNonlinIters = %i\n",flag);
    return(1);
  }

  /* If (dense_order == -1), tell integrator to use tstop */
  if (dense_order == -1) {
    idense = 0;
  } else {    /* otherwise tell integrator to use dense output */
    idense = 1;
  }


  /* Call ARKodeSetMaxNumSteps to increase default (for testing) */
  flag = ARKodeSetMaxNumSteps(arkode_mem, 10000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return(1);

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

  /* Write all solver parameters to stdout */
  flag = ARKodeWriteParameters(arkode_mem, stdout);
  if (check_flag(&flag, "ARKodeWriteParameters", 1)) return(1);

  /* In loop, call ARKode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t = T0;
  realtype tout = dTout;
  realtype u, uerr, errI=0.0, err2=0.0;
  int Nt=0;
  printf("         t             u          error\n");
  printf("   ----------------------------------------\n");
  while (Tf - t > 1.0e-8) {
    if (!idense) 
      flag = ARKodeSetStopTime(arkode_mem, tout);
    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_flag(&flag, "ARKode", 1)) break;
    if (flag >= 0) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
    u = NV_Ith_S(y,0);
    uerr = fabs(u - atan(t));
    errI = (errI > uerr) ? errI : uerr;
    err2 += uerr*uerr;
    Nt++;
    printf("  %12.8f  %12.8f  %12.5e\n", t, u, uerr);
  }
  err2 = sqrt(err2 / Nt);
  printf("   ----------------------------------------\n");


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
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Error: max = %g, rms = %g\n", errI, err2);
  printf("   Oversolve = %g\n\n", reltol/err2);

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
  realtype lamda = rdata[0];
  realtype u = NV_Ith_S(y,0);

  NV_Ith_S(ydot,0) = lamda*u + 1.0/(1.0+t*t) - lamda*atan(t);

  return(0);
}

/* Jacobian routine to compute J(t,y) = df/dy. */

static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *rdata = (realtype *) user_data;
  realtype lamda = rdata[0];
  DENSE_ELEM(J,0,0) = lamda;

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
