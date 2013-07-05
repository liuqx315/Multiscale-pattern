/*---------------------------------------------------------------
 $Revision: $
 $Date: $
-----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Example problem:
 
 The following test simulates a simple 1D heat equation,
    u_t = k*u_xx + f
 for t in [0, 10], x in [0, 1], with initial conditions
    u(0,x) =  0
 Dirichlet boundary conditions, i.e. 
    u_t(t,0) = u_t(t,1) = 0,
 and a point-source heating term,
    f = 1 for x=0.5.
 
 The spatial derivatives are computed using second-order 
 centered differences, with the data distributed over N points 
 on a uniform spatial grid.

 The number of spatial points N and the parameter k, as well as 
 the desired relative and absolute solver tolerances, are 
 provided in the input file input_heat1D.txt.
 
 This program solves the problem with either an ERK or DIRK
 method.  For the DIRK method, we use a Newton iteration with 
 the PCG linear solver, and a user-supplied Jacobian-vector 
 product routine.

 100 outputs are printed at equal intervals, and run statistics 
 are printed at the end.
---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode.h>
#include <nvector/nvector_serial.h>
#include <arkode/arkode_pcg.h>
#include <arkode/arkode_spgmr.h>
#include <arkode/arkode_spbcgs.h>
#include <arkode/arkode_sptfqmr.h>
#include <sundials/sundials_types.h>


/* user data structure */
typedef struct {  
  long int N;    /* number of intervals     */
  realtype dx;   /* mesh spacing            */
  realtype k;    /* diffusion coefficient   */
} *UserData;


/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
	       N_Vector fy, void *user_data, N_Vector tmp);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);

/* Parameter input helper function */
int init_from_file(void *ark_mem, char *fname, ARKRhsFn f, 
		   ARKRhsFn fe, ARKRhsFn fi, realtype T0, 
		   N_Vector y0, int *ImEx, int *dorder, 
		   realtype *RTol, realtype *ATol);


/* Main Program */
int main() {

  /* general problem parameters */
  realtype T0 = RCONST(0.0);
  realtype Tf = RCONST(1.0);
  int Nt = 10;
  UserData udata = NULL;
  realtype *data;
  long int N, i;

  /* declare solver parameters */
  int flag, dense_order, imex;

  /* general problem variables */
  int idense;
  N_Vector y = NULL;
  N_Vector ytrue = NULL;
  N_Vector yerr  = NULL;
  void *arkode_mem = NULL;
  void *arktrue_mem = NULL;

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  if (check_flag((void *) udata, "malloc", 2)) return 1;

  /* read problem parameter and tolerances from input file:
     N - number of spatial discretization points
     k - diffusion coefficient */
  double k;
  FILE *FID;
  FID = fopen("input_heat1D.txt","r");
  flag = fscanf(FID,"  N = %li\n", &N);
  flag = fscanf(FID,"  k = %lf\n", &k);
  fclose(FID);

  /* store the inputs in the UserData structure */
  udata->N = N;
  udata->k = k;

  /* open solver diagnostics output file for writing */
  FILE *DFID;
  DFID=fopen("diags_ark_heat1D.txt","w");
  
  /* Initial problem output */
  printf("\n1D Heat PDE test problem:\n");
  printf("  N = %li\n", udata->N);
  printf("  diffusion coefficient:  k = %g\n", udata->k);

  /* Create serial vector of length N for initial condition */
  y = N_VNew_Serial(N);
  if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;
  ytrue = N_VNew_Serial(N);
  if (check_flag((void *) ytrue, "N_VNew_Serial", 0)) return 1;
  yerr = N_VNew_Serial(N);
  if (check_flag((void *) yerr, "N_VNew_Serial", 0)) return 1;

  /* set spatial mesh spacing */
  udata->dx = RCONST(1.0)/(1.0*N-1.0);

  /* output mesh to disk */
  FID=fopen("heat_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->dx*i);
  fclose(FID);
  
  /* Set initial conditions into y, ytrue */
  N_VConst(0.0, y);
  N_VConst(0.0, ytrue);


  /* Call ARKodeCreate to create the solver memory */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *) arkode_mem, "ARKodeCreate", 0)) return 1;
  arktrue_mem = ARKodeCreate();
  if (check_flag((void *) arktrue_mem, "ARKodeCreate", 0)) return 1;

  /* Call init_from_file helper routine to read and set solver parameters */
  realtype rtol, atol;
  flag = init_from_file(arkode_mem, "solve_params.txt", f, NULL, NULL,
			T0, y, &imex, &dense_order, &rtol, &atol);
  if (check_flag(&flag, "init_from_file", 1)) return 1;
  if (rtol <= 0.0)  rtol = 1.e-6;
  if (atol <= 0.0)  atol = 1.e-10;
  realtype reltol = rtol;
  realtype abstol = atol;
  realtype reltol2 = reltol*1.0e-2;
  realtype abstol2 = abstol*1.0e-2;

  /* If (dense_order == -1), tell integrator to use tstop */
  if (dense_order == -1) 
    idense = 0;
  else          /* otherwise tell integrator to use dense output */
    idense = 1;
  
  /* Reference solution will be computed with default implicit method */
  flag = ARKodeInit(arktrue_mem, NULL, f, T0, ytrue);
  if (check_flag(&flag, "ARKodeInit", 1)) return 1;
  
  /* Call ARKodeSetUserData to pass rdata to user functions */
  flag = ARKodeSetUserData(arkode_mem, (void *) udata);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;
  flag = ARKodeSetUserData(arktrue_mem, (void *) udata);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;

  /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
  flag = ARKodeSetDiagnostics(arkode_mem, DFID);
  if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return 1;

  /* Call ARKodeSetMaxNumSteps to increase default (for testing) */
  flag = ARKodeSetMaxNumSteps(arkode_mem, 10000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;
  flag = ARKodeSetMaxNumSteps(arktrue_mem, 10000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;

  /* Call ARKodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;
  flag = ARKodeSStolerances(arktrue_mem, reltol2, abstol2);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Specify the linear solver */
  flag = ARKPcg(arkode_mem, 0, N);
  if (check_flag(&flag, "ARKPcg", 1)) return 1;
  flag = ARKPcg(arktrue_mem, 0, N);
  if (check_flag(&flag, "ARKPcg", 1)) return 1;

  /* Set the Jacobian routine to Jac (user-supplied) */
  if (imex == 0) { /* purely implicit */
    flag = ARKSpilsSetJacTimesVecFn(arkode_mem, Jac);
    if (check_flag(&flag, "ARKSpilsSetJacTimesVecFn", 1)) return 1;
  }
  flag = ARKSpilsSetJacTimesVecFn(arktrue_mem, Jac);
  if (check_flag(&flag, "ARKSpilsSetJacTimesVecFn", 1)) return 1;

  /* Open output stream for results, access data arrays */
  FILE *UFID=fopen("heat1D.txt","w");
  data = N_VGetArrayPointer(y);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;

  /* output initial condition to disk */
  for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[i]);
  fprintf(UFID,"\n");

  /* Write all solver parameters to stdout */
  printf("\n");
  flag = ARKodeWriteParameters(arkode_mem, stdout);
  if (check_flag(&flag, "ARKodeWriteParameters", 1)) return 1;

  /* In loop, call ARKode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t  = T0;
  realtype t2 = T0;
  realtype dTout = Tf/Nt;
  realtype tout = T0+dTout;
  realtype u, uerr, errI=0.0, err2=0.0;
  printf("        t      ||u||_rms    ||uerr||\n");
  printf("   ------------------------------------\n");
  int iout;
  for (iout=0; iout<Nt; iout++) {

    if (!idense)
      flag = ARKodeSetStopTime(arktrue_mem, tout);
    flag = ARKode(arktrue_mem, tout, ytrue, &t2, ARK_NORMAL);
    if (!idense)
      flag = ARKodeSetStopTime(arkode_mem, tout);
    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_flag(&flag, "ARKode", 1)) break;
    if (flag >= 0) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
    u = N_VDotProd(y,y);
    u = sqrt(u/N);

    N_VLinearSum( 1.0, ytrue, -1.0, y, yerr );
    uerr = N_VDotProd(yerr,yerr);
    errI = (errI > N_VMaxNorm(yerr)) ? errI : N_VMaxNorm(yerr);
    err2 += uerr;
    printf("  %10.6f  %10.6f  %12.5e\n", t, u, sqrt(uerr/N));

    /* output results to disk */
    for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[i]);
    fprintf(UFID,"\n");
  }
  err2 = sqrt(err2 / N / Nt);
  printf("   ------------------------------------\n");
  fclose(UFID);
    

  /* Print some final statistics */
  long int nst, nst_a, nfe, nfi, nsetups, nli, nJv, nlcf, nni, ncfn, netf;
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
  flag = ARKSpilsGetNumLinIters(arkode_mem, &nli);
  check_flag(&flag, "ARKSpilsGetNumLinIters", 1);
  flag = ARKSpilsGetNumJtimesEvals(arkode_mem, &nJv);
  check_flag(&flag, "ARKSpilsGetNumJtimesEvals", 1);
  flag = ARKSpilsGetNumConvFails(arkode_mem, &nlcf);
  check_flag(&flag, "ARKSpilsGetNumConvFails", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", 
	 nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total linear iterations = %li\n", nli);
  printf("   Total number of Jacobian-vector products = %li\n", nJv);
  printf("   Total number of linear solver convergence failures = %li\n", nlcf);
  printf("   Total number of nonlinear iterations = %li\n", nni);
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Error: max = %g, rms = %g\n", errI, err2);
  printf("   Oversolve = %g\n\n", reltol/(err2+1.e-10*reltol));

  /* Free vectors */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(ytrue);
  N_VDestroy_Serial(yerr);

  /* Free user data */
  free(udata);

  /* Free integrator memory */
  ARKodeFree(&arkode_mem);
  ARKodeFree(&arktrue_mem);

  /* close solver diagnostics output file */
  fclose(DFID);

  return 0;
}


/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  /* clear out ydot (to be careful) */
  N_VConst(0.0, ydot);

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of intervals, background values */
  long int N  = udata->N;
  realtype k  = udata->k;
  realtype dx = udata->dx;

  /* access data arrays */
  realtype *Y = N_VGetArrayPointer(y);
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return 1;
  realtype *Ydot = N_VGetArrayPointer(ydot);
  if (check_flag((void *) Ydot, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over domain, computing all equations */
  realtype c1 = k/dx/dx;
  realtype c2 = -RCONST(2.0)*k/dx/dx;
  long int i;
  long int isource = N/2;
  Ydot[0] = 0.0;                 /* left boundary condition */
  for (i=1; i<N-1; i++)
    Ydot[i] = c1*Y[i-1] + c2*Y[i] + c1*Y[i+1];
  Ydot[N-1] = 0.0;               /* right boundary condition */
  Ydot[isource] += 1.0;          /* f */

  return 0;
}



/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
	       N_Vector fy, void *user_data, N_Vector tmp)
{
  /* clear out result (to be careful) */
  N_VConst(0.0, Jv);

  /* shortcuts to number of intervals, background values */
  UserData udata = (UserData) user_data;
  long int N  = udata->N;
  realtype k  = udata->k;
  realtype dx = udata->dx;

  /* access data arrays */
  realtype *V = N_VGetArrayPointer(v);
  if (check_flag((void *) V, "N_VGetArrayPointer", 0)) return 1;
  realtype *JV = N_VGetArrayPointer(Jv);
  if (check_flag((void *) JV, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over domain, computing all Jacobian-vector products */
  realtype c1 = k/dx/dx;
  realtype c2 = -RCONST(2.0)*k/dx/dx;
  long int i;
  JV[0] = 0.0;
  for (i=1; i<N-1; i++)
    JV[i] = c1*V[i-1] + c2*V[i] + c1*V[i+1];
  JV[N-1] = 0.0;

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
