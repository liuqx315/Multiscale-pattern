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
 * This program solves the problem with the BDF method,
 * Newton iteration with the CVDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * Output is printed every 1.0 units of time (10 total).
 * Run statistics (optional outputs) are printed at the end.
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
  realtype T0 = 0.0;
  realtype Tf = 10.0;
  realtype dTout = 1.0;
  long int NEQ = 1;

  /* general problem variables */
  int flag;
  N_Vector y = NULL;
  void *cvode_mem = NULL;

  /* read problem parameter and tolerances from input file:
     lamda  - problem stiffness parameter
     reltol - desired relative tolerance
     abstol - desired absolute tolerance */
  double reltol_, abstol_, lamda_;
  FILE *FID;
  FID=fopen("input_analytic.txt","r");
  fscanf(FID,"  lamda = %lf\n", &lamda_);
  fscanf(FID,"  reltol = %lf\n", &reltol_);
  fscanf(FID,"  abstol = %lf\n", &abstol_);
  fclose(FID);

  /* convert the inputs to 'realtype' format */
  realtype reltol = reltol_;
  realtype abstol = abstol_;
  realtype lamda  = lamda_;

  /* Initial problem output */
  printf("\nAnalytical ODE test problem\n");
  printf("    lamda = %g\n",lamda);
  printf("   reltol = %.1e\n",reltol);
  printf("   abstol = %.1e\n\n",abstol);


  /* Create serial vector of length NEQ for initial condition */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  /* Initialize y to 0 */
  NV_Ith_S(y,0) = 0.0;

  /* Call CVodeCreate to create the solver memory and specify the 
     Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
  
  /* Call CVodeInit to initialize the integrator memory and specify the
     user's right hand side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y */
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  /* Call CVodeSetUserData to pass lamda to user functions */
  flag = CVodeSetUserData(cvode_mem, (void *) &lamda);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)) return(1);

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t = T0;
  realtype tout = dTout;
  realtype u;
  while (Tf - t > 1.0e-15) {
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    u = NV_Ith_S(y,0);
    printf("  t = %4.1f,   y = %.6f,   error = %.5e\n", t, u, atan(t)-u);

    if (check_flag(&flag, "CVode", 1)) break;
    if (flag == CV_SUCCESS) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
  }

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
  printf("   Total number of error test failures = %li\n\n", netf);

  /* Free y vector */
  N_VDestroy_Serial(y);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

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
