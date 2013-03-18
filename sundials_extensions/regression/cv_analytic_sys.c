/*---------------------------------------------------------------
 $Revision: $
 $Date: $
-----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Example problem:
 
 The following is a simple example problem with analytical 
 solution,
    dy/dt = A*y
 where A = V*D*Vinv, 
      V = [1 -1 1; -1 2 1; 0 -1 2];
      Vinv = 0.25*[5 1 -3; 2 2 -2; 1 1 1];
      D = [-0.5 0 0; 0 -0.1 0; 0 0 lam];
 where lam is a large negative number. The analytical solution to
 this problem is 
   Y(t) = V*exp(D*t)*Vinv*Y0
 for t in the interval [0.0, 0.05], with initial condition: 
 y(0) = [1,1,1]'.
 
 The stiffness of the problem is directly proportional to the 
 value of "lamda", which is specified through an input file.  The 
 value of lamda should be negative to result in a well-posed ODE; 
 for values with magnitude larger than 100 the problem becomes 
 quite stiff.

 In the example input file, we choose lamda = -100.
 
 This program solves the problem with the BDF method,
 Newton iteration with the CVDENSE dense linear solver, and a
 user-supplied Jacobian routine.
 Output is printed every 1.0 units of time (10 total).
 Run statistics (optional outputs) are printed at the end.
---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <math.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);
static int sol(realtype t, realtype lam, N_Vector y);
static int dense_MM(DlsMat A, DlsMat B, DlsMat C);


/* Main Program */
int main() {

  /* general problem parameters */
  realtype T0 = RCONST(0.0);
  realtype Tf = RCONST(0.05);
  realtype dTout = RCONST(0.005);
  long int NEQ = 3;

  /* general problem variables */
  int flag;
  N_Vector y = NULL;
  N_Vector ytrue = NULL;
  void *cvode_mem = NULL;

  /* read problem parameter and tolerances from input file:
     lamda  - problem stiffness parameter */
  double lamda_;
  FILE *FID;
  FID=fopen("input_analytic_sys.txt","r");
  flag = fscanf(FID,"  lamda = %lf\n", &lamda_);
  fclose(FID);

  /* set the tolerances, and convert the input to 'realtype' format */
  realtype reltol = 1.0e-6;
  realtype abstol = 1.0e-10;
  realtype lamda  = lamda_;

  /* Initial problem output */
  printf("\nAnalytical ODE test problem:\n");
  printf("    lamda = %g\n",lamda);
  printf("   reltol = %.1e\n",reltol);
  printf("   abstol = %.1e\n\n",abstol);


  /* Create serial vector of length NEQ for initial condition */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  ytrue = N_VNew_Serial(NEQ);
  if (check_flag((void *)ytrue, "N_VNew_Serial", 0)) return 1;

  /* Initialize y to 0 */
  NV_Ith_S(y,0) = 1.0;
  NV_Ith_S(y,1) = 1.0;
  NV_Ith_S(y,2) = 1.0;

  /* Call CVodeCreate to create the solver memory and specify the 
     Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return 1;
  
  /* Call CVodeInit to initialize the integrator memory and specify the
     user's right hand side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y */
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)) return 1;

  /* Call CVodeSetUserData to pass lamda to user functions */
  flag = CVodeSetUserData(cvode_mem, (void *) &lamda);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return 1;

  /* Call CVodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSStolerances", 1)) return 1;

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)) return 1;

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return 1;

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t = T0;
  realtype tout = dTout;
  realtype y0, y1, y2, yt0, yt1, yt2;
  realtype y0err, y1err, y2err, errI=0.0, err2=0.0;
  int Nt=0;
  printf("      t        y0        y1        y2        err0        err1        err2\n");
  printf("   --------------------------------------------------------------------------\n");
  while (Tf - t > 1.0e-15) {

    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    y0 = NV_Ith_S(y,0);
    y1 = NV_Ith_S(y,1);
    y2 = NV_Ith_S(y,2);

    flag = sol(t, lamda, ytrue);
    yt0 = NV_Ith_S(ytrue,0);
    yt1 = NV_Ith_S(ytrue,1);
    yt2 = NV_Ith_S(ytrue,2);

    y0err = fabs(y0-yt0);
    y1err = fabs(y1-yt1);
    y2err = fabs(y2-yt2);

    errI = (errI > y0err) ? errI : y0err;
    errI = (errI > y1err) ? errI : y1err;
    errI = (errI > y2err) ? errI : y2err;
    err2 += y0err*y0err + y1err*y1err + y2err*y2err;
    Nt++;

    printf("  %8.4f  %8.5f  %8.5f  %8.5f  %10.3e  %10.3e  %10.3e\n", 
	   t, y0, y1, y2, y0err, y1err, y2err);

    if (check_flag(&flag, "CVode", 1)) break;
    if (flag == CV_SUCCESS) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
  }
  err2 = sqrt(err2 / 3.0 / Nt);
  printf("   --------------------------------------------------------------------------\n");

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
  printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Error: max = %g, rms = %g\n", errI, err2);
  printf("   Oversolve = %g\n\n", reltol/err2);

  /* Free y vector */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(ytrue);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  return 0;
}


/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  realtype *rdata = (realtype *) user_data;
  realtype lam = rdata[0];
  realtype y0 = NV_Ith_S(y,0);
  realtype y1 = NV_Ith_S(y,1);
  realtype y2 = NV_Ith_S(y,2);
  realtype yd0, yd1, yd2;
  
  /* f(t,y) = V*D*Vi*y, where 
        V = [1 -1 1; -1 2 1; 0 -1 2] 
        Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1]
        D = [-0.5 0 0; 0 -0.1 0; 0 0 lam] */

  /*   yd = Vi*y */
  yd0 = 0.25*(5.0*y0 + 1.0*y1 - 3.0*y2);
  yd1 = 0.25*(2.0*y0 + 2.0*y1 - 2.0*y2);
  yd2 = 0.25*(1.0*y0 + 1.0*y1 + 1.0*y2);

  /*   y = D*yd */
  y0  = -0.5*yd0;
  y1  = -0.1*yd1;
  y2  =  lam*yd2;

  /*   yd = V*y */
  yd0 =  1.0*y0 - 1.0*y1 + 1.0*y2;
  yd1 = -1.0*y0 + 2.0*y1 + 1.0*y2;
  yd2 =  0.0*y0 - 1.0*y1 + 2.0*y2;

  NV_Ith_S(ydot,0) = yd0;
  NV_Ith_S(ydot,1) = yd1;
  NV_Ith_S(ydot,2) = yd2;

  return 0;
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

  realtype *rdata = (realtype *) user_data;
  realtype lam = rdata[0];
  DlsMat V  = NewDenseMat(3,3);
  DlsMat D  = NewDenseMat(3,3);
  DlsMat Vi = NewDenseMat(3,3);

  /* initialize temporary matrices to zero */
  DenseScale(0.0, V);
  DenseScale(0.0, D);
  DenseScale(0.0, Vi);

  /* J = V*D*Vi, where 
        V = [1 -1 1; -1 2 1; 0 -1 2] 
        Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1]
        D = [-0.5 0 0; 0 -0.1 0; 0 0 lam] */
  DENSE_ELEM(V,0,0) =  1.0;
  DENSE_ELEM(V,0,1) = -1.0;
  DENSE_ELEM(V,0,2) =  1.0;
  DENSE_ELEM(V,1,0) = -1.0;
  DENSE_ELEM(V,1,1) =  2.0;
  DENSE_ELEM(V,1,2) =  1.0;
  DENSE_ELEM(V,2,0) =  0.0;
  DENSE_ELEM(V,2,1) = -1.0;
  DENSE_ELEM(V,2,2) =  2.0;

  DENSE_ELEM(Vi,0,0) =  0.25*5.0;
  DENSE_ELEM(Vi,0,1) =  0.25*1.0;
  DENSE_ELEM(Vi,0,2) = -0.25*3.0;
  DENSE_ELEM(Vi,1,0) =  0.25*2.0;
  DENSE_ELEM(Vi,1,1) =  0.25*2.0;
  DENSE_ELEM(Vi,1,2) = -0.25*2.0;
  DENSE_ELEM(Vi,2,0) =  0.25*1.0;
  DENSE_ELEM(Vi,2,1) =  0.25*1.0;
  DENSE_ELEM(Vi,2,2) =  0.25*1.0;

  DENSE_ELEM(D,0,0) = -0.5;
  DENSE_ELEM(D,1,1) = -0.1;
  DENSE_ELEM(D,2,2) = lam;

  /* J = D*Vi */
  if (dense_MM(D,Vi,J) != 0) {
    fprintf(stderr, "matmul error\n");
    return 1;
  }

  /* D = V*J [= V*D*Vi] */
  if (dense_MM(V,J,D) != 0) {
    fprintf(stderr, "matmul error\n");
    return 1;
  }

  /* J = D [= V*D*Vi] */
  DenseCopy(D, J);

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
             NULL pointer  */
static int check_flag(void *flagvalue, char *funcname, int opt) {

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


/* sol routine to compute the ODE solution y(t). */
static int sol(realtype t, realtype lam, N_Vector y) {

  realtype y0, y1, y2;
  realtype x0 = 1.0, x1 = 1.0, x2 = 1.0;

  /* y = V*exp(D*t)*Vi*y0, where 
        V = [1 -1 1; -1 2 1; 0 -1 2] 
        Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1]
        D = [-0.5 0 0; 0 -0.1 0; 0 0 lam]
        y0 = [1, 1, 1]' */

  /*   y = Vi*x [= Vi*y0] */
  y0 = 0.25*(5.0*x0 + 1.0*x1 - 3.0*x2);
  y1 = 0.25*(2.0*x0 + 2.0*x1 - 2.0*x2);
  y2 = 0.25*(1.0*x0 + 1.0*x1 + 1.0*x2);

  /*   x = exp(D*t)*y */
  x0  = exp(-0.5*t)*y0;
  x1  = exp(-0.1*t)*y1;
  x2  = exp( lam*t)*y2;

  /*   y = V*x */
  y0 =  1.0*x0 - 1.0*x1 + 1.0*x2;
  y1 = -1.0*x0 + 2.0*x1 + 1.0*x2;
  y2 =  0.0*x0 - 1.0*x1 + 2.0*x2;

  NV_Ith_S(y,0) = y0;
  NV_Ith_S(y,1) = y1;
  NV_Ith_S(y,2) = y2;

  return 0;
}


/* DlsMat matrix-multiply utility routine: C = A*B. */
static int dense_MM(DlsMat A, DlsMat B, DlsMat C) {

  /* check for legal dimensions */
  if ((A->N != B->M) || (C->M != A->M) || (C->N != B->N)) {
    fprintf(stderr, "\n matmul error: dimension mismatch\n\n");
    return 1;
  }
    
  /* access data and extents */
  realtype **adata = A->cols;
  realtype **bdata = B->cols;
  realtype **cdata = C->cols;
  long int m = C->M;
  long int n = C->N;
  long int l = A->N;
  int i, j, k;

  /* initialize output */
  DenseScale(0.0, C);

  /* perform multiply (not optimal, but fine for 3x3 matrices) */
  for (i=0; i<m; i++) 
    for (j=0; j<n; j++) 
      for (k=0; k<l; k++) 
	cdata[i][j] += adata[i][k] * bdata[k][j];

  return 0;
}

/*---- end of file ----*/
