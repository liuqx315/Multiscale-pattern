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
 *    u_t = d1*u_xx + a - (w+1)*u + v*u^2
 *    v_t = d2*v_xx + w*u - v*u^2
 *    w_t = d3*w_xx + (b-w)/ep - w*u
 * for t in [0, 80], x in [0, 1], with initial conditions 
 * Y0 = [u0,v0,w0], where
 *    u0(x) = a + 0.1*sin(pi*x)
 *    v0(x) = b/a + 0.1*sin(pi*x)
 *    w0(x) = b + 0.1*sin(pi*x),
 * and with stationary boundary conditions, i.e. 
 *    u_t(t,0) = u_t(t,1) = v_t(t,0) = v_t(t,1) = w_t(t,0) = w_t(t,1) = 0.
 * Note: these can also be implemented as Dirichlet boundary conditions 
 * with values identical to the initial conditions.
 * 
 * The spatial derivatives are computed using second-order centered 
 * differences, with the data distributed over N points on a uniform 
 * spatial grid.
 *
 * The number of spatial points N, the parameters a, b, d1, d2, d3 
 * and ep, as well as the desired relative and absolute solver 
 * tolerances, are provided in the input file input_brusselator1D.txt.
 * 
 * This program solves the problem with the BDF method, using a
 * Newton iteration with the CVBAND band linear solver, and a
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
#include <cvode/cvode_band.h>        /* prototype for CVBand */
#include <sundials/sundials_band.h>  /* definitions of type DlsMat and macros */
#include <sundials/sundials_types.h> /* definition of type realtype */


/* accessor macros between (x,v) location and 1D NVector array */
#define IDX(x,v) (3*(x)+v)


/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(long int N, long int mu, long int ml,
               realtype t, N_Vector y, N_Vector fy, 
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);



/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = 0.0;
  realtype Tf = 10.0;
  int Nt = 100;
  int Nvar = 3;
  realtype a, b, d1, d2, d3, ep, reltol, abstol, dx, pi;
  realtype *Ydata;
  long int N, NEQ, i;

  /* general problem variables */
  int flag;
  N_Vector y = NULL;
  N_Vector umask = NULL;
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  void *cvode_mem = NULL;

  /* read problem parameter and tolerances from input file:
     N - number of spatial discretization points
     a - constant forcing on u
     b - steady-state value of w
     d1 - diffusion coefficient for u
     d2 - diffusion coefficient for v
     d3 - diffusion coefficient for w
     ep - stiffness parameter
     reltol - desired relative tolerance
     abstol - desired absolute tolerance */
  double a_, b_, d1_, d2_, d3_, ep_, reltol_, abstol_;
  FILE *FID;
  FID=fopen("input_brusselator1D.txt","r");
  fscanf(FID,"  N = %li\n", &N);
  fscanf(FID,"  a = %lf\n", &a_);
  fscanf(FID,"  b = %lf\n", &b_);
  fscanf(FID,"  d1 = %lf\n", &d1_);
  fscanf(FID,"  d2 = %lf\n", &d2_);
  fscanf(FID,"  d3 = %lf\n", &d3_);
  fscanf(FID,"  ep = %lf\n", &ep_);
  fscanf(FID,"  reltol = %lf\n", &reltol_);
  fscanf(FID,"  abstol = %lf\n", &abstol_);
  fclose(FID);

  /* convert the inputs to 'realtype' format */
  a  = a_;
  b  = b_;
  d1 = d1_;
  d2 = d2_;
  d3 = d3_;
  ep = ep_;
  reltol = reltol_;
  abstol = abstol_;

  /* set total allocated vector length */
  NEQ = Nvar*N;

  /* Initial problem output */
  printf("\n1D Brusselator PDE test problem:\n");
  printf("    N = %li,  NEQ = %li\n",N,NEQ);
  printf("    problem parameters:  a = %g,  b = %g,  ep = %g\n",a,b,ep);
  printf("    diffusion coefficients:  d1 = %g,  d2 = %g,  d3 = %g\n",d1,d2,d3);
  printf("    reltol = %.1e,  abstol = %.1e\n\n",reltol,abstol);


  /* Create serial vector of length NEQ for initial condition */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  /* Access data array for new NVector y, and set shortcuts for each component */
  Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *)y, "N_VGetArrayPointer", 0)) return(1);

  /* Set initial conditions into y, grouped according to location */
  pi = 4.0*atan(1.0);
  dx = 1.0/(N-1);
  for (i=0; i<N; i++) {
    Ydata[IDX(i,0)] =  a  + 0.1*sin(pi*dx*i);  /* u */
    Ydata[IDX(i,1)] = b/a + 0.1*sin(pi*dx*i);  /* v */
    Ydata[IDX(i,2)] =  b  + 0.1*sin(pi*dx*i);  /* w */
  }

  /* Create serial vector masks for each solution component */
  umask = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  vmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  wmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  /* Set mask array values for each solution component */
  N_VConst(0.0, umask);
  Ydata = N_VGetArrayPointer(umask);
  if (check_flag((void *)y, "N_VGetArrayPointer", 0)) return(1);
  for (i=0; i<N; i++)  Ydata[IDX(i,0)] = 1.0;

  N_VConst(0.0, vmask);
  Ydata = N_VGetArrayPointer(vmask);
  if (check_flag((void *)y, "N_VGetArrayPointer", 0)) return(1);
  for (i=0; i<N; i++)  Ydata[IDX(i,1)] = 1.0;

  N_VConst(0.0, wmask);
  Ydata = N_VGetArrayPointer(wmask);
  if (check_flag((void *)y, "N_VGetArrayPointer", 0)) return(1);
  for (i=0; i<N; i++)  Ydata[IDX(i,2)] = 1.0;


  /* set user data to contain problem-defining parameters */
  realtype rdata[7] = {a, b, d1, d2, d3, ep, dx};

  
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
  flag = CVodeSetUserData(cvode_mem, (void *) rdata);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);

  /* Call CVBand to specify the CVBAND band linear solver */
  flag = CVBand(cvode_mem, NEQ, 4, 4);
  if (check_flag(&flag, "CVBand", 1)) return(1);

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = CVDlsSetBandJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDlsSetBandJacFn", 1)) return(1);

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t = T0;
  realtype dTout = Tf/Nt;
  realtype tout = dTout;
  realtype u, v, w;
  printf("        t      ||u||_rms   ||v||_rms   ||w||_rms\n");
  printf("   ----------------------------------------------\n");
  int iout;
  for (iout=0; iout<Nt; iout++) {
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    u = N_VWL2Norm(y,umask);
    v = N_VWL2Norm(y,vmask);
    w = N_VWL2Norm(y,wmask);
    printf("  %10.6f  %10.6f  %10.6f  %10.6f\n", t, u, v, w);

    if (check_flag(&flag, "CVode", 1)) break;
    if (flag == CV_SUCCESS) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
  }
  printf("   ----------------------------------------------\n");

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

  /* Free vectors */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(umask);
  N_VDestroy_Serial(vmask);
  N_VDestroy_Serial(wmask);

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
  /* clear out ydot (to be careful) */
  N_VConst(0.0, ydot);

  /* access parameters from problem data */
  realtype *rdata = (realtype *) user_data;
  realtype a  = rdata[0];
  realtype b  = rdata[1];
  realtype d1 = rdata[2];
  realtype d2 = rdata[3];
  realtype d3 = rdata[4];
  realtype ep = rdata[5];
  realtype dx = rdata[6];

  /* determine problem size */
  long int N, M, i;
  N_VSpace(y, &M, &N);  /* M is total vector length */
  N = M/3;              /* N is the number of spatial mesh points */

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *)y, "N_VGetArrayPointer", 0)) return(1);
  realtype *dYdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *)y, "N_VGetArrayPointer", 0)) return(1);

  /* iterate over domain, computing all equations */
  realtype uconst = d1/dx/dx;
  realtype vconst = d2/dx/dx;
  realtype wconst = d3/dx/dx;
  realtype u, ul, ur, v, vl, vr, w, wl, wr;
  for (i=1; i<N-1; i++) {

    /* set shortcuts */
    u = Ydata[IDX(i,0)];  ul = Ydata[IDX(i-1,0)];  ur = Ydata[IDX(i+1,0)];
    v = Ydata[IDX(i,1)];  vl = Ydata[IDX(i-1,1)];  vr = Ydata[IDX(i+1,1)];
    w = Ydata[IDX(i,2)];  wl = Ydata[IDX(i-1,2)];  wr = Ydata[IDX(i+1,2)];

    /* u_t = d1*u_xx + a - (w+1)*u + v*u^2 */
    dYdata[IDX(i,0)] = (ul - 2.0*u + ur)*uconst + a - (w+1.0)*u + v*u*u;

    /* v_t = d2*v_xx + w*u - v*u^2 */
    dYdata[IDX(i,1)] = (vl - 2.0*v + vr)*vconst + w*u - v*u*u;

    /* w_t = d3*w_xx + (b-w)/ep - w*u */
    dYdata[IDX(i,2)] = (wl - 2.0*w + wr)*wconst + (b-w)/ep - w*u;

  }

  /* enforce stationary boundaries */
  dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
  dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;

  return(0);
}


/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(long int M, long int mu, long int ml,
               realtype t, N_Vector y, N_Vector fy, 
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  /* clear out Jacobian (to be careful) */
  SetToZero(J);

  /* access parameters from problem data */
  realtype *rdata = (realtype *) user_data;
  realtype d1 = rdata[2];
  realtype d2 = rdata[3];
  realtype d3 = rdata[4];
  realtype ep = rdata[5];
  realtype dx = rdata[6];

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *)y, "N_VGetArrayPointer", 0)) return(1);

  /* set scaled diffusion differencing parameters */
  realtype uconst = d1/dx/dx;
  realtype vconst = d2/dx/dx;
  realtype wconst = d3/dx/dx;

  /* iterate over space, setting Jacobian entries */
  realtype u, v, w;
  long int i, N = M/3;   /* M is the total data length, so M/3 is each var. */
  for (i=1; i<N-1; i++) {

    /* set shortcuts */
    u = Ydata[IDX(i,0)];
    v = Ydata[IDX(i,1)];
    w = Ydata[IDX(i,2)];

    /* set diffusion components (by row, even though that's less efficient) */
    /*     d1*u_xx  */
    BAND_ELEM(J,IDX(i,0),IDX(i-1,0)) = uconst;
    BAND_ELEM(J,IDX(i,0),IDX(i+1,0)) = uconst;
    BAND_ELEM(J,IDX(i,0),IDX(i,0))   = -2.0*uconst;
    /*     d2*v_xx  */
    BAND_ELEM(J,IDX(i,1),IDX(i-1,1)) = vconst;
    BAND_ELEM(J,IDX(i,1),IDX(i+1,1)) = vconst;
    BAND_ELEM(J,IDX(i,1),IDX(i,1))   = -2.0*vconst;
    /*     d3*w_xx  */
    BAND_ELEM(J,IDX(i,2),IDX(i-1,2)) = wconst;
    BAND_ELEM(J,IDX(i,2),IDX(i+1,2)) = wconst;
    BAND_ELEM(J,IDX(i,2),IDX(i,2))   = -2.0*wconst;
    
    /* set reaction components*/
    /*     all vars wrt u */
    BAND_ELEM(J,IDX(i,0),IDX(i,0)) += 2.0*u*v - (w+1.0);
    BAND_ELEM(J,IDX(i,1),IDX(i,0)) += w - 2.0*u*v;
    BAND_ELEM(J,IDX(i,2),IDX(i,0)) -= w;
    /*     all vars wrt v */
    BAND_ELEM(J,IDX(i,0),IDX(i,1)) += u*u;
    BAND_ELEM(J,IDX(i,1),IDX(i,1)) -= u*u;
    /*     all vars wrt w */
    BAND_ELEM(J,IDX(i,0),IDX(i,2)) -= u;
    BAND_ELEM(J,IDX(i,1),IDX(i,2)) += u;
    BAND_ELEM(J,IDX(i,2),IDX(i,2)) += (-1.0/ep - u);

  }

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
