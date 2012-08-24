/* -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following test simulates a brusselator problem from chemical 
 * kinetics.  This is n PDE system with 3 components, Y = [u,v,w], 
 * satisfying the equations,
 *    u_t = du*u_xx + a - (w+1)*u + v*u^2
 *    v_t = dv*v_xx + w*u - v*u^2
 *    w_t = dw*w_xx + (b-w)/ep - w*u
 * for t in [0, 80], x in [0, 1], with initial conditions
 *    u(0,x) =  a  + 0.1*sin(pi*x)
 *    v(0,x) = b/a + 0.1*sin(pi*x)
 *    w(0,x) =  b  + 0.1*sin(pi*x),
 * and with stationary boundary conditions, i.e. 
 *    u_t(t,0) = u_t(t,1) = v_t(t,0) = v_t(t,1) = w_t(t,0) = w_t(t,1) = 0.
 * Note: these can also be implemented as Dirichlet boundary conditions 
 * with values identical to the initial conditions.
 * 
 * The spatial derivatives are computed using second-order centered 
 * differences, with the data distributed over N points on a uniform 
 * spatial grid.
 *
 * The number of spatial points N, the parameters a, b, du, dv, dw 
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
#include <stdlib.h>
#include <math.h>

/* Header files with a description of contents used */

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_band.h>        /* prototype for CVBand */
#include <sundials/sundials_band.h>  /* definitions of type DlsMat and macros */
#include <sundials/sundials_types.h> /* definition of type realtype */


/* accessor macros between (x,v) location and 1D NVector array */
#define IDX(x,v) (3*(x)+v)

/* constants */
#define ONE (RCONST(1.0))
#define TWO (RCONST(2.0))


/* user data structure */
typedef struct {  
  long int N;    /* number of intervals     */
  realtype dx;   /* mesh spacing            */
  realtype a;    /* constant forcing on u   */
  realtype b;    /* steady-state value of w */
  realtype du;   /* diffusion coeff for u   */
  realtype dv;   /* diffusion coeff for v   */
  realtype dw;   /* diffusion coeff for w   */
  realtype ep;   /* stiffness parameter     */
} *UserData;



/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(long int N, long int mu, long int ml,
               realtype t, N_Vector y, N_Vector fy, 
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);
static int LaplaceMatrix(realtype c, DlsMat Jac, UserData udata);
static int ReactionJac(realtype c, N_Vector y, DlsMat Jac, UserData udata);



/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);
  realtype Tf = RCONST(10.0);
  int Nt = 100;
  int Nvar = 3;
  UserData udata = NULL;
  realtype *data;
  long int N, NEQ, i;

  /* general problem variables */
  int flag;
  N_Vector y = NULL;
  N_Vector ytrue = NULL;
  N_Vector yerr  = NULL;
  N_Vector umask = NULL;
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  void *cvode_mem = NULL;
  void *cvtrue_mem = NULL;

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  if (check_flag((void *)udata, "malloc", 2)) return(1);

  /* read problem parameter and tolerances from input file:
     N - number of spatial discretization points
     a - constant forcing on u
     b - steady-state value of w
     du - diffusion coefficient for u
     dv - diffusion coefficient for v
     dw - diffusion coefficient for w
     ep - stiffness parameter
     reltol - desired relative tolerance
     abstol - desired absolute tolerance */
  double a, b, du, dv, dw, ep, reltol, abstol;
  FILE *FID;
  FID=fopen("input_brusselator1D.txt","r");
  fscanf(FID,"  N = %li\n", &N);
  fscanf(FID,"  a = %lf\n", &a);
  fscanf(FID,"  b = %lf\n", &b);
  fscanf(FID,"  du = %lf\n", &du);
  fscanf(FID,"  dv = %lf\n", &dv);
  fscanf(FID,"  dw = %lf\n", &dw);
  fscanf(FID,"  ep = %lf\n", &ep);
  fscanf(FID,"  reltol = %lf\n", &reltol);
  fscanf(FID,"  abstol = %lf\n", &abstol);
  fclose(FID);

  /* store the inputs in the UserData structure */
  udata->N  = N;
  udata->a  = a;
  udata->b  = b;
  udata->du = du;
  udata->dv = dv;
  udata->dw = dw;
  udata->ep = ep;

  /* set total allocated vector length */
  NEQ = Nvar*udata->N;

  /* Initial problem output */
  printf("\n1D Brusselator PDE test problem:\n");
  printf("    N = %li,  NEQ = %li\n", udata->N, NEQ);
  printf("    problem parameters:  a = %g,  b = %g,  ep = %g\n",
	 udata->a, udata->b, udata->ep);
  printf("    diffusion coefficients:  du = %g,  dv = %g,  dw = %g\n", 
	 udata->du, udata->dv, udata->dw);
  printf("    reltol = %.1e,  abstol = %.1e\n\n", reltol, abstol);
  realtype reltol2 = reltol*1.0e-3;
  realtype abstol2 = abstol*1.0e-3;

  /* Create serial vector of length NEQ for initial condition */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  ytrue = N_VNew_Serial(NEQ);
  if (check_flag((void *)ytrue, "N_VNew_Serial", 0)) return(1);
  yerr = N_VNew_Serial(NEQ);
  if (check_flag((void *)yerr, "N_VNew_Serial", 0)) return(1);

  /* set spatial mesh spacing */
  udata->dx = ONE/(N-1);

  /* output mesh to disk */
  FID=fopen("bruss_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->dx*i);
  fclose(FID);
  

  /* Access data array for new NVector y */
  data = N_VGetArrayPointer(y);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return(1);

  /* Set initial conditions into y, ytrue */
  realtype pi = RCONST(4.0)*atan(ONE);
  for (i=0; i<N; i++) {
    data[IDX(i,0)] =  a  + RCONST(0.1)*sin(pi*i*udata->dx);  /* u */
    data[IDX(i,1)] = b/a + RCONST(0.1)*sin(pi*i*udata->dx);  /* v */
    data[IDX(i,2)] =  b  + RCONST(0.1)*sin(pi*i*udata->dx);  /* w */
  }
  N_VScale(1.0, y, ytrue);

  /* Create serial vector masks for each solution component */
  umask = N_VNew_Serial(NEQ);
  if (check_flag((void *)umask, "N_VNew_Serial", 0)) return(1);
  vmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)vmask, "N_VNew_Serial", 0)) return(1);
  wmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)wmask, "N_VNew_Serial", 0)) return(1);

  /* Set mask array values for each solution component */
  N_VConst(0.0, umask);
  data = N_VGetArrayPointer(umask);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return(1);
  for (i=0; i<N; i++)  data[IDX(i,0)] = ONE;

  N_VConst(0.0, vmask);
  data = N_VGetArrayPointer(vmask);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return(1);
  for (i=0; i<N; i++)  data[IDX(i,1)] = ONE;

  N_VConst(0.0, wmask);
  data = N_VGetArrayPointer(wmask);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return(1);
  for (i=0; i<N; i++)  data[IDX(i,2)] = ONE;


  /* Call CVodeCreate to create the solver memory and specify the 
     Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
  cvtrue_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvtrue_mem, "CVodeCreate", 0)) return(1);
  
  /* Call CVodeInit to initialize the integrator memory and specify the
     user's right hand side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y */
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);
  flag = CVodeInit(cvtrue_mem, f, T0, ytrue);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);
  
  /* Call CVodeSetUserData to pass rdata to user functions */
  flag = CVodeSetUserData(cvode_mem, (void *) udata);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);
  flag = CVodeSetUserData(cvtrue_mem, (void *) udata);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* Call CVodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = CVodeSStolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);
  flag = CVodeSStolerances(cvtrue_mem, reltol2, abstol2);
  if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);

  /* Call CVBand to specify the CVBAND band linear solver */
  flag = CVBand(cvode_mem, NEQ, 4, 4);
  if (check_flag(&flag, "CVBand", 1)) return(1);
  flag = CVBand(cvtrue_mem, NEQ, 4, 4);
  if (check_flag(&flag, "CVBand", 1)) return(1);

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = CVDlsSetBandJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDlsSetBandJacFn", 1)) return(1);
  flag = CVDlsSetBandJacFn(cvtrue_mem, Jac);
  if (check_flag(&flag, "CVDlsSetBandJacFn", 1)) return(1);

  /* Open output stream for results, access data arrays */
  FILE *UFID=fopen("bruss_u.txt","w");
  FILE *VFID=fopen("bruss_v.txt","w");
  FILE *WFID=fopen("bruss_w.txt","w");
  data = N_VGetArrayPointer(y);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return(1);

  /* output initial condition to disk */
  for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[IDX(i,0)]);
  for (i=0; i<N; i++)  fprintf(VFID," %.16e", data[IDX(i,1)]);
  for (i=0; i<N; i++)  fprintf(WFID," %.16e", data[IDX(i,2)]);
  fprintf(UFID,"\n");
  fprintf(VFID,"\n");
  fprintf(WFID,"\n");

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t  = T0;
  realtype t2 = T0;
  realtype dTout = Tf/Nt;
  realtype tout = dTout;
  realtype u, v, w, uerr, verr, werr, errI=0.0, err2=0.0;
  printf("        t      ||u||_rms   ||v||_rms   ||w||_rms     ||uerr||      ||verr||      ||werr||\n");
  printf("   ---------------------------------------------------------------------------------------\n");
  int iout;
  for (iout=0; iout<Nt; iout++) {

    flag = CVode(cvtrue_mem, tout, ytrue, &t2, CV_NORMAL);
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    u = N_VWL2Norm(y,umask);
    u = sqrt(u*u/N);
    v = N_VWL2Norm(y,vmask);
    v = sqrt(v*v/N);
    w = N_VWL2Norm(y,wmask);
    w = sqrt(w*w/N);

    N_VLinearSum( 1.0, ytrue, -1.0, y, yerr );
    uerr = N_VWL2Norm(yerr,umask);  
    uerr = sqrt(uerr*uerr/N);
    verr = N_VWL2Norm(yerr,vmask);  
    verr = sqrt(verr*verr/N);
    werr = N_VWL2Norm(yerr,wmask);  
    werr = sqrt(werr*werr/N);
    errI = (errI > N_VMaxNorm(yerr)) ? errI : N_VMaxNorm(yerr);
    err2 += uerr*uerr + verr*verr + werr*werr;
    printf("  %10.6f  %10.6f  %10.6f  %10.6f  %12.5e  %12.5e  %12.5e\n", 
	   t, u, v, w, uerr, verr, werr);

    /* output results to disk */
    for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[IDX(i,0)]);
    for (i=0; i<N; i++)  fprintf(VFID," %.16e", data[IDX(i,1)]);
    for (i=0; i<N; i++)  fprintf(WFID," %.16e", data[IDX(i,2)]);
    fprintf(UFID,"\n");
    fprintf(VFID,"\n");
    fprintf(WFID,"\n");

    if (check_flag(&flag, "CVode", 1)) break;
    if (flag == CV_SUCCESS) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
  }
  err2 = sqrt(err2 / 3.0 / N / Nt);
  printf("   ---------------------------------------------------------------------------------------\n");
  fclose(UFID);
  fclose(VFID);
  fclose(WFID);
    

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
  printf("   Total number of error test failures = %li\n", netf);
  printf("   Error: max = %g, rms = %g\n", errI, err2);
  printf("   Oversolve = %g\n\n", reltol/err2);

  /* Free vectors */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(ytrue);
  N_VDestroy_Serial(yerr);
  N_VDestroy_Serial(umask);
  N_VDestroy_Serial(vmask);
  N_VDestroy_Serial(wmask);

  /* Free user data */
  free(udata);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);
  CVodeFree(&cvtrue_mem);

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

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of intervals, background values */
  long int N  = udata->N;
  realtype a  = udata->a;
  realtype b  = udata->b;
  realtype ep = udata->ep;
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;
  realtype dx = udata->dx;

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *)Ydata, "N_VGetArrayPointer", 0)) return(1);
  realtype *dYdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *)dYdata, "N_VGetArrayPointer", 0)) return(1);

  /* iterate over domain, computing all equations */
  realtype uconst = du/dx/dx;
  realtype vconst = dv/dx/dx;
  realtype wconst = dw/dx/dx;
  realtype u, ul, ur, v, vl, vr, w, wl, wr;
  long int i;
  for (i=1; i<N-1; i++) {

    /* set shortcuts */
    u = Ydata[IDX(i,0)];  ul = Ydata[IDX(i-1,0)];  ur = Ydata[IDX(i+1,0)];
    v = Ydata[IDX(i,1)];  vl = Ydata[IDX(i-1,1)];  vr = Ydata[IDX(i+1,1)];
    w = Ydata[IDX(i,2)];  wl = Ydata[IDX(i-1,2)];  wr = Ydata[IDX(i+1,2)];

    /* u_t = du*u_xx + a - (w+1)*u + v*u^2 */
    dYdata[IDX(i,0)] = (ul - TWO*u + ur)*uconst + a - (w+ONE)*u + v*u*u;

    /* v_t = dv*v_xx + w*u - v*u^2 */
    dYdata[IDX(i,1)] = (vl - TWO*v + vr)*vconst + w*u - v*u*u;

    /* w_t = dw*w_xx + (b-w)/ep - w*u */
    dYdata[IDX(i,2)] = (wl - TWO*w + wr)*wconst + (b-w)/ep - w*u;

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

  /* problem data */
  UserData udata = (UserData) user_data;

  /* Fill in the Laplace matrix */
  if (LaplaceMatrix(ONE, J, udata)) {
    printf("Jacobian calculation error in calling LaplaceMatrix!\n");
    return(1);
  }

  /* Add in the Jacobian of the reaction terms matrix */
  if (ReactionJac(ONE, y, J, udata)) {
    printf("Jacobian calculation error in calling ReactionJac!\n");
    return(1);
  }

  return(0);
}



/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Routine to compute the stiffness matrix from (L*y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
static int LaplaceMatrix(realtype c, DlsMat Jac, UserData udata)
{
  /* shortcut to number of intervals */
  long int N = udata->N;

  /* set shortcuts */
  long int i;
  realtype dx = udata->dx;
  
  /* iterate over intervals, filling in Jacobian entries */
  for (i=1; i<N-1; i++) {

    /* Jacobian of (L*y) at this node */
    BAND_ELEM(Jac,IDX(i,0),IDX(i-1,0)) += c*udata->du/dx/dx;
    BAND_ELEM(Jac,IDX(i,1),IDX(i-1,1)) += c*udata->dv/dx/dx;
    BAND_ELEM(Jac,IDX(i,2),IDX(i-1,2)) += c*udata->dw/dx/dx;
    BAND_ELEM(Jac,IDX(i,0),IDX(i,0)) += -c*TWO*udata->du/dx/dx;
    BAND_ELEM(Jac,IDX(i,1),IDX(i,1)) += -c*TWO*udata->dv/dx/dx;
    BAND_ELEM(Jac,IDX(i,2),IDX(i,2)) += -c*TWO*udata->dw/dx/dx;
    BAND_ELEM(Jac,IDX(i,0),IDX(i+1,0)) += c*udata->du/dx/dx;
    BAND_ELEM(Jac,IDX(i,1),IDX(i+1,1)) += c*udata->dv/dx/dx;
    BAND_ELEM(Jac,IDX(i,2),IDX(i+1,2)) += c*udata->dw/dx/dx;
  }

  return(0);
}



/* Routine to compute the Jacobian matrix from R(y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
static int ReactionJac(realtype c, N_Vector y, DlsMat Jac, UserData udata)
{

  /* shortcuts to number of intervals, background values */
  long int N  = udata->N;

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *)Ydata, "N_VGetArrayPointer", 0)) return(1);

  /* set shortcuts */
  long int i;
  realtype u, v, w;
  realtype ep = udata->ep;
  
  /* iterate over nodes, filling in Jacobian entries */
  for (i=1; i<N-1; i++) {

    /* set nodal value shortcuts (shifted index due to start at first interior node) */
    u = Ydata[IDX(i,0)];
    v = Ydata[IDX(i,1)];
    w = Ydata[IDX(i,2)];

    /* all vars wrt u */
    BAND_ELEM(Jac,IDX(i,0),IDX(i,0)) += c*(TWO*u*v-(w+ONE));
    BAND_ELEM(Jac,IDX(i,1),IDX(i,0)) += c*(w - TWO*u*v);
    BAND_ELEM(Jac,IDX(i,2),IDX(i,0)) += c*(-w);

    /* all vars wrt v */
    BAND_ELEM(Jac,IDX(i,0),IDX(i,1)) += c*(u*u);
    BAND_ELEM(Jac,IDX(i,1),IDX(i,1)) += c*(-u*u);

    /* all vars wrt w */
    BAND_ELEM(Jac,IDX(i,0),IDX(i,2)) += c*(-u);
    BAND_ELEM(Jac,IDX(i,1),IDX(i,2)) += c*(u);
    BAND_ELEM(Jac,IDX(i,2),IDX(i,2)) += c*(-ONE/ep - u);

  }

  return(0);
}



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
