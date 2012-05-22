/*---------------------------------------------------------------
 $Revision: $
 $Date: $
-----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Example problem:
 
 The following test simulates a brusselator problem from chemical 
 kinetics.  This is n PDE system with 3 components, Y = [u,v,w], 
 satisfying the equations,
    u_t = du*u_xx + a - (w+1)*u + v*u^2
    v_t = dv*v_xx + w*u - v*u^2
    w_t = dw*w_xx + (b-w)/ep - w*u
 for t in [0, 80], x in [0, 1], with initial conditions
    u(0,x) =  a  + 0.1*sin(pi*x)
    v(0,x) = b/a + 0.1*sin(pi*x)
    w(0,x) =  b  + 0.1*sin(pi*x),
 and with stationary boundary conditions, i.e. 
    u_t(t,0) = u_t(t,1) = 0,
    v_t(t,0) = v_t(t,1) = 0,
    w_t(t,0) = w_t(t,1) = 0.
 Note: these can also be implemented as Dirichlet boundary 
 conditions with values identical to the initial conditions.
 
 The spatial derivatives are computed using second-order 
 centered differences, with the data distributed over N points 
 on a uniform spatial grid.

 The number of spatial points N, the parameters a, b, du, dv, 
 dw and ep, as well as the desired relative and absolute solver 
 tolerances, are provided in the input file 
 input_brusselator1D.txt.
 
 This program solves the problem with the BDF method, using a
 Newton iteration with the CVBAND band linear solver, and a
 user-supplied Jacobian routine.

 100 outputs are printed at equal intervals, and run statistics 
 are printed at the end.
---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Header files with a description of contents used */

#include <arkode/arkode.h>           /* prototypes for ARKODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <arkode/arkode_band.h>        /* prototype for ARKBand */
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
  N_Vector umask = NULL;
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  void *arkode_mem = NULL;

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


  /* Create serial vector of length NEQ for initial condition */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

  /* set spatial mesh spacing */
  udata->dx = ONE/(N-1);

  /* output mesh to disk */
  FID=fopen("bruss_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->dx*i);
  fclose(FID);
  

  /* Access data array for new NVector y */
  data = N_VGetArrayPointer(y);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return(1);

  /* Set initial conditions into y */
  realtype pi = RCONST(4.0)*atan(ONE);
  for (i=0; i<N; i++) {
    data[IDX(i,0)] =  a  + RCONST(0.1)*sin(pi*i*udata->dx);  /* u */
    data[IDX(i,1)] = b/a + RCONST(0.1)*sin(pi*i*udata->dx);  /* v */
    data[IDX(i,2)] =  b  + RCONST(0.1)*sin(pi*i*udata->dx);  /* w */
  }

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


  /* Call ARKodeCreate to create the solver memory and specify the 
     Backward Differentiation Formula and the use of a Newton iteration */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return(1);
  
  /* Call ARKodeInit to initialize the integrator memory and specify the
     user's right hand side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y */
  flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return(1);
  
  /* Call ARKodeSetUserData to pass rdata to user functions */
  flag = ARKodeSetUserData(arkode_mem, (void *) udata);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return(1);

  /* Call ARKodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return(1);

  /* Call ARKBand to specify the ARKBAND band linear solver */
  flag = ARKBand(arkode_mem, NEQ, 4, 4);
  if (check_flag(&flag, "ARKBand", 1)) return(1);

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = ARKDlsSetBandJacFn(arkode_mem, Jac);
  if (check_flag(&flag, "ARKDlsSetBandJacFn", 1)) return(1);

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

  /* In loop, call ARKode, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t = T0;
  realtype dTout = Tf/Nt;
  realtype tout = dTout;
  realtype u, v, w;
  printf("        t      ||u||_rms   ||v||_rms   ||w||_rms\n");
  printf("   ----------------------------------------------\n");
  int iout;
  for (iout=0; iout<Nt; iout++) {
    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
    u = N_VWL2Norm(y,umask);
    v = N_VWL2Norm(y,vmask);
    w = N_VWL2Norm(y,wmask);
    printf("  %10.6f  %10.6f  %10.6f  %10.6f\n", t, u, v, w);

    /* output results to disk */
    for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[IDX(i,0)]);
    for (i=0; i<N; i++)  fprintf(VFID," %.16e", data[IDX(i,1)]);
    for (i=0; i<N; i++)  fprintf(WFID," %.16e", data[IDX(i,2)]);
    fprintf(UFID,"\n");
    fprintf(VFID,"\n");
    fprintf(WFID,"\n");

    if (check_flag(&flag, "ARKode", 1)) break;
    if (flag == ARK_SUCCESS) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
  }
  printf("   ----------------------------------------------\n");
  fclose(UFID);
  fclose(VFID);
  fclose(WFID);
    

  /* Print some final statistics */
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf;
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKodeGetNumSteps", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe);
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

  /* Free user data */
  free(udata);

  /* Free integrator memory */
  ARKodeFree(&arkode_mem);

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
