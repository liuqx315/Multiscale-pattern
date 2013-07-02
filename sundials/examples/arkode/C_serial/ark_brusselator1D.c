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

 This program solves the problem with the BDF method, using a
 Newton iteration with the CVBAND band linear solver, and a
 user-supplied Jacobian routine.

 100 outputs are printed at equal intervals, and run statistics 
 are printed at the end.
---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode.h>            /* prototypes for ARKode fcts., consts. */
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., macros */
#include <arkode/arkode_band.h>       /* prototype for ARKBand solver */
#include <sundials/sundials_band.h>   /* defs. of DlsMat and BAND_ELEM */
#include <sundials/sundials_types.h>  /* def. of type 'realtype' */

/* accessor macros between (x,v) location and 1D NVector array */
#define IDX(x,v) (3*(x)+v)

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

/* Private helper functions  */
static int LaplaceMatrix(realtype c, DlsMat Jac, UserData udata);
static int ReactionJac(realtype c, N_Vector y, DlsMat Jac, UserData udata);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);

/* Main Program */
int main()
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);    /* initial time */
  realtype Tf = RCONST(10.0);   /* final time */
  int Nt = 100;                 /* total number of output times */
  int Nvar = 3;                 /* number of solution fields */
  UserData udata = NULL;
  realtype *data;
  long int N = 201;             /* spatial mesh size */
  realtype a = 0.6;             /* problem parameters */
  realtype b = 2.0;
  realtype du = 0.025;
  realtype dv = 0.025;
  realtype dw = 0.025;
  realtype ep = 1.0e-5;         /* stiffness parameter */
  realtype reltol = 1.0e-6;     /* tolerances */
  realtype abstol = 1.0e-10;
  long int NEQ, i;

  /* general problem variables */
  int flag;                     /* reusable error-checking flag */
  N_Vector y = NULL;            /* empty vector for storing solution */
  N_Vector umask = NULL;        /* empty mask vectors for viewing solution components */
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  void *arkode_mem = NULL;      /* empty ARKode memory structure */

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));

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

  /* Initialize data structures */
  y = N_VNew_Serial(NEQ);           /* Create serial vector for solution */
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  udata->dx = RCONST(1.0)/(N-1);    /* set spatial mesh spacing */
  data = N_VGetArrayPointer(y);     /* Access data array for new NVector y */
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  umask = N_VNew_Serial(NEQ);       /* Create serial vector masks */
  if (check_flag((void *)umask, "N_VNew_Serial", 0)) return 1;
  vmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)vmask, "N_VNew_Serial", 0)) return 1;
  wmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)wmask, "N_VNew_Serial", 0)) return 1;

  /* Set initial conditions into y */
  realtype pi = RCONST(4.0)*atan(RCONST(1.0));
  for (i=0; i<N; i++) {
    data[IDX(i,0)] =  a  + RCONST(0.1)*sin(pi*i*udata->dx);  /* u */
    data[IDX(i,1)] = b/a + RCONST(0.1)*sin(pi*i*udata->dx);  /* v */
    data[IDX(i,2)] =  b  + RCONST(0.1)*sin(pi*i*udata->dx);  /* w */
  }

  /* Set mask array values for each solution component */
  N_VConst(0.0, umask);
  data = N_VGetArrayPointer(umask);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,0)] = RCONST(1.0);

  N_VConst(0.0, vmask);
  data = N_VGetArrayPointer(vmask);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,1)] = RCONST(1.0);

  N_VConst(0.0, wmask);
  data = N_VGetArrayPointer(wmask);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,2)] = RCONST(1.0);

  /* Create the solver memory */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return 1;

  /* Call ARKodeInit to initialize the integrator memory and specify the
     hand-side side function in y'=f(t,y), the inital time T0, and
     the initial dependent variable vector y.  Note: since this
     problem is fully implicit, we set f_E to NULL and f_I to f. */
  flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return 1;

  /* Set routines */
  flag = ARKodeSetUserData(arkode_mem, (void *) udata);     /* Pass udata to user functions */
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);    /* Specify tolerances */
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Linear solver specification */
  flag = ARKBand(arkode_mem, NEQ, 4, 4);          /* Specify the band linear solver */
  if (check_flag(&flag, "ARKBand", 1)) return 1;
  flag = ARKDlsSetBandJacFn(arkode_mem, Jac);     /* Set the Jacobian routine */
  if (check_flag(&flag, "ARKDlsSetBandJacFn", 1)) return 1;

  /* output spatial mesh to disk */
  FILE *FID = fopen("bruss_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->dx*i);
  fclose(FID);

  /* Open output streams for results, access data array */
  FILE *UFID=fopen("bruss_u.txt","w");
  FILE *VFID=fopen("bruss_v.txt","w");
  FILE *WFID=fopen("bruss_w.txt","w");
  data = N_VGetArrayPointer(y);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;

  /* output initial condition to disk */
  for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[IDX(i,0)]);
  for (i=0; i<N; i++)  fprintf(VFID," %.16e", data[IDX(i,1)]);
  for (i=0; i<N; i++)  fprintf(WFID," %.16e", data[IDX(i,2)]);
  fprintf(UFID,"\n");
  fprintf(VFID,"\n");
  fprintf(WFID,"\n");

  /* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached */
  realtype t = T0;
  realtype dTout = (Tf-T0)/Nt;
  realtype tout = T0+dTout;
  realtype u, v, w;
  printf("        t      ||u||_rms   ||v||_rms   ||w||_rms\n");
  printf("   ----------------------------------------------\n");
  int iout;
  for (iout=0; iout<Nt; iout++) {

    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);    /* call integrator */
    if (check_flag(&flag, "ARKode", 1)) break;
    u = N_VWL2Norm(y,umask);                               /* access/print solution statistics */
    u = sqrt(u*u/N);
    v = N_VWL2Norm(y,vmask);
    v = sqrt(v*v/N);
    w = N_VWL2Norm(y,wmask);
    w = sqrt(w*w/N);
    printf("  %10.6f  %10.6f  %10.6f  %10.6f\n", t, u, v, w);
    if (flag >= 0) {                                       /* successful solve: update output time */
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                               /* unsuccessful solve: break */
      fprintf(stderr,"Solver failure, stopping integration\n");
      break;
    }

    /* output results to disk */
    for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[IDX(i,0)]);
    for (i=0; i<N; i++)  fprintf(VFID," %.16e", data[IDX(i,1)]);
    for (i=0; i<N; i++)  fprintf(WFID," %.16e", data[IDX(i,2)]);
    fprintf(UFID,"\n");
    fprintf(VFID,"\n");
    fprintf(WFID,"\n");
  }
  printf("   ----------------------------------------------\n");
  fclose(UFID);
  fclose(VFID);
  fclose(WFID);

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
  printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total linear solver setups = %li\n", nsetups);
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of linear solver convergence failures = %li\n", ncfn);
  printf("   Total number of error test failures = %li\n\n", netf);

  /* Clean up and return with successful completion */
  N_VDestroy_Serial(y);         /* Free vectors */
  N_VDestroy_Serial(umask);
  N_VDestroy_Serial(vmask);
  N_VDestroy_Serial(wmask);
  free(udata);                  /* Free user data */
  ARKodeFree(&arkode_mem);      /* Free integrator memory */
  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  N_VConst(0.0, ydot);                        /* initialize ydot to zero */
  UserData udata = (UserData) user_data;      /* access problem data */
  long int N  = udata->N;                     /* set variable shortcuts */
  realtype a  = udata->a;
  realtype b  = udata->b;
  realtype ep = udata->ep;
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;
  realtype dx = udata->dx;
  realtype *Ydata = N_VGetArrayPointer(y);     /* access data arrays */
  if (check_flag((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;
  realtype *dYdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *)dYdata, "N_VGetArrayPointer", 0)) return 1;

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

    /* Fill in ODE RHS for u */
    dYdata[IDX(i,0)] = (ul - RCONST(2.0)*u + ur)*uconst + a - (w+RCONST(1.0))*u + v*u*u;

    /* Fill in ODE RHS for v */
    dYdata[IDX(i,1)] = (vl - RCONST(2.0)*v + vr)*vconst + w*u - v*u*u;

    /* Fill in ODE RHS for w */
    dYdata[IDX(i,2)] = (wl - RCONST(2.0)*w + wr)*wconst + (b-w)/ep - w*u;
  }

  /* enforce stationary boundaries */
  dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
  dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;

  return 0;     /* Return with success */
}

/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(long int M, long int mu, long int ml,
               realtype t, N_Vector y, N_Vector fy, 
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  SetToZero(J);                              /* Initialize Jacobian to zero */
  UserData udata = (UserData) user_data;     /* access problem data */

  /* Fill in the Laplace matrix */
  LaplaceMatrix(RCONST(1.0), J, udata);

  /* Add in the Jacobian of the reaction terms matrix */
  ReactionJac(RCONST(1.0), y, J, udata);

  return 0;                                  /* Return with success */
}

/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Routine to compute the stiffness matrix from (L*y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
static int LaplaceMatrix(realtype c, DlsMat Jac, UserData udata)
{
  long int i;                /* set shortcuts */
  long int N = udata->N;
  realtype dx = udata->dx;

  /* iterate over intervals, filling in Jacobian of (L*y) */
  for (i=1; i<N-1; i++) {
    BAND_ELEM(Jac,IDX(i,0),IDX(i-1,0)) += c*udata->du/dx/dx;
    BAND_ELEM(Jac,IDX(i,1),IDX(i-1,1)) += c*udata->dv/dx/dx;
    BAND_ELEM(Jac,IDX(i,2),IDX(i-1,2)) += c*udata->dw/dx/dx;
    BAND_ELEM(Jac,IDX(i,0),IDX(i,0)) += -c*RCONST(2.0)*udata->du/dx/dx;
    BAND_ELEM(Jac,IDX(i,1),IDX(i,1)) += -c*RCONST(2.0)*udata->dv/dx/dx;
    BAND_ELEM(Jac,IDX(i,2),IDX(i,2)) += -c*RCONST(2.0)*udata->dw/dx/dx;
    BAND_ELEM(Jac,IDX(i,0),IDX(i+1,0)) += c*udata->du/dx/dx;
    BAND_ELEM(Jac,IDX(i,1),IDX(i+1,1)) += c*udata->dv/dx/dx;
    BAND_ELEM(Jac,IDX(i,2),IDX(i+1,2)) += c*udata->dw/dx/dx;
  }

  return 0;                  /* Return with success */
}

/* Routine to compute the Jacobian matrix from R(y), scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
static int ReactionJac(realtype c, N_Vector y, DlsMat Jac, UserData udata)
{
  long int N  = udata->N;                      /* set shortcuts */
  long int i;
  realtype u, v, w;
  realtype ep = udata->ep;
  realtype *Ydata = N_VGetArrayPointer(y);     /* access solution array */
  if (check_flag((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over nodes, filling in Jacobian of reaction terms */
  for (i=1; i<N-1; i++) {

    u = Ydata[IDX(i,0)];                       /* set nodal value shortcuts */
    v = Ydata[IDX(i,1)];
    w = Ydata[IDX(i,2)];

    /* all vars wrt u */
    BAND_ELEM(Jac,IDX(i,0),IDX(i,0)) += c*(RCONST(2.0)*u*v-(w+RCONST(1.0)));
    BAND_ELEM(Jac,IDX(i,1),IDX(i,0)) += c*(w - RCONST(2.0)*u*v);
    BAND_ELEM(Jac,IDX(i,2),IDX(i,0)) += c*(-w);

    /* all vars wrt v */
    BAND_ELEM(Jac,IDX(i,0),IDX(i,1)) += c*(u*u);
    BAND_ELEM(Jac,IDX(i,1),IDX(i,1)) += c*(-u*u);

    /* all vars wrt w */
    BAND_ELEM(Jac,IDX(i,0),IDX(i,2)) += c*(-u);
    BAND_ELEM(Jac,IDX(i,1),IDX(i,2)) += c*(u);
    BAND_ELEM(Jac,IDX(i,2),IDX(i,2)) += c*(-RCONST(1.0)/ep - u);

  }

  return 0;                                   /* Return with success */
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
