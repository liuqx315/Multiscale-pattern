/*---------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2014, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 ----------------------------------------------------------------
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
 
 This program solves the problem with the DIRK method, using a
 Newton iteration.  The inner linear systems are solved using 
 the ARKSUPERLUMT linear solver.

 100 outputs are printed at equal intervals, and run statistics 
 are printed at the end.
---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode.h>
#include <arkode/arkode_superlumt.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>


/* accessor macros between (x,v) location and 1D NVector array */
#define IDX(x,v) (3*(x)+v)

/* constants */
#define ONE (RCONST(1.0))
#define TWO (RCONST(2.0))

/* user data structure */
typedef struct {  
  int N;         /* number of intervals     */
  realtype dx;   /* mesh spacing            */
  realtype a;    /* constant forcing on u   */
  realtype b;    /* steady-state value of w */
  realtype du;   /* diffusion coeff for u   */
  realtype dv;   /* diffusion coeff for v   */
  realtype dw;   /* diffusion coeff for w   */
  realtype ep;   /* stiffness parameter     */
  SlsMat R;      /* temporary storage       */
} *UserData;


/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(realtype t, N_Vector y, N_Vector fy, 
	       SlsMat J, void *user_data, 
	       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int JacI(realtype t, N_Vector y, N_Vector fy, 
		SlsMat J, void *user_data, 
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);
static int LaplaceMatrix(SlsMat Jac, UserData udata);
static int ReactionJac(N_Vector y, SlsMat Jac, UserData udata);


/* Parameter input helper function */
int init_from_file(void *ark_mem, char *fname, ARKRhsFn f, 
		   ARKRhsFn fe, ARKRhsFn fi, realtype T0, 
		   N_Vector y0, int *ImEx, int *dorder, 
		   int *fxpt, realtype *RTol, realtype *ATol);


/* Main Program */
int main(int argc, char *argv[])
{
  /* general problem parameters */
  realtype T0 = RCONST(0.0);
  realtype Tf = RCONST(10.0);
  int Nt = 10;
  int Nvar = 3;
  UserData udata = NULL;
  realtype *data;
  int i, N, num_threads;
  long int NEQ, NNZ;

  /* declare solver parameters */
  int flag, dense_order, imex, fixedpt;

  /* general problem variables */
  int idense;
  N_Vector y = NULL;
  N_Vector ytrue = NULL;
  N_Vector yerr  = NULL;
  N_Vector umask = NULL;
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  void *arkode_mem = NULL;
  void *arktrue_mem = NULL;

  /* if an argument supplied, set num_threads */
  num_threads = 1;
  if (argc > 1)
    num_threads = strtol(argv[1], NULL, 0);

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  if (check_flag((void *) udata, "malloc", 2)) return 1;

  /* read problem parameter and tolerances from input file:
     N - number of spatial discretization points
     a - constant forcing on u
     b - steady-state value of w
     du - diffusion coefficient for u
     dv - diffusion coefficient for v
     dw - diffusion coefficient for w
     ep - stiffness parameter */
  double a, b, du, dv, dw, ep;
  FILE *FID;
  FID=fopen("input_brusselator1D.txt","r");
  flag = fscanf(FID,"  N = %i\n", &N);
  flag = fscanf(FID,"  a = %lf\n", &a);
  flag = fscanf(FID,"  b = %lf\n", &b);
  flag = fscanf(FID,"  du = %lf\n", &du);
  flag = fscanf(FID,"  dv = %lf\n", &dv);
  flag = fscanf(FID,"  dw = %lf\n", &dw);
  flag = fscanf(FID,"  ep = %lf\n", &ep);
  fclose(FID);

  /* store the inputs in the UserData structure */
  udata->N  = N;
  udata->a  = a;
  udata->b  = b;
  udata->du = du;
  udata->dv = dv;
  udata->dw = dw;
  udata->ep = ep;
  udata->R  = NULL;

  /* open solver diagnostics output file for writing */
  FILE *DFID;
  DFID=fopen("diags_ark_bruss1D_slu.txt","w");
  
  /* set total allocated vector length */
  NEQ = Nvar*udata->N;

  /* Initial problem output */
  printf("\n1D Brusselator PDE test problem:\n");
  printf("    N = %i,  NEQ = %li\n", udata->N, NEQ);
  printf("    num_threads = %i\n", num_threads);
  printf("    problem parameters:  a = %g,  b = %g,  ep = %g\n",
	 udata->a, udata->b, udata->ep);
  printf("    diffusion coefficients:  du = %g,  dv = %g,  dw = %g\n", 
	 udata->du, udata->dv, udata->dw);

  /* Create serial vector of length NEQ for initial condition */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;
  ytrue = N_VNew_Serial(NEQ);
  if (check_flag((void *) ytrue, "N_VNew_Serial", 0)) return 1;
  yerr = N_VNew_Serial(NEQ);
  if (check_flag((void *) yerr, "N_VNew_Serial", 0)) return 1;

  /* set spatial mesh spacing */
  udata->dx = ONE/(N-1);

  /* output mesh to disk */
  FID=fopen("bruss_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->dx*i);
  fclose(FID);
  
  /* Access data array for new NVector y */
  data = N_VGetArrayPointer(y);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;

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
  if (check_flag((void *) umask, "N_VNew_Serial", 0)) return 1;
  vmask = N_VNew_Serial(NEQ);
  if (check_flag((void *) vmask, "N_VNew_Serial", 0)) return 1;
  wmask = N_VNew_Serial(NEQ);
  if (check_flag((void *) wmask, "N_VNew_Serial", 0)) return 1;

  /* Set mask array values for each solution component */
  N_VConst(0.0, umask);
  data = N_VGetArrayPointer(umask);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,0)] = ONE;

  N_VConst(0.0, vmask);
  data = N_VGetArrayPointer(vmask);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,1)] = ONE;

  N_VConst(0.0, wmask);
  data = N_VGetArrayPointer(wmask);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,2)] = ONE;


  /* Call ARKodeCreate to create the solver memory */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *) arkode_mem, "ARKodeCreate", 0)) return 1;
  arktrue_mem = ARKodeCreate();
  if (check_flag((void *) arktrue_mem, "ARKodeCreate", 0)) return 1;
  
  /* Call init_from_file helper routine to read and set solver parameters */
  realtype rtol, atol;
  flag = init_from_file(arkode_mem, "solve_params.txt", f, fe, fi, T0,
			y, &imex, &dense_order, &fixedpt, &rtol, &atol);
  if (check_flag(&flag, "init_from_file", 1)) return 1;
  if (rtol <= 0.0)  rtol = 1.e-6;
  if (atol <= 0.0)  atol = 1.e-10;
  realtype reltol  = rtol;
  realtype abstol  = atol;
  realtype reltol2 = rtol*1.0e-2;
  realtype abstol2 = atol*1.0e-2;

  /* If (dense_order == -1), tell integrator to use tstop */
  if (dense_order == -1) {
    idense = 0;
  } else {    /* otherwise tell integrator to use dense output */
    idense = 1;
  }

  /* Reference solution uses default implicit method */
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

  /* Call ARKodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;
  flag = ARKodeSStolerances(arktrue_mem, reltol2, abstol2);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Specify the SUPERLUMT sparse linear solver */
  NNZ = 5*NEQ;
  flag = ARKSuperLUMT(arkode_mem, num_threads, NEQ, NNZ);
  if (check_flag(&flag, "ARKSuperLUMT", 1)) return 1;
  flag = ARKSuperLUMT(arktrue_mem, num_threads, NEQ, NNZ);
  if (check_flag(&flag, "ARKSuperLUMT", 1)) return 1;

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
  if (check_flag(&flag, "ARKSlsSetSparseJacFn", 1)) return 1;

  /* Open output stream for results, access data arrays */
  FILE *UFID=fopen("bruss_u.txt","w");
  FILE *VFID=fopen("bruss_v.txt","w");
  FILE *WFID=fopen("bruss_w.txt","w");
  data = N_VGetArrayPointer(y);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;

  /* output initial condition to disk */
  for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[IDX(i,0)]);
  for (i=0; i<N; i++)  fprintf(VFID," %.16e", data[IDX(i,1)]);
  for (i=0; i<N; i++)  fprintf(WFID," %.16e", data[IDX(i,2)]);
  fprintf(UFID,"\n");
  fprintf(VFID,"\n");
  fprintf(WFID,"\n");

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
  realtype u, v, w, uerr, verr, werr, errI=0.0, err2=0.0;
  printf("        t      ||u||_rms   ||v||_rms   ||w||_rms     ||uerr||      ||verr||      ||werr||\n");
  printf("   ---------------------------------------------------------------------------------------\n");
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
      return(1);
    }
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
  }
  err2 = sqrt(err2 / 3.0 / N / Nt);
  printf("   ---------------------------------------------------------------------------------------\n");
  fclose(UFID);
  fclose(VFID);
  fclose(WFID);
    

  /* Print some final statistics */
  long int nst, nst_a, nfe, nfi, nsetups, nje, nni, ncfn, netf, nli, nlcf;
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

  /* Free vectors */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(ytrue);
  N_VDestroy_Serial(yerr);
  N_VDestroy_Serial(umask);
  N_VDestroy_Serial(vmask);
  N_VDestroy_Serial(wmask);

  /* Free user data */
  DestroySparseMat(udata->R);
  free(udata);

  /* Free integrator memory */
  ARKodeFree(&arkode_mem);
  ARKodeFree(&arktrue_mem);

  /* close solver diagnostics output file */
  fclose(DFID);

  return 0;
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
  int N       = udata->N;
  realtype a  = udata->a;
  realtype b  = udata->b;
  realtype ep = udata->ep;
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;
  realtype dx = udata->dx;

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;
  realtype *dYdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *) dYdata, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over domain, computing all equations */
  realtype uconst = du/dx/dx;
  realtype vconst = dv/dx/dx;
  realtype wconst = dw/dx/dx;
  realtype u, ul, ur, v, vl, vr, w, wl, wr;
  int i;
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

  return 0;
}


/* fe routine to compute the diffusion portion of f(t,y). */
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  /* clear out ydot (to be careful) */
  N_VConst(0.0, ydot);

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of intervals, background values */
  int N       = udata->N;
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;
  realtype dx = udata->dx;

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;
  realtype *dYdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *) dYdata, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over domain, computing all equations */
  realtype uconst = du/dx/dx;
  realtype vconst = dv/dx/dx;
  realtype wconst = dw/dx/dx;
  realtype u, ul, ur, v, vl, vr, w, wl, wr;
  int i;
  for (i=1; i<N-1; i++) {

    /* set shortcuts */
    u = Ydata[IDX(i,0)];  ul = Ydata[IDX(i-1,0)];  ur = Ydata[IDX(i+1,0)];
    v = Ydata[IDX(i,1)];  vl = Ydata[IDX(i-1,1)];  vr = Ydata[IDX(i+1,1)];
    w = Ydata[IDX(i,2)];  wl = Ydata[IDX(i-1,2)];  wr = Ydata[IDX(i+1,2)];

    /* u_t = du*u_xx */
    dYdata[IDX(i,0)] = (ul - TWO*u + ur)*uconst;

    /* v_t = dv*v_xx */
    dYdata[IDX(i,1)] = (vl - TWO*v + vr)*vconst;

    /* w_t = dw*w_xx */
    dYdata[IDX(i,2)] = (wl - TWO*w + wr)*wconst;

  }

  /* enforce stationary boundaries */
  dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
  dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;

  return 0;
}


/* fi routine to compute the reaction portion of f(t,y). */
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  /* clear out ydot (to be careful) */
  N_VConst(0.0, ydot);

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of intervals, background values */
  int N       = udata->N;
  realtype a  = udata->a;
  realtype b  = udata->b;
  realtype ep = udata->ep;

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;
  realtype *dYdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *) dYdata, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over domain, computing all equations */
  realtype u, v, w;
  int i;
  for (i=1; i<N-1; i++) {

    /* set shortcuts */
    u = Ydata[IDX(i,0)];
    v = Ydata[IDX(i,1)];
    w = Ydata[IDX(i,2)];

    /* u_t = a - (w+1)*u + v*u^2 */
    dYdata[IDX(i,0)] = a - (w+ONE)*u + v*u*u;

    /* v_t = w*u - v*u^2 */
    dYdata[IDX(i,1)] = w*u - v*u*u;

    /* w_t = (b-w)/ep - w*u */
    dYdata[IDX(i,2)] = (b-w)/ep - w*u;

  }

  /* enforce stationary boundaries */
  dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
  dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;

  return 0;
}


/* Jacobian routine to compute J(t,y) = df/dy. */
static int Jac(realtype t, N_Vector y, N_Vector fy, 
	       SlsMat J, void *user_data, 
	       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  /* problem data */
  UserData udata = (UserData) user_data;
  int N = udata->N;

  /* ensure that Jac is the correct size */
  if ((J->M != N*3) || (J->N != N*3)) {
    printf("Jacobian calculation error: matrix is the wrong size!\n");
    return 1;
  }
  
  /* Fill in the Laplace matrix */
  if (LaplaceMatrix(J, udata)) {
    printf("Jacobian calculation error in calling LaplaceMatrix!\n");
    return 1;
  }

  /* Create empty reaction Jacobian matrix (if not done already) */
  if (udata->R == NULL) {
    udata->R = NewSparseMat(J->M, J->N, J->NNZ);
    if (udata->R == NULL) {
      printf("Jacobian calculation error in allocating R matrix!\n");
      return 1;
    }
  }
      
  /* Add in the Jacobian of the reaction terms matrix */
  if (ReactionJac(y, udata->R, udata)) {
    printf("Jacobian calculation error in calling ReactionJac!\n");
    return 1;
  }

  /* printf("Laplace matrix:\n"); */
  /* PrintSparseMat(J); */

  /* printf("Reaction matrix:\n"); */
  /* PrintSparseMat(udata->R); */

  /* Add R to J */
  if (SlsAddMat(J,udata->R) != 0) {
    printf("Jacobian calculation error in adding sparse matrices!\n");
    return 1;
  }

  return 0;
}


/* Jacobian routine to compute J(t,y) = dfi/dy. */
static int JacI(realtype t, N_Vector y, N_Vector fy, 
		SlsMat J, void *user_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  /* problem data */
  UserData udata = (UserData) user_data;

  /* Add in the Jacobian of the reaction terms matrix */
  if (ReactionJac(y, J, udata)) {
    printf("Jacobian calculation error in calling ReactionJac!\n");
    return 1;
  }

  return 0;
}




/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Routine to compute the stiffness matrix from (L*y) */
static int LaplaceMatrix(SlsMat Lap, UserData udata)
{
  /* shortcut to number of intervals */
  int N = udata->N;

  /* set shortcuts */
  int i, nz=0;
  realtype uconst = (udata->du)/(udata->dx)/(udata->dx);
  realtype uconst2 = -TWO*uconst;
  realtype vconst = (udata->dv)/(udata->dx)/(udata->dx);
  realtype vconst2 = -TWO*vconst;
  realtype wconst = (udata->dw)/(udata->dx)/(udata->dx);
  realtype wconst2 = -TWO*wconst;

  /* clear out matrix */
  SlsSetToZero(Lap);

  /* set first column to zero */
  Lap->colptrs[IDX(0,0)] = nz;
  Lap->colptrs[IDX(0,1)] = nz;
  Lap->colptrs[IDX(0,2)] = nz;
  
  /* iterate over nodes, filling in Laplacian entries depending on these */
  for (i=1; i<N-1; i++) {

    /* dependence on u at this node */
    Lap->colptrs[IDX(i,0)] = nz;
    if (i>1) {                /* node to left */
      Lap->data[nz] = uconst;
      Lap->rowvals[nz++] = IDX(i-1,0);
    }

    Lap->data[nz] = uconst2;  /* self */
    Lap->rowvals[nz++] = IDX(i,0);

    if (i<N-2) {              /* node to right */
      Lap->data[nz] = uconst;
      Lap->rowvals[nz++] = IDX(i+1,0);
    }

    /* dependence on v at this node */
    Lap->colptrs[IDX(i,1)] = nz;
    if (i>1) {                /* node to left */
      Lap->data[nz] = vconst;
      Lap->rowvals[nz++] = IDX(i-1,1);
    }

    Lap->data[nz] = vconst2;  /* self */
    Lap->rowvals[nz++] = IDX(i,1);

    if (i<N-2) {              /* node to right */
      Lap->data[nz] = vconst;
      Lap->rowvals[nz++] = IDX(i+1,1);
    }

    /* dependence on w at this node */
    Lap->colptrs[IDX(i,2)] = nz;
    if (i>1) {                /* node to left */
      Lap->data[nz] = wconst;
      Lap->rowvals[nz++] = IDX(i-1,2);
    }

    Lap->data[nz] = wconst2;  /* self */
    Lap->rowvals[nz++] = IDX(i,2);

    if (i<N-2) {              /* node to right */
      Lap->data[nz] = wconst;
      Lap->rowvals[nz++] = IDX(i+1,2);
    }

  }

  /* set last column to zero */
  Lap->colptrs[IDX(N-1,0)] = nz;
  Lap->colptrs[IDX(N-1,1)] = nz;
  Lap->colptrs[IDX(N-1,2)] = nz;
  
  /* end of data */
  Lap->colptrs[IDX(N-1,2)+1] = nz;

  return 0;
}



/* Routine to compute the Jacobian matrix from R(y) */
static int ReactionJac(N_Vector y, SlsMat Jac, UserData udata)
{

  /* shortcuts to number of intervals, background values */
  int N = udata->N;

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;

  /* set shortcuts */
  int i, nz=0;
  realtype u, v, w;
  realtype ep=udata->ep;

  /* clear out matrix */
  SlsSetToZero(Jac);

  /* set first matrix column to zero */
  Jac->colptrs[IDX(0,0)] = 0;
  Jac->colptrs[IDX(0,1)] = 0;
  Jac->colptrs[IDX(0,2)] = 0;
  
  /* iterate over interior nodes, filling in Jacobian entries */
  for (i=1; i<N-1; i++) {

    /* set nodal value shortcuts */
    u = Ydata[IDX(i,0)];
    v = Ydata[IDX(i,1)];
    w = Ydata[IDX(i,2)];

    /* dependence on u at this node */
    Jac->colptrs[IDX(i,0)] = nz;

    Jac->rowvals[nz] = IDX(i,0);        /* fu wrt u */
    Jac->data[nz++] = TWO*u*v - w - ONE;

    Jac->rowvals[nz] = IDX(i,1);        /* fv wrt u */
    Jac->data[nz++] = w - TWO*u*v;

    Jac->rowvals[nz] = IDX(i,2);        /* fw wrt u */
    Jac->data[nz++] = -w;

    /* dependence on v at this node */
    Jac->colptrs[IDX(i,1)] = nz;

    Jac->rowvals[nz] = IDX(i,0);        /* fu wrt v */
    Jac->data[nz++] = u*u;

    Jac->rowvals[nz] = IDX(i,1);        /* fv wrt v */
    Jac->data[nz++] = -u*u;

    /* dependence on w at this node */
    Jac->colptrs[IDX(i,2)] = nz;

    Jac->rowvals[nz] = IDX(i,0);        /* fu wrt w */
    Jac->data[nz++] = -u;

    Jac->rowvals[nz] = IDX(i,1);        /* fv wrt w */
    Jac->data[nz++] = u;

    Jac->rowvals[nz] = IDX(i,2);        /* fw wrt w */
    Jac->data[nz++] = -ONE/ep - u;

  }

  /* set last matrix column to zero */
  Jac->colptrs[IDX(N-1,0)] = nz;
  Jac->colptrs[IDX(N-1,1)] = nz;
  Jac->colptrs[IDX(N-1,2)] = nz;

  /* end of data */
  Jac->colptrs[IDX(N-1,2)+1] = nz;

  return 0;
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
