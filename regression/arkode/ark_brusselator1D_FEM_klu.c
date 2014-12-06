/*---------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2014, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 ----------------------------------------------------------------
 Example problem:
 
 The following test simulates a brusselator problem from chemical 
 kinetics.  This is a PDE system with 3 components, Y = [u,v,w], 
 satisfying the equations,
    u_t = du*u_xx + a - (w+1)*u + v*u^2
    v_t = dv*v_xx + w*u - v*u^2
    w_t = dw*w_xx + (b-w)/ep - w*u
 for t in [0, 80], x in [0, 1], with initial conditions
    u(0,x) =  a  + 0.1*sin(pi*x)
    v(0,x) = b/a + 0.1*sin(pi*x)
    w(0,x) =  b  + 0.1*sin(pi*x),
 and with stationary boundary conditions, i.e. 
    u_t(t,0) = u_t(t,1) = 0
    v_t(t,0) = v_t(t,1) = 0
    w_t(t,0) = w_t(t,1) = 0.
 
 Here, we use a piecewise linear Galerkin finite element 
 discretization in space, where all element-wise integrals are 
 computed using 3-node Gaussian quadrature (since we will have 
 quartic polynomials in the reaction terms for the u_t and v_t 
 equations (including the test function)).  The time derivative 
 terms for this system will include a mass matrix, giving rise 
 to an ODE system of the form
      M y_t = L y + R(y),
 where M is the block mass matrix for each component, L is 
 the block Laplace operator for each component, and R(y) is 
 a 3x3 block comprised of the nonlinear reaction terms for 
 each component.  Since it it highly inefficient to rewrite 
 this system as
      y_t = M^{-1}(L y + R(y)),
 we solve this system using ARKode, with a user-supplied mass
 matrix.  We therefore provide functions to evaluate the ODE RHS 
    f(t,y) = L y + R(y),
 its Jacobian
    J(t,y) = L + dR/dy,
 and the mass matrix, M.

 The number of spatial intervals N, the parameters a, b, du, dv, 
 dw and ep are provided in the input file 
 input_brusselator1D.txt.

 This program solves the problem with the DIRK method, using a
 Newton iteration with the ARKKLU sparse linear solver.

 100 outputs are printed at equal time intervals, and run 
 statistics are printed at the end.
---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode.h>
#include <nvector/nvector_serial.h>
#include <arkode/arkode_klu.h>
#include <sundials/sundials_types.h>

#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>


/* accessor macros between (x,v) location and 1D NVector array */
/* [variables are grouped according to spatial location] */
#define IDX(x,v) (3*(x)+v)

/* constants */
#define ZERO (RCONST(0.0))
#define ONE  (RCONST(1.0))
#define TWO  (RCONST(2.0))
#define HALF (RCONST(0.5))

/* Gaussian quadrature nodes, weights and formula (3 node, 7th-order accurate) */
#define X1(xl,xr)   (HALF*(xl+xr) - HALF*(xr-xl)*RCONST(0.774596669241483377035853079956))
#define X2(xl,xr)   (HALF*(xl+xr))
#define X3(xl,xr)   (HALF*(xl+xr) + HALF*(xr-xl)*RCONST(0.774596669241483377035853079956))
#define W1          (RCONST(0.55555555555555555555555555555556))
#define W2          (RCONST(0.88888888888888888888888888888889))
#define W3          (RCONST(0.55555555555555555555555555555556))
#define Quad(f1,f2,f3,xl,xr) (HALF*(xr-xl)*(W1*f1 + W2*f2 + W3*f3))

/* evaluation macros for variables, basis functions and basis derivatives */
#define ChiL(xl,xr,x) ((xr-x)/(xr-xl))
#define ChiR(xl,xr,x) ((x-xl)/(xr-xl))
#define ChiL_x(xl,xr) (ONE/(xl-xr))
#define ChiR_x(xl,xr) (ONE/(xr-xl))
#define Eval(ul,ur,xl,xr,x) (ul*ChiL(xl,xr,x) + ur*ChiR(xl,xr,x))
#define Eval_x(ul,ur,xl,xr) (ul*ChiL_x(xl,xr) + ur*ChiR_x(xl,xr))


/* user data structure */
typedef struct {  
  int N;         /* number of intervals     */
  realtype *x;   /* mesh node locations     */
  realtype a;    /* constant forcing on u   */
  realtype b;    /* steady-state value of w */
  realtype du;   /* diffusion coeff for u   */
  realtype dv;   /* diffusion coeff for v   */
  realtype dw;   /* diffusion coeff for w   */
  realtype ep;   /* stiffness parameter     */
  N_Vector tmp;  /* temporary vector        */
  SlsMat R;      /* temporary storage       */
} *UserData;


/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f_diff(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int f_rx(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int MassMatrix(realtype t, SlsMat M, void *user_data, 
		      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
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
int main() {

  /* general problem parameters */
  realtype T0 = RCONST(0.0);
  realtype Tf = RCONST(10.0);
  int Nt = 10;
  int Nvar = 3;
  UserData udata = NULL;
  realtype *data;
  int i, N;
  long int NEQ, NNZ;

  /* declare solver parameters */
  int flag, dense_order, imex, fixedpt;

  /* general problem variables */
  int idense;
  N_Vector y = NULL;
  N_Vector ytrue = NULL;
  N_Vector yerr = NULL;
  N_Vector umask = NULL;
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  void *arkode_mem = NULL;
  void *arktrue_mem = NULL;

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  udata->x = NULL;
  udata->tmp = NULL;
  if (check_flag((void *)udata, "malloc", 2)) return 1;

  /* read problem parameter and tolerances from input file:
     N - number of spatial intervals
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

  /* open solver diagnostics output file for writing */
  FILE *DFID;
  DFID=fopen("diags_ark_bruss1D_FEM_klu.txt","w");
  
  /* set total allocated vector length (N-1 intervals, Dirichlet end points) */
  NEQ = Nvar*udata->N;

  /* Initial problem output */
  printf("\n1D FEM Brusselator PDE test problem:\n");
  printf("    N = %i,  NEQ = %li\n", udata->N, NEQ);
  printf("    problem parameters:  a = %g,  b = %g,  ep = %g\n",
	 udata->a, udata->b, udata->ep);
  printf("    diffusion coefficients:  du = %g,  dv = %g,  dw = %g\n", 
	 udata->du, udata->dv, udata->dw);

  /* Create serial vectors of length NEQ for initial condition, etc. */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;
  ytrue = N_VNew_Serial(NEQ);
  if (check_flag((void *) ytrue, "N_VNew_Serial", 0)) return 1;
  yerr = N_VNew_Serial(NEQ);
  if (check_flag((void *) yerr, "N_VNew_Serial", 0)) return 1;


  /* allocate and set up spatial mesh; this [arbitrarily] clusters 
     more intervals near the end points of the interval */
  udata->x = (realtype *) malloc(N*sizeof(realtype));
  if (check_flag((void *)udata->x, "malloc", 2)) return 1;
  realtype h=10.0/(N-1);
  realtype z;
  for (i=0; i<N; i++) {
    z = -5.0 + h*i;
    udata->x[i] = 0.5/atan(5.0)*atan(z) + 0.5;
  }

  /* output mesh to disk */
  FID=fopen("bruss_FEM_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->x[i]);
  fclose(FID);

  /* allocate space for temporary N_Vector inside user data structure */
  udata->tmp = N_VNew_Serial(NEQ);
  if (check_flag((void *) udata->tmp, "N_VNew_Serial", 0)) return 1;


  /* Access data array for new NVector y */
  data = N_VGetArrayPointer(y);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;

  /* Set initial conditions into y, ytrue */
  realtype pi=RCONST(4.0)*atan(ONE);
  for (i=0; i<N; i++) {
    data[IDX(i,0)] =  a  + RCONST(0.1)*sin(pi*udata->x[i]);  /* u */
    data[IDX(i,1)] = b/a + RCONST(0.1)*sin(pi*udata->x[i]);  /* v */
    data[IDX(i,2)] =  b  + RCONST(0.1)*sin(pi*udata->x[i]);  /* w */
  }
  N_VScale(1.0, y, ytrue);

  
  /* Create serial vector masks for each solution component */
  umask = N_VNew_Serial(NEQ);
  if (check_flag((void *)umask, "N_VNew_Serial", 0)) return 1;
  vmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)vmask, "N_VNew_Serial", 0)) return 1;
  wmask = N_VNew_Serial(NEQ);
  if (check_flag((void *)wmask, "N_VNew_Serial", 0)) return 1;

  /* Set mask array values for each solution component */
  N_VConst(0.0, umask);
  data = N_VGetArrayPointer(umask);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,0)] = ONE;

  N_VConst(0.0, vmask);
  data = N_VGetArrayPointer(vmask);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
  for (i=0; i<N; i++)  data[IDX(i,1)] = ONE;

  N_VConst(0.0, wmask);
  data = N_VGetArrayPointer(wmask);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;
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
  flag = ARKodeSetMaxNumSteps(arktrue_mem, 10000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;

  /* Call ARKodeSStolerances to specify the scalar relative and 
     absolute tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;
  flag = ARKodeSStolerances(arktrue_mem, reltol2, abstol2);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Call ARKodeResStolerance to specify the scalar 
     absolute residual tolerance */
  flag = ARKodeResStolerance(arkode_mem, abstol);
  if (check_flag(&flag, "ARKodeResStolerance", 1)) return 1;
  flag = ARKodeResStolerance(arktrue_mem, abstol2);
  if (check_flag(&flag, "ARKodeResStolerance", 1)) return 1;

  /* Specify the KLU sparse linear solver */
  NNZ = 15*NEQ;
  flag = ARKKLU(arkode_mem, NEQ, NNZ);
  if (check_flag(&flag, "ARKKLU", 1)) return 1;
  flag = ARKKLU(arktrue_mem, NEQ, NNZ);
  if (check_flag(&flag, "ARKKLU", 1)) return 1;

  /* Set the Jacobian routine (user-supplied) */
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


  /* Specify the mass matrix linear solver */
  flag = ARKMassKLU(arkode_mem, NEQ, NNZ, MassMatrix);
  if (check_flag(&flag, "ARKMassKLU", 1)) return 1;
  flag = ARKMassKLU(arktrue_mem, NEQ, NNZ, MassMatrix);
  if (check_flag(&flag, "ARKMassKLU", 1)) return 1;


  /* Open output stream for results, access data arrays */
  FILE *UFID=fopen("bruss_FEM_u.txt","w");
  FILE *VFID=fopen("bruss_FEM_v.txt","w");
  FILE *WFID=fopen("bruss_FEM_w.txt","w");
  data = N_VGetArrayPointer(y);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;

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
  long int nst, nst_a, nfe, nfi, nsetups, nje, nni, ncfn;
  long int netf, nli, nlcf, nms, nMv;
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
  flag = ARKodeGetNumMassSolves(arkode_mem, &nms);
  check_flag(&flag, "ARKodeGetNumMassSolves", 1);
  flag = ARKodeGetNumMassMultiplies(arkode_mem, &nMv);
  check_flag(&flag, "ARKodeGetNumMassMultiplies", 1);
  flag = ARKSlsGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKSlsGetNumJacEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", 
	 nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total mass matrix solves = %li\n", nms);
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
  N_VDestroy_Serial(udata->tmp);
  free(udata->x);
  free(udata);

  /* Free integrator memory */
  ARKodeFree(&arkode_mem);
  ARKodeFree(&arktrue_mem);

  /* close solver diagnostics output file */
  fclose(DFID);

  return 0;
}


/*------------------------------
  Functions called by the solver
 *------------------------------*/


/* Routine to compute the ODE RHS function f(t,y), where system is of the form
        M y_t = f(t,y) := Ly + R(y) 
   This routine only computes the f(t,y), leaving (M y_t) alone. */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* local data */
  int ier;

  /* problem data */
  UserData udata = (UserData) user_data;

  /* clear out RHS (to be careful) */
  N_VConst(0.0, ydot);

  /* add reaction terms to RHS */
  ier = f_rx(t, y, ydot, user_data);
  if (ier != 0)  return ier;
  
  /* add diffusion terms to RHS */
  ier = f_diff(t, y, ydot, user_data);
  if (ier != 0)  return ier;
  
  return 0;
}


/* Routine to compute the explicit components of the ODE RHS function f(t,y). */
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* local data */
  int ier;

  /* problem data */
  UserData udata = (UserData) user_data;

  /* clear out RHS (to be careful) */
  N_VConst(0.0, ydot);

  /* add diffusion terms to RHS */
  ier = f_diff(t, y, ydot, user_data);
  if (ier != 0)  return ier;
  
  return 0;
}


/* Routine to compute the implicit components of the ODE RHS function f(t,y). */
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* local data */
  int ier;

  /* problem data */
  UserData udata = (UserData) user_data;

  /* clear out RHS (to be careful) */
  N_VConst(0.0, ydot);

  /* add reaction terms to RHS */
  ier = f_rx(t, y, ydot, user_data);
  if (ier != 0)  return ier;
  
  return 0;
}


/* Routine to compute the diffusion portion of the ODE RHS function f(t,y). */
static int f_diff(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of intervals, background values */
  int N = udata->N;
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;
  realtype *RHSdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *)RHSdata, "N_VGetArrayPointer", 0)) return 1;

  /* set shortcuts */
  long int i;
  realtype ul, ur, vl, vr, wl, wr;
  realtype u, v, w, xl, xr, f1;
  booleantype left, right;
  
  /* iterate over intervals, filling in residual function */
  for (i=0; i<N-1; i++) {

    /* set booleans to determine whether equations exist on the left/right */
    left  = (i==0)     ? FALSE : TRUE;
    right = (i==(N-2)) ? FALSE : TRUE;

    /* set nodal value shortcuts (interval index aligns with left node) */
    ul = Ydata[IDX(i,0)];
    vl = Ydata[IDX(i,1)];
    wl = Ydata[IDX(i,2)];
    ur = Ydata[IDX(i+1,0)];
    vr = Ydata[IDX(i+1,1)];
    wr = Ydata[IDX(i+1,2)];

    /* set mesh shortcuts */
    xl = udata->x[i];
    xr = udata->x[i+1];

    /* evaluate L*y on this subinterval
       NOTE: all f values are the same since constant on interval */
    /*    left test function */
    if (left) {
      /*  u */
      f1 = -du * Eval_x(ul,ur,xl,xr) * ChiL_x(xl,xr);
      RHSdata[IDX(i,0)] += Quad(f1,f1,f1,xl,xr);

      /*  v */
      f1 = -dv * Eval_x(vl,vr,xl,xr) * ChiL_x(xl,xr);
      RHSdata[IDX(i,1)] += Quad(f1,f1,f1,xl,xr);
      
      /*  w */
      f1 = -dw * Eval_x(wl,wr,xl,xr) * ChiL_x(xl,xr);
      RHSdata[IDX(i,2)] += Quad(f1,f1,f1,xl,xr);
    }
    /*    right test function */
    if (right) {
      /*  u */
      f1 = -du * Eval_x(ul,ur,xl,xr) * ChiR_x(xl,xr);
      RHSdata[IDX(i+1,0)] += Quad(f1,f1,f1,xl,xr);

      /*  v */
      f1 = -dv * Eval_x(vl,vr,xl,xr) * ChiR_x(xl,xr);
      RHSdata[IDX(i+1,1)] += Quad(f1,f1,f1,xl,xr);

      /*  w */
      f1 = -dw * Eval_x(wl,wr,xl,xr) * ChiR_x(xl,xr);
      RHSdata[IDX(i+1,2)] += Quad(f1,f1,f1,xl,xr);
    }
  }

  return 0;
}



/* Routine to compute the reaction portion of the ODE RHS function f(t,y). */
static int f_rx(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of intervals, background values */
  int N = udata->N;
  realtype a  = udata->a;
  realtype b  = udata->b;
  realtype ep = udata->ep;

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *)Ydata, "N_VGetArrayPointer", 0)) return 1;
  realtype *RHSdata = N_VGetArrayPointer(ydot);
  if (check_flag((void *)RHSdata, "N_VGetArrayPointer", 0)) return 1;

  /* set shortcuts */
  long int i;
  realtype ul, ur, vl, vr, wl, wr;
  realtype u, v, w, xl, xr, f1, f2, f3;
  booleantype left, right;
  
  /* iterate over intervals, filling in residual function */
  for (i=0; i<N-1; i++) {

    /* set booleans to determine whether equations exist on the left/right */
    left  = (i==0)     ? FALSE : TRUE;
    right = (i==(N-2)) ? FALSE : TRUE;

    /* set nodal value shortcuts (interval index aligns with left node) */
    ul = Ydata[IDX(i,0)];
    vl = Ydata[IDX(i,1)];
    wl = Ydata[IDX(i,2)];
    ur = Ydata[IDX(i+1,0)];
    vr = Ydata[IDX(i+1,1)];
    wr = Ydata[IDX(i+1,2)];

    /* set mesh shortcuts */
    xl = udata->x[i];
    xr = udata->x[i+1];

    /* evaluate R(y) on this subinterval */
    /*    left test function */
    if (left) {
      /*  u */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = (a - (w+ONE)*u + v*u*u) * ChiL(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = (a - (w+ONE)*u + v*u*u) * ChiL(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = (a - (w+ONE)*u + v*u*u) * ChiL(xl,xr,X3(xl,xr));
      RHSdata[IDX(i,0)] += Quad(f1,f2,f3,xl,xr);
    
      /*  v */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = (w*u - v*u*u) * ChiL(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = (w*u - v*u*u) * ChiL(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = (w*u - v*u*u) * ChiL(xl,xr,X3(xl,xr));
      RHSdata[IDX(i,1)] += Quad(f1,f2,f3,xl,xr);
    
      /*  w */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = ((b-w)/ep - w*u) * ChiL(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = ((b-w)/ep - w*u) * ChiL(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = ((b-w)/ep - w*u) * ChiL(xl,xr,X3(xl,xr));
      RHSdata[IDX(i,2)] += Quad(f1,f2,f3,xl,xr);
    }
    /*    right test function */
    if (right) {
      /*  u */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = (a - (w+ONE)*u + v*u*u) * ChiR(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = (a - (w+ONE)*u + v*u*u) * ChiR(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = (a - (w+ONE)*u + v*u*u) * ChiR(xl,xr,X3(xl,xr));
      RHSdata[IDX(i+1,0)] += Quad(f1,f2,f3,xl,xr);
    
      /*  v */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = (w*u - v*u*u) * ChiR(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = (w*u - v*u*u) * ChiR(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = (w*u - v*u*u) * ChiR(xl,xr,X3(xl,xr));
      RHSdata[IDX(i+1,1)] += Quad(f1,f2,f3,xl,xr);
    
      /*  w */
      u = Eval(ul,ur,xl,xr,X1(xl,xr));
      v = Eval(vl,vr,xl,xr,X1(xl,xr));
      w = Eval(wl,wr,xl,xr,X1(xl,xr));
      f1 = ((b-w)/ep - w*u) * ChiR(xl,xr,X1(xl,xr));
      u = Eval(ul,ur,xl,xr,X2(xl,xr));
      v = Eval(vl,vr,xl,xr,X2(xl,xr));
      w = Eval(wl,wr,xl,xr,X2(xl,xr));
      f2 = ((b-w)/ep - w*u) * ChiR(xl,xr,X2(xl,xr));
      u = Eval(ul,ur,xl,xr,X3(xl,xr));
      v = Eval(vl,vr,xl,xr,X3(xl,xr));
      w = Eval(wl,wr,xl,xr,X3(xl,xr));
      f3 = ((b-w)/ep - w*u) * ChiR(xl,xr,X3(xl,xr));
      RHSdata[IDX(i+1,2)] += Quad(f1,f2,f3,xl,xr);
    }
  }

  return 0;
}



/* Interface routine to compute the Jacobian of the full RHS function, f(y) */
static int Jac(realtype t, N_Vector y, N_Vector fy, 
               SlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

  /* temporary variables */
  int ier;
  UserData udata = (UserData) user_data;
  int N = udata->N;

  /* ensure that Jac is the correct size */
  if ((J->M != N*3) || (J->N != N*3)) {
    printf("Jacobian calculation error: matrix is the wrong size!\n");
    return 1;
  }
  
  /* Fill in the Laplace matrix */
  ier = LaplaceMatrix(J, udata);
  if (ier != 0) {
    fprintf(stderr,"Jac: error in filling Laplace matrix = %i\n",ier);
    return 1;
  }

  /* Create empty reaction Jacobian matrix (if not done already) */
  if (udata->R == NULL) {
    udata->R = NewSparseMat(J->M, J->N, J->NNZ);
    if (udata->R == NULL) {
      printf("Jac: error in allocating R matrix!\n");
      return 1;
    }
  }
      
  /* Add in the Jacobian of the reaction terms matrix */
  ier = ReactionJac(y, udata->R, udata);
  if (ier != 0) {
    fprintf(stderr,"Jac: error in filling reaction Jacobian = %i\n",ier);
    return 1;
  }

  /* printf("Laplace matrix:\n"); */
  /* PrintSparseMat(J); */

  /* printf("Reaction matrix:\n"); */
  /* PrintSparseMat(udata->R); */

  /* Add R to J */
  ier = SlsAddMat(J,udata->R);
  if (ier != 0) {
    printf("Jac: error in adding sparse matrices = %i!\n",ier);
    return 1;
  }

  return 0;
}



/* Interface routine to compute the Jacobian of the implicit RHS function, fi(y) */
static int JacI(realtype t, N_Vector y, N_Vector fy, 
		SlsMat J, void *user_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)  {

  /* temporary variables */
  int ier;
  UserData udata = (UserData) user_data;

  /* Add in the Jacobian of the reaction terms matrix */
  ier = ReactionJac(y, J, udata);
  if (ier != 0) {
    fprintf(stderr,"Jac: error in filling reaction Jacobian = %i\n",ier);
    return 1;
  }

  return 0;
}



/* Routine to compute the mass matrix multiplying y_t. */
static int MassMatrix(realtype t, SlsMat M, void *user_data, 
		      N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

  /* clear out mass matrix */
  SlsSetToZero(M);

  /* user data structure */
  UserData udata = (UserData) user_data;

  /* set shortcuts */
  int N = udata->N;
  int i, nz=0;

  /* iterate over columns, filling in matrix entries */
  realtype xl, xr, f1, f2, f3, dtmp;
  for (i=0; i<N; i++) {

    /* dependence on u at this node */
    M->colptrs[IDX(i,0)] = nz;

    /*    left u trial function */
    if (i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      M->data[nz] = Quad(f1,f2,f3,xl,xr);
      M->rowvals[nz++] = IDX(i-1,0);
    }
    /*    this u trial function */
    dtmp = ZERO;
    if (i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiL(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiL(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiL(xl,xr,X3(xl,xr));
      dtmp += Quad(f1,f2,f3,xl,xr);
    }
    if (i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      f1 = ChiR(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiR(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiR(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      dtmp += Quad(f1,f2,f3,xl,xr);
    }
    M->data[nz] = dtmp;
    M->rowvals[nz++] = IDX(i,0);
    /*    right u trial function */
    if (i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      M->data[nz] = Quad(f1,f2,f3,xl,xr);
      M->rowvals[nz++] = IDX(i+1,0);
    }


    /* dependence on v at this node */
    M->colptrs[IDX(i,1)] = nz;

    /*    left v trial function */
    if (i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      M->data[nz] = Quad(f1,f2,f3,xl,xr);
      M->rowvals[nz++] = IDX(i-1,1);
    }
    /*    this v trial function */
    dtmp = ZERO;
    if (i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiL(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiL(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiL(xl,xr,X3(xl,xr));
      dtmp += Quad(f1,f2,f3,xl,xr);
    }
    if (i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      f1 = ChiR(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiR(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiR(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      dtmp += Quad(f1,f2,f3,xl,xr);
    }
    M->data[nz] = dtmp;
    M->rowvals[nz++] = IDX(i,1);
    /*    right v trial function */
    if (i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      M->data[nz] = Quad(f1,f2,f3,xl,xr);
      M->rowvals[nz++] = IDX(i+1,1);
    }


    /* dependence on w at this node */
    M->colptrs[IDX(i,2)] = nz;

    /*    left w trial function */
    if (i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      M->data[nz] = Quad(f1,f2,f3,xl,xr);
      M->rowvals[nz++] = IDX(i-1,2);
    }
    /*    this w trial function */
    dtmp = ZERO;
    if (i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiL(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiL(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiL(xl,xr,X3(xl,xr));
      dtmp += Quad(f1,f2,f3,xl,xr);
    }
    if (i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      f1 = ChiR(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiR(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiR(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      dtmp += Quad(f1,f2,f3,xl,xr);
    }
    M->data[nz] = dtmp;
    M->rowvals[nz++] = IDX(i,2);
    /*    right w trial function */
    if (i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
      f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
      f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
      M->data[nz] = Quad(f1,f2,f3,xl,xr);
      M->rowvals[nz++] = IDX(i+1,2);
    }

  }

  /* signal end of data */
  M->colptrs[IDX(N-1,2)+1] = nz;

  return 0;
}





/*-------------------------------
 * Private helper functions
 *-------------------------------*/



/* Routine to compute the Laplace matrix */
static int LaplaceMatrix(SlsMat L, UserData udata)
{

  /* set shortcuts, local variables */
  int N = udata->N;
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;
  int i, nz=0;
  realtype xl, xr, dtmp;
  
  /* clear out matrix */
  SlsSetToZero(L);

  /* iterate over columns, filling in Laplace matrix entries */
  for (i=0; i<N; i++) {

    /* dependence on u at this node */
    L->colptrs[IDX(i,0)] = nz;

    if (i>1) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      L->data[nz] = (-du) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i-1,0);
    }
    if (i<N-1 && i>0) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      L->data[nz] = (-du) * Quad(ONE,ONE,ONE,xl,xr) * ChiR_x(xl,xr) * ChiR_x(xl,xr);
      xl = udata->x[i];
      xr = udata->x[i+1];
      L->data[nz] += (-du) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiL_x(xl,xr);
      L->rowvals[nz++] = IDX(i,0);
    }
    if (i<N-2) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      L->data[nz] = (-du) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i+1,0);
    }

    /* dependence on v at this node */
    L->colptrs[IDX(i,1)] = nz;

    if (i>1) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      L->data[nz] = (-dv) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i-1,1);
    }
    if (i>0 && i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      L->data[nz] = (-dv) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiL_x(xl,xr);
      xl = udata->x[i-1];
      xr = udata->x[i];
      L->data[nz] += (-dv) * Quad(ONE,ONE,ONE,xl,xr) * ChiR_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i,1);
    }
    if (i<N-2) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      L->data[nz] = (-dv) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i+1,1);
    }

    /* dependence on w at this node */
    L->colptrs[IDX(i,2)] = nz;

    if (i>1) {
      xl = udata->x[i-1];
      xr = udata->x[i];
      L->data[nz] = (-dw) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i-1,2);
    }
    if (i>0 && i<N-1) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      L->data[nz] = (-dw) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiL_x(xl,xr);
      xl = udata->x[i-1];
      xr = udata->x[i];
      L->data[nz] += (-dw) * Quad(ONE,ONE,ONE,xl,xr) * ChiR_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i,2);
    }
    if (i<N-2) {
      xl = udata->x[i];
      xr = udata->x[i+1];
      L->data[nz] = (-dw) * Quad(ONE,ONE,ONE,xl,xr) * ChiL_x(xl,xr) * ChiR_x(xl,xr);
      L->rowvals[nz++] = IDX(i+1,2);
    }

  }

  /* signal end of data */
  L->colptrs[IDX(N-1,2)+1] = nz;

  return 0;
}



/* Routine to compute the Jacobian matrix from R(y) */
static int ReactionJac(N_Vector y, SlsMat Jac, UserData udata)
{
  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;

  /* set shortcuts, local variables */
  int N = udata->N;
  int i, nz=0;
  realtype ep = udata->ep;
  realtype ul, uc, ur, vl, vc, vr, wl, wc, wr;
  realtype u1l, u2l, u3l, v1l, v2l, v3l, w1l, w2l, w3l;
  realtype u1r, u2r, u3r, v1r, v2r, v3r, w1r, w2r, w3r;
  realtype xl, xc, xr, df1, df2, df3;
  realtype dQdf1l, dQdf2l, dQdf3l, ChiL1l, ChiL2l, ChiL3l, ChiR1l, ChiR2l, ChiR3l;
  realtype dQdf1r, dQdf2r, dQdf3r, ChiL1r, ChiL2r, ChiL3r, ChiR1r, ChiR2r, ChiR3r;
  
  /* clear out matrix */
  SlsSetToZero(Jac);

  /* iterate over columns, filling in reaction Jacobian */
  for (i=0; i<N; i++) {

    /* set mesh shortcuts */
    if (i>0)
      xl = udata->x[i-1];
    xc = udata->x[i];
    if (i<N-1)
      xr = udata->x[i+1];

    /* set nodal value shortcuts */
    if (i>0) {
      ul = Ydata[IDX(i-1,0)];
      vl = Ydata[IDX(i-1,1)];
      wl = Ydata[IDX(i-1,2)];
    }
    uc = Ydata[IDX(i,0)];
    vc = Ydata[IDX(i,1)];
    wc = Ydata[IDX(i,2)];
    if (i<N-1) {
      ur = Ydata[IDX(i+1,0)];
      vr = Ydata[IDX(i+1,1)];
      wr = Ydata[IDX(i+1,2)];
    }
    if (i>0) {
      u1l = Eval(ul,uc,xl,xc,X1(xl,xc));
      v1l = Eval(vl,vc,xl,xc,X1(xl,xc));
      w1l = Eval(wl,wc,xl,xc,X1(xl,xc));
      u2l = Eval(ul,uc,xl,xc,X2(xl,xc));
      v2l = Eval(vl,vc,xl,xc,X2(xl,xc));
      w2l = Eval(wl,wc,xl,xc,X2(xl,xc));
      u3l = Eval(ul,uc,xl,xc,X3(xl,xc));
      v3l = Eval(vl,vc,xl,xc,X3(xl,xc));
      w3l = Eval(wl,wc,xl,xc,X3(xl,xc));
    }
    if (i<N-1) {
      u1r = Eval(uc,ur,xc,xr,X1(xc,xr));
      v1r = Eval(vc,vr,xc,xr,X1(xc,xr));
      w1r = Eval(wc,wr,xc,xr,X1(xc,xr));
      u2r = Eval(uc,ur,xc,xr,X2(xc,xr));
      v2r = Eval(vc,vr,xc,xr,X2(xc,xr));
      w2r = Eval(wc,wr,xc,xr,X2(xc,xr));
      u3r = Eval(uc,ur,xc,xr,X3(xc,xr));
      v3r = Eval(vc,vr,xc,xr,X3(xc,xr));
      w3r = Eval(wc,wr,xc,xr,X3(xc,xr));
    }

    /* set partial derivative shortcuts */
    if (i>0) {
      dQdf1l = Quad(ONE, ZERO, ZERO, xl, xc);
      dQdf2l = Quad(ZERO, ONE, ZERO, xl, xc);
      dQdf3l = Quad(ZERO, ZERO, ONE, xl, xc);
      ChiL1l = ChiL(xl,xc,X1(xl,xc));
      ChiL2l = ChiL(xl,xc,X2(xl,xc));
      ChiL3l = ChiL(xl,xc,X3(xl,xc));
      ChiR1l = ChiR(xl,xc,X1(xl,xc));
      ChiR2l = ChiR(xl,xc,X2(xl,xc));
      ChiR3l = ChiR(xl,xc,X3(xl,xc));
    }
    if (i<N-1) {
      dQdf1r = Quad(ONE, ZERO, ZERO, xc, xr);
      dQdf2r = Quad(ZERO, ONE, ZERO, xc, xr);
      dQdf3r = Quad(ZERO, ZERO, ONE, xc, xr);
      ChiL1r = ChiL(xc,xr,X1(xc,xr));
      ChiL2r = ChiL(xc,xr,X2(xc,xr));
      ChiL3r = ChiL(xc,xr,X3(xc,xr));
      ChiR1r = ChiR(xc,xr,X1(xc,xr));
      ChiR2r = ChiR(xc,xr,X2(xc,xr));
      ChiR3r = ChiR(xc,xr,X3(xc,xr));
    }


    /*** evaluate dR/dy at this node ***/


    /* dependence on u at this node */
    Jac->colptrs[IDX(i,0)] = nz;

    if (i>1) {
      /*  dR_ul/duc */
      df1 = (-(w1l+ONE) + TWO*v1l*u1l) * ChiL1l * ChiR1l;
      df2 = (-(w2l+ONE) + TWO*v2l*u2l) * ChiL2l * ChiR2l;
      df3 = (-(w3l+ONE) + TWO*v3l*u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,0);

      /*  dR_vl/duc */
      df1 = (w1l - TWO*v1l*u1l) * ChiL1l * ChiR1l;
      df2 = (w2l - TWO*v2l*u2l) * ChiL2l * ChiR2l;
      df3 = (w3l - TWO*v3l*u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,1);

      /*  dR_wl/duc */
      df1 = (-w1l) * ChiL1l * ChiR1l;
      df2 = (-w2l) * ChiL2l * ChiR2l;
      df3 = (-w3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,2);
    }
    if (i>0 && i<N-1) {
      /*  dR_uc/duc */
      df1 = (-(w1r+ONE) + TWO*v1r*u1r) * ChiL1r * ChiL1r;
      df2 = (-(w2r+ONE) + TWO*v2r*u2r) * ChiL2r * ChiL2r;
      df3 = (-(w3r+ONE) + TWO*v3r*u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;

      df1 = (-(w1l+ONE) + TWO*v1l*u1l) * ChiR1l * ChiR1l;
      df2 = (-(w2l+ONE) + TWO*v2l*u2l) * ChiR2l * ChiR2l;
      df3 = (-(w3l+ONE) + TWO*v3l*u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] += dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i,0);

      /*  dR_vc/duc */
      df1 = (w1l - TWO*v1l*u1l) * ChiR1l * ChiR1l;
      df2 = (w2l - TWO*v2l*u2l) * ChiR2l * ChiR2l;
      df3 = (w3l - TWO*v3l*u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;

      df1 = (w1r - TWO*v1r*u1r) * ChiL1r * ChiL1r;
      df2 = (w2r - TWO*v2r*u2r) * ChiL2r * ChiL2r;
      df3 = (w3r - TWO*v3r*u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] += dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i,1);

      /*  dR_wc/duc */
      df1 = (-w1r) * ChiL1r * ChiL1r;
      df2 = (-w2r) * ChiL2r * ChiL2r;
      df3 = (-w3r) * ChiL3r * ChiL3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;

      df1 = (-w1l) * ChiR1l * ChiR1l;
      df2 = (-w2l) * ChiR2l * ChiR2l;
      df3 = (-w3l) * ChiR3l * ChiR3l;
      Jac->data[nz] += dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i,2);
    }
    if (i<N-2) {
      /*  dR_ur/duc */
      df1 = (-(w1r+ONE) + TWO*v1r*u1r) * ChiL1r * ChiR1r;
      df2 = (-(w2r+ONE) + TWO*v2r*u2r) * ChiL2r * ChiR2r;
      df3 = (-(w3r+ONE) + TWO*v3r*u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,0);

      /*  dR_vr/duc */
      df1 = (w1r - TWO*v1r*u1r) * ChiL1r * ChiR1r;
      df2 = (w2r - TWO*v2r*u2r) * ChiL2r * ChiR2r;
      df3 = (w3r - TWO*v3r*u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,1);

      /*  dR_wr/duc */
      df1 = (-w1r) * ChiL1r * ChiR1r;
      df2 = (-w2r) * ChiL2r * ChiR2r;
      df3 = (-w3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,2);
    }


    /* dependence on v at this node */
    Jac->colptrs[IDX(i,1)] = nz;

    if (i>1) {
      /*  dR_ul/dvc */
      df1 = (u1l*u1l) * ChiL1l * ChiR1l;
      df2 = (u2l*u2l) * ChiL2l * ChiR2l;
      df3 = (u3l*u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,0);

      /*  dR_vl/dvc */
      df1 = (-u1l*u1l) * ChiL1l * ChiR1l;
      df2 = (-u2l*u2l) * ChiL2l * ChiR2l;
      df3 = (-u3l*u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,1);
    }
    if (i>0 && i<N-1) {
      /*  dR_uc/dvc */
      df1 = (u1l*u1l) * ChiR1l * ChiR1l;
      df2 = (u2l*u2l) * ChiR2l * ChiR2l;
      df3 = (u3l*u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;

      df1 = (u1r*u1r) * ChiL1r * ChiL1r;
      df2 = (u2r*u2r) * ChiL2r * ChiL2r;
      df3 = (u3r*u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] += dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i,0);

      /*  dR_vc/dvc */
      df1 = (-u1l*u1l) * ChiR1l * ChiR1l;
      df2 = (-u2l*u2l) * ChiR2l * ChiR2l;
      df3 = (-u3l*u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;

      df1 = (-u1r*u1r) * ChiL1r * ChiL1r;
      df2 = (-u2r*u2r) * ChiL2r * ChiL2r;
      df3 = (-u3r*u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] += dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i,1);
    }
    if (i<N-2) {
      /*  dR_ur/dvc */
      df1 = (u1r*u1r) * ChiL1r * ChiR1r;
      df2 = (u2r*u2r) * ChiL2r * ChiR2r;
      df3 = (u3r*u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,0);

      /*  dR_vr/dvc */
      df1 = (-u1r*u1r) * ChiL1r * ChiR1r;
      df2 = (-u2r*u2r) * ChiL2r * ChiR2r;
      df3 = (-u3r*u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,1);
    }


    /* dependence on w at this node */
    Jac->colptrs[IDX(i,2)] = nz;

    if (i>1) {
      /*  dR_ul/dwc */
      df1 = (-u1l) * ChiL1l * ChiR1l;
      df2 = (-u2l) * ChiL2l * ChiR2l;
      df3 = (-u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,0);

      /*  dR_vl/dwc */
      df1 = (u1l) * ChiL1l * ChiR1l;
      df2 = (u2l) * ChiL2l * ChiR2l;
      df3 = (u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,1);

      /*  dR_wl/dwc */
      df1 = (-ONE/ep - u1l) * ChiL1l * ChiR1l;
      df2 = (-ONE/ep - u2l) * ChiL2l * ChiR2l;
      df3 = (-ONE/ep - u3l) * ChiL3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;
      Jac->rowvals[nz++] = IDX(i-1,2);
    }
    if (i>0 && i<N-1) {
      /*  dR_uc/dwc */
      df1 = (-u1l) * ChiR1l * ChiR1l;
      df2 = (-u2l) * ChiR2l * ChiR2l;
      df3 = (-u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;

      df1 = (-u1r) * ChiL1r * ChiL1r;
      df2 = (-u2r) * ChiL2r * ChiL2r;
      df3 = (-u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] += dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i,0);

      /*  dR_vc/dwc */
      df1 = (u1l) * ChiR1l * ChiR1l;
      df2 = (u2l) * ChiR2l * ChiR2l;
      df3 = (u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;

      df1 = (u1r) * ChiL1r * ChiL1r;
      df2 = (u2r) * ChiL2r * ChiL2r;
      df3 = (u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] += dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i,1);

      /*  dR_wc/dwc */
      df1 = (-ONE/ep - u1l) * ChiR1l * ChiR1l;
      df2 = (-ONE/ep - u2l) * ChiR2l * ChiR2l;
      df3 = (-ONE/ep - u3l) * ChiR3l * ChiR3l;
      Jac->data[nz] = dQdf1l*df1 + dQdf2l*df2 + dQdf3l*df3;

      df1 = (-ONE/ep - u1r) * ChiL1r * ChiL1r;
      df2 = (-ONE/ep - u2r) * ChiL2r * ChiL2r;
      df3 = (-ONE/ep - u3r) * ChiL3r * ChiL3r;
      Jac->data[nz] += dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i,2);
    }
    if (i<N-2) {
      /*  dR_ur/dwc */
      df1 = (-u1r) * ChiL1r * ChiR1r;
      df2 = (-u2r) * ChiL2r * ChiR2r;
      df3 = (-u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,0);

      /*  dR_vr/dwc */
      df1 = (u1r) * ChiL1r * ChiR1r;
      df2 = (u2r) * ChiL2r * ChiR2r;
      df3 = (u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,1);

      /*  dR_wr/dwc */
      df1 = (-ONE/ep - u1r) * ChiL1r * ChiR1r;
      df2 = (-ONE/ep - u2r) * ChiL2r * ChiR2r;
      df3 = (-ONE/ep - u3r) * ChiL3r * ChiR3r;
      Jac->data[nz] = dQdf1r*df1 + dQdf2r*df2 + dQdf3r*df3;
      Jac->rowvals[nz++] = IDX(i+1,2);
    }

  }

  /* signal end of data */
  Jac->colptrs[IDX(N-1,2)+1] = nz;

  return 0;
}



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


/*---- end of file ----*/
