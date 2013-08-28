/*---------------------------------------------------------------
 $Revision: $
 $Date: $
-----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
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
 equations (including the test function)).  While this system of 
 equations does not have any constraint equations, the time 
 derivative terms will include a mass matrix, giving rise to an 
 ODE system of the form
      M y_t = L y + R(y),
 where M is the 3x3 block mass matrix for each component, L is 
 the 3x3 block Laplace operator for each component, and R(y) is 
 comprised of the nonlinear reaction terms for each component.  
 Since it it highly inefficient to rewrite this system as
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

 We use a vector-valued absolute tolerance, where the values are 
 set as the input scalar value multiplied by the width of the 
 support for the corresponding basis function.  On a uniform mesh 
 this would result in a constant set of values, but on a 
 non-uniform mesh this spreads these weights in an integral 
 sense.
 
 This program solves the problem with the DIRK method, using a
 Newton iteration with the ARKBAND band linear solver.

 100 outputs are printed at equal time intervals, and run 
 statistics are printed at the end.
---------------------------------------------------------------*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode.h>
#include <nvector/nvector_serial.h>
#include <arkode/arkode_band.h>
#include <arkode/arkode_spgmr.h>
#include <arkode/arkode_spbcgs.h>
#include <arkode/arkode_sptfqmr.h>
#include <arkode/arkode_pcg.h>
#include <sundials/sundials_types.h>

#define USE_ITERATIVE
#define USE_SPGMR
/* #define USE_SPBCG */
/* #define USE_SPTFQMR */

#define MASS_USE_ITERATIVE
#define MASS_USE_SPGMR
/* #define MASS_USE_PCG */
/* #define MASS_USE_SPBCG */
/* #define MASS_USE_SPTFQMR */

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
  long int N;    /* number of intervals     */
  realtype *x;   /* mesh node locations     */
  realtype a;    /* constant forcing on u   */
  realtype b;    /* steady-state value of w */
  realtype du;   /* diffusion coeff for u   */
  realtype dv;   /* diffusion coeff for v   */
  realtype dw;   /* diffusion coeff for w   */
  realtype ep;   /* stiffness parameter     */
  N_Vector tmp;  /* temporary vector        */
} *UserData;


/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int MassMatrix(long int N, long int mu, long int ml, realtype t, 
		      DlsMat M, void *user_data, N_Vector tmp1, 
		      N_Vector tmp2, N_Vector tmp3);
static int Jac(long int N, long int mu, long int ml,
               realtype t, N_Vector y, N_Vector fy, 
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int JacI(long int N, long int mu, long int ml,
               realtype t, N_Vector y, N_Vector fy, 
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int MassTimes(N_Vector v, N_Vector Mv, 
		     realtype t, void *user_data);
static int JacV(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
		N_Vector fy, void *user_data, N_Vector tmp1);
static int JacVI(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
		 N_Vector fy, void *user_data, N_Vector tmp1);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);
static int LaplaceMatrix(DlsMat Jac, UserData udata);
static int ReactionJac(N_Vector y, DlsMat Jac, UserData udata);
static int LaplaceProduct(N_Vector v, N_Vector Lv, UserData udata);
static int ReactionProduct(N_Vector v, N_Vector Jv, 
			   N_Vector y, UserData udata);

/* Parameter input helper function */
int init_from_file(void *ark_mem, char *fname, ARKRhsFn f, 
		   ARKRhsFn fe, ARKRhsFn fi, realtype T0, 
		   N_Vector y0, int *ImEx, int *dorder, 
		   int *fxpt, realtype *RTol, realtype *ATol);


/* Main Program */
int main() {

  /* general problem parameters */
  realtype T0 = RCONST(0.0);
  /* realtype Tf = RCONST(10.0); */
  realtype Tf = RCONST(0.0001);
  int Nt = 10;
  int Nvar = 3;
  UserData udata = NULL;
  realtype *data;
  long int N, NEQ, i;

  /* declare solver parameters */
  int flag, dense_order, imex, fixedpt;

  /* general problem variables */
  int idense;
  N_Vector y = NULL;
  N_Vector ytrue = NULL;
  N_Vector yerr = NULL;
  N_Vector avtol = NULL;
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
  flag = fscanf(FID,"  N = %li\n", &N);
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
  DFID=fopen("diags_ark_bruss1D_FEM.txt","w");
  
  /* set total allocated vector length (N-1 intervals, Dirichlet end points) */
  NEQ = Nvar*udata->N;

  /* Initial problem output */
  printf("\n1D FEM Brusselator PDE test problem:\n");
  printf("    N = %li,  NEQ = %li\n", udata->N, NEQ);
  printf("    problem parameters:  a = %g,  b = %g,  ep = %g\n",
	 udata->a, udata->b, udata->ep);
  printf("    diffusion coefficients:  du = %g,  dv = %g,  dw = %g\n", 
	 udata->du, udata->dv, udata->dw);

  /* Create serial vectors of length NEQ for initial condition & abs tol */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;
  ytrue = N_VNew_Serial(NEQ);
  if (check_flag((void *) ytrue, "N_VNew_Serial", 0)) return 1;
  yerr = N_VNew_Serial(NEQ);
  if (check_flag((void *) yerr, "N_VNew_Serial", 0)) return 1;
  avtol = N_VNew_Serial(NEQ);
  if (check_flag((void *) avtol, "N_VNew_Serial", 0)) return 1;


  /* allocate and set up spatial mesh; this [arbitrarily] clusters 
     more intervals near the endpoints of the interval */
  udata->x = (realtype *) malloc(N*sizeof(realtype));
  if (check_flag((void *)udata->x, "malloc", 2)) return 1;
  realtype z, h = ONE/(N-1);
  for (i=0; i<N; i++) {
    z = h*i - 0.5;
    udata->x[i] = RCONST(4.0)*z*z*z + RCONST(0.5);
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
  realtype pi = RCONST(4.0)*atan(ONE);
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

  /* Access data array for absolute tolerance NVector avtol */
  data = N_VGetArrayPointer(avtol);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;

  /* Set support widths into avtol */
  i = 0; {
    data[IDX(i,0)] = abstol * (udata->x[i+1] - udata->x[i]);  /* u */
    data[IDX(i,1)] = abstol * (udata->x[i+1] - udata->x[i]);  /* v */
    data[IDX(i,2)] = abstol * (udata->x[i+1] - udata->x[i]);  /* w */
  }
  for (i=1; i<N-1; i++) {
    data[IDX(i,0)] = abstol * (udata->x[i+1] - udata->x[i-1]);  /* u */
    data[IDX(i,1)] = abstol * (udata->x[i+1] - udata->x[i-1]);  /* v */
    data[IDX(i,2)] = abstol * (udata->x[i+1] - udata->x[i-1]);  /* w */
  }
  i=N-1; {
    data[IDX(i,0)] = abstol * (udata->x[i] - udata->x[i-1]);  /* u */
    data[IDX(i,1)] = abstol * (udata->x[i] - udata->x[i-1]);  /* v */
    data[IDX(i,2)] = abstol * (udata->x[i] - udata->x[i-1]);  /* w */
  }


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

  /* Call ARKodeSVtolerances to specify the scalar relative and 
     vector absolute tolerances */
  flag = ARKodeSVtolerances(arkode_mem, reltol, avtol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;
  flag = ARKodeSVtolerances(arktrue_mem, reltol2, avtol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Free avtol */
  N_VDestroy_Serial(avtol);

  /* Specify the linear solver */
#ifdef USE_ITERATIVE
#ifdef USE_SPGMR
  flag = ARKSpgmr(arkode_mem, 0, 500);
  if (check_flag(&flag, "ARKSpgmr", 1)) return 1;
  flag = ARKSpgmr(arktrue_mem, 0, 500);
  if (check_flag(&flag, "ARKSpgmr", 1)) return 1;
#endif
#ifdef USE_SPBCG
  flag = ARKSpbcg(arkode_mem, 0, 500);
  if (check_flag(&flag, "ARKSpbcg", 1)) return 1;
  flag = ARKSpbcg(arktrue_mem, 0, 500);
  if (check_flag(&flag, "ARKSpbcg", 1)) return 1;
#endif
#ifdef USE_SPTFQMR
  flag = ARKSptfqmr(arkode_mem, 0, 500);
  if (check_flag(&flag, "ARKSptfqmr", 1)) return 1;
  flag = ARKSptfqmr(arktrue_mem, 0, 500);
  if (check_flag(&flag, "ARKSptfqmr", 1)) return 1;
#endif
#else
  flag = ARKBand(arkode_mem, NEQ, 5, 5);
  if (check_flag(&flag, "ARKBand", 1)) return 1;
  flag = ARKBand(arktrue_mem, NEQ, 5, 5);
  if (check_flag(&flag, "ARKBand", 1)) return 1;
#endif

  /* Set the Jacobian routine to Jac (user-supplied) */
#ifdef USE_ITERATIVE
  switch (imex) {
  case 0:         /* purely implicit */
    flag = ARKSpilsSetJacTimesVecFn(arkode_mem, JacV);   break;
  case 1:         /* purely explicit */
    break;
  default:        /* imex */
    flag = ARKSpilsSetJacTimesVecFn(arkode_mem, JacVI);  break;
  }
  if (check_flag(&flag, "ARKSpilsSetJacTimesVecFn", 1)) return 1;
  flag = ARKSpilsSetJacTimesVecFn(arktrue_mem, JacV);
  if (check_flag(&flag, "ARKSpilsSetJacTimesVecFn", 1)) return 1;
#else
  switch (imex) {
  case 0:         /* purely implicit */
    flag = ARKDlsSetBandJacFn(arkode_mem, Jac);   break;
  case 1:         /* purely explicit */
    break;
  default:        /* imex */
    flag = ARKDlsSetBandJacFn(arkode_mem, JacI);  break;
  }
  if (check_flag(&flag, "ARKDlsSetBandJacFn", 1)) return 1;
  flag = ARKDlsSetBandJacFn(arktrue_mem, Jac);
  if (check_flag(&flag, "ARKDlsSetBandJacFn", 1)) return 1;
#endif


  /* Specify the mass matrix linear solver */
#ifdef MASS_USE_ITERATIVE
#ifdef MASS_USE_SPGMR
  flag = ARKMassSpgmr(arkode_mem, 0, 500, MassTimes, (void *) udata);
  if (check_flag(&flag, "ARKMassSpgmr", 1)) return 1;
  flag = ARKMassSpgmr(arktrue_mem, 0, 500, MassTimes, (void *) udata);
  if (check_flag(&flag, "ARKMassSpgmr", 1)) return 1;
#endif
#ifdef MASS_USE_PCG
  flag = ARKMassPcg(arkode_mem, 0, 500, MassTimes, (void *) udata);
  if (check_flag(&flag, "ARKMassPcg", 1)) return 1;
  flag = ARKMassPcg(arktrue_mem, 0, 500, MassTimes, (void *) udata);
  if (check_flag(&flag, "ARKMassPcg", 1)) return 1;
#endif
#ifdef MASS_USE_SPBCG
  flag = ARKMassSpbcg(arkode_mem, 0, 500, MassTimes, (void *) udata);
  if (check_flag(&flag, "ARKMassSpbcg", 1)) return 1;
  flag = ARKMassSpbcg(arktrue_mem, 0, 500, MassTimes, (void *) udata);
  if (check_flag(&flag, "ARKMassSpbcg", 1)) return 1;
#endif
#ifdef MASS_USE_SPTFQMR
  flag = ARKMassSptfqmr(arkode_mem, 0, 500, MassTimes, (void *) udata);
  if (check_flag(&flag, "ARKMassSptfqmr", 1)) return 1;
  flag = ARKMassSptfqmr(arktrue_mem, 0, 500, MassTimes, (void *) udata);
  if (check_flag(&flag, "ARKMassSptfqmr", 1)) return 1;
#endif
#else
  flag = ARKMassBand(arkode_mem, NEQ, 5, 5, MassMatrix);
  if (check_flag(&flag, "ARKMassBand", 1)) return 1;
  flag = ARKMassBand(arktrue_mem, NEQ, 5, 5, MassMatrix);
  if (check_flag(&flag, "ARKMassBand", 1)) return 1;
#endif


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
  realtype t = T0;
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
  long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn;
  long int netf, nli, nlcf, nJv, nms, nmi, nmcf, nMv;
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
#ifdef USE_ITERATIVE
  flag = ARKSpilsGetNumLinIters(arkode_mem, &nli);
  check_flag(&flag, "ARKSpilsGetNumLinIters", 1);
  flag = ARKSpilsGetNumConvFails(arkode_mem, &nlcf);
  check_flag(&flag, "ARKSpilsGetNumConvFails", 1);
  flag = ARKSpilsGetNumJtimesEvals(arkode_mem, &nJv);
  check_flag(&flag, "ARKSpilsGetNumJtimesEvals", 1);
#else
  flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKDlsGetNumJacEvals", 1);
  flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKDlsGetNumRhsEvals", 1);
#endif
#ifdef MASS_USE_ITERATIVE
  flag = ARKSpilsGetNumMassIters(arkode_mem, &nmi);
  check_flag(&flag, "ARKSpilsGetNumMassIters", 1);
  flag = ARKSpilsGetNumMassConvFails(arkode_mem, &nmcf);
  check_flag(&flag, "ARKSpilsGetNumMassConvFails", 1);
#endif

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %li (attempted = %li)\n", 
	 nst, nst_a);
  printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
  printf("   Total mass matrix solves = %li\n", nms);
  printf("   Total linear solver setups = %li\n", nsetups);
#ifdef USE_ITERATIVE
  printf("   Total linear iterations = %li\n", nli);
  printf("   Total linear convergence failures = %li\n", nlcf);
  printf("   Total J*v evaluations = %li\n", nJv);
  printf("   Total M*v evaluations = %li\n", nMv);
#else
  printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
#endif
#ifdef MASS_USE_ITERATIVE
  printf("   Total mass matrix solver iters = %li\n", nmi);
  printf("   Total mass matrix solver sonvergence failures = %li\n", nmcf);
#endif
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

  /* fill reaction portion of RHS array */
  ier = fi(t, y, ydot, user_data);
  if (ier != 0)  return ier;
  
  /* fill diffusion portion of RHS array */
  ier = fe(t, y, udata->tmp, user_data);
  if (ier != 0)  return ier;
  
  /* add components together */
  N_VLinearSum( 1.0, udata->tmp, 1.0, ydot, ydot );
  
  return 0;
}


/* Routine to compute the diffusion portion of the ODE RHS function f(t,y). */
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* clear out RHS (to be careful) */
  N_VConst(0.0, ydot);

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of intervals, background values */
  long int N  = udata->N;
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
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

  /* clear out RHS (to be careful) */
  N_VConst(0.0, ydot);

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of intervals, background values */
  long int N  = udata->N;
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
static int Jac(long int N, long int mu, long int ml,
               realtype t, N_Vector y, N_Vector fy, 
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

  /* temporary variables */
  int ier;
  UserData udata = (UserData) user_data;

  /* ensure that matrix and user data structure match */
  if (N != 3*udata->N) {
    fprintf(stderr,"Jac error: dimension mismatch, %li != %li\n",
	    N, 3*udata->N);
    return 1;
  }

  /* clear out Jacobian matrix */
  SetToZero(J);

  /* add Laplace matrix to J */
  ier = LaplaceMatrix(J, udata);
  if (ier != 0) {
    fprintf(stderr,"Jac: error in filling Laplace matrix = %i\n",ier);
    return 1;
  }

  /* add reaction Jacobian to J */
  ier = ReactionJac(y, J, udata);
  if (ier != 0) {
    fprintf(stderr,"Jac: error in filling reaction Jacobian = %i\n",ier);
    return 1;
  }

  return 0;
}



/* Interface routine to compute the Jacobian of the implicit RHS function, fi(y) */
static int JacI(long int N, long int mu, long int ml,
               realtype t, N_Vector y, N_Vector fy, 
               DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)  {

  /* temporary variables */
  int ier;
  UserData udata = (UserData) user_data;

  /* ensure that matrix and user data structure match */
  if (N != 3*udata->N) {
    fprintf(stderr,"JacI error: dimension mismatch, %li != %li\n",
	    N, 3*udata->N);
    return 1;
  }

  /* clear out Jacobian matrix */
  SetToZero(J);

  /* add reaction Jacobian to J */
  ier = ReactionJac(y, J, udata);
  if (ier != 0) {
    fprintf(stderr,"Jac: error in filling reaction Jacobian = %i\n",ier);
    return 1;
  }

  return 0;
}



/* Routine to compute the mass matrix multiplying y_t. */
static int MassMatrix(long int N, long int mu, long int ml, realtype t, 
		      DlsMat M, void *user_data, N_Vector tmp1, 
		      N_Vector tmp2, N_Vector tmp3) {

  /* clear out mass matrix */
  SetToZero(M);

  /* user data structure */
  UserData udata = (UserData) user_data;

  /* ensure that matrix and user data structure match */
  if (N != 3*udata->N) {
    fprintf(stderr,"MassMatrix error: dimension mismatch, %li != %li\n",
	    N, 3*udata->N);
    return 1;
  }

  /* iterate over intervals, filling in matrix entries */
  long int i;
  realtype xl, xr, f1, f2, f3, qd;
  for (i=0; i<udata->N-1; i++) {

    /* set mesh shortcuts */
    xl = udata->x[i];
    xr = udata->x[i+1];

    /*    left basis and test functions */
    f1 = ChiL(xl,xr,X1(xl,xr)) * ChiL(xl,xr,X1(xl,xr));
    f2 = ChiL(xl,xr,X2(xl,xr)) * ChiL(xl,xr,X2(xl,xr));
    f3 = ChiL(xl,xr,X3(xl,xr)) * ChiL(xl,xr,X3(xl,xr));
    qd = Quad(f1,f2,f3,xl,xr);
    BAND_ELEM(M,IDX(i,0),IDX(i,0)) += qd;
    BAND_ELEM(M,IDX(i,1),IDX(i,1)) += qd;
    BAND_ELEM(M,IDX(i,2),IDX(i,2)) += qd;

    /*    right basis and test functions */
    f1 = ChiR(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
    f2 = ChiR(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
    f3 = ChiR(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
    qd = Quad(f1,f2,f3,xl,xr);
    BAND_ELEM(M,IDX(i+1,0),IDX(i+1,0)) += qd;
    BAND_ELEM(M,IDX(i+1,1),IDX(i+1,1)) += qd;
    BAND_ELEM(M,IDX(i+1,2),IDX(i+1,2)) += qd;

    /*    left and right basis and test functions */
    f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
    f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
    f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
    qd = Quad(f1,f2,f3,xl,xr);
    BAND_ELEM(M,IDX(i,0),IDX(i+1,0)) += qd;
    BAND_ELEM(M,IDX(i,1),IDX(i+1,1)) += qd;
    BAND_ELEM(M,IDX(i,2),IDX(i+1,2)) += qd;
    BAND_ELEM(M,IDX(i+1,0),IDX(i,0)) += qd;
    BAND_ELEM(M,IDX(i+1,1),IDX(i,1)) += qd;
    BAND_ELEM(M,IDX(i+1,2),IDX(i,2)) += qd;

  }

  return 0;
}



/* Interface routine to compute the Jacobian-vector product 
   for the full RHS function, f(y) */
static int JacV(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
		N_Vector fy, void *user_data, N_Vector tmp1) {

  /* temporary variables */
  int ier;
  UserData udata = (UserData) user_data;
  
  /* clear out result */
  N_VConst(0.0, Jv);

  /* add Laplace matrix product to Jv */
  ier = LaplaceProduct(v, Jv, udata);
  if (ier != 0) {
    fprintf(stderr,"JacV: error in filling Laplace matrix product = %i\n",ier);
    return 1;
  }

  /* add reaction Jacobian to J */
  ier = ReactionProduct(v, Jv, y, udata);
  if (ier != 0) {
    fprintf(stderr,"JacV: error in filling reaction Jacobian product = %i\n",ier);
    return 1;
  }

  return 0;
}



/* Interface routine to compute the Jacobian-vector product 
   for the implicit RHS function, fi(y) */
static int JacVI(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
		 N_Vector fy, void *user_data, N_Vector tmp1) {

  /* temporary variables */
  int ier;
  UserData udata = (UserData) user_data;
  
  /* clear out result */
  N_VConst(0.0, Jv);

  /* add reaction Jacobian to J */
  ier = ReactionProduct(v, Jv, y, udata);
  if (ier != 0) {
    fprintf(stderr,"JacVI: error in filling reaction Jacobian product = %i\n",ier);
    return 1;
  }

  return 0;
}



/* Interface routine to compute the mass-matrix-vector product */
static int MassTimes(N_Vector v, N_Vector Mv, 
		     realtype t, void *user_data) {

  /* clear out result */
  N_VConst(0.0, Mv);

  /* user data structure */
  UserData udata = (UserData) user_data;

  /* shortcut to number of nodes */
  long int N = udata->N;

  /* access v and Mv data arrays */
  realtype *Vdata = N_VGetArrayPointer(v);
  if (check_flag((void *) Vdata, "N_VGetArrayPointer", 0)) return 1;
  realtype *MVdata = N_VGetArrayPointer(Mv);
  if (check_flag((void *) MVdata, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over intervals, filling in matrix entries */
  long int i;
  realtype xl, xr, f1, f2, f3, qd;
  for (i=0; i<N-1; i++) {

    /* set mesh shortcuts */
    xl = udata->x[i];
    xr = udata->x[i+1];

    /*    left basis and test functions */
    f1 = ChiL(xl,xr,X1(xl,xr)) * ChiL(xl,xr,X1(xl,xr));
    f2 = ChiL(xl,xr,X2(xl,xr)) * ChiL(xl,xr,X2(xl,xr));
    f3 = ChiL(xl,xr,X3(xl,xr)) * ChiL(xl,xr,X3(xl,xr));
    qd = Quad(f1,f2,f3,xl,xr);
    MVdata[IDX(i,0)] += qd * Vdata[IDX(i,0)];
    MVdata[IDX(i,1)] += qd * Vdata[IDX(i,1)];
    MVdata[IDX(i,2)] += qd * Vdata[IDX(i,2)];

    /*    right basis and test functions */
    f1 = ChiR(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
    f2 = ChiR(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
    f3 = ChiR(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
    qd = Quad(f1,f2,f3,xl,xr);
    MVdata[IDX(i+1,0)] += qd * Vdata[IDX(i+1,0)];
    MVdata[IDX(i+1,1)] += qd * Vdata[IDX(i+1,1)];
    MVdata[IDX(i+1,2)] += qd * Vdata[IDX(i+1,2)];

    /*    left and right basis and test functions */
    f1 = ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
    f2 = ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
    f3 = ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
    qd = Quad(f1,f2,f3,xl,xr);
    MVdata[IDX(i,0)] += qd * Vdata[IDX(i+1,0)];
    MVdata[IDX(i,1)] += qd * Vdata[IDX(i+1,1)];
    MVdata[IDX(i,2)] += qd * Vdata[IDX(i+1,2)];
    MVdata[IDX(i+1,0)] += qd * Vdata[IDX(i,0)];
    MVdata[IDX(i+1,1)] += qd * Vdata[IDX(i,1)];
    MVdata[IDX(i+1,2)] += qd * Vdata[IDX(i,2)];

  }

  return 0;
}



/*-------------------------------
 * Private helper functions
 *-------------------------------*/



/* Routine to compute the Laplace matrix */
static int LaplaceMatrix(DlsMat L, UserData udata)
{
  /* set shortcuts, local variables */
  long int N = udata->N;
  long int i;
  realtype xl, xr;
  booleantype left, right;
  
  /* iterate over intervals, filling in Laplace matrix entries */
  for (i=0; i<N-1; i++) {

    /* set mesh shortcuts */
    xl = udata->x[i];
    xr = udata->x[i+1];

    /* set booleans to determine whether equations exist on the left/right */
    left  = (i==0)     ? FALSE : TRUE;
    right = (i==(N-2)) ? FALSE : TRUE;

    /*    left basis and test functions */
    if (left) {
      BAND_ELEM(L,IDX(i,0),IDX(i,0)) += ChiL_x(xl,xr) * ChiL_x(xl,xr);
      BAND_ELEM(L,IDX(i,1),IDX(i,1)) += ChiL_x(xl,xr) * ChiL_x(xl,xr);
      BAND_ELEM(L,IDX(i,2),IDX(i,2)) += ChiL_x(xl,xr) * ChiL_x(xl,xr);
    }

    /*    right basis and test functions */
    if (right) {
      BAND_ELEM(L,IDX(i+1,0),IDX(i+1,0)) += ChiR_x(xl,xr) * ChiR_x(xl,xr);
      BAND_ELEM(L,IDX(i+1,1),IDX(i+1,1)) += ChiR_x(xl,xr) * ChiR_x(xl,xr);
      BAND_ELEM(L,IDX(i+1,2),IDX(i+1,2)) += ChiR_x(xl,xr) * ChiR_x(xl,xr);
    }

    /*    left and right basis and test functions */
    if (left && right) {
      BAND_ELEM(L,IDX(i,0),IDX(i+1,0)) += ChiL_x(xl,xr) * ChiR_x(xl,xr);
      BAND_ELEM(L,IDX(i,1),IDX(i+1,1)) += ChiL_x(xl,xr) * ChiR_x(xl,xr);
      BAND_ELEM(L,IDX(i,2),IDX(i+1,2)) += ChiL_x(xl,xr) * ChiR_x(xl,xr);
      BAND_ELEM(L,IDX(i+1,0),IDX(i,0)) += ChiL_x(xl,xr) * ChiR_x(xl,xr);
      BAND_ELEM(L,IDX(i+1,1),IDX(i,1)) += ChiL_x(xl,xr) * ChiR_x(xl,xr);
      BAND_ELEM(L,IDX(i+1,2),IDX(i,2)) += ChiL_x(xl,xr) * ChiR_x(xl,xr);
    }

  }

  return 0;
}



/* Routine to compute the Jacobian matrix from R(y) */
static int ReactionJac(N_Vector y, DlsMat Jac, UserData udata)
{

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;

  /* set shortcuts, local variables */
  long int N = udata->N;
  long int i;
  realtype ep = udata->ep;
  realtype ul, ur, vl, vr, wl, wr;
  realtype u1, u2, u3, v1, v2, v3, w1, w2, w3, chi1, chi2, chi3;
  realtype u, v, w, xl, xr, f1, f2, f3;
  booleantype left, right;
  
  /* iterate over intervals, filling in reaction Jacobian */
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

    /* evaluate dR/dy on this subinterval */

    /* left basis function, evaluated at all quadrature nodes */
    if (left) {
      u1 = ul * ChiL(xl,xr,X1(xl,xr));
      v1 = vl * ChiL(xl,xr,X1(xl,xr));
      w1 = wl * ChiL(xl,xr,X1(xl,xr));
      u2 = ul * ChiL(xl,xr,X2(xl,xr));
      v2 = vl * ChiL(xl,xr,X2(xl,xr));
      w2 = wl * ChiL(xl,xr,X2(xl,xr));
      u3 = ul * ChiL(xl,xr,X3(xl,xr));
      v3 = vl * ChiL(xl,xr,X3(xl,xr));
      w3 = wl * ChiL(xl,xr,X3(xl,xr));
    } else {
      u1 = u2 = u3 = v1 = v2 = v3 = w1 = w2 = w3 = 0.0;
    }

    /* left trial function, evaluated at all quadrature nodes */
    if (left) {
      chi1 = ChiL(xl,xr,X1(xl,xr));
      chi2 = ChiL(xl,xr,X2(xl,xr));
      chi3 = ChiL(xl,xr,X3(xl,xr));
    } else {
      chi1 = chi2 = chi3 = 0.0;
    }

    /* left basis, left trial Jacobian entries */
    /*   dR_u/du */
    f1 = (TWO*u1*v1-(w1+ONE)) * chi1;
    f2 = (TWO*u2*v2-(w2+ONE)) * chi2;
    f3 = (TWO*u3*v3-(w3+ONE)) * chi3;
    BAND_ELEM(Jac,IDX(i,0),IDX(i,0)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_u/dv */
    f1 = u1 * u1 * chi1;
    f2 = u2 * u2 * chi2;
    f3 = u3 * u3 * chi3;
    BAND_ELEM(Jac,IDX(i,0),IDX(i,1)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_u/dw */
    f1 = -u1 * chi1;
    f2 = -u2 * chi2;
    f3 = -u3 * chi3;
    BAND_ELEM(Jac,IDX(i,0),IDX(i,2)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_v/du */
    f1 = (w1 - TWO*u1*v1) * chi1;
    f2 = (w2 - TWO*u2*v2) * chi2;
    f3 = (w3 - TWO*u3*v3) * chi3;
    BAND_ELEM(Jac,IDX(i,1),IDX(i,0)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_v/dv */
    f1 = -u1 * u1 * chi1;
    f2 = -u2 * u2 * chi2;
    f3 = -u3 * u3 * chi3;
    BAND_ELEM(Jac,IDX(i,1),IDX(i,1)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_v/dw */
    f1 = u1 * chi1;
    f2 = u2 * chi2;
    f3 = u3 * chi3;
    BAND_ELEM(Jac,IDX(i,1),IDX(i,2)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_w/du */
    f1 = -w1 * chi1;
    f2 = -w2 * chi2;
    f3 = -w3 * chi3;
    BAND_ELEM(Jac,IDX(i,2),IDX(i,0)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_w/dw */
    f1 = (-ONE/ep - u1) * chi1;
    f2 = (-ONE/ep - u2) * chi2;
    f3 = (-ONE/ep - u3) * chi3;
    BAND_ELEM(Jac,IDX(i,2),IDX(i,2)) += Quad(f1,f2,f3,xl,xr);


    /* right trial function, evaluated at all quadrature nodes */
    if (right) {
      chi1 = ChiR(xl,xr,X1(xl,xr));
      chi2 = ChiR(xl,xr,X2(xl,xr));
      chi3 = ChiR(xl,xr,X3(xl,xr));
    } else {
      chi1 = chi2 = chi3 = 0.0;
    }

    /* left basis, right trial Jacobian entries */
    /*   dR_u/du */
    f1 = (TWO*u1*v1-(w1+ONE)) * chi1;
    f2 = (TWO*u2*v2-(w2+ONE)) * chi2;
    f3 = (TWO*u3*v3-(w3+ONE)) * chi3;
    BAND_ELEM(Jac,IDX(i+1,0),IDX(i,0)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_u/dv */
    f1 = u1 * u1 * chi1;
    f2 = u2 * u2 * chi2;
    f3 = u3 * u3 * chi3;
    BAND_ELEM(Jac,IDX(i+1,0),IDX(i,1)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_u/dw */
    f1 = -u1 * chi1;
    f2 = -u2 * chi2;
    f3 = -u3 * chi3;
    BAND_ELEM(Jac,IDX(i+1,0),IDX(i,2)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_v/du */
    f1 = (w1 - TWO*u1*v1) * chi1;
    f2 = (w2 - TWO*u2*v2) * chi2;
    f3 = (w3 - TWO*u3*v3) * chi3;
    BAND_ELEM(Jac,IDX(i+1,1),IDX(i,0)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_v/dv */
    f1 = -u1 * u1 * chi1;
    f2 = -u2 * u2 * chi2;
    f3 = -u3 * u3 * chi3;
    BAND_ELEM(Jac,IDX(i+1,1),IDX(i,1)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_v/dw */
    f1 = u1 * chi1;
    f2 = u2 * chi2;
    f3 = u3 * chi3;
    BAND_ELEM(Jac,IDX(i+1,1),IDX(i,2)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_w/du */
    f1 = -w1 * chi1;
    f2 = -w2 * chi2;
    f3 = -w3 * chi3;
    BAND_ELEM(Jac,IDX(i+1,2),IDX(i,0)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_w/dw */
    f1 = (-ONE/ep - u1) * chi1;
    f2 = (-ONE/ep - u2) * chi2;
    f3 = (-ONE/ep - u3) * chi3;
    BAND_ELEM(Jac,IDX(i+1,2),IDX(i,2)) += Quad(f1,f2,f3,xl,xr);


    /* right basis function, evaluated at all quadrature nodes */
    if (left) {
      u1 = ul * ChiR(xl,xr,X1(xl,xr));
      v1 = vl * ChiR(xl,xr,X1(xl,xr));
      w1 = wl * ChiR(xl,xr,X1(xl,xr));
      u2 = ul * ChiR(xl,xr,X2(xl,xr));
      v2 = vl * ChiR(xl,xr,X2(xl,xr));
      w2 = wl * ChiR(xl,xr,X2(xl,xr));
      u3 = ul * ChiR(xl,xr,X3(xl,xr));
      v3 = vl * ChiR(xl,xr,X3(xl,xr));
      w3 = wl * ChiR(xl,xr,X3(xl,xr));
    } else {
      u1 = u2 = u3 = v1 = v2 = v3 = w1 = w2 = w3 = 0.0;
    }

    /* right basis, right trial Jacobian entries */
    /*   dR_u/du */
    f1 = (TWO*u1*v1-(w1+ONE)) * chi1;
    f2 = (TWO*u2*v2-(w2+ONE)) * chi2;
    f3 = (TWO*u3*v3-(w3+ONE)) * chi3;
    BAND_ELEM(Jac,IDX(i+1,0),IDX(i+1,0)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_u/dv */
    f1 = u1 * u1 * chi1;
    f2 = u2 * u2 * chi2;
    f3 = u3 * u3 * chi3;
    BAND_ELEM(Jac,IDX(i+1,0),IDX(i+1,1)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_u/dw */
    f1 = -u1 * chi1;
    f2 = -u2 * chi2;
    f3 = -u3 * chi3;
    BAND_ELEM(Jac,IDX(i+1,0),IDX(i+1,2)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_v/du */
    f1 = (w1 - TWO*u1*v1) * chi1;
    f2 = (w2 - TWO*u2*v2) * chi2;
    f3 = (w3 - TWO*u3*v3) * chi3;
    BAND_ELEM(Jac,IDX(i+1,1),IDX(i+1,0)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_v/dv */
    f1 = -u1 * u1 * chi1;
    f2 = -u2 * u2 * chi2;
    f3 = -u3 * u3 * chi3;
    BAND_ELEM(Jac,IDX(i+1,1),IDX(i+1,1)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_v/dw */
    f1 = u1 * chi1;
    f2 = u2 * chi2;
    f3 = u3 * chi3;
    BAND_ELEM(Jac,IDX(i+1,1),IDX(i+1,2)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_w/du */
    f1 = -w1 * chi1;
    f2 = -w2 * chi2;
    f3 = -w3 * chi3;
    BAND_ELEM(Jac,IDX(i+1,2),IDX(i+1,0)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_w/dw */
    f1 = (-ONE/ep - u1) * chi1;
    f2 = (-ONE/ep - u2) * chi2;
    f3 = (-ONE/ep - u3) * chi3;
    BAND_ELEM(Jac,IDX(i+1,2),IDX(i+1,2)) += Quad(f1,f2,f3,xl,xr);


    /* left trial function, evaluated at all quadrature nodes */
    if (left) {
      chi1 = ChiL(xl,xr,X1(xl,xr));
      chi2 = ChiL(xl,xr,X2(xl,xr));
      chi3 = ChiL(xl,xr,X3(xl,xr));
    } else {
      chi1 = chi2 = chi3 = 0.0;
    }

    /* right basis, left trial Jacobian entries */
    /*   dR_u/du */
    f1 = (TWO*u1*v1-(w1+ONE)) * chi1;
    f2 = (TWO*u2*v2-(w2+ONE)) * chi2;
    f3 = (TWO*u3*v3-(w3+ONE)) * chi3;
    BAND_ELEM(Jac,IDX(i,0),IDX(i+1,0)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_u/dv */
    f1 = u1 * u1 * chi1;
    f2 = u2 * u2 * chi2;
    f3 = u3 * u3 * chi3;
    BAND_ELEM(Jac,IDX(i,0),IDX(i+1,1)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_u/dw */
    f1 = -u1 * chi1;
    f2 = -u2 * chi2;
    f3 = -u3 * chi3;
    BAND_ELEM(Jac,IDX(i,0),IDX(i+1,2)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_v/du */
    f1 = (w1 - TWO*u1*v1) * chi1;
    f2 = (w2 - TWO*u2*v2) * chi2;
    f3 = (w3 - TWO*u3*v3) * chi3;
    BAND_ELEM(Jac,IDX(i,1),IDX(i+1,0)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_v/dv */
    f1 = -u1 * u1 * chi1;
    f2 = -u2 * u2 * chi2;
    f3 = -u3 * u3 * chi3;
    BAND_ELEM(Jac,IDX(i,1),IDX(i+1,1)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_v/dw */
    f1 = u1 * chi1;
    f2 = u2 * chi2;
    f3 = u3 * chi3;
    BAND_ELEM(Jac,IDX(i,1),IDX(i+1,2)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_w/du */
    f1 = -w1 * chi1;
    f2 = -w2 * chi2;
    f3 = -w3 * chi3;
    BAND_ELEM(Jac,IDX(i,2),IDX(i+1,0)) += Quad(f1,f2,f3,xl,xr);
    /*   dR_w/dw */
    f1 = (-ONE/ep - u1) * chi1;
    f2 = (-ONE/ep - u2) * chi2;
    f3 = (-ONE/ep - u3) * chi3;
    BAND_ELEM(Jac,IDX(i,2),IDX(i+1,2)) += Quad(f1,f2,f3,xl,xr);

  }

  return 0;
}



/* Routine to compute the product of the Laplace matrix and a given vector */
static int LaplaceProduct(N_Vector v, N_Vector Lv, UserData udata)
{
  /* set shortcuts, local variables */
  long int N = udata->N;
  long int i;
  realtype xl, xr;
  booleantype left, right;
  
  /* access v and Lv data arrays */
  realtype *Vdata = N_VGetArrayPointer(v);
  if (check_flag((void *) Vdata, "N_VGetArrayPointer", 0)) return 1;
  realtype *LVdata = N_VGetArrayPointer(Lv);
  if (check_flag((void *) LVdata, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over intervals, filling in result */
  for (i=0; i<N-1; i++) {

    /* set mesh shortcuts */
    xl = udata->x[i];
    xr = udata->x[i+1];

    /* set booleans to determine whether equations exist on the left/right */
    left  = (i==0)     ? FALSE : TRUE;
    right = (i==(N-2)) ? FALSE : TRUE;

    /*    left basis and test functions */
    if (left) {
      LVdata[IDX(i,0)] += ChiL_x(xl,xr) * ChiL_x(xl,xr) * Vdata[IDX(i,0)];
      LVdata[IDX(i,1)] += ChiL_x(xl,xr) * ChiL_x(xl,xr) * Vdata[IDX(i,1)];
      LVdata[IDX(i,2)] += ChiL_x(xl,xr) * ChiL_x(xl,xr) * Vdata[IDX(i,2)];
    }

    /*    right basis and test functions */
    if (right) {
      LVdata[IDX(i+1,0)] += ChiR_x(xl,xr) * ChiR_x(xl,xr) * Vdata[IDX(i+1,0)];
      LVdata[IDX(i+1,1)] += ChiR_x(xl,xr) * ChiR_x(xl,xr) * Vdata[IDX(i+1,1)];
      LVdata[IDX(i+1,2)] += ChiR_x(xl,xr) * ChiR_x(xl,xr) * Vdata[IDX(i+1,2)];
    }

    /*    left and right basis and test functions */
    if (left && right) {
      LVdata[IDX(i,0)] += ChiL_x(xl,xr) * ChiR_x(xl,xr) * Vdata[IDX(i+1,0)];
      LVdata[IDX(i,1)] += ChiL_x(xl,xr) * ChiR_x(xl,xr) * Vdata[IDX(i+1,1)];
      LVdata[IDX(i,2)] += ChiL_x(xl,xr) * ChiR_x(xl,xr) * Vdata[IDX(i+1,2)];
      LVdata[IDX(i+1,0)] += ChiL_x(xl,xr) * ChiR_x(xl,xr) * Vdata[IDX(i,0)];
      LVdata[IDX(i+1,1)] += ChiL_x(xl,xr) * ChiR_x(xl,xr) * Vdata[IDX(i,1)];
      LVdata[IDX(i+1,2)] += ChiL_x(xl,xr) * ChiR_x(xl,xr) * Vdata[IDX(i,2)];
    }

  }

  return 0;
}



/* Routine to compute the product of the Jacobian matrix from R(y), times a given vector */
static int ReactionProduct(N_Vector V, N_Vector Jv, N_Vector y, UserData udata)
{

  /* set shortcuts, local variables */
  long int N = udata->N;
  long int i;
  realtype ep = udata->ep;
  realtype ul, ur, vl, vr, wl, wr;
  realtype u1, u2, u3, v1, v2, v3, w1, w2, w3, chi1, chi2, chi3;
  realtype u, v, w, xl, xr, f1, f2, f3;
  booleantype left, right;

  /* access v and Lv data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;
  realtype *Vdata = N_VGetArrayPointer(V);
  if (check_flag((void *) Vdata, "N_VGetArrayPointer", 0)) return 1;
  realtype *JVdata = N_VGetArrayPointer(Jv);
  if (check_flag((void *) JVdata, "N_VGetArrayPointer", 0)) return 1;

  /* iterate over intervals, filling in reaction Jacobian */
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

    /* evaluate dR/dy on this subinterval */

    /* left basis function, evaluated at all quadrature nodes */
    if (left) {
      u1 = ul * ChiL(xl,xr,X1(xl,xr));
      v1 = vl * ChiL(xl,xr,X1(xl,xr));
      w1 = wl * ChiL(xl,xr,X1(xl,xr));
      u2 = ul * ChiL(xl,xr,X2(xl,xr));
      v2 = vl * ChiL(xl,xr,X2(xl,xr));
      w2 = wl * ChiL(xl,xr,X2(xl,xr));
      u3 = ul * ChiL(xl,xr,X3(xl,xr));
      v3 = vl * ChiL(xl,xr,X3(xl,xr));
      w3 = wl * ChiL(xl,xr,X3(xl,xr));
    } else {
      u1 = u2 = u3 = v1 = v2 = v3 = w1 = w2 = w3 = 0.0;
    }

    /* left trial function, evaluated at all quadrature nodes */
    if (left) {
      chi1 = ChiL(xl,xr,X1(xl,xr));
      chi2 = ChiL(xl,xr,X2(xl,xr));
      chi3 = ChiL(xl,xr,X3(xl,xr));
    } else {
      chi1 = chi2 = chi3 = 0.0;
    }

    /* left basis, left trial Jacobian entries */
    /*   dR_u/du */
    f1 = (TWO*u1*v1-(w1+ONE)) * chi1;
    f2 = (TWO*u2*v2-(w2+ONE)) * chi2;
    f3 = (TWO*u3*v3-(w3+ONE)) * chi3;
    JVdata[IDX(i,0)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,0)];
    /*   dR_u/dv */
    f1 = u1 * u1 * chi1;
    f2 = u2 * u2 * chi2;
    f3 = u3 * u3 * chi3;
    JVdata[IDX(i,0)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,1)];
    /*   dR_u/dw */
    f1 = -u1 * chi1;
    f2 = -u2 * chi2;
    f3 = -u3 * chi3;
    JVdata[IDX(i,0)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,2)];
    /*   dR_v/du */
    f1 = (w1 - TWO*u1*v1) * chi1;
    f2 = (w2 - TWO*u2*v2) * chi2;
    f3 = (w3 - TWO*u3*v3) * chi3;
    JVdata[IDX(i,1)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,0)];
    /*   dR_v/dv */
    f1 = -u1 * u1 * chi1;
    f2 = -u2 * u2 * chi2;
    f3 = -u3 * u3 * chi3;
    JVdata[IDX(i,1)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,1)];
    /*   dR_v/dw */
    f1 = u1 * chi1;
    f2 = u2 * chi2;
    f3 = u3 * chi3;
    JVdata[IDX(i,1)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,2)];
    /*   dR_w/du */
    f1 = -w1 * chi1;
    f2 = -w2 * chi2;
    f3 = -w3 * chi3;
    JVdata[IDX(i,2)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,0)];
    /*   dR_w/dw */
    f1 = (-ONE/ep - u1) * chi1;
    f2 = (-ONE/ep - u2) * chi2;
    f3 = (-ONE/ep - u3) * chi3;
    JVdata[IDX(i,2)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,2)];


    /* right trial function, evaluated at all quadrature nodes */
    if (right) {
      chi1 = ChiR(xl,xr,X1(xl,xr));
      chi2 = ChiR(xl,xr,X2(xl,xr));
      chi3 = ChiR(xl,xr,X3(xl,xr));
    } else {
      chi1 = chi2 = chi3 = 0.0;
    }

    /* left basis, right trial Jacobian entries */
    /*   dR_u/du */
    f1 = (TWO*u1*v1-(w1+ONE)) * chi1;
    f2 = (TWO*u2*v2-(w2+ONE)) * chi2;
    f3 = (TWO*u3*v3-(w3+ONE)) * chi3;
    JVdata[IDX(i+1,0)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,0)];
    /*   dR_u/dv */
    f1 = u1 * u1 * chi1;
    f2 = u2 * u2 * chi2;
    f3 = u3 * u3 * chi3;
    JVdata[IDX(i+1,0)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,1)];
    /*   dR_u/dw */
    f1 = -u1 * chi1;
    f2 = -u2 * chi2;
    f3 = -u3 * chi3;
    JVdata[IDX(i+1,0)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,2)];
    /*   dR_v/du */
    f1 = (w1 - TWO*u1*v1) * chi1;
    f2 = (w2 - TWO*u2*v2) * chi2;
    f3 = (w3 - TWO*u3*v3) * chi3;
    JVdata[IDX(i+1,1)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,0)];
    /*   dR_v/dv */
    f1 = -u1 * u1 * chi1;
    f2 = -u2 * u2 * chi2;
    f3 = -u3 * u3 * chi3;
    JVdata[IDX(i+1,1)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,1)];
    /*   dR_v/dw */
    f1 = u1 * chi1;
    f2 = u2 * chi2;
    f3 = u3 * chi3;
    JVdata[IDX(i+1,1)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,2)];
    /*   dR_w/du */
    f1 = -w1 * chi1;
    f2 = -w2 * chi2;
    f3 = -w3 * chi3;
    JVdata[IDX(i+1,2)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,0)];
    /*   dR_w/dw */
    f1 = (-ONE/ep - u1) * chi1;
    f2 = (-ONE/ep - u2) * chi2;
    f3 = (-ONE/ep - u3) * chi3;
    JVdata[IDX(i+1,2)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i,2)];


    /* right basis function, evaluated at all quadrature nodes */
    if (left) {
      u1 = ul * ChiR(xl,xr,X1(xl,xr));
      v1 = vl * ChiR(xl,xr,X1(xl,xr));
      w1 = wl * ChiR(xl,xr,X1(xl,xr));
      u2 = ul * ChiR(xl,xr,X2(xl,xr));
      v2 = vl * ChiR(xl,xr,X2(xl,xr));
      w2 = wl * ChiR(xl,xr,X2(xl,xr));
      u3 = ul * ChiR(xl,xr,X3(xl,xr));
      v3 = vl * ChiR(xl,xr,X3(xl,xr));
      w3 = wl * ChiR(xl,xr,X3(xl,xr));
    } else {
      u1 = u2 = u3 = v1 = v2 = v3 = w1 = w2 = w3 = 0.0;
    }

    /* right basis, right trial Jacobian entries */
    /*   dR_u/du */
    f1 = (TWO*u1*v1-(w1+ONE)) * chi1;
    f2 = (TWO*u2*v2-(w2+ONE)) * chi2;
    f3 = (TWO*u3*v3-(w3+ONE)) * chi3;
    JVdata[IDX(i+1,0)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,0)];
    /*   dR_u/dv */
    f1 = u1 * u1 * chi1;
    f2 = u2 * u2 * chi2;
    f3 = u3 * u3 * chi3;
    JVdata[IDX(i+1,0)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,1)];
    /*   dR_u/dw */
    f1 = -u1 * chi1;
    f2 = -u2 * chi2;
    f3 = -u3 * chi3;
    JVdata[IDX(i+1,0)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,2)];
    /*   dR_v/du */
    f1 = (w1 - TWO*u1*v1) * chi1;
    f2 = (w2 - TWO*u2*v2) * chi2;
    f3 = (w3 - TWO*u3*v3) * chi3;
    JVdata[IDX(i+1,1)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,0)];
    /*   dR_v/dv */
    f1 = -u1 * u1 * chi1;
    f2 = -u2 * u2 * chi2;
    f3 = -u3 * u3 * chi3;
    JVdata[IDX(i+1,1)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,1)];
    /*   dR_v/dw */
    f1 = u1 * chi1;
    f2 = u2 * chi2;
    f3 = u3 * chi3;
    JVdata[IDX(i+1,1)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,2)];
    /*   dR_w/du */
    f1 = -w1 * chi1;
    f2 = -w2 * chi2;
    f3 = -w3 * chi3;
    JVdata[IDX(i+1,2)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,0)];
    /*   dR_w/dw */
    f1 = (-ONE/ep - u1) * chi1;
    f2 = (-ONE/ep - u2) * chi2;
    f3 = (-ONE/ep - u3) * chi3;
    JVdata[IDX(i+1,2)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,2)];


    /* left trial function, evaluated at all quadrature nodes */
    if (left) {
      chi1 = ChiL(xl,xr,X1(xl,xr));
      chi2 = ChiL(xl,xr,X2(xl,xr));
      chi3 = ChiL(xl,xr,X3(xl,xr));
    } else {
      chi1 = chi2 = chi3 = 0.0;
    }

    /* right basis, left trial Jacobian entries */
    /*   dR_u/du */
    f1 = (TWO*u1*v1-(w1+ONE)) * chi1;
    f2 = (TWO*u2*v2-(w2+ONE)) * chi2;
    f3 = (TWO*u3*v3-(w3+ONE)) * chi3;
    JVdata[IDX(i,0)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,0)];
    /*   dR_u/dv */
    f1 = u1 * u1 * chi1;
    f2 = u2 * u2 * chi2;
    f3 = u3 * u3 * chi3;
    JVdata[IDX(i,0)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,1)];
    /*   dR_u/dw */
    f1 = -u1 * chi1;
    f2 = -u2 * chi2;
    f3 = -u3 * chi3;
    JVdata[IDX(i,0)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,2)];
    /*   dR_v/du */
    f1 = (w1 - TWO*u1*v1) * chi1;
    f2 = (w2 - TWO*u2*v2) * chi2;
    f3 = (w3 - TWO*u3*v3) * chi3;
    JVdata[IDX(i,1)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,0)];
    /*   dR_v/dv */
    f1 = -u1 * u1 * chi1;
    f2 = -u2 * u2 * chi2;
    f3 = -u3 * u3 * chi3;
    JVdata[IDX(i,1)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,1)];
    /*   dR_v/dw */
    f1 = u1 * chi1;
    f2 = u2 * chi2;
    f3 = u3 * chi3;
    JVdata[IDX(i,1)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,2)];
    /*   dR_w/du */
    f1 = -w1 * chi1;
    f2 = -w2 * chi2;
    f3 = -w3 * chi3;
    JVdata[IDX(i,2)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,0)];
    /*   dR_w/dw */
    f1 = (-ONE/ep - u1) * chi1;
    f2 = (-ONE/ep - u2) * chi2;
    f3 = (-ONE/ep - u3) * chi3;
    JVdata[IDX(i,2)] += Quad(f1,f2,f3,xl,xr) * Vdata[IDX(i+1,2)];

  }

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
