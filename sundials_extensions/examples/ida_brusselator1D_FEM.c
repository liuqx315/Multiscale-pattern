/* -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following test simulates a brusselator problem from chemical 
 * kinetics.  This is a PDE system with 3 components, Y = [u,v,w], 
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
 * 
 * Here, we use a piecewise linear Galerkin finite element 
 * discretization in space, where all element-wise integrals are 
 * computed using 3-node Gaussian quadrature (since we will have 
 * quartic polynomials in the reaction terms for the u_t and v_t 
 * equations (including the test function)).  While this system of 
 * equations does not have any constraint equations, the time 
 * derivative terms will include a mass matrix, giving rise to an 
 * ODE system of the form
 *      M y_t = L y + R(y),
 * where M is the 3x3 block mass matrix for each component, L is 
 * the 3x3 block Laplace operator for each component, and R(y) is 
 * comprised of the nonlinear reaction terms for each component.  
 * Since it it highly inefficient to rewrite this system as
 *      y_t = M^{-1}(L y + R(y)),
 * we do not wish to use CVODE for time integration, instead using 
 * IDA to write our ODE system in DAE form
 *      0 = F(t,y,y_t) = M y_t - L y - R(y).
 * We therefore provide functions to evaluate the residual F(t,y,y_t) 
 * and its Jacobian,  J = dF/dy + c*dF/dy_t.  In addition, since we 
 * only know initial conditions for y, we must call IDACalcIC to 
 * compute compatible initial conditions for y_t.
 *
 * The number of spatial intervals N, the parameters a, b, du, dv, dw 
 * and ep, and the desired relative and absolute solver tolerances,
 * are provided in the input file input_brusselator1D.txt.
 *
 * We use a vector-valued absolute tolerance, where the values are 
 * set as the input scalar value multiplied by the width of the 
 * support for the corresponding basis function.  On a uniform mesh 
 * this would result in a constant set of values, but on a 
 * non-uniform mesh this spreads these weights in an integral sense.
 * 
 * This program solves the problem with the BDF method, using a
 * Newton iteration with the IDABAND band linear solver, and the 
 * built-in difference-quotient Jacobian.
 *
 * 100 outputs are printed at equal time intervals, and run 
 * statistics are printed at the end.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Header files with a description of contents used */

#include <ida/ida.h>                 /* IDA functions & constants */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <ida/ida_band.h>            /* prototype for IDABand */
#include <sundials/sundials_band.h>  /* definitions of type DlsMat and macros */
#include <sundials/sundials_types.h> /* definition of type realtype */


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
  realtype fsc;  /* residual scaling factor */
} *UserData;



/* User-supplied Functions Called by the Solver */
static int F(realtype t, N_Vector y, N_Vector y_t, N_Vector r, void *user_data);

/* Private function to check function return values */
static int check_flag(void *flagvalue, char *funcname, int opt);
static int Calc_RHS(realtype t, N_Vector y, N_Vector RHS, UserData udata);
static int MassMatrix(realtype c, DlsMat Jac, UserData udata);
static int Calc_y0dot(realtype t0, N_Vector y, N_Vector y_t, 
		      N_Vector tmp, UserData udata);


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
  N_Vector y_t = NULL;
  N_Vector avtol = NULL;
  N_Vector umask = NULL;
  N_Vector vmask = NULL;
  N_Vector wmask = NULL;
  void *ida_mem = NULL;

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  udata->x = NULL;
  if (check_flag((void *)udata, "malloc", 2)) return(1);

  /* read problem parameter and tolerances from input file:
     N - number of spatial intervals
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

  /* set total allocated vector length (N-1 intervals, end points are homogeneous Dirichlet) */
  NEQ = Nvar*udata->N;

  /* Initial problem output */
  printf("\n1D FEM Brusselator PDE test problem:\n");
  printf("    N = %li,  NEQ = %li\n", udata->N, NEQ);
  printf("    problem parameters:  a = %g,  b = %g,  ep = %g\n",
	 udata->a, udata->b, udata->ep);
  printf("    diffusion coefficients:  du = %g,  dv = %g,  dw = %g\n", 
	 udata->du, udata->dv, udata->dw);
  printf("    reltol = %.1e,  abstol = %.1e\n\n", reltol, abstol);


  /* Create serial vectors of length NEQ for initial condition & abs tol */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  y_t = N_VNew_Serial(NEQ);
  if (check_flag((void *)y_t, "N_VNew_Serial", 0)) return(1);
  avtol = N_VNew_Serial(NEQ);
  if (check_flag((void *)avtol, "N_VNew_Serial", 0)) return(1);


  /* allocate and set up spatial mesh; this clusters more intervals 
     near the endpoints of the interval */
  udata->x = (realtype *) malloc(N*sizeof(realtype));
  if (check_flag((void *)udata->x, "malloc", 2)) return(1);
  realtype z, h = ONE/(N-1);
  for (i=0; i<N; i++) {
    z = h*i - 0.5;
    udata->x[i] = RCONST(4.0)*z*z*z + RCONST(0.5);
  }

  /* output mesh to disk */
  FID=fopen("bruss_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->x[i]);
  fclose(FID);


  /* set up residual scaling factor as the reciprocal of the average cell width */
  realtype hsum=0.0;
  for (i=0; i<N-1; i++)  hsum += udata->x[i+1] - udata->x[i];
  udata->fsc = hsum/(N-1);

  
  /* Access data array for new NVector y */
  data = N_VGetArrayPointer(y);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return(1);

  /* Set initial conditions into y */
  realtype pi = RCONST(4.0)*atan(ONE);
  for (i=0; i<N; i++) {
    data[IDX(i,0)] =  a  + RCONST(0.1)*sin(pi*udata->x[i]);  /* u */
    data[IDX(i,1)] = b/a + RCONST(0.1)*sin(pi*udata->x[i]);  /* v */
    data[IDX(i,2)] =  b  + RCONST(0.1)*sin(pi*udata->x[i]);  /* w */
  }


  /* Access data array for absolute tolerance NVector avtol */
  data = N_VGetArrayPointer(avtol);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return(1);

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

  
  /* Call IDACreate to create the solver memory */
  ida_mem = IDACreate();
  if (check_flag((void *)ida_mem, "IDACreate", 0)) return(1);
  
  /* Call IDASetUserData to pass rdata to user functions */
  flag = IDASetUserData(ida_mem, (void *) udata);
  if (check_flag(&flag, "IDASetUserData", 1)) return(1);

  /* Call IDAInit to initialize the integrator memory, specify the
     residual function, the inital time T0, and the initial dependent
     variable vectors y and y_t */
  flag = IDAInit(ida_mem, F, T0, y, y_t);
  if (check_flag(&flag, "IDAInit", 1)) return(1);
  
  /* Call IDASVtolerances to set tolerances */
  flag = IDASVtolerances(ida_mem, reltol, avtol);
  if(check_flag(&flag, "IDASVtolerances", 1)) return(1);

  /* Compute the initial condition for y_t to satisfy residual,
     using avtol for temporary storage since it's no longer used */
  if (Calc_y0dot(T0, y, y_t, avtol, udata)) {
    printf("Error in calculating initial condition for y_t!\n");
    return(1);
  }

  /* Free avtol */
  N_VDestroy_Serial(avtol);

  /* Call IDABand to specify the linear solver */
  flag = IDABand(ida_mem, NEQ, 5, 5);
  if (check_flag(&flag, "IDABand", 1)) return(1);

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

  /* In loop, call IDASolve, print results, and test for error.
     Break out of loop when the final output time has been reached */
  realtype t = T0;
  realtype dTout = Tf/Nt;
  realtype tout = dTout;
  realtype u, v, w;
  printf("        t      ||u||_rms   ||v||_rms   ||w||_rms\n");
  printf("   ----------------------------------------------\n");
  int iout;
  for (iout=0; iout<Nt; iout++) {

    flag = IDASolve(ida_mem, tout, &t, y, y_t, IDA_NORMAL);
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

    if (check_flag(&flag, "IDASolve", 1)) break;
    if (flag == IDA_SUCCESS) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
  }
  printf("   ----------------------------------------------\n");
  fclose(UFID);
  fclose(VFID);
  fclose(WFID);

  /* Print some final statistics */
  long int nst, nre, nje, nreLS, nni, ncfn, netf;
  flag = IDAGetNumSteps(ida_mem, &nst);
  check_flag(&flag, "IDAGetNumSteps", 1);
  flag = IDAGetNumResEvals(ida_mem, &nre);
  check_flag(&flag, "IDAGetNumResEvals", 1);
  flag = IDADlsGetNumJacEvals(ida_mem, &nje);
  check_flag(&flag, "IDADlsGetNumJacEvals", 1);
  flag = IDAGetNumNonlinSolvIters(ida_mem, &nni);
  check_flag(&flag, "IDAGetNumNonlinSolvIters", 1);
  flag = IDAGetNumErrTestFails(ida_mem, &netf);
  check_flag(&flag, "IDAGetNumErrTestFails", 1);
  flag = IDAGetNumNonlinSolvConvFails(ida_mem, &ncfn);
  check_flag(&flag, "IDAGetNumNonlinSolvConvFails", 1);
  flag = IDADlsGetNumResEvals(ida_mem, &nreLS);
  check_flag(&flag, "IDADlsGetNumResEvals", 1);

  printf("\nFinal Solver Statistics:\n");
  printf("   Total internal solver steps = %li\n", nst);
  printf("   Total Res evals = %li\n", nre+nreLS);
  printf("   Total number of Jacobian evaluations = %li\n", nje);
  printf("   Total number of Newton iterations = %li\n", nni);
  printf("   Total number of error test failures = %li\n\n", netf);
  printf("   Total number of nonlinear convergence failures = %li\n", ncfn);

  /* Free vectors */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(y_t);
  N_VDestroy_Serial(umask);
  N_VDestroy_Serial(vmask);
  N_VDestroy_Serial(wmask);

  /* Free user data */
  free(udata->x);
  free(udata);

  /* Free integrator memory */
  IDAFree(&ida_mem);

  return(0);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* Routine to compute DAE residual function, (M y_t - f(t,y)) */
static int F(realtype t, N_Vector y, N_Vector y_t, N_Vector r, void *user_data)
{
  /* clear out r (to be careful) */
  N_VConst(0.0, r);

  /* problem data */
  UserData udata = (UserData) user_data;

  /* shortcuts to number of nodes */
  long int N = udata->N;

  /* access data arrays */
  realtype *dYdata = N_VGetArrayPointer(y_t);
  if (check_flag((void *)dYdata, "N_VGetArrayPointer", 0)) return(1);
  realtype *ResData = N_VGetArrayPointer(r);
  if (check_flag((void *)ResData, "N_VGetArrayPointer", 0)) return(1);

  /* first, fill in RHS portion, f(t,y) */
  if (Calc_RHS(t, y, r, udata)) {
    printf("Residual calculation error in calling Calc_RHS!\n");
    return(1);
  }
  
  /* set r = -fsc * r */
  N_VScale(-(udata->fsc), r, r);

  /* iterate over intervals, adding in (M y_t) term */
  long int i;
  realtype udotl, udotr, vdotl, vdotr, wdotl, wdotr;
  realtype xl, xr, f1, f2, f3;
  for (i=0; i<N-1; i++) {

    /* set nodal value shortcuts (interval index aligns with left node) */
    udotl = dYdata[IDX(i,0)];
    vdotl = dYdata[IDX(i,1)];
    wdotl = dYdata[IDX(i,2)];
    udotr = dYdata[IDX(i+1,0)];
    vdotr = dYdata[IDX(i+1,1)];
    wdotr = dYdata[IDX(i+1,2)];

    /* set mesh shortcuts */
    xl = udata->x[i];
    xr = udata->x[i+1];

    /* evaluate M*y_t on this subinterval */
    /*    left test function */
    /*  u */
    f1 = Eval(udotl,udotr,xl,xr,X1(xl,xr)) * ChiL(xl,xr,X1(xl,xr));
    f2 = Eval(udotl,udotr,xl,xr,X2(xl,xr)) * ChiL(xl,xr,X2(xl,xr));
    f3 = Eval(udotl,udotr,xl,xr,X3(xl,xr)) * ChiL(xl,xr,X3(xl,xr));
    ResData[IDX(i,0)] += udata->fsc * Quad(f1,f2,f3,xl,xr);
    
    /*  v */
    f1 = Eval(vdotl,vdotr,xl,xr,X1(xl,xr)) * ChiL(xl,xr,X1(xl,xr));
    f2 = Eval(vdotl,vdotr,xl,xr,X2(xl,xr)) * ChiL(xl,xr,X2(xl,xr));
    f3 = Eval(vdotl,vdotr,xl,xr,X3(xl,xr)) * ChiL(xl,xr,X3(xl,xr));
    ResData[IDX(i,1)] += udata->fsc * Quad(f1,f2,f3,xl,xr);
    
    /*  w */
    f1 = Eval(wdotl,wdotr,xl,xr,X1(xl,xr)) * ChiL(xl,xr,X1(xl,xr));
    f2 = Eval(wdotl,wdotr,xl,xr,X2(xl,xr)) * ChiL(xl,xr,X2(xl,xr));
    f3 = Eval(wdotl,wdotr,xl,xr,X3(xl,xr)) * ChiL(xl,xr,X3(xl,xr));
    ResData[IDX(i,2)] += udata->fsc * Quad(f1,f2,f3,xl,xr);

    /*    right test function */
    /*  u */
    f1 = Eval(udotl,udotr,xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
    f2 = Eval(udotl,udotr,xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
    f3 = Eval(udotl,udotr,xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
    ResData[IDX(i+1,0)] += udata->fsc * Quad(f1,f2,f3,xl,xr);
    
    /*  v */
    f1 = Eval(vdotl,vdotr,xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
    f2 = Eval(vdotl,vdotr,xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
    f3 = Eval(vdotl,vdotr,xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
    ResData[IDX(i+1,1)] += udata->fsc * Quad(f1,f2,f3,xl,xr);
    
    /*  w */
    f1 = Eval(wdotl,wdotr,xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
    f2 = Eval(wdotl,wdotr,xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
    f3 = Eval(wdotl,wdotr,xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
    ResData[IDX(i+1,2)] += udata->fsc * Quad(f1,f2,f3,xl,xr);
  }

  return(0);
}


/*-------------------------------
 * Private helper functions
 *-------------------------------*/


/* Routine to compute the ODE RHS function f(t,y), where system is of the form
        M y_t = f(t,y) := Ly + R(y) 
   This routine only computes the f(t,y), leaving (M y_t) alone. */
static int Calc_RHS(realtype t, N_Vector y, N_Vector RHS, UserData udata)
{
  /* clear out RHS (to be careful) */
  N_VConst(0.0, RHS);

  /* shortcuts to number of intervals, background values */
  long int N  = udata->N;
  realtype a  = udata->a;
  realtype b  = udata->b;
  realtype ep = udata->ep;
  realtype du = udata->du;
  realtype dv = udata->dv;
  realtype dw = udata->dw;

  /* access data arrays */
  realtype *Ydata = N_VGetArrayPointer(y);
  if (check_flag((void *)Ydata, "N_VGetArrayPointer", 0)) return(1);
  realtype *RHSdata = N_VGetArrayPointer(RHS);
  if (check_flag((void *)RHSdata, "N_VGetArrayPointer", 0)) return(1);

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

  return(0);
}



/* Routine to compute the mass matrix multiplying y_t, scaled by the factor c.
   We add the result into Jac and do not erase what was already there */
static int MassMatrix(realtype c, DlsMat Jac, UserData udata)
{
  /* shortcut to number of nodes */
  long int N = udata->N;

  /* iterate over intervals, filling in matrix entries */
  long int i;
  realtype xl, xr, f1, f2, f3;
  for (i=0; i<N-1; i++) {

    /* set mesh shortcuts */
    xl = udata->x[i];
    xr = udata->x[i+1];

    /* Jacobian of (M*y_t) on this subinterval, all vars have identical Jacobians, 
       and all require multiplication by c to account for the fact that these are 
       derivatives wrt y_t */

    /*    left basis and test functions */
    f1 = c * ChiL(xl,xr,X1(xl,xr)) * ChiL(xl,xr,X1(xl,xr));
    f2 = c * ChiL(xl,xr,X2(xl,xr)) * ChiL(xl,xr,X2(xl,xr));
    f3 = c * ChiL(xl,xr,X3(xl,xr)) * ChiL(xl,xr,X3(xl,xr));
    BAND_ELEM(Jac,IDX(i,0),IDX(i,0)) += Quad(f1,f2,f3,xl,xr);
    BAND_ELEM(Jac,IDX(i,1),IDX(i,1)) += Quad(f1,f2,f3,xl,xr);
    BAND_ELEM(Jac,IDX(i,2),IDX(i,2)) += Quad(f1,f2,f3,xl,xr);

    /*    right basis and test functions */
    f1 = c * ChiR(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
    f2 = c * ChiR(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
    f3 = c * ChiR(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
    BAND_ELEM(Jac,IDX(i+1,0),IDX(i+1,0)) += Quad(f1,f2,f3,xl,xr);
    BAND_ELEM(Jac,IDX(i+1,1),IDX(i+1,1)) += Quad(f1,f2,f3,xl,xr);
    BAND_ELEM(Jac,IDX(i+1,2),IDX(i+1,2)) += Quad(f1,f2,f3,xl,xr);

    /*    left and right basis and test functions */
    f1 = c * ChiL(xl,xr,X1(xl,xr)) * ChiR(xl,xr,X1(xl,xr));
    f2 = c * ChiL(xl,xr,X2(xl,xr)) * ChiR(xl,xr,X2(xl,xr));
    f3 = c * ChiL(xl,xr,X3(xl,xr)) * ChiR(xl,xr,X3(xl,xr));
    BAND_ELEM(Jac,IDX(i,0),IDX(i+1,0)) += Quad(f1,f2,f3,xl,xr);
    BAND_ELEM(Jac,IDX(i,1),IDX(i+1,1)) += Quad(f1,f2,f3,xl,xr);
    BAND_ELEM(Jac,IDX(i,2),IDX(i+1,2)) += Quad(f1,f2,f3,xl,xr);
    BAND_ELEM(Jac,IDX(i+1,0),IDX(i,0)) += Quad(f1,f2,f3,xl,xr);
    BAND_ELEM(Jac,IDX(i+1,1),IDX(i,1)) += Quad(f1,f2,f3,xl,xr);
    BAND_ELEM(Jac,IDX(i+1,2),IDX(i,2)) += Quad(f1,f2,f3,xl,xr);

  }

  return(0);
}



/* Routine to compute initial conditions for y_t
   by creating/solving the linear system M*y_t = RHS */
static int Calc_y0dot(realtype t0, N_Vector y, N_Vector y_t, 
		      N_Vector tmp, UserData udata)
{
  /* shortcuts to number of nodes, etc */
  long int N = udata->N;
  long int nrows = 3*N;
  long int bw = 5;
  long int smu = 10;
  realtype c = RCONST(1.0);
  long int *p = (long int *) malloc(nrows*sizeof(long int));

  /* create/fill Mass matrix M */
  DlsMat M = NewBandMat(nrows, bw, bw, smu);
  SetToZero(M);
  if (MassMatrix(c, M, udata)) {
    printf("Calc_y0dot error in calling MassMatrix!\n");
    return(1);
  }

  /* store (RHS) in y_t N_Vector */
  if (Calc_RHS(t0, y, y_t, udata)) {
    printf("Calc_y0dot error in calling CalcRHS!\n");
    return(1);
  }

  /* factor mass matrix */
  if (BandGBTRF(M, p)) {
    printf("Calc_y0dot error in factoring mass matrix!\n");
    return(1);
  }

  /* solve with factored matrix, using storage from y_t for rhs and sol */
  realtype *b = N_VGetArrayPointer(y_t);
  if (check_flag((void *)b, "N_VGetArrayPointer", 0)) return(1);
  BandGBTRS(M, p, b);

  /* check whether y_t satisfies the initial residual, storing resid in tmp */
  realtype resid;
  if (F(t0, y, y_t, tmp, (void *) udata)) {
    printf("Calc_y0dot error in calling residual routine!\n");
    return(1);
  }
  resid = N_VMaxNorm(tmp);
  printf(" Calc_y0dot:  || F(t0,y0,y0_t) ||_max = %g\n\n", resid);

  /* if the initial residual is too large, return with an error */
  if (resid > 1.0e-7)  return(1);

  /* clean up */
  free(p);
  DestroyMat(M);

  return(0);
}




/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer  */
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
