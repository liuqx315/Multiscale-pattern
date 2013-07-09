/* Header files */
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode.h>           /* prototypes for ARKode fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* def. of type 'realtype' */

using namespace std;

/* constants */
#define ONE (RCONST(1.0))
#define pi RCONST(3.1415926535897932)

/* user data structure */
typedef struct {  
  long int N;     /* grids                   */
  realtype dx;    /* mesh spacing            */
  realtype delta; /* ratio of dt/dx          */
} *UserData;

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const string funcname, int opt);

int main(int argc, const char * argv[])
{

   /* general problem parameters */
  realtype T0 = RCONST(0.0);
  realtype Tf = RCONST(150.0);
  int Nt = 30000;
  UserData udata = NULL;
  realtype *data;
  long int N, i;

  /* general problem variables */
  int flag;
  N_Vector y = NULL;
  void *arkode_mem = NULL;

  /* allocate udata structure */
  udata = (UserData) malloc(sizeof(*udata));
  if (check_flag((void *) udata, "malloc", 2)) return 1;

  /* read problem parameter and tolerances from input file:
     N - number of spatial discretization points
     delta - ratio of dt/dx */
  double delta;
  FILE *FID;
  FID = fopen("input_WENO1D.txt","r");
  flag = fscanf(FID,"  N = %li\n", &N);
  flag = fscanf(FID,"  delta = %lf\n", &delta);
  fclose(FID);

  /* store the inputs in the UserData structure */
  udata->N = N;
  udata->delta = delta;

  /* open solver diagnostics output file for writing */
  //FILE *DFID;
  //DFID=fopen("diags_ark_WENO1D.txt","w");
  
  /* Initial problem output */
  printf("\n1D WENO ODE test problem:\n");
  printf("  N = %li\n", udata->N);
  printf("  ratio of dt/dx:  delta = %g\n", udata->delta);

  /* Create serial vector of length N for initial condition */
  y = N_VNew_Serial(N);
  if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;

  /* set spatial mesh spacing */
  udata->dx = RCONST(2.0)/(1.0*N-1.0);

  /* output mesh to disk */
  FID=fopen("WENO1D_mesh.txt","w");
  for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->dx*i);
  fclose(FID);

  /* Access data array for new NVector y */
  data = N_VGetArrayPointer(y);
  if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;

  /* Set initial conditions into y */
  //realtype pi = RCONST(4.0)*atan(ONE);
  for (i=0; i<N; i++) {
    /* discontinuous initial condition */
    // if (i*udata->dx<0.5||(i*udata->dx<=2&&i*udata->dx>1.5))
    //  data[i]=1.0;
    //else
    //  data[i]=0.0;
    /* continuous initial condition */
    data[i] = sin(pi*i*udata->dx)*sin(pi*i*udata->dx)*sin(pi*i*udata->dx)*sin(pi*i*udata->dx)*sin(pi*i*udata->dx)*sin(pi*i*udata->dx)*sin(pi*i*udata->dx)*sin(pi*i*udata->dx)*sin(pi*i*udata->dx);
  }

  /* Call ARKodeCreate to create the solver memory */
  arkode_mem = ARKodeCreate();
  if (check_flag((void *) arkode_mem, "ARKodeCreate", 0)) return 1;
  
  /* Set relative tolerance and absolute tolerance  */
  //realtype reltol = 1.e+12;
  //realtype abstol = 1.e+12;
  realtype reltol = 1.e-3;
  realtype abstol = 1.e-6;
  
  /* Solution will be computed with default explicit method */
  flag = ARKodeInit(arkode_mem, f, NULL, T0, y);
  if (check_flag(&flag, "ARKodeInit", 1)) return 1;

  /* Call ARKodeSetUserData to pass rdata to user functions */
  flag = ARKodeSetUserData(arkode_mem, (void *) udata);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;

  /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
  //  flag = ARKodeSetDiagnostics(arkode_mem, DFID);
  //if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return 1;

  /* Call ARKodeSetInitStep to initialize time step */
  flag = ARKodeSetInitStep(arkode_mem, (udata->delta*udata->dx));
  if (check_flag(&flag, "ARKodeSetInitStep", 1)) return 1;

  /* Call ARKodeSetMinStep to set min time step */
  //flag = ARKodeSetMinStep(arkode_mem, (udata->delta*udata->dx));
  //if (check_flag(&flag, "ARKodeSetMinStep", 1)) return 1;

  /* Call ARKodeSetMaxStep to set max time step */
  //flag = ARKodeSetMaxStep(arkode_mem, (udata->delta*udata->dx));
  //if (check_flag(&flag, "ARKodeSetMaxStep", 1)) return 1;

  /* Call ARKodeSetMaxNumSteps to increase default (for testing) */
  flag = ARKodeSetMaxNumSteps(arkode_mem, 100000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;

  /* Call ARKodeSStolerances to specify the scalar relative and absolute
     tolerances */
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  /* Open output stream for results, access data arrays */
  FILE *UFID=fopen("WENO1D.txt","w");
  data = N_VGetArrayPointer(y);
  if (check_flag((void *)data, "N_VGetArrayPointer", 0)) return 1;

  /* output initial condition to disk */
   for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[i]);
  fprintf(UFID,"\n");


  /* In loop, call ARKode, print results.
     Break out of loop when the final output time has been reached */
  realtype t  = T0;
  //  realtype dTout = udata->delta*udata->dx;
  realtype dTout = Tf/Nt;
  realtype tout = T0+dTout;
  int iout;
  for(iout=0;iout<Nt;iout++){

    // stop exactly at this time
    //flag = ARKodeSetStopTime(arkode_mem, tout);   
    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
    if (check_flag(&flag, "ARKode", 1)) break;
    
    if (flag >= 0) {
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }
  
    if (iout>29990){
    /* output last nine results to disk */
     for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[i]);
     fprintf(UFID,"\n");
    }
  }
  
  fclose(UFID);

  /* Free vectors */
  N_VDestroy_Serial(y);

  /* Free user data */
  free(udata);

  /* Free integrator memory */
  ARKodeFree(&arkode_mem);

  /* close solver diagnostics output file */
  //fclose(DFID);

  return 0;
}


/*--------------------------------
 * Functions called by the solver
 *--------------------------------*/

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  /* clear out ydot (to be careful) */
  N_VConst(0.0, ydot);
  
  /* declare variables */
  long int i; 
  realtype IS0_p,IS1_p,IS2_p,IS0_n,IS1_n,IS2_n,Epsilon;
  realtype alpha_0p,alpha_1p,alpha_2p,alpha_0n,alpha_1n,alpha_2n;
  realtype w0_p,w1_p,w2_p,w0_n,w1_n,w2_n,u_tpph,u_tpnh,u_tnph,u_tnnh;

  /* problem data */
  UserData udata = (UserData) user_data;
  
  /* create y_positive and y_negative */
  N_Vector yp = NULL;
  N_Vector yn = NULL;

  /* shortcuts to number of intervals, background values */
  long int N  = udata->N;
  realtype delta  = udata->delta;
  realtype dx = udata->dx;

  /* Create serial vector of length N for initial condition */
  yp = N_VNew_Serial(N);
  if (check_flag((void *) yp, "N_VNew_Serial", 0)) return 1;
  yn = N_VNew_Serial(N);
  if (check_flag((void *) yn, "N_VNew_Serial", 0)) return 1;

  /* set positive y and negative y */
  N_VLinearSum( 0.5, y, 0.5, y, yp );
  N_VLinearSum( 0.5, y, -0.5, y, yn );

  /* access data arrays */
  realtype *Y = N_VGetArrayPointer(y);
  if (check_flag((void *) Y, "N_VGetArrayPointer", 0)) return 1;
  realtype *Ydot = N_VGetArrayPointer(ydot);
  if (check_flag((void *) Ydot, "N_VGetArrayPointer", 0)) return 1;
  realtype *YP = N_VGetArrayPointer(yp);
  if (check_flag((void *) YP, "N_VGetArrayPointer", 0)) return 1;
  realtype *YN = N_VGetArrayPointer(yn);
  if (check_flag((void *) YN, "N_VGetArrayPointer", 0)) return 1;
  
  /* iterate over domain, computing all equations */
  /* i=0, i=1, i=2, i=N-3, i=N-2, i=N-1, should be considered particularly */
    i=0;
    IS0_p=(13.0/12.0)*(YP[N-3]-2.0*YP[N-2]+YP[0])*(YP[N-3]-2.0*YP[N-2]+YP[0])+(1.0/4.0)*(YP[N-3]-4.0*YP[N-2]+3.0*YP[0])*(YP[N-3]-4.0*YP[N-2]+3.0*YP[0]);

    IS1_p=(13.0/12.0)*(YP[N-2]-2.0*YP[0]+YP[1])*(YP[N-2]-2.0*YP[0]+YP[1])+(1.0/4.0)*(YP[N-2]-YP[1])*(YP[N-2]-YP[1]);

    IS2_p=(13.0/12.0)*(YP[0]-2.0*YP[1]+YP[2])*(YP[0]-2.0*YP[1]+YP[2])+(1.0/4.0)*(3.0*YP[0]-4.0*YP[1]+YP[2])*(3.0*YP[0]-4.0*YP[1]+YP[2]);

    IS0_n=(13.0/12.0)*(YN[1]-2.0*YN[2]+YN[3])*(YN[1]-2.0*YN[2]+YN[3])+(1.0/4.0)*(3.0*YN[1]-4.0*YN[2]+YN[3])*(3.0*YN[1]-4.0*YN[2]+YN[3]);

    IS1_n=(13.0/12.0)*(YN[0]-2.0*YN[1]+YN[2])*(YN[0]-2.0*YN[1]+YN[2])+(1.0/4.0)*(YN[0]-YN[2])*(YN[0]-YN[2]);

    IS2_n=(13.0/12.0)*(YN[N-2]-2.0*YN[0]+YN[1])*(YN[N-2]-2.0*YN[0]+YN[1])+(1.0/4.0)*(YN[N-2]-4.0*YN[0]+3.0*YN[1])*(YN[N-2]-4.0*YN[0]+3.0*YN[1]);

    Epsilon=0.000001;
    alpha_0p=(1.0/10.0)*(1.0/(Epsilon+IS0_p))*(1.0/(Epsilon+IS0_p));
    alpha_1p=(6.0/10.0)*(1.0/(Epsilon+IS1_p))*(1.0/(Epsilon+IS1_p));
    alpha_2p=(3.0/10.0)*(1.0/(Epsilon+IS2_p))*(1.0/(Epsilon+IS2_p));
    alpha_0n=(1.0/10.0)*(1.0/(Epsilon+IS0_n))*(1.0/(Epsilon+IS0_n));
    alpha_1n=(6.0/10.0)*(1.0/(Epsilon+IS1_n))*(1.0/(Epsilon+IS1_n));
    alpha_2n=(3.0/10.0)*(1.0/(Epsilon+IS2_n))*(1.0/(Epsilon+IS2_n));

    w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
    w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
    w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
    w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
    w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
    w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);

    u_tpph=w0_p*((2.0/6.0)*YP[N-3]-(7.0/6.0)*YP[N-2]+(11.0/6.0)*YP[0])+w1_p*((-1.0/6.0)*YP[N-2]+(5.0/6.0)*YP[0]+(2.0/6.0)*YP[1])+w2_p*((2.0/6.0)*YP[0]+(5.0/6.0)*YP[1]-(1.0/6.0)*YP[2]);

    u_tnph=w2_n*((-1.0/6.0)*YN[N-2]+(5.0/6.0)*YN[0]+(2.0/6.0)*YN[1])+w1_n*((2.0/6.0)*YN[0]+(5.0/6.0)*YN[1]-(1.0/6.0)*YN[2])+w0_n*((11.0/6.0)*YN[1]-(7.0/6.0)*YN[2]+(2.0/6.0)*YN[3]);

    u_tpnh=w0_p*((2.0/6.0)*YP[N-4]-(7.0/6.0)*YP[N-3]+(11.0/6.0)*YP[N-2])+w1_p*((-1.0/6.0)*YP[N-3]+(5.0/6.0)*YP[N-2]+(2.0/6.0)*YP[0])+w2_p*((2.0/6.0)*YP[N-2]+(5.0/6.0)*YP[0]-(1.0/6.0)*YP[1]);

    u_tnnh=w2_n*((-1.0/6.0)*YN[N-3]+(5.0/6.0)*YN[N-2]+(2.0/6.0)*YN[0])+w1_n*((2.0/6.0)*YN[N-2]+(5.0/6.0)*YN[0]-(1.0/6.0)*YN[1])+w0_n*((11.0/6.0)*YN[0]-(7.0/6.0)*YN[1]+(2.0/6.0)*YN[2]);

    Ydot[0]=-(1.0/udata->dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
   
    i=1;
    IS0_p=(13.0/12.0)*(YP[N-2]-2.0*YP[i-1]+YP[i])*(YP[N-2]-2.0*YP[i-1]+YP[i])+(1.0/4.0)*(YP[N-2]-4.0*YP[i-1]+3.0*YP[i])*(YP[N-2]-4.0*YP[i-1]+3.0*YP[i]);

    IS1_p=(13.0/12.0)*(YP[i-1]-2.0*YP[i]+YP[i+1])*(YP[i-1]-2.0*YP[i]+YP[i+1])+(1.0/4.0)*(YP[i-1]-YP[i+1])*(YP[i-1]-YP[i+1]);

    IS2_p=(13.0/12.0)*(YP[i]-2.0*YP[i+1]+YP[i+2])*(YP[i]-2.0*YP[i+1]+YP[i+2])+(1.0/4.0)*(3.0*YP[i]-4.0*YP[i+1]+YP[i+2])*(3.0*YP[i]-4.0*YP[i+1]+YP[i+2]);

    IS0_n=(13.0/12.0)*(YN[i+1]-2.0*YN[i+2]+YN[i+3])*(YN[i+1]-2.0*YN[i+2]+YN[i+3])+(1.0/4.0)*(3.0*YN[i+1]-4.0*YN[i+2]+YN[i+3])*(3.0*YN[i+1]-4.0*YN[i+2]+YN[i+3]);

    IS1_n=(13.0/12.0)*(YN[i]-2.0*YN[i+1]+YN[i+2])*(YN[i]-2.0*YN[i+1]+YN[i+2])+(1.0/4.0)*(YN[i]-YN[i+2])*(YN[i]-YN[i+2]);

    IS2_n=(13.0/12.0)*(YN[i-1]-2.0*YN[i]+YN[i+1])*(YN[i-1]-2.0*YN[i]+YN[i+1])+(1.0/4.0)*(YN[i-1]-4.0*YN[i]+3.0*YN[i+1])*(YN[i-1]-4.0*YN[i]+3.0*YN[i+1]);

    Epsilon=0.000001;
    alpha_0p=(1.0/10.0)*(1.0/(Epsilon+IS0_p))*(1.0/(Epsilon+IS0_p));
    alpha_1p=(6.0/10.0)*(1.0/(Epsilon+IS1_p))*(1.0/(Epsilon+IS1_p));
    alpha_2p=(3.0/10.0)*(1.0/(Epsilon+IS2_p))*(1.0/(Epsilon+IS2_p));
    alpha_0n=(1.0/10.0)*(1.0/(Epsilon+IS0_n))*(1.0/(Epsilon+IS0_n));
    alpha_1n=(6.0/10.0)*(1.0/(Epsilon+IS1_n))*(1.0/(Epsilon+IS1_n));
    alpha_2n=(3.0/10.0)*(1.0/(Epsilon+IS2_n))*(1.0/(Epsilon+IS2_n));

    w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
    w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
    w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
    w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
    w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
    w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);

    u_tpph=w0_p*((2.0/6.0)*YP[N-2]-(7.0/6.0)*YP[i-1]+(11.0/6.0)*YP[i])+w1_p*((-1.0/6.0)*YP[i-1]+(5.0/6.0)*YP[i]+(2.0/6.0)*YP[i+1])+w2_p*((2.0/6.0)*YP[i]+(5.0/6.0)*YP[i+1]-(1.0/6.0)*YP[i+2]);

    u_tnph=w2_n*((-1.0/6.0)*YN[i-1]+(5.0/6.0)*YN[i]+(2.0/6.0)*YN[i+1])+w1_n*((2.0/6.0)*YN[i]+(5.0/6.0)*YN[i+1]-(1.0/6.0)*YN[i+2])+w0_n*((11.0/6.0)*YN[i+1]-(7.0/6.0)*YN[i+2]+(2.0/6.0)*YN[i+3]);

    u_tpnh=w0_p*((2.0/6.0)*YP[N-3]-(7.0/6.0)*YP[N-2]+(11.0/6.0)*YP[i-1])+w1_p*((-1.0/6.0)*YP[N-2]+(5.0/6.0)*YP[i-1]+(2.0/6.0)*YP[i])+w2_p*((2.0/6.0)*YP[i-1]+(5.0/6.0)*YP[i]-(1.0/6.0)*YP[i+1]);

    u_tnnh=w2_n*((-1.0/6.0)*YN[N-2]+(5.0/6.0)*YN[i-1]+(2.0/6.0)*YN[i])+w1_n*((2.0/6.0)*YN[i-1]+(5.0/6.0)*YN[i]-(1.0/6.0)*YN[i+1])+w0_n*((11.0/6.0)*YN[i]-(7.0/6.0)*YN[i+1]+(2.0/6.0)*YN[i+2]);

    Ydot[i]=-(1.0/udata->dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
   
      i=2;
      IS0_p=(13.0/12.0)*(YP[i-2]-2.0*YP[i-1]+YP[i])*(YP[i-2]-2.0*YP[i-1]+YP[i])+(1.0/4.0)*(YP[i-2]-4.0*YP[i-1]+3.0*YP[i])*(YP[i-2]-4.0*YP[i-1]+3.0*YP[i]);

    IS1_p=(13.0/12.0)*(YP[i-1]-2.0*YP[i]+YP[i+1])*(YP[i-1]-2.0*YP[i]+YP[i+1])+(1.0/4.0)*(YP[i-1]-YP[i+1])*(YP[i-1]-YP[i+1]);

    IS2_p=(13.0/12.0)*(YP[i]-2.0*YP[i+1]+YP[i+2])*(YP[i]-2.0*YP[i+1]+YP[i+2])+(1.0/4.0)*(3.0*YP[i]-4.0*YP[i+1]+YP[i+2])*(3.0*YP[i]-4.0*YP[i+1]+YP[i+2]);

    IS0_n=(13.0/12.0)*(YN[i+1]-2.0*YN[i+2]+YN[i+3])*(YN[i+1]-2.0*YN[i+2]+YN[i+3])+(1.0/4.0)*(3.0*YN[i+1]-4.0*YN[i+2]+YN[i+3])*(3.0*YN[i+1]-4.0*YN[i+2]+YN[i+3]);

    IS1_n=(13.0/12.0)*(YN[i]-2.0*YN[i+1]+YN[i+2])*(YN[i]-2.0*YN[i+1]+YN[i+2])+(1.0/4.0)*(YN[i]-YN[i+2])*(YN[i]-YN[i+2]);

    IS2_n=(13.0/12.0)*(YN[i-1]-2.0*YN[i]+YN[i+1])*(YN[i-1]-2.0*YN[i]+YN[i+1])+(1.0/4.0)*(YN[i-1]-4.0*YN[i]+3.0*YN[i+1])*(YN[i-1]-4.0*YN[i]+3.0*YN[i+1]);

    Epsilon=0.000001;
    alpha_0p=(1.0/10.0)*(1.0/(Epsilon+IS0_p))*(1.0/(Epsilon+IS0_p));
    alpha_1p=(6.0/10.0)*(1.0/(Epsilon+IS1_p))*(1.0/(Epsilon+IS1_p));
    alpha_2p=(3.0/10.0)*(1.0/(Epsilon+IS2_p))*(1.0/(Epsilon+IS2_p));
    alpha_0n=(1.0/10.0)*(1.0/(Epsilon+IS0_n))*(1.0/(Epsilon+IS0_n));
    alpha_1n=(6.0/10.0)*(1.0/(Epsilon+IS1_n))*(1.0/(Epsilon+IS1_n));
    alpha_2n=(3.0/10.0)*(1.0/(Epsilon+IS2_n))*(1.0/(Epsilon+IS2_n));

    w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
    w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
    w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
    w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
    w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
    w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);

    u_tpph=w0_p*((2.0/6.0)*YP[i-2]-(7.0/6.0)*YP[i-1]+(11.0/6.0)*YP[i])+w1_p*((-1.0/6.0)*YP[i-1]+(5.0/6.0)*YP[i]+(2.0/6.0)*YP[i+1])+w2_p*((2.0/6.0)*YP[i]+(5.0/6.0)*YP[i+1]-(1.0/6.0)*YP[i+2]);

    u_tnph=w2_n*((-1.0/6.0)*YN[i-1]+(5.0/6.0)*YN[i]+(2.0/6.0)*YN[i+1])+w1_n*((2.0/6.0)*YN[i]+(5.0/6.0)*YN[i+1]-(1.0/6.0)*YN[i+2])+w0_n*((11.0/6.0)*YN[i+1]-(7.0/6.0)*YN[i+2]+(2.0/6.0)*YN[i+3]);

    u_tpnh=w0_p*((2.0/6.0)*YP[N-2]-(7.0/6.0)*YP[i-2]+(11.0/6.0)*YP[i-1])+w1_p*((-1.0/6.0)*YP[i-2]+(5.0/6.0)*YP[i-1]+(2.0/6.0)*YP[i])+w2_p*((2.0/6.0)*YP[i-1]+(5.0/6.0)*YP[i]-(1.0/6.0)*YP[i+1]);

    u_tnnh=w2_n*((-1.0/6.0)*YN[i-2]+(5.0/6.0)*YN[i-1]+(2.0/6.0)*YN[i])+w1_n*((2.0/6.0)*YN[i-1]+(5.0/6.0)*YN[i]-(1.0/6.0)*YN[i+1])+w0_n*((11.0/6.0)*YN[i]-(7.0/6.0)*YN[i+1]+(2.0/6.0)*YN[i+2]);

    Ydot[i]=-(1.0/udata->dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));

  for (i=3; i<N-3; i++){

    IS0_p=(13.0/12.0)*(YP[i-2]-2.0*YP[i-1]+YP[i])*(YP[i-2]-2.0*YP[i-1]+YP[i])+(1.0/4.0)*(YP[i-2]-4.0*YP[i-1]+3.0*YP[i])*(YP[i-2]-4.0*YP[i-1]+3.0*YP[i]);

    IS1_p=(13.0/12.0)*(YP[i-1]-2.0*YP[i]+YP[i+1])*(YP[i-1]-2.0*YP[i]+YP[i+1])+(1.0/4.0)*(YP[i-1]-YP[i+1])*(YP[i-1]-YP[i+1]);

    IS2_p=(13.0/12.0)*(YP[i]-2.0*YP[i+1]+YP[i+2])*(YP[i]-2.0*YP[i+1]+YP[i+2])+(1.0/4.0)*(3.0*YP[i]-4.0*YP[i+1]+YP[i+2])*(3.0*YP[i]-4.0*YP[i+1]+YP[i+2]);

    IS0_n=(13.0/12.0)*(YN[i+1]-2.0*YN[i+2]+YN[i+3])*(YN[i+1]-2.0*YN[i+2]+YN[i+3])+(1.0/4.0)*(3.0*YN[i+1]-4.0*YN[i+2]+YN[i+3])*(3.0*YN[i+1]-4.0*YN[i+2]+YN[i+3]);

    IS1_n=(13.0/12.0)*(YN[i]-2.0*YN[i+1]+YN[i+2])*(YN[i]-2.0*YN[i+1]+YN[i+2])+(1.0/4.0)*(YN[i]-YN[i+2])*(YN[i]-YN[i+2]);

    IS2_n=(13.0/12.0)*(YN[i-1]-2.0*YN[i]+YN[i+1])*(YN[i-1]-2.0*YN[i]+YN[i+1])+(1.0/4.0)*(YN[i-1]-4.0*YN[i]+3.0*YN[i+1])*(YN[i-1]-4.0*YN[i]+3.0*YN[i+1]);

    Epsilon=0.000001;
    alpha_0p=(1.0/10.0)*(1.0/(Epsilon+IS0_p))*(1.0/(Epsilon+IS0_p));
    alpha_1p=(6.0/10.0)*(1.0/(Epsilon+IS1_p))*(1.0/(Epsilon+IS1_p));
    alpha_2p=(3.0/10.0)*(1.0/(Epsilon+IS2_p))*(1.0/(Epsilon+IS2_p));
    alpha_0n=(1.0/10.0)*(1.0/(Epsilon+IS0_n))*(1.0/(Epsilon+IS0_n));
    alpha_1n=(6.0/10.0)*(1.0/(Epsilon+IS1_n))*(1.0/(Epsilon+IS1_n));
    alpha_2n=(3.0/10.0)*(1.0/(Epsilon+IS2_n))*(1.0/(Epsilon+IS2_n));

    w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
    w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
    w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
    w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
    w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
    w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);

    u_tpph=w0_p*((2.0/6.0)*YP[i-2]-(7.0/6.0)*YP[i-1]+(11.0/6.0)*YP[i])+w1_p*((-1.0/6.0)*YP[i-1]+(5.0/6.0)*YP[i]+(2.0/6.0)*YP[i+1])+w2_p*((2.0/6.0)*YP[i]+(5.0/6.0)*YP[i+1]-(1.0/6.0)*YP[i+2]);

    u_tnph=w2_n*((-1.0/6.0)*YN[i-1]+(5.0/6.0)*YN[i]+(2.0/6.0)*YN[i+1])+w1_n*((2.0/6.0)*YN[i]+(5.0/6.0)*YN[i+1]-(1.0/6.0)*YN[i+2])+w0_n*((11.0/6.0)*YN[i+1]-(7.0/6.0)*YN[i+2]+(2.0/6.0)*YN[i+3]);

    u_tpnh=w0_p*((2.0/6.0)*YP[i-3]-(7.0/6.0)*YP[i-2]+(11.0/6.0)*YP[i-1])+w1_p*((-1.0/6.0)*YP[i-2]+(5.0/6.0)*YP[i-1]+(2.0/6.0)*YP[i])+w2_p*((2.0/6.0)*YP[i-1]+(5.0/6.0)*YP[i]-(1.0/6.0)*YP[i+1]);

    u_tnnh=w2_n*((-1.0/6.0)*YN[i-2]+(5.0/6.0)*YN[i-1]+(2.0/6.0)*YN[i])+w1_n*((2.0/6.0)*YN[i-1]+(5.0/6.0)*YN[i]-(1.0/6.0)*YN[i+1])+w0_n*((11.0/6.0)*YN[i]-(7.0/6.0)*YN[i+1]+(2.0/6.0)*YN[i+2]);

    Ydot[i]=-(1.0/udata->dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
  }
    
    i=N-3;  
    IS0_p=(13.0/12.0)*(YP[i-2]-2.0*YP[i-1]+YP[i])*(YP[i-2]-2.0*YP[i-1]+YP[i])+(1.0/4.0)*(YP[i-2]-4.0*YP[i-1]+3.0*YP[i])*(YP[i-2]-4.0*YP[i-1]+3.0*YP[i]);

    IS1_p=(13.0/12.0)*(YP[i-1]-2.0*YP[i]+YP[i+1])*(YP[i-1]-2.0*YP[i]+YP[i+1])+(1.0/4.0)*(YP[i-1]-YP[i+1])*(YP[i-1]-YP[i+1]);

    IS2_p=(13.0/12.0)*(YP[i]-2.0*YP[i+1]+YP[i+2])*(YP[i]-2.0*YP[i+1]+YP[i+2])+(1.0/4.0)*(3.0*YP[i]-4.0*YP[i+1]+YP[i+2])*(3.0*YP[i]-4.0*YP[i+1]+YP[i+2]);

    IS0_n=(13.0/12.0)*(YN[i+1]-2.0*YN[i+2]+YN[1])*(YN[i+1]-2.0*YN[i+2]+YN[1])+(1.0/4.0)*(3.0*YN[i+1]-4.0*YN[i+2]+YN[1])*(3.0*YN[i+1]-4.0*YN[i+2]+YN[1]);

    IS1_n=(13.0/12.0)*(YN[i]-2.0*YN[i+1]+YN[i+2])*(YN[i]-2.0*YN[i+1]+YN[i+2])+(1.0/4.0)*(YN[i]-YN[i+2])*(YN[i]-YN[i+2]);

    IS2_n=(13.0/12.0)*(YN[i-1]-2.0*YN[i]+YN[i+1])*(YN[i-1]-2.0*YN[i]+YN[i+1])+(1.0/4.0)*(YN[i-1]-4.0*YN[i]+3.0*YN[i+1])*(YN[i-1]-4.0*YN[i]+3.0*YN[i+1]);

    Epsilon=0.000001;
    alpha_0p=(1.0/10.0)*(1.0/(Epsilon+IS0_p))*(1.0/(Epsilon+IS0_p));
    alpha_1p=(6.0/10.0)*(1.0/(Epsilon+IS1_p))*(1.0/(Epsilon+IS1_p));
    alpha_2p=(3.0/10.0)*(1.0/(Epsilon+IS2_p))*(1.0/(Epsilon+IS2_p));
    alpha_0n=(1.0/10.0)*(1.0/(Epsilon+IS0_n))*(1.0/(Epsilon+IS0_n));
    alpha_1n=(6.0/10.0)*(1.0/(Epsilon+IS1_n))*(1.0/(Epsilon+IS1_n));
    alpha_2n=(3.0/10.0)*(1.0/(Epsilon+IS2_n))*(1.0/(Epsilon+IS2_n));

    w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
    w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
    w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
    w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
    w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
    w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);

    u_tpph=w0_p*((2.0/6.0)*YP[i-2]-(7.0/6.0)*YP[i-1]+(11.0/6.0)*YP[i])+w1_p*((-1.0/6.0)*YP[i-1]+(5.0/6.0)*YP[i]+(2.0/6.0)*YP[i+1])+w2_p*((2.0/6.0)*YP[i]+(5.0/6.0)*YP[i+1]-(1.0/6.0)*YP[i+2]);

    u_tnph=w2_n*((-1.0/6.0)*YN[i-1]+(5.0/6.0)*YN[i]+(2.0/6.0)*YN[i+1])+w1_n*((2.0/6.0)*YN[i]+(5.0/6.0)*YN[i+1]-(1.0/6.0)*YN[i+2])+w0_n*((11.0/6.0)*YN[i+1]-(7.0/6.0)*YN[i+2]+(2.0/6.0)*YN[1]);

    u_tpnh=w0_p*((2.0/6.0)*YP[i-3]-(7.0/6.0)*YP[i-2]+(11.0/6.0)*YP[i-1])+w1_p*((-1.0/6.0)*YP[i-2]+(5.0/6.0)*YP[i-1]+(2.0/6.0)*YP[i])+w2_p*((2.0/6.0)*YP[i-1]+(5.0/6.0)*YP[i]-(1.0/6.0)*YP[i+1]);

    u_tnnh=w2_n*((-1.0/6.0)*YN[i-2]+(5.0/6.0)*YN[i-1]+(2.0/6.0)*YN[i])+w1_n*((2.0/6.0)*YN[i-1]+(5.0/6.0)*YN[i]-(1.0/6.0)*YN[i+1])+w0_n*((11.0/6.0)*YN[i]-(7.0/6.0)*YN[i+1]+(2.0/6.0)*YN[i+2]);

    Ydot[i]=-(1.0/udata->dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));

    i=N-2;    
    IS0_p=(13.0/12.0)*(YP[i-2]-2.0*YP[i-1]+YP[i])*(YP[i-2]-2.0*YP[i-1]+YP[i])+(1.0/4.0)*(YP[i-2]-4.0*YP[i-1]+3.0*YP[i])*(YP[i-2]-4.0*YP[i-1]+3.0*YP[i]);

    IS1_p=(13.0/12.0)*(YP[i-1]-2.0*YP[i]+YP[i+1])*(YP[i-1]-2.0*YP[i]+YP[i+1])+(1.0/4.0)*(YP[i-1]-YP[i+1])*(YP[i-1]-YP[i+1]);

    IS2_p=(13.0/12.0)*(YP[i]-2.0*YP[i+1]+YP[1])*(YP[i]-2.0*YP[i+1]+YP[1])+(1.0/4.0)*(3.0*YP[i]-4.0*YP[i+1]+YP[1])*(3.0*YP[i]-4.0*YP[i+1]+YP[1]);

    IS0_n=(13.0/12.0)*(YN[i+1]-2.0*YN[1]+YN[2])*(YN[i+1]-2.0*YN[1]+YN[2])+(1.0/4.0)*(3.0*YN[i+1]-4.0*YN[1]+YN[2])*(3.0*YN[i+1]-4.0*YN[1]+YN[2]);

    IS1_n=(13.0/12.0)*(YN[i]-2.0*YN[i+1]+YN[1])*(YN[i]-2.0*YN[i+1]+YN[1])+(1.0/4.0)*(YN[i]-YN[1])*(YN[i]-YN[1]);

    IS2_n=(13.0/12.0)*(YN[i-1]-2.0*YN[i]+YN[i+1])*(YN[i-1]-2.0*YN[i]+YN[i+1])+(1.0/4.0)*(YN[i-1]-4.0*YN[i]+3.0*YN[i+1])*(YN[i-1]-4.0*YN[i]+3.0*YN[i+1]);

    Epsilon=0.000001;
    alpha_0p=(1.0/10.0)*(1.0/(Epsilon+IS0_p))*(1.0/(Epsilon+IS0_p));
    alpha_1p=(6.0/10.0)*(1.0/(Epsilon+IS1_p))*(1.0/(Epsilon+IS1_p));
    alpha_2p=(3.0/10.0)*(1.0/(Epsilon+IS2_p))*(1.0/(Epsilon+IS2_p));
    alpha_0n=(1.0/10.0)*(1.0/(Epsilon+IS0_n))*(1.0/(Epsilon+IS0_n));
    alpha_1n=(6.0/10.0)*(1.0/(Epsilon+IS1_n))*(1.0/(Epsilon+IS1_n));
    alpha_2n=(3.0/10.0)*(1.0/(Epsilon+IS2_n))*(1.0/(Epsilon+IS2_n));

    w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
    w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
    w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
    w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
    w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
    w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);

    u_tpph=w0_p*((2.0/6.0)*YP[i-2]-(7.0/6.0)*YP[i-1]+(11.0/6.0)*YP[i])+w1_p*((-1.0/6.0)*YP[i-1]+(5.0/6.0)*YP[i]+(2.0/6.0)*YP[i+1])+w2_p*((2.0/6.0)*YP[i]+(5.0/6.0)*YP[i+1]-(1.0/6.0)*YP[1]);

    u_tnph=w2_n*((-1.0/6.0)*YN[i-1]+(5.0/6.0)*YN[i]+(2.0/6.0)*YN[i+1])+w1_n*((2.0/6.0)*YN[i]+(5.0/6.0)*YN[i+1]-(1.0/6.0)*YN[1])+w0_n*((11.0/6.0)*YN[i+1]-(7.0/6.0)*YN[1]+(2.0/6.0)*YN[2]);

    u_tpnh=w0_p*((2.0/6.0)*YP[i-3]-(7.0/6.0)*YP[i-2]+(11.0/6.0)*YP[i-1])+w1_p*((-1.0/6.0)*YP[i-2]+(5.0/6.0)*YP[i-1]+(2.0/6.0)*YP[i])+w2_p*((2.0/6.0)*YP[i-1]+(5.0/6.0)*YP[i]-(1.0/6.0)*YP[i+1]);

    u_tnnh=w2_n*((-1.0/6.0)*YN[i-2]+(5.0/6.0)*YN[i-1]+(2.0/6.0)*YN[i])+w1_n*((2.0/6.0)*YN[i-1]+(5.0/6.0)*YN[i]-(1.0/6.0)*YN[i+1])+w0_n*((11.0/6.0)*YN[i]-(7.0/6.0)*YN[i+1]+(2.0/6.0)*YN[1]);

    Ydot[i]=-(1.0/udata->dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));

    i=N-1;
    IS0_p=(13.0/12.0)*(YP[i-2]-2.0*YP[i-1]+YP[i])*(YP[i-2]-2.0*YP[i-1]+YP[i])+(1.0/4.0)*(YP[i-2]-4.0*YP[i-1]+3.0*YP[i])*(YP[i-2]-4.0*YP[i-1]+3.0*YP[i]);

    IS1_p=(13.0/12.0)*(YP[i-1]-2.0*YP[i]+YP[1])*(YP[i-1]-2.0*YP[i]+YP[1])+(1.0/4.0)*(YP[i-1]-YP[1])*(YP[i-1]-YP[1]);

    IS2_p=(13.0/12.0)*(YP[i]-2.0*YP[1]+YP[2])*(YP[i]-2.0*YP[1]+YP[2])+(1.0/4.0)*(3.0*YP[i]-4.0*YP[1]+YP[2])*(3.0*YP[i]-4.0*YP[1]+YP[2]);

    IS0_n=(13.0/12.0)*(YN[1]-2.0*YN[2]+YN[3])*(YN[1]-2.0*YN[2]+YN[3])+(1.0/4.0)*(3.0*YN[1]-4.0*YN[2]+YN[3])*(3.0*YN[1]-4.0*YN[2]+YN[3]);

    IS1_n=(13.0/12.0)*(YN[i]-2.0*YN[1]+YN[2])*(YN[i]-2.0*YN[1]+YN[2])+(1.0/4.0)*(YN[i]-YN[2])*(YN[i]-YN[2]);

    IS2_n=(13.0/12.0)*(YN[i-1]-2.0*YN[i]+YN[1])*(YN[i-1]-2.0*YN[i]+YN[1])+(1.0/4.0)*(YN[i-1]-4.0*YN[i]+3.0*YN[1])*(YN[i-1]-4.0*YN[i]+3.0*YN[1]);

    Epsilon=0.000001;
    alpha_0p=(1.0/10.0)*(1.0/(Epsilon+IS0_p))*(1.0/(Epsilon+IS0_p));
    alpha_1p=(6.0/10.0)*(1.0/(Epsilon+IS1_p))*(1.0/(Epsilon+IS1_p));
    alpha_2p=(3.0/10.0)*(1.0/(Epsilon+IS2_p))*(1.0/(Epsilon+IS2_p));
    alpha_0n=(1.0/10.0)*(1.0/(Epsilon+IS0_n))*(1.0/(Epsilon+IS0_n));
    alpha_1n=(6.0/10.0)*(1.0/(Epsilon+IS1_n))*(1.0/(Epsilon+IS1_n));
    alpha_2n=(3.0/10.0)*(1.0/(Epsilon+IS2_n))*(1.0/(Epsilon+IS2_n));

    w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
    w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
    w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
    w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
    w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
    w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);

    u_tpph=w0_p*((2.0/6.0)*YP[i-2]-(7.0/6.0)*YP[i-1]+(11.0/6.0)*YP[i])+w1_p*((-1.0/6.0)*YP[i-1]+(5.0/6.0)*YP[i]+(2.0/6.0)*YP[1])+w2_p*((2.0/6.0)*YP[i]+(5.0/6.0)*YP[1]-(1.0/6.0)*YP[2]);

    u_tnph=w2_n*((-1.0/6.0)*YN[i-1]+(5.0/6.0)*YN[i]+(2.0/6.0)*YN[1])+w1_n*((2.0/6.0)*YN[i]+(5.0/6.0)*YN[1]-(1.0/6.0)*YN[2])+w0_n*((11.0/6.0)*YN[1]-(7.0/6.0)*YN[2]+(2.0/6.0)*YN[3]);

    u_tpnh=w0_p*((2.0/6.0)*YP[i-3]-(7.0/6.0)*YP[i-2]+(11.0/6.0)*YP[i-1])+w1_p*((-1.0/6.0)*YP[i-2]+(5.0/6.0)*YP[i-1]+(2.0/6.0)*YP[i])+w2_p*((2.0/6.0)*YP[i-1]+(5.0/6.0)*YP[i]-(1.0/6.0)*YP[1]);

    u_tnnh=w2_n*((-1.0/6.0)*YN[i-2]+(5.0/6.0)*YN[i-1]+(2.0/6.0)*YN[i])+w1_n*((2.0/6.0)*YN[i-1]+(5.0/6.0)*YN[i]-(1.0/6.0)*YN[1])+w0_n*((11.0/6.0)*YN[i]-(7.0/6.0)*YN[1]+(2.0/6.0)*YN[2]);

    Ydot[i]=-(1.0/udata->dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
    
 /* Free vectors */
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(yn);

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
             NULL pointer  
*/

static int check_flag(void *flagvalue, const string funcname, int opt)
{
  int *errflag;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == NULL) {
    cerr << "\nSUNDIALS_ERROR: " << funcname << " failed - returned NULL pointer\n\n";
    return 1; }

  // Check if flag < 0
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      cerr << "\nSUNDIALS_ERROR: " << funcname << " failed with flag = " << *errflag << "\n\n";
      return 1; 
    }
  }
  
  // Check if function returned NULL pointer - no memory allocated
  else if (opt == 2 && flagvalue == NULL) {
    cerr << "\nMEMORY_ERROR: " << funcname << " failed - returned NULL pointer\n\n";
    return 1; }

  return 0;
}

/*---- end of file ----*/

