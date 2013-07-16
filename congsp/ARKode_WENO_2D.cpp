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

/* accessor macros between (i,j,Nx,Ny,k) location and 1D NVector array */
#define idx(i,j,Nx,Ny,k) ((k)*(Nx)*(Ny)+(j)*(Nx)+i)
/* accessor macros between (i,j,Nx) location and 1D NVector array */
#define idx_v(i,j,Nx) ((j)*(Nx)+i)

/* constants */
#define PI RCONST(3.1415926535897932)
#define ONE RCONST(1.0)
#define TWO RCONST(2.0)

/* user data structure */
typedef struct {
    long int Nx;   /* number of x grids        */
    long int Ny;   /* number of x grids        */
    realtype dx;   /* x direction mesh spacing */
    realtype dy;   /* y direction mesh spacing */
    realtype Lx;   /* max value in x direction */
    realtype Ly;   /* max value in y direction */
    realtype gama; /* gas constant             */
} *UserData;

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const string funcname, int opt);

/* Set value of tao in whole domain*/
//static int Gettao(N_Vector yxup, N_Vector yxdown, N_Vector yyup, N_Vector yydown, N_Vector taoup, N_Vector taodown, long int Nx, long int Ny, realtype gama);
static int Gettao(realtype *yxupdata, realtype *yxdowndata, realtype *yyupdata, realtype *yydowndata, N_Vector taoup, N_Vector taodown, long int Nx, long int Ny, realtype gama);

/* Set value of J in whole domain*/
//static int GetCj(N_Vector yxup, N_Vector yxdown, N_Vector yyup, N_Vector yydown, N_Vector Cjup, N_Vector Cjdown, long int Nx, long int Ny, realtype gama);
static int GetCj(realtype *yxupdata, realtype *yxdowndata, realtype *yyupdata, realtype *yydowndata, N_Vector Cjup, N_Vector Cjdown, long int Nx, long int Ny, realtype gama);

/* Set left eigenvectors in x direction */
//static int Setlfxegm(realtype *Ydata, realtype **lfxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype egv1x, realtype egv2x, realtype egv3x, realtype egv4x, int flag);
static int Setlfxegm(realtype *Ydata, realtype **lfxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype *egvx, int flag);

/* Set right eigenvectors in x direction */
static int Setrhxegm(realtype *Ydata, realtype **rhxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag);

/* Set left eigenvectors in y direction */
//static int Setlfyegm(realtype *Ydata, realtype **lfyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype egv1y, realtype egv2y, realtype egv3y, realtype egv4y, int flag);
static int Setlfyegm(realtype *Ydata, realtype **lfyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype *egvy, int flag);

/* Set right eigenvectors in y direction */
static int Setrhyegm(realtype *Ydata, realtype **rhyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag);

/* Fill in the indicators of smoothness on x direction */
static int SetISX(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx,realtype *yxpdata, realtype *yxndata, long int i, long int j, long int Nx, long int Ny, int flag);

/* Fill in the indicators of smoothness on y direction */
static int SetISY(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny,realtype *yypdata, realtype *yyndata, long int i, long int j, long int Nx, long int Ny, int flag);

/* Fill in the stencil weights on x direction */
static int Setalphawx(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx, realtype *alpha_0px, realtype *alpha_1px, realtype *alpha_2px, realtype *alpha_0nx, realtype *alpha_1nx, realtype *alpha_2nx, realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype Epsilon);

/* Fill in the stencil weights on y direction */
static int Setalphawy(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny, realtype *alpha_0py, realtype *alpha_1py, realtype *alpha_2py, realtype *alpha_0ny, realtype *alpha_1ny, realtype *alpha_2ny, realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype Epsilon);

/* Get the derivative on x direction */
static int SetUx(realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype *yxpdata, realtype *yxndata, realtype *u_tpphx, realtype *u_tnphx, realtype *u_tpnhx, realtype *u_tnnhx, long int i, long int j, long int Nx, long int Ny, int flag);

/* Get the derivative on y direction */
static int SetUy(realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype *yypdata, realtype *yyndata, realtype *u_tpphy, realtype *u_tnphy, realtype *u_tpnhy, realtype *u_tnnhy, long int i, long int j, long int Nx, long int Ny, int flag);

int main(int argc, const char * argv[])
{
    /* general problem parameters */
    realtype T0 = RCONST(0.0);
    realtype Tf = RCONST(0.1);
    int Nt = 1;
    int Nvar = 4;
    UserData udata = NULL;
    realtype *data;
    long int Nx, Ny, NEQ, i, j;
    realtype Lx, Ly, gama, p, r;
    
    /* declare solver parameters */
    int flag;
    
    /* general problem variables */
    N_Vector y = NULL;
    void *arkode_mem = NULL;
    
    /* allocate udata structure */
    udata = (UserData) malloc(sizeof(*udata));
    if (check_flag((void *) udata, "malloc", 2)) return 1;
    
    /* read problem parameter and tolerances from input file:*/
    double dx, dy;
    FILE *FID;
    FID = fopen("input_WENO2D.txt","r");
    flag = fscanf(FID," Nx = %li\n", &Nx);
    flag = fscanf(FID," Ny = %li\n", &Ny);
    flag = fscanf(FID," Lx = %lf\n", &Lx);
    flag = fscanf(FID," Ly = %lf\n", &Ly);
    flag = fscanf(FID," gama = %lf\n", &gama);
    fclose(FID);
    
    /* store the inputs in the UserData structure */
    udata->Nx = Nx;
    udata->Ny = Ny;
    udata->Lx = Lx;
    udata->Ly = Ly;
    udata->gama = gama;
    
    /* open solver diagnostics output file for writing */
    //FILE *DFID;
    //DFID=fopen("diags_ark_WENO2D.txt","w");
    
    /* set total allocated vector length */
    NEQ = Nvar*udata->Nx*udata->Ny;
    
    /* Initial problem output */
    printf("\n2D gas dynamic test problem:\n");
    printf("    Nx = %li,  Ny = %li, NEQ = %li\n", udata->Nx, udata->Ny, NEQ);
    printf("    problem parameters:  Lx = %g,  Ly = %g\n, gama = %g\n", udata->Lx, udata->Ly, udata->gama);
    
    /* Create serial vector of length NEQ for initial condition */
    y = N_VNew_Serial(NEQ);
    if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;
    
    /* set spatial mesh spacing */
    udata->dx = Lx/Nx;
    udata->dy = Ly/Ny;
    
    /* output x_mesh to disk */
    FID=fopen("WENO2D_meshx.txt","w");
    for (i=0; i<Nx; i++)  fprintf(FID,"  %.16e\n", udata->dx*(i+0.5));
    fclose(FID);
    
    /* output y_mesh to disk */
    FID=fopen("WENO2D_meshy.txt","w");
    for (i=0; i<Ny; i++)  fprintf(FID,"  %.16e\n", udata->dy*(i+0.5));
    fclose(FID);
    
    /* Access data array for new NVector y */
    data = N_VGetArrayPointer(y);
    if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
    
    for(j=0;j<Ny;j++){
        for (i=0; i<Nx; i++) {
             if (udata->dx*(i+0.5)<0.5&&udata->dy*(j+0.5)>0.5){
            data[idx(i,j,Nx,Ny,0)] =  0.5323;  /* rou */
            data[idx(i,j,Nx,Ny,1)] =  1.206*0.5323;  /* qx */
            data[idx(i,j,Nx,Ny,2)] =  0.5323*0.0;  /* qy */
            data[idx(i,j,Nx,Ny,3)] =  0.5323*(0.3/((udata->gama-1.0)*0.5323)+0.5*(1.206*1.206));  /* E */
 }
 if (udata->dx*(i+0.5)>0.5&&udata->dy*(j+0.5)>0.5){
            data[idx(i,j,Nx,Ny,0)] =  1.5;  /* rou */
            data[idx(i,j,Nx,Ny,1)] =  1.5*0.0;  /* qx */
            data[idx(i,j,Nx,Ny,2)] =  1.5*0.0;  /* qy */
            data[idx(i,j,Nx,Ny,3)] =  1.5*(1.5/((udata->gama-1.0)*1.5));  /* E */
 }
 if (udata->dx*(i+0.5)<0.5&&udata->dy*(j+0.5)<0.5){
            data[idx(i,j,Nx,Ny,0)] =  0.138;  /* rou */
            data[idx(i,j,Nx,Ny,1)] =  1.206*0.138;  /* qx */
            data[idx(i,j,Nx,Ny,2)] =  1.206*0.138;  /* qy */
            data[idx(i,j,Nx,Ny,3)] =  0.138*(0.029/((udata->gama-1.0)*0.138)+0.5*(2*1.206*1.206));  /* E */
 }
 if (udata->dx*(i+0.5)>0.5&&udata->dy*(j+0.5)<0.5){
            data[idx(i,j,Nx,Ny,0)] =  0.5323;  /* rou */
            data[idx(i,j,Nx,Ny,1)] =  1.206*0.0;  /* qx */
            data[idx(i,j,Nx,Ny,2)] =  1.206*0.5323;  /* qy */
            data[idx(i,j,Nx,Ny,3)] =  0.5323*(0.3/((udata->gama-1.0)*0.5323)+0.5*(1.206*1.206));  /* E */
 }

        }
    }
    
    /* Call ARKodeCreate to create the solver memory */
    arkode_mem = ARKodeCreate();
    if (check_flag((void *) arkode_mem, "ARKodeCreate", 0)) return 1;
    
    /* Set solver parameters */
    realtype reltol = 1.0e+12;
    realtype abstol = 1.0e+12;
    //realtype reltol  = 1.e-3;
    //realtype abstol  = 1.e-6;
    
    /* Solution uses default explicit method */
    flag = ARKodeInit(arkode_mem, f, NULL, T0, y);
    if (check_flag(&flag, "ARKodeInit", 1)) return 1;
    
    /* Call ARKodeSetUserData to pass rdata to user functions */
    flag = ARKodeSetUserData(arkode_mem, (void *) udata);
    if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;
    
    /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
    //flag = ARKodeSetDiagnostics(arkode_mem, DFID);
    //if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return 1;
    
    /* Call ARKodeSetInitStep to initialize time step */
    flag = ARKodeSetInitStep(arkode_mem, 0.1);
    if (check_flag(&flag, "ARKodeSetInitStep", 1)) return 1;
    
    /* Call ARKodeSetMinStep to set min time step */
    flag = ARKodeSetMinStep(arkode_mem, 0.1);
    if (check_flag(&flag, "ARKodeSetMinStep", 1)) return 1;
    
    /* Call ARKodeSetMaxStep to set max time step */
    flag = ARKodeSetMaxStep(arkode_mem, 0.1);
    if (check_flag(&flag, "ARKodeSetMaxStep", 1)) return 1;
    
    /* Call ARKodeSetMaxNumSteps to increase default (for testing) */
    flag = ARKodeSetMaxNumSteps(arkode_mem, 100000);
    if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;
    
    /* Call ARKodeSStolerances to specify the scalar relative and absolute
     tolerances */
    flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
    if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;
    
    /* Open output stream for results, access data arrays */
    FILE *RFID=fopen("2DWENO_rou.txt","w");
    FILE *QXFID=fopen("2DWENO_qx.txt","w");
    FILE *QYFID=fopen("2DWENO_qy.txt","w");
    FILE *EFID=fopen("2DWENO_E.txt","w");
    data = N_VGetArrayPointer(y);
    if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
    
    /* output initial condition to disk */
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            fprintf(RFID," %.16e", data[idx(i,j,Nx,Ny,0)]);
            fprintf(QXFID," %.16e", data[idx(i,j,Nx,Ny,1)]);
            fprintf(QYFID," %.16e", data[idx(i,j,Nx,Ny,2)]);
            fprintf(EFID," %.16e", data[idx(i,j,Nx,Ny,3)]);
        }
    }
    fprintf(RFID,"\n");
    fprintf(QXFID,"\n");
    fprintf(QYFID,"\n");
    fprintf(EFID,"\n");
    
    /* Write all solver parameters to stdout */
    // printf("\n");
    //flag = ARKodeWriteParameters(arkode_mem, stdout);
    //if (check_flag(&flag, "ARKodeWriteParameters", 1)) return 1;
    
    /* In loop, call ARKode, print results, and test for error.
     Break out of loop when the final output time has been reached */
    realtype t  = T0;
    realtype dTout = Tf/Nt;
    realtype tout = T0+dTout;
    int iout;
    for (iout=0; iout<Nt; iout++) {
        
        /* stop exactly at this time */
        //flag = ARKodeSetStopTime(arkode_mem, tout);
        flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
        if (check_flag(&flag, "ARKode", 1)) break;
        if (flag >= 0) {
            tout += dTout;
            tout = (tout > Tf) ? Tf : tout;
        }
        /* output results to disk */
        for(j=0;j<Ny;j++){
            for(i=0;i<Nx;i++){
                fprintf(RFID," %.16e", data[idx(i,j,Nx,Ny,0)]);
                fprintf(QXFID," %.16e", data[idx(i,j,Nx,Ny,1)]);
                fprintf(QYFID," %.16e", data[idx(i,j,Nx,Ny,2)]);
                fprintf(EFID," %.16e", data[idx(i,j,Nx,Ny,3)]);
            }
        }
        fprintf(RFID,"\n");
        fprintf(QXFID,"\n");
        fprintf(QYFID,"\n");
        fprintf(EFID,"\n");
    }
    fclose(RFID);
    fclose(QXFID);
    fclose(QYFID);
    fclose(EFID);
    
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


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* f routine to compute the ODE RHS function f(t,y). */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    /* declare variables */
    long int NEQ, NEQS, i, j, k;
    //realtype egv1x, egv2x, egv3x, egv4x, egvmaxtempx, egvmaxx, egv1xtmx, egv3xtmx;
    //realtype egv1y, egv2y, egv3y, egv4y, egvmaxtempy, egvmaxy, egv1ytmy, egv3ytmy;
    int flag;
    realtype p, Epsilon;
    Epsilon = 0.000001;
    
    /* create relative arrays */
    realtype *IS0_px = new realtype [4];
    realtype *IS1_px = new realtype [4];
    realtype *IS2_px = new realtype [4];
    realtype *IS0_nx = new realtype [4];
    realtype *IS1_nx = new realtype [4];
    realtype *IS2_nx = new realtype [4];
    realtype *IS0_py = new realtype [4];
    realtype *IS1_py = new realtype [4];
    realtype *IS2_py = new realtype [4];
    realtype *IS0_ny = new realtype [4];
    realtype *IS1_ny = new realtype [4];
    realtype *IS2_ny = new realtype [4];
    
    realtype *alpha_0px = new realtype [4];
    realtype *alpha_1px = new realtype [4];
    realtype *alpha_2px = new realtype [4];
    realtype *alpha_0nx = new realtype [4];
    realtype *alpha_1nx = new realtype [4];
    realtype *alpha_2nx = new realtype [4];
    realtype *alpha_0py = new realtype [4];
    realtype *alpha_1py = new realtype [4];
    realtype *alpha_2py = new realtype [4];
    realtype *alpha_0ny = new realtype [4];
    realtype *alpha_1ny = new realtype [4];
    realtype *alpha_2ny = new realtype [4];
    
    realtype *w0_px = new realtype [4];
    realtype *w1_px = new realtype [4];
    realtype *w2_px = new realtype [4];
    realtype *w0_nx = new realtype [4];
    realtype *w1_nx = new realtype [4];
    realtype *w2_nx = new realtype [4];
    realtype *w0_py = new realtype [4];
    realtype *w1_py = new realtype [4];
    realtype *w2_py = new realtype [4];
    realtype *w0_ny = new realtype [4];
    realtype *w1_ny = new realtype [4];
    realtype *w2_ny = new realtype [4];
    
    realtype *u_tpphx = new realtype [4];
    realtype *u_tpnhx = new realtype [4];
    realtype *u_tnphx = new realtype [4];
    realtype *u_tnnhx = new realtype [4];
    realtype *u_tpphy = new realtype [4];
    realtype *u_tpnhy = new realtype [4];
    realtype *u_tnphy = new realtype [4];
    realtype *u_tnnhy = new realtype [4];
    
    realtype **lfxegm = new realtype *[4];
    for (i=0;i<4;i++){
        lfxegm[i] = new realtype [4];
    }
    
    realtype **rhxegm = new realtype *[4];
    for (i=0;i<4;i++){
        rhxegm[i] = new realtype [4];
    }
    
    realtype **lfyegm = new realtype *[4];
    for (i=0;i<4;i++){
        lfyegm[i] = new realtype [4];
    }
    
    realtype **rhyegm = new realtype *[4];
    for (i=0;i<4;i++){
        rhyegm[i] = new realtype [4];
    }
    
    /* clear out ydot (to be careful) */
    N_VConst(0.0, ydot);
    
    /* problem data */
    UserData udata = (UserData) user_data;
    
    /* shortcuts to number of grads, background values */
    long int Nx = udata->Nx;
    long int Ny = udata->Ny;
    realtype dx = udata->dx;
    realtype dy = udata->dy;
    realtype Lx = udata->Lx;
    realtype Ly = udata->Ly;
    realtype gama = udata->gama;
    
    /* fill in the value of NEQ and NEQS */
    NEQ = 4*Nx*Ny;
    NEQS = 2*Nx*Ny;
    
    /* create relative vectors */
    N_Vector taodown = NULL;
    N_Vector taoup = NULL;
    N_Vector Cjdown = NULL;
    N_Vector Cjup = NULL;
    /*
    N_Vector yxup = NULL;
    N_Vector yxdown = NULL;
    N_Vector yyup = NULL;
    N_Vector yydown = NULL;
    N_Vector yxp = NULL;
    N_Vector yxn = NULL;
    N_Vector yyp = NULL;
    N_Vector yyn = NULL;
    N_Vector yxnew = NULL;
    N_Vector yynew = NULL;
    N_Vector yxupback = NULL;
    N_Vector yyupback = NULL;
    N_Vector yxdownback = NULL;
    N_Vector yydownback = NULL;
    */
    
    /* Create serial vector of length NEQ and NEQS */
    //  yxup = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yxup, "N_VNew_Serial", 0)) return 1;
    //yyup = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yyup, "N_VNew_Serial", 0)) return 1;
    taoup = N_VNew_Serial(NEQ);
    if (check_flag((void *) taoup, "N_VNew_Serial", 0)) return 1;
    Cjup = N_VNew_Serial(NEQS);
    if (check_flag((void *) Cjup, "N_VNew_Serial", 0)) return 1;
    //yxdown = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yxdown, "N_VNew_Serial", 0)) return 1;
    //yydown = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yydown, "N_VNew_Serial", 0)) return 1;
    taodown = N_VNew_Serial(NEQ);
    if (check_flag((void *) taodown, "N_VNew_Serial", 0)) return 1;
    Cjdown = N_VNew_Serial(NEQS);
    if (check_flag((void *) Cjdown, "N_VNew_Serial", 0)) return 1;
    //yxp = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yxp, "N_VNew_Serial", 0)) return 1;
    //yxn = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yxn, "N_VNew_Serial", 0)) return 1;
    //yyp = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yyp, "N_VNew_Serial", 0)) return 1;
    //yyn = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yyn, "N_VNew_Serial", 0)) return 1;
    //yxnew = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yxnew, "N_VNew_Serial", 0)) return 1;
    //yynew = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yynew, "N_VNew_Serial", 0)) return 1;
    //yxupback = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yxupback, "N_VNew_Serial", 0)) return 1;
    //yyupback = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yyupback, "N_VNew_Serial", 0)) return 1;
    //yxdownback = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yxdownback, "N_VNew_Serial", 0)) return 1;
    //yydownback = N_VNew_Serial(NEQ);
    //if (check_flag((void *) yydownback, "N_VNew_Serial", 0)) return 1;
    
    /* fill in the value of tao and J in the whole domain */
  //  flag=Gettao(y, tao, Nx, Ny, gama);
  //  if (flag!=0) printf("error in Gettao function \n");
  //  flag=GetCj(y, Cj, Nx, Ny, gama);
  //  if (flag!=0) printf("error in GetCj function \n");
    
    realtype *yxupdata = new realtype [4*Nx*Ny];
    realtype *yyupdata = new realtype [4*Nx*Ny];
    realtype *yxdowndata = new realtype [4*Nx*Ny];
    realtype *yydowndata = new realtype [4*Nx*Ny];
    realtype *yxpdata = new realtype [4*Nx*Ny];
    realtype *yxndata = new realtype [4*Nx*Ny];
    realtype *yypdata = new realtype [4*Nx*Ny];
    realtype *yyndata = new realtype [4*Nx*Ny];
    realtype *yxnewdata = new realtype [4*Nx*Ny];
    realtype *yynewdata = new realtype [4*Nx*Ny];
    realtype *yxupbackdata = new realtype [4*Nx*Ny];
    realtype *yyupbackdata = new realtype [4*Nx*Ny];
    realtype *yxdownbackdata = new realtype [4*Nx*Ny];
    realtype *yydownbackdata = new realtype [4*Nx*Ny];
    realtype *egvx = new realtype [4*Nx*Ny];
    realtype *egvy = new realtype [4*Nx*Ny];
    realtype *egvxmax = new realtype [Nx*Ny];
    realtype *egvymax = new realtype [Nx*Ny];

    /* access data arrays */
    realtype *Ydata = N_VGetArrayPointer(y);
    if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;
    realtype *dYdata = N_VGetArrayPointer(ydot);
    if (check_flag((void *) dYdata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yxupdata = N_VGetArrayPointer(yxup);
    //if (check_flag((void *) yxupdata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yyupdata = N_VGetArrayPointer(yyup);
    //if (check_flag((void *) yyupdata, "N_VGetArrayPointer", 0)) return 1;
    realtype *taoupdata = N_VGetArrayPointer(taoup);
    if (check_flag((void *) taoupdata, "N_VGetArrayPointer", 0)) return 1;
    realtype *Cjupdata = N_VGetArrayPointer(Cjup);
    if (check_flag((void *) Cjupdata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yxdowndata = N_VGetArrayPointer(yxdown);
    //if (check_flag((void *) yxdowndata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yydowndata = N_VGetArrayPointer(yydown);
    //if (check_flag((void *) yydowndata, "N_VGetArrayPointer", 0)) return 1;
    realtype *taodowndata = N_VGetArrayPointer(taodown);
    if (check_flag((void *) taodowndata, "N_VGetArrayPointer", 0)) return 1;
    realtype *Cjdowndata = N_VGetArrayPointer(Cjdown);
    if (check_flag((void *) Cjdowndata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yxpdata = N_VGetArrayPointer(yxp);
    //if (check_flag((void *) yxpdata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yxndata = N_VGetArrayPointer(yxn);
    //if (check_flag((void *) yxndata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yypdata = N_VGetArrayPointer(yyp);
    //if (check_flag((void *) yypdata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yyndata = N_VGetArrayPointer(yyn);
    //if (check_flag((void *) yyndata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yxnewdata = N_VGetArrayPointer(yxnew);
    //if (check_flag((void *) yxnewdata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yynewdata = N_VGetArrayPointer(yynew);
    //if (check_flag((void *) yynewdata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yxupbackdata = N_VGetArrayPointer(yxupback);
    //if (check_flag((void *) yxupbackdata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yyupbackdata = N_VGetArrayPointer(yyupback);
    //if (check_flag((void *) yyupbackdata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yxdownbackdata = N_VGetArrayPointer(yxdownback);
    //if (check_flag((void *) yxdownbackdata, "N_VGetArrayPointer", 0)) return 1;
    //realtype *yydownbackdata = N_VGetArrayPointer(yydownback);
    //if (check_flag((void *) yydownbackdata, "N_VGetArrayPointer", 0)) return 1;
    
    /* compute max absolue eigenvalue and fill in ypdata, yndata, taopdata, taondata, Cjpdata, Cjpdata */
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
	    // flag = Setlfxegm(Ydata, lfxegm, i, j, Nx, Ny, gama, egv1x, egv2x, egv3x, egv4x, 2);
	    flag = Setlfxegm(Ydata, lfxegm, i, j, Nx, Ny, gama, egvx, 2);
            if (flag!=0) printf("error in Setlfxegm function \n");
	   
	    //if (i==0&&j==0)
	    //printf(" Ydata : i = %li, j = %li, Ydata[idx(i, j, Nx, Ny, 0)] = %f, Ydata[idx(i, j, Nx, Ny, 1)] = %f, Ydata[idx(i, j, Nx, Ny, 2)] = %f, Ydata[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, Ydata[idx(i, j, Nx, Ny, 0)], Ydata[idx(i, j, Nx, Ny, 1)],Ydata[idx(i, j, Nx, Ny, 2)],Ydata[idx(i, j, Nx, Ny, 3)]);
		      //printf(" maxeigenccccccccccccccc : i = %li, j = %li, egv1 = %f, egv2 = %f, egv3 = %f, egv4 = %f, egvmaxx = %f\n", i, j, egvx[idx(i,j,Nx,Ny,0)], egvx[idx(i,j,Nx,Ny,1)],egvx[idx(i,j,Nx,Ny,2)],egvx[idx(i,j,Nx,Ny,3)], egvxmax[idx_v(i,j,Nx)]);
            for (k=0;k<4;k++){
	      //if (i==0&&j==0){
	      //printf(" lfxegm : i = %li, j = %li, lfxegm[k][0] = %f, lfxegm[k][1] = %f, lfxegm[k][2] = %f, lfxegm[k][3] = %f\n", i, j, lfxegm[k][0], lfxegm[k][1],lfxegm[k][2],lfxegm[k][3]);
	      //}
                yxnewdata[idx(i, j, Nx, Ny, k)] = lfxegm[k][0]*Ydata[idx(i, j, Nx, Ny, 0)]+lfxegm[k][1]*Ydata[idx(i, j, Nx, Ny, 1)]+lfxegm[k][2]*Ydata[idx(i, j, Nx, Ny, 2)]+lfxegm[k][3]*Ydata[idx(i, j, Nx, Ny, 3)];
            }
            //printf("   yxnew problem parameters: i = %li, j = %li, yxnew0 = %g,  yxnew1 = %g, yxnew2 = %g,  yxnew3 = %g\n", i, j, yxnewdata[idx(i, j, Nx, Ny, 0)], yxnewdata[idx(i, j, Nx, Ny, 1)], yxnewdata[idx(i, j, Nx, Ny, 2)], yxnewdata[idx(i, j, Nx, Ny, 3)]);
           
            egvxmax[idx_v(i,j,Nx)] = (fabs(egvx[idx(i,j,Nx,Ny,0)])>fabs(egvx[idx(i,j,Nx,Ny,2)])) ? fabs(egvx[idx(i,j,Nx,Ny,0)]) : fabs(egvx[idx(i,j,Nx,Ny,2)]);
	 
	    //if (i==0&&j==0){    
	    //printf(" maxeigen : i = %li, j = %li, egv1 = %f, egv2 = %f, egv3 = %f, egv4 = %f, egvmaxx = %f\n", i, j, egvx[idx(i,j,Nx,Ny,0)], egvx[idx(i,j,Nx,Ny,1)],egvx[idx(i,j,Nx,Ny,2)],egvx[idx(i,j,Nx,Ny,3)], egvxmax[idx_v(i,j,Nx)]);
	    //}
            yxpdata[idx(i, j, Nx, Ny, 0)]=0.5*(egvx[idx(i,j,Nx,Ny,0)]*yxnewdata[idx(i, j, Nx, Ny, 0)]+egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 0)]);
            yxndata[idx(i, j, Nx, Ny, 0)]=0.5*(egvx[idx(i,j,Nx,Ny,0)]*yxnewdata[idx(i, j, Nx, Ny, 0)]-egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 0)]);
            yxpdata[idx(i, j, Nx, Ny, 1)]=0.5*(egvx[idx(i,j,Nx,Ny,1)]*yxnewdata[idx(i, j, Nx, Ny, 1)]+egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 1)]);
            yxndata[idx(i, j, Nx, Ny, 1)]=0.5*(egvx[idx(i,j,Nx,Ny,1)]*yxnewdata[idx(i, j, Nx, Ny, 1)]-egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 1)]);
            yxpdata[idx(i, j, Nx, Ny, 2)]=0.5*(egvx[idx(i,j,Nx,Ny,2)]*yxnewdata[idx(i, j, Nx, Ny, 2)]+egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 2)]);
            yxndata[idx(i, j, Nx, Ny, 2)]=0.5*(egvx[idx(i,j,Nx,Ny,2)]*yxnewdata[idx(i, j, Nx, Ny, 2)]-egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 2)]);
            yxpdata[idx(i, j, Nx, Ny, 3)]=0.5*(egvx[idx(i,j,Nx,Ny,3)]*yxnewdata[idx(i, j, Nx, Ny, 3)]+egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 3)]);
            yxndata[idx(i, j, Nx, Ny, 3)]=0.5*(egvx[idx(i,j,Nx,Ny,3)]*yxnewdata[idx(i, j, Nx, Ny, 3)]-egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 3)]);
        }
    }
    
    for(i=0; i<Nx; i++){
        for (j=0; j<Ny; j++){
	  //flag = Setlfyegm(Ydata, lfyegm, i, j, Nx, Ny, gama, egv1y, egv2y, egv3y, egv4y, 2);
	  flag = Setlfyegm(Ydata, lfyegm, i, j, Nx, Ny, gama, egvy, 2);
            if (flag!=0) printf("error in Setlfyegm function \n");
            for (k=0;k<4;k++){
                yynewdata[idx(i, j, Nx, Ny, k)] = lfyegm[k][0]*Ydata[idx(i, j, Nx, Ny, 0)]+lfyegm[k][1]*Ydata[idx(i, j, Nx, Ny, 1)]+lfyegm[k][2]*Ydata[idx(i, j, Nx, Ny, 2)]+lfyegm[k][3]*Ydata[idx(i, j, Nx, Ny, 3)];
            }
            
	    egvymax[idx_v(i,j,Nx)] = (fabs(egvy[idx(i,j,Nx,Ny,0)])>fabs(egvy[idx(i,j,Nx,Ny,2)])) ? fabs(egvy[idx(i,j,Nx,Ny,0)]) : fabs(egvy[idx(i,j,Nx,Ny,2)]);
            
            yypdata[idx(i, j, Nx, Ny, 0)]=0.5*(egvy[idx(i,j,Nx,Ny,0)]*yynewdata[idx(i, j, Nx, Ny, 0)]+egvymax[idx_v(i,j,Nx)]*yynewdata[idx(i, j, Nx, Ny, 0)]);
            yyndata[idx(i, j, Nx, Ny, 0)]=0.5*(egvy[idx(i,j,Nx,Ny,0)]*yynewdata[idx(i, j, Nx, Ny, 0)]-egvymax[idx_v(i,j,Nx)]*yynewdata[idx(i, j, Nx, Ny, 0)]);
            yypdata[idx(i, j, Nx, Ny, 1)]=0.5*(egvy[idx(i,j,Nx,Ny,1)]*yynewdata[idx(i, j, Nx, Ny, 1)]+egvymax[idx_v(i,j,Nx)]*yynewdata[idx(i, j, Nx, Ny, 1)]);
            yyndata[idx(i, j, Nx, Ny, 1)]=0.5*(egvy[idx(i,j,Nx,Ny,1)]*yynewdata[idx(i, j, Nx, Ny, 1)]-egvymax[idx_v(i,j,Nx)]*yynewdata[idx(i, j, Nx, Ny, 1)]);
            yypdata[idx(i, j, Nx, Ny, 2)]=0.5*(egvy[idx(i,j,Nx,Ny,2)]*yynewdata[idx(i, j, Nx, Ny, 2)]+egvymax[idx_v(i,j,Nx)]*yynewdata[idx(i, j, Nx, Ny, 2)]);
            yyndata[idx(i, j, Nx, Ny, 2)]=0.5*(egvy[idx(i,j,Nx,Ny,2)]*yynewdata[idx(i, j, Nx, Ny, 2)]-egvymax[idx_v(i,j,Nx)]*yynewdata[idx(i, j, Nx, Ny, 2)]);
            yypdata[idx(i, j, Nx, Ny, 3)]=0.5*(egvy[idx(i,j,Nx,Ny,3)]*yynewdata[idx(i, j, Nx, Ny, 3)]+egvymax[idx_v(i,j,Nx)]*yynewdata[idx(i, j, Nx, Ny, 3)]);
            yyndata[idx(i, j, Nx, Ny, 3)]=0.5*(egvy[idx(i,j,Nx,Ny,3)]*yynewdata[idx(i, j, Nx, Ny, 3)]-egvymax[idx_v(i,j,Nx)]*yynewdata[idx(i, j, Nx, Ny, 3)]);
	    // printf("    problem parameters:  i = %li, j=%li, yypdata0 = %g,  yypdata1 = %g, yypdata2 = %g,  yypdata3 = %g\n", i, j, yypdata[idx(i, j, Nx, Ny, 0)], yypdata[idx(i, j, Nx, Ny, 1)], yypdata[idx(i, j, Nx, Ny, 2)], yypdata[idx(i, j, Nx, Ny, 3)]);
        }
    }
    
    /* iterate over domain, computing all equations */
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
            
            /* get derivative on x direction */
            flag = SetISX(IS0_px, IS1_px, IS2_px, IS0_nx, IS1_nx, IS2_nx, yxpdata, yxndata, i, j, Nx, Ny, 2);
            if (flag!=0) printf("error in SetISX function \n");
            
            flag = Setalphawx(IS0_px, IS1_px, IS2_px, IS0_nx, IS1_nx, IS2_nx, alpha_0px, alpha_1px, alpha_2px, alpha_0nx, alpha_1nx, alpha_2nx, w0_px, w1_px, w2_px, w0_nx, w1_nx, w2_nx, Epsilon);
            if (flag!=0) printf("error in Setalphawx function \n");
            
            flag = SetUx(w0_px, w1_px, w2_px, w0_nx, w1_nx, w2_nx, yxpdata, yxndata, u_tpphx, u_tnphx, u_tpnhx, u_tnnhx, i, j, Nx, Ny, 2);
            if (flag!=0) printf("error in SetUx function \n");
            
	    //if (i==0&j==0){
            for(k=0;k<4;k++){
               // yxdata[idx(i, j, Nx, Ny, k)]=-(1.0/dx)*((u_tpphx[k]-u_tpnhx[k])+(u_tnphx[k]-u_tnnhx[k]));
	      yxupdata[idx(i, j, Nx, Ny, k)]=(u_tpphx[k]+u_tnphx[k]);
	      yxdowndata[idx(i, j, Nx, Ny, k)]=(u_tpnhx[k]+u_tnnhx[k]);
	      //printf("yx: i=%li, j=%li, k=%li, u_tpphx[k]=%f, u_tnphx[k]=%f,u_tpnhx[k]=%f,u_tnnhx[k]=%f\n",i,j, k,u_tpphx[k], u_tnphx[k], u_tpnhx[k],u_tnnhx[k]);
            }
	       
	    //	    }
        }
    }
    
    for (i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            
            /* get derivative on y direction */
            flag = SetISY(IS0_py, IS1_py, IS2_py, IS0_ny, IS1_ny, IS2_ny, yypdata, yyndata, i, j, Nx, Ny, 2);
            if (flag!=0) printf("error in SetISY function \n");
            
            flag = Setalphawy(IS0_py, IS1_py, IS2_py, IS0_ny, IS1_ny, IS2_ny, alpha_0py, alpha_1py, alpha_2py, alpha_0ny, alpha_1ny, alpha_2ny, w0_py, w1_py, w2_py, w0_ny, w1_ny, w2_ny, Epsilon);
            if (flag!=0) printf("error in Setalphawy function \n");
            
            flag = SetUy(w0_py, w1_py, w2_py, w0_ny, w1_ny, w2_ny, yypdata, yyndata, u_tpphy, u_tnphy, u_tpnhy, u_tnnhy, i, j, Nx, Ny, 2);
            if (flag!=0) printf("error in SetUy function \n");
            
            for(k=0;k<4;k++){
                //yydata[idx(i, j, Nx, Ny, k)]=-(1.0/dy)*((u_tpphy[k]-u_tpnhy[k])+(u_tnphy[k]-u_tnnhy[k]));
                yyupdata[idx(i, j, Nx, Ny, k)]=u_tpphy[k]+u_tnphy[k];
                yydowndata[idx(i, j, Nx, Ny, k)]=u_tpnhy[k]+u_tnnhy[k];
		//printf("yy: i=%li, j=%li, k=%li, u_tpphy[k]=%f, u_tnphy[k]=%f,u_tpnhy[k]=%f,u_tnnhy[k]=%f\n",i,j, k,u_tpphy[k], u_tnphy[k], u_tpnhy[k],u_tnnhy[k]);
            }
        }
    }
    
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
            flag = Setrhxegm(Ydata, rhxegm, i, j, Nx, Ny, gama, 2);
            if (flag!=0) printf("error in Setrhxegm function \n");
            for (k=0;k<4;k++){
	      //if (i==0&&j==0){
	      //printf(" rhxegm : i = %li, j = %li, rhxegm[k][0] = %f, rhxegm[k][1] = %f, rhxegm[k][2] = %f, rhxegm[k][3] = %f\n", i, j, rhxegm[k][0], rhxegm[k][1],rhxegm[k][2],rhxegm[k][3]);
	      //}
                yxupbackdata[idx(i, j, Nx, Ny, k)] = rhxegm[k][0]*yxupdata[idx(i, j, Nx, Ny, 0)]+rhxegm[k][1]*yxupdata[idx(i, j, Nx, Ny, 1)]+rhxegm[k][2]*yxupdata[idx(i, j, Nx, Ny, 2)]+rhxegm[k][3]*yxupdata[idx(i, j, Nx, Ny, 3)];
                yxdownbackdata[idx(i, j, Nx, Ny, k)] = rhxegm[k][0]*yxdowndata[idx(i, j, Nx, Ny, 0)]+rhxegm[k][1]*yxdowndata[idx(i, j, Nx, Ny, 1)]+rhxegm[k][2]*yxdowndata[idx(i, j, Nx, Ny, 2)]+rhxegm[k][3]*yxdowndata[idx(i, j, Nx, Ny, 3)];
            }
       
	    //printf(" yxupback : i = %li, j = %li, yxupbackdata[idx(i, j, Nx, Ny, 0)]=%f, yxupbackdata[idx(i, j, Nx, Ny, 1)]=%f, yxupbackdata[idx(i, j, Nx, Ny, 2)]=%f, yxupbackdata[idx(i, j, Nx, Ny, 3)]=%f\n",i,j,yxupbackdata[idx(i, j, Nx, Ny, 0)],yxupbackdata[idx(i, j, Nx, Ny, 1)],yxupbackdata[idx(i, j, Nx, Ny, 2)],yxupbackdata[idx(i, j, Nx, Ny, 3)]);
        }
    }
    
    for(i=0; i<Nx; i++){
        for (j=0; j<Ny; j++){
            flag = Setrhyegm(Ydata, rhyegm, i, j, Nx, Ny, gama, 2);
            if (flag!=0) printf("error in Setrhyegm function \n");
            for (k=0;k<4;k++){
	      //  if (i==0&&j==0){
	      //printf(" rhyegm : i = %li, j = %li, rhyegm[k][0] = %f, rhyegm[k][1] = %f, rhyegm[k][2] = %f, rhyegm[k][3] = %f\n", i, j, rhyegm[k][0], rhyegm[k][1],rhyegm[k][2],rhyegm[k][3]);
	      //  }
                yyupbackdata[idx(i, j, Nx, Ny, k)] = rhyegm[k][0]*yyupdata[idx(i, j, Nx, Ny, 0)]+rhyegm[k][1]*yyupdata[idx(i, j, Nx, Ny, 1)]+rhyegm[k][2]*yyupdata[idx(i, j, Nx, Ny, 2)]+rhyegm[k][3]*yyupdata[idx(i, j, Nx, Ny, 3)];
                yydownbackdata[idx(i, j, Nx, Ny, k)] = rhyegm[k][0]*yydowndata[idx(i, j, Nx, Ny, 0)]+rhyegm[k][1]*yydowndata[idx(i, j, Nx, Ny, 1)]+rhyegm[k][2]*yydowndata[idx(i, j, Nx, Ny, 2)]+rhyegm[k][3]*yydowndata[idx(i, j, Nx, Ny, 3)];
            }
	    
	    //if(i==0&&j==0){
	    //	printf(" yydownbackdata : i = %li, j = %li, yydownbackdata[idx(i, j, Nx, Ny, 0)] = %f, yydownbackdata[idx(i, j, Nx, Ny, 1)] = %f, yydownbackdata[idx(i, j, Nx, Ny, 2)] = %f, yydownbackdata[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, yydownbackdata[idx(i, j, Nx, Ny, 0)], yydownbackdata[idx(i, j, Nx, Ny, 1)],yydownbackdata[idx(i, j, Nx, Ny, 2)],yydownbackdata[idx(i, j, Nx, Ny, 3)]);
	    //	}
        }
    }
    
    flag = Gettao(yxupbackdata, yxdownbackdata, yyupbackdata, yydownbackdata, taoup, taodown, Nx, Ny, gama);
    if (flag!=0) printf("error in Gettao function \n");
    flag = GetCj(yxupbackdata, yxdownbackdata, yyupbackdata, yydownbackdata, Cjup, Cjdown, Nx, Ny, gama);
    if (flag!=0) printf("error in GetCj function \n");
    
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
            /* get derivative both on x and y direction */
            //for (k=0; k<4; k++){
                //  if(k==2)
                //        dYdata[idx(i, j, Nx, Ny, k)]=yxdata[idx(i, j, Nx, Ny, k)]+yydata[idx(i, j, Nx, Ny, k)]-0.1;
                // else
	      //if (i==0&&j==0){
	      //printf(" dYdata : i = %li, j = %li, yxupbackdata[idx(i, j, Nx, Ny, 1)] = %f, yxdownbackdata[idx(i, j, Nx, Ny, 1)] = %f, yyupbackdata[idx(i, j, Nx, Ny, 2)] = %f, yydownbackdata[idx(i, j, Nx, Ny, 2)] = %f\n", i, j, yxupbackdata[idx(i, j, Nx, Ny, 1)], yxdownbackdata[idx(i, j, Nx, Ny, 1)],yyupbackdata[idx(i, j, Nx, Ny, 2)],yydownbackdata[idx(i, j, Nx, Ny, 2)]);
	      //printf(" taodata1 : i = %li, j = %li, taoupdata[idx(i, j, Nx, Ny, 0)] = %f, taodowndata[idx(i, j, Nx, Ny, 0)] = %f, taoupdata[idx(i, j, Nx, Ny, 1)] = %f, taodowndata[idx(i, j, Nx, Ny, 1)] = %f\n", i, j, taoupdata[idx(i, j, Nx, Ny, 0)], taodowndata[idx(i, j, Nx, Ny, 0)], taoupdata[idx(i, j, Nx, Ny, 1)], taodowndata[idx(i, j, Nx, Ny, 1)]);
	      //}
                dYdata[idx(i, j, Nx, Ny, 0)]=(-1.0/dx)*(yxupbackdata[idx(i, j, Nx, Ny, 1)]-yxdownbackdata[idx(i, j, Nx, Ny, 1)])+(-1.0/dy)*(yyupbackdata[idx(i, j, Nx, Ny, 2)]-yydownbackdata[idx(i, j, Nx, Ny, 2)]);
                dYdata[idx(i, j, Nx, Ny, 1)]=(-1.0/dx)*(taoupdata[idx(i, j, Nx, Ny, 0)]-taodowndata[idx(i, j, Nx, Ny, 0)])+(-1.0/dy)*(taoupdata[idx(i, j, Nx, Ny, 1)]-taodowndata[idx(i, j, Nx, Ny, 1)]);
                dYdata[idx(i, j, Nx, Ny, 2)]=(-1.0/dx)*(taoupdata[idx(i, j, Nx, Ny, 2)]-taodowndata[idx(i, j, Nx, Ny, 2)])+(-1.0/dy)*(taoupdata[idx(i, j, Nx, Ny, 3)]-taodowndata[idx(i, j, Nx, Ny, 3)]);
                dYdata[idx(i, j, Nx, Ny, 4)]=(-1.0/dx)*(Cjupdata[idx(i, j, Nx, Ny, 0)]-Cjdowndata[idx(i, j, Nx, Ny, 0)])+(-1.0/dy)*(Cjupdata[idx(i, j, Nx, Ny, 1)]-Cjdowndata[idx(i, j, Nx, Ny, 1)]);
		//if (i==0&&j==0){
		//printf(" dYdata : i = %li, j = %li, dYdata[idx(i, j, Nx, Ny, 0)] = %f, dYdata[idx(i, j, Nx, Ny, 1)] = %f, dYdata[idx(i, j, Nx, Ny, 2)] = %f, dYdata[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, dYdata[idx(i, j, Nx, Ny, 0)], dYdata[idx(i, j, Nx, Ny, 1)],dYdata[idx(i, j, Nx, Ny, 2)],dYdata[idx(i, j, Nx, Ny, 3)]);
		//}
		//  }
		//printf(" dYdata : i = %li, j = %li, dYdata[idx(i, j, Nx, Ny, 0)] = %f, dYdata[idx(i, j, Nx, Ny, 1)] = %f, dYdata[idx(i, j, Nx, Ny, 2)] = %f, dYdata[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, dYdata[idx(i, j, Nx, Ny, 0)], dYdata[idx(i, j, Nx, Ny, 1)],dYdata[idx(i, j, Nx, Ny, 2)],dYdata[idx(i, j, Nx, Ny, 3)]);
        }
    }
    
    /* delete relative arrays */
    delete []IS0_px;
    delete []IS1_px;
    delete []IS2_px;
    delete []IS0_nx;
    delete []IS1_nx;
    delete []IS2_nx;
    delete []IS0_py;
    delete []IS1_py;
    delete []IS2_py;
    delete []IS0_ny;
    delete []IS1_ny;
    delete []IS2_ny;
    
    delete []alpha_0px;
    delete []alpha_1px;
    delete []alpha_2px;
    delete []alpha_0nx;
    delete []alpha_1nx;
    delete []alpha_2nx;
    delete []alpha_0py;
    delete []alpha_1py;
    delete []alpha_2py;
    delete []alpha_0ny;
    delete []alpha_1ny;
    delete []alpha_2ny;
    
    delete []w0_px;
    delete []w1_px;
    delete []w2_px;
    delete []w0_nx;
    delete []w1_nx;
    delete []w2_nx;
    delete []w0_py;
    delete []w1_py;
    delete []w2_py;
    delete []w0_ny;
    delete []w1_ny;
    delete []w2_ny;
    
    delete []u_tpphx;
    delete []u_tpnhx;
    delete []u_tnphx;
    delete []u_tnnhx;
    delete []u_tpphy;
    delete []u_tpnhy;
    delete []u_tnphy;
    delete []u_tnnhy;
    
    delete []yxupdata;
    delete []yyupdata;
    delete []yxdowndata;
    delete []yydowndata;
    delete []yxpdata;
    delete []yxndata;
    delete []yypdata;
    delete []yyndata;
    delete []yxnewdata;
    delete []yynewdata;
    delete []yxupbackdata;
    delete []yyupbackdata;
    delete []yxdownbackdata;
    delete []yydownbackdata;

    for (i=0;i<4;i++){
        delete[] lfxegm[i];
    }
    delete []lfxegm;
    for (i=0;i<4;i++){
        delete[] rhxegm[i];
    }
    delete []rhxegm;
    for (i=0;i<4;i++){
        delete[] lfyegm[i];
    }
    delete []lfyegm;
    for (i=0;i<4;i++){
        delete[] rhyegm[i];
    }
    delete []rhyegm;

    delete []egvx;
    delete []egvy;
    delete []egvxmax;
    delete []egvymax;

    /* Free vectors */
    //N_VDestroy_Serial(yxp);
    //N_VDestroy_Serial(yxn);
    //N_VDestroy_Serial(yyp);
    //N_VDestroy_Serial(yyn);
    N_VDestroy_Serial(taoup);
    N_VDestroy_Serial(Cjup);
    //N_VDestroy_Serial(yxup);
    //N_VDestroy_Serial(yyup);
    N_VDestroy_Serial(taodown);
    N_VDestroy_Serial(Cjdown);
    //N_VDestroy_Serial(yxdown);
    //N_VDestroy_Serial(yydown);
    //N_VDestroy_Serial(yxnew);
    //N_VDestroy_Serial(yynew);
    //N_VDestroy_Serial(yxupback);
    //N_VDestroy_Serial(yyupback);
    //N_VDestroy_Serial(yxdownback);
    //N_VDestroy_Serial(yydownback);
    
    return 0;
}

/* fill in left eigenvectors for x component*/
//static int Setlfxegm(realtype *Ydata, realtype **lfxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype egv1x, realtype egv2x, realtype egv3x, realtype egv4x, int flag)
static int Setlfxegm(realtype *Ydata, realtype **lfxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype *egvx, int flag)
{
    realtype rou, vx, vy, p, a, vxnext, vynext, pnext, vxcur, vycur, pcur, h, hcur, hnext;
    if (i!=Nx-1){
      vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i+1, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i+1, j, Nx, Ny, 0)]);
      vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i+1, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i+1, j, Nx, Ny, 0)]);
      p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i+1, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i+1, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
      h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i+1, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i+1, j, Nx, Ny, 0)]));
      //vxnext = Ydata[idx(i+1, j, Nx, Ny, 1)]/Ydata[idx(i+1, j, Nx, Ny, 0)];
      //vynext = Ydata[idx(i+1, j, Nx, Ny, 2)]/Ydata[idx(i+1, j, Nx, Ny, 0)];
      //pnext = (gama-1.0)*Ydata[idx(i+1, j, Nx, Ny, 3)]-0.5*Ydata[idx(i+1, j, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
      //hnext = (Ydata[idx(i+1, j, Nx, Ny, 3)]+pnext)/Ydata[idx(i+1, j, Nx, Ny, 0)];
      rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i+1, j, Nx, Ny, 0)]);
        }
    /* periodic boundary condition */
    if (flag == 0)
    {
        if (i==Nx-1){
      vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(0, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)]);
      vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(0, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)]);
      p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(0, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
      h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(0, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)]));
      //vxnext = Ydata[idx(0, j, Nx, Ny, 1)]/Ydata[idx(0, j, Nx, Ny, 0)];
      //    vynext = Ydata[idx(0, j, Nx, Ny, 2)]/Ydata[idx(0, j, Nx, Ny, 0)];
      //    pnext = (gama-1.0)*Ydata[idx(0, j, Nx, Ny, 3)]-0.5*Ydata[idx(0, j, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
      //    hnext = (Ydata[idx(0, j, Nx, Ny, 3)]+pnext)/Ydata[idx(0, j, Nx, Ny, 0)];
          rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)]);
        }
    }
    /* reflecting boundary condition */
    if (flag == 1)
    {
        if (i==Nx-1){
      vx = -(Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
      vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
      p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
      h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]));
      //     vxnext = -Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
      //    vynext = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
      //    pnext = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
      //    hnext = (Ydata[idx(i, j, Nx, Ny, 3)]+pnext)/Ydata[idx(i, j, Nx, Ny, 0)];
          rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    /* nature or transmissive boundary condition */
    if (flag == 2)
    {
        if (i==Nx-1){
      vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
      vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
      p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
      h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]));
      //       vxnext = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
      //    vynext = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
      //    pnext = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
      //    hnext = (Ydata[idx(i, j, Nx, Ny, 3)]+pnext)/Ydata[idx(i, j, Nx, Ny, 0)];
          rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    
    //  vxcur = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
    //vycur = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
    //pcur = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vxcur*vxcur+vycur*vycur);
    //hcur = (Ydata[idx(i, j, Nx, Ny, 3)]+pcur)/Ydata[idx(i, j, Nx, Ny, 0)];
    
    //vx = 0.5*(vxcur+vxnext);
    //vy = 0.5*(vycur+vynext);
    //p = 0.5*(pcur+pnext);
    a = sqrt(gama*p/rou);
    //h = 0.5*(hcur+hnext);
    // egv1x = vx-a;
    //egv2x = vx;
    //egv3x = vx+a;
    //egv4x = vx;
    egvx[idx(i,j,Nx,Ny,0)]=vx-a;
    egvx[idx(i,j,Nx,Ny,1)]=vx;
    egvx[idx(i,j,Nx,Ny,2)]=vx+a;
    egvx[idx(i,j,Nx,Ny,3)]=vx;
    
    lfxegm[0][0] = ((gama-1.0)*h)/(2.0*a*a)+vx/(2.0*a)-0.5;
    lfxegm[1][0] = ((1.0-gama)*h)/(a*a)+2.0;
    lfxegm[2][0] = ((gama-1.0)*h)/(2.0*a*a)-vx/(2.0*a)-0.5;
    lfxegm[3][0] = ((1.0-gama)*h*vy)/(a*a)+vy;
    lfxegm[0][1] = ((1.0-gama)*vx)/(2*a*a)-1.0/(2.0*a);
    lfxegm[1][1] = (gama-1.0)*vx/(a*a);
    lfxegm[2][1] = ((1.0-gama)*vx)/(2*a*a)+1.0/(2.0*a);
    lfxegm[3][1] = (gama-1.0)*vx*vy/(a*a);
    lfxegm[0][2] = ((1.0-gama)*vy)/(2*a*a);
    lfxegm[1][2] = (gama-1.0)*vy/(a*a);
    lfxegm[2][2] = ((1.0-gama)*vy)/(2*a*a);
    lfxegm[3][2] = (gama-1.0)*vy*vy/(a*a)+1.0;
    lfxegm[0][3] = (gama-1.0)/(2*a*a);
    lfxegm[1][3] = (1.0-gama)/(a*a);
    lfxegm[2][3] = (gama-1.0)/(2*a*a);
    lfxegm[3][3] = (1.0-gama)*vy/(a*a);
    //printf("lfxegm : i=%li, j=%li, lfxegm[0][0]=%f, lfxegm[1][0]=%f, lfxegm[2][0]=%f, lfxegm[3][0]=%f, lfxegm[0][1]=%f, lfxegm[1][1]=%f, lfxegm[2][1]=%f, lfxegm[3][1]=%f, lfxegm[0][2]=%f, lfxegm[1][2]=%f, lfxegm[2][2]=%f, lfxegm[3][2]=%f, lfxegm[0][3]=%f, lfxegm[1][3]=%f, lfxegm[2][3]=%f, lfxegm[3][3]=%f\n", i, j, lfxegm[0][0], lfxegm[1][0], lfxegm[2][0], lfxegm[3][0], lfxegm[0][1], lfxegm[1][1], lfxegm[2][1], lfxegm[3][1], lfxegm[0][2], lfxegm[1][2], lfxegm[2][2], lfxegm[3][2], lfxegm[0][3], lfxegm[1][3], lfxegm[2][3], lfxegm[3][3]);
    
    return 0;
}

/* fill in right eigenvectors for x component*/
static int Setrhxegm(realtype *Ydata, realtype **rhxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag)
{
    realtype rou, vx, vy, p, a, vxnext, vynext, pnext, vxcur, vycur, pcur, h, hcur, hnext;
    if (i!=Nx-1){
      vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i+1, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i+1, j, Nx, Ny, 0)]);
      vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i+1, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i+1, j, Nx, Ny, 0)]);
      p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i+1, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i+1, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
      h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i+1, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i+1, j, Nx, Ny, 0)]));
      //vxnext = Ydata[idx(i+1, j, Nx, Ny, 1)]/Ydata[idx(i+1, j, Nx, Ny, 0)];
      //vynext = Ydata[idx(i+1, j, Nx, Ny, 2)]/Ydata[idx(i+1, j, Nx, Ny, 0)];
      //pnext = (gama-1.0)*Ydata[idx(i+1, j, Nx, Ny, 3)]-0.5*Ydata[idx(i+1, j, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
      //hnext = (Ydata[idx(i+1, j, Nx, Ny, 3)]+pnext)/Ydata[idx(i+1, j, Nx, Ny, 0)];
        rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i+1, j, Nx, Ny, 0)]);
    }
    /* periodic boundary condition */
    if (flag == 0)
    {
        if (i==Nx-1){
      vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(0, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)]);
      vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(0, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)]);
      p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(0, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
      h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(0, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)]));
      //      vxnext = Ydata[idx(0, j, Nx, Ny, 1)]/Ydata[idx(0, j, Nx, Ny, 0)];
      //    vynext = Ydata[idx(0, j, Nx, Ny, 2)]/Ydata[idx(0, j, Nx, Ny, 0)];
      //    pnext = (gama-1.0)*Ydata[idx(0, j, Nx, Ny, 3)]-0.5*Ydata[idx(0, j, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
      //    hnext = (Ydata[idx(0, j, Nx, Ny, 3)]+pnext)/Ydata[idx(0, j, Nx, Ny, 0)];
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)]);
        }
    }
    /* reflecting boundary condition */
    if (flag == 1)
    {
        if (i==Nx-1){
      vx = -(Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
      vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
      p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
      h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]));
      //        vxnext = -Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
      //    vynext = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
      //    pnext = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
      //    hnext = (Ydata[idx(i, j, Nx, Ny, 3)]+pnext)/Ydata[idx(i, j, Nx, Ny, 0)];
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    /* nature or transmissive boundary condition */
    if (flag == 2)
    {
        if (i==Nx-1){
      vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
      vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
      p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
      h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]));
      //     vxnext = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
      //    vynext = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
      //    pnext = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
      //    hnext = (Ydata[idx(i, j, Nx, Ny, 3)]+pnext)/Ydata[idx(i, j, Nx, Ny, 0)];
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    
    //vxcur = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
    //vycur = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
    //pcur = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vxcur*vxcur+vycur*vycur);
    //hcur = (Ydata[idx(i, j, Nx, Ny, 3)]+pcur)/Ydata[idx(i, j, Nx, Ny, 0)];
    
    //vx = 0.5*(vxcur+vxnext);
    //vy = 0.5*(vycur+vynext);
    //p = 0.5*(pcur+pnext);
    a = sqrt(gama*p/rou);
    //h = 0.5*(hcur+hnext);
    
    rhxegm[0][0] = 1.0;
    rhxegm[1][0] = vx-a;
    rhxegm[2][0] = vy;
    rhxegm[3][0] = vx*vx+vy*vy-h-a*vx+(2*a*a)/(gama-1.0);
    rhxegm[0][1] = 1.0;
    rhxegm[1][1] = vx;
    rhxegm[2][1] = 0.0;
    rhxegm[3][1] = vx*vx+a*a/(gama-1.0)-h;
    rhxegm[0][2] = 1.0;
    rhxegm[1][2] = vx+a;
    rhxegm[2][2] = vy;
    rhxegm[3][2] = vx*vx+vy*vy-h+a*vx+(2*a*a)/(gama-1.0);
    rhxegm[0][3] = 0.0;
    rhxegm[1][3] = 0.0;
    rhxegm[2][3] = 1.0;
    rhxegm[3][3] = vy;
    
    return 0;
}

/* fill in left eigenvectors for y component*/
//static int Setlfyegm(realtype *Ydata, realtype **lfyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype egv1y, realtype egv2y, realtype egv3y, realtype egv4y, int flag)
static int Setlfyegm(realtype *Ydata, realtype **lfyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype *egvy, int flag)
{
    realtype rou, vx, vy, p, a, vxnext, vynext, pnext, vxcur, vycur, pcur, h, hcur, hnext;
    if (j!=Ny-1){
        vxnext = Ydata[idx(i, j+1, Nx, Ny, 1)]/Ydata[idx(i, j+1, Nx, Ny, 0)];
        vynext = Ydata[idx(i, j+1, Nx, Ny, 2)]/Ydata[idx(i, j+1, Nx, Ny, 0)];
        pnext = (gama-1.0)*Ydata[idx(i, j+1, Nx, Ny, 3)]-0.5*Ydata[idx(i, j+1, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
        hnext = (Ydata[idx(i, j+1, Nx, Ny, 3)]+pnext)/Ydata[idx(i, j+1, Nx, Ny, 0)];
        rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j+1, Nx, Ny, 0)]);
    }
    /* periodic boundary condition */
    if (flag == 0)
    {
        if (j==Ny-1){
            vxnext = Ydata[idx(i, 0, Nx, Ny, 1)]/Ydata[idx(i, 0, Nx, Ny, 0)];
            vynext = Ydata[idx(i, 0, Nx, Ny, 2)]/Ydata[idx(i, 0, Nx, Ny, 0)];
            pnext = (gama-1.0)*Ydata[idx(i, 0, Nx, Ny, 3)]-0.5*Ydata[idx(i, 0, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
            hnext = (Ydata[idx(i, 0, Nx, Ny, 3)]+pnext)/Ydata[idx(i, 0, Nx, Ny, 0)];
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, 0, Nx, Ny, 0)]);
        }
    }
    /* reflecting boundary condition */
    if (flag == 1)
    {
        if (j==Ny-1){
            vxnext = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
            vynext = -Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
            pnext = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
            hnext = (Ydata[idx(i, j, Nx, Ny, 3)]+pnext)/Ydata[idx(i, j, Nx, Ny, 0)];
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    /* nature or transmissive boundary condition */
    if (flag == 2)
    {
        if (j==Ny-1){
            vxnext = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
            vynext = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
            pnext = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
            hnext = (Ydata[idx(i, j, Nx, Ny, 3)]+pnext)/Ydata[idx(i, j, Nx, Ny, 0)];
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    
    vxcur = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
    vycur = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
    pcur = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vxcur*vxcur+vycur*vycur);
    hcur = (Ydata[idx(i, j, Nx, Ny, 3)]+pcur)/Ydata[idx(i, j, Nx, Ny, 0)];
    
    vx = 0.5*(vxcur+vxnext);
    vy = 0.5*(vycur+vynext);
    p = 0.5*(pcur+pnext);
    a = sqrt(gama*p/rou);
    h = 0.5*(hcur+hnext);
    //egv1y = vy-a;
    //egv2y = vy;
    //egv3y = vy+a;
    //egv4y = vy;
    egvy[idx(i,j,Nx,Ny,0)]=vy-a;
    egvy[idx(i,j,Nx,Ny,1)]=vy;
    egvy[idx(i,j,Nx,Ny,2)]=vy+a;
    egvy[idx(i,j,Nx,Ny,3)]=vy;

    lfyegm[0][0] = ((gama-1.0)*h)/(2.0*a*a)+vy/(2.0*a)-0.5;
    lfyegm[1][0] = ((1.0-gama)*h)/(a*a)+2.0;
    lfyegm[2][0] = ((gama-1.0)*h)/(2.0*a*a)-vy/(2.0*a)-0.5;
    lfyegm[3][0] = ((1.0-gama)*h*vx)/(a*a)+vx;
    lfyegm[0][1] = ((1.0-gama)*vx)/(2*a*a);
    lfyegm[1][1] = (gama-1.0)*vx/(a*a);
    lfyegm[2][1] = ((1.0-gama)*vx)/(2*a*a);
    lfyegm[3][1] = (gama-1.0)*vx*vx/(a*a)+1.0;
    lfyegm[0][2] = ((1.0-gama)*vy)/(2*a*a)-1.0/(2.0*a);
    lfyegm[1][2] = (gama-1.0)*vy/(a*a);
    lfyegm[2][2] = ((1.0-gama)*vy)/(2*a*a)+1.0/(2.0*a);
    lfyegm[3][2] = (gama-1.0)*vx*vy/(a*a);
    lfyegm[0][3] = (gama-1.0)/(2*a*a);
    lfyegm[1][3] = (1.0-gama)/(a*a);
    lfyegm[2][3] = (gama-1.0)/(2*a*a);
    lfyegm[3][3] = (1.0-gama)*vx/(a*a);
    
    return 0;
}

/* fill in right eigenvectors for y component*/
static int Setrhyegm(realtype *Ydata, realtype **rhyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag)
{
    realtype rou, vx, vy, p, a, vxnext, vynext, pnext, vxcur, vycur, pcur, h, hcur, hnext;
    if (j!=Ny-1){
        vxnext = Ydata[idx(i, j+1, Nx, Ny, 1)]/Ydata[idx(i, j+1, Nx, Ny, 0)];
        vynext = Ydata[idx(i, j+1, Nx, Ny, 2)]/Ydata[idx(i, j+1, Nx, Ny, 0)];
        pnext = (gama-1.0)*Ydata[idx(i, j+1, Nx, Ny, 3)]-0.5*Ydata[idx(i, j+1, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
        hnext = (Ydata[idx(i, j+1, Nx, Ny, 3)]+pnext)/Ydata[idx(i, j+1, Nx, Ny, 0)];
        rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j+1, Nx, Ny, 0)]);
    }
    /* periodic boundary condition */
    if (flag == 0)
    {
        if (j==Ny-1){
            vxnext = Ydata[idx(i, 0, Nx, Ny, 1)]/Ydata[idx(i, 0, Nx, Ny, 0)];
            vynext = Ydata[idx(i, 0, Nx, Ny, 2)]/Ydata[idx(i, 0, Nx, Ny, 0)];
            pnext = (gama-1.0)*Ydata[idx(i, 0, Nx, Ny, 3)]-0.5*Ydata[idx(i, 0, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
            hnext = (Ydata[idx(i, 0, Nx, Ny, 3)]+pnext)/Ydata[idx(i, 0, Nx, Ny, 0)];
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, 0, Nx, Ny, 0)]);
        }
    }
    /* reflecting boundary condition */
    if (flag == 1)
    {
        if (j==Ny-1){
            vxnext = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
            vynext = -Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
            pnext = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
            hnext = (Ydata[idx(i, j, Nx, Ny, 3)]+pnext)/Ydata[idx(i, j, Nx, Ny, 0)];
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    /* nature or transmissive boundary condition */
    if (flag == 2)
    {
        if (j==Ny-1){
            vxnext = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
            vynext = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
            pnext = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vxnext*vxnext+vynext*vynext);
            hnext = (Ydata[idx(i, j, Nx, Ny, 3)]+pnext)/Ydata[idx(i, j, Nx, Ny, 0)];
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    
    vxcur = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
    vycur = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
    pcur = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vxcur*vxcur+vycur*vycur);
    hcur = (Ydata[idx(i, j, Nx, Ny, 3)]+pcur)/Ydata[idx(i, j, Nx, Ny, 0)];
    
    vx = 0.5*(vxcur+vxnext);
    vy = 0.5*(vycur+vynext);
    p = 0.5*(pcur+pnext);
    a = sqrt(gama*p/rou);
    h = 0.5*(hcur+hnext);
    
    rhyegm[0][0] = 1.0;
    rhyegm[1][0] = vx;
    rhyegm[2][0] = vy-a;
    rhyegm[3][0] = vx*vx+vy*vy-h-vy*a+(2.0*a*a)/(gama-1.0);
    rhyegm[0][1] = 1.0;
    rhyegm[1][1] = 0.0;
    rhyegm[2][1] = vy;
    rhyegm[3][1] = vy*vy+a*a/(gama-1.0)-h;
    rhyegm[0][2] = 1.0;
    rhyegm[1][2] = vx;
    rhyegm[2][2] = vy+a;
    rhyegm[3][2] = vx*vx+vy*vy-h+vy*a+(2.0*a*a)/(gama-1.0);
    rhyegm[0][3] = 0.0;
    rhyegm[1][3] = 1.0;
    rhyegm[2][3] = 0.0;
    rhyegm[3][3] = vx;
    
    return 0;
}

/* Fill in the indicators of smoothness on x direction */
static int SetISX(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx,realtype *yxpdata, realtype *yxndata, long int i, long int j, long int Nx, long int Ny, int flag)
{
    /* consider Nx */
    if (Nx<4){
        cerr << "\nNx should be more than 3!\n";
        return 1;
    }
    else {
        /* declaration */
        long int k, n;
        realtype *yp = new realtype [28];
        realtype *yn = new realtype [28];
        
        if (i>2&&i<Nx-3){
            for (n=0;n<4;n++){
                for (k=0;k<7;k++){
                    yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                    yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                }
            }
        }
        
        if (flag == 0){
            if (i<3){
                for (n=0;n<4;n++){
                    for (k=3;k<7;k++){
                        yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==0){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(Nx-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-3,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(Nx-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-3,j,Nx,Ny,n)];
                    }
                }
                
                if (i==1){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-2,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-2,j,Nx,Ny,n)];
                    }
                }
                
                if (i==2){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                    }
                }
            }
            
            if(i>Nx-4){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-3){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(0,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(0,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(0,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(0,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(1,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(1,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(0,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(0,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(2,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(2,j,Nx,Ny,n)];
                    }
                }
            }
        }
        
        
        if (flag == 1){
            /* consider the situations of different i */
            if (i<3){
                for (n=0;n<4;n++){
                    for (k=3;k<7;k++){
                        yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                    }
                }
                if (i==0){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                    }
		    //        yp[2+3*7] = -yxpdata[idx(i,j,Nx,Ny,3)];
		    // yp[1+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yp[0+3*7] = -yxpdata[idx(i+2,j,Nx,Ny,3)];
                    //yn[2+3*7] = -yxndata[idx(i,j,Nx,Ny,3)];
                    //yn[1+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yn[0+3*7] = -yxndata[idx(i+2,j,Nx,Ny,3)];
                }
                if (i==1){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yxpdata[idx(i-1,j,Nx,Ny,3)];
                    //yp[1+3*7] = -yxpdata[idx(i-1,j,Nx,Ny,3)];
                    //yp[0+3*7] = -yxpdata[idx(i,j,Nx,Ny,3)];
                    //yn[2+3*7] = -yxndata[idx(i-1,j,Nx,Ny,3)];
                    //yn[1+3*7] = -yxndata[idx(i-1,j,Nx,Ny,3)];
		    // yn[0+3*7] = -yxndata[idx(i,j,Nx,Ny,3)];
                }
                if (i==2){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yxpdata[idx(i-1,j,Nx,Ny,3)];
                    //yp[1+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yp[0+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yn[2+3*7] = -yxndata[idx(i-1,j,Nx,Ny,3)];
                    //yn[1+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                    //yn[0+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                }
            }
            
            if(i>Nx-4){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-3){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yp[5+3*7] = -yxpdata[idx(i+2,j,Nx,Ny,3)];
                    //yp[6+3*7] = -yxpdata[idx(i+2,j,Nx,Ny,3)];
                    //yn[4+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yn[5+3*7] = -yxndata[idx(i+2,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yxndata[idx(i+2,j,Nx,Ny,3)];
                }
                if (i==Nx-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yn[4+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yp[5+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yn[5+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yp[6+3*7] = -yxpdata[idx(i,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yxndata[idx(i,j,Nx,Ny,3)];
                }
                if (i==Nx-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yxpdata[idx(i,j,Nx,Ny,3)];
                    //yn[4+3*7] = -yxndata[idx(i,j,Nx,Ny,3)];
                    //yp[5+3*7] = -yxpdata[idx(i-1,j,Nx,Ny,3)];
                    //yn[5+3*7] = -yxndata[idx(i-1,j,Nx,Ny,3)];
                    //yp[6+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                }
            }
        }
        
        
        if (flag == 2){
            if (i<3){
                for (n=0;n<4;n++){
                    for (k=3;k<7;k++){
                        yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==0){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                    }
                }
                
                if (i==1){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                    }
                }
                
                if (i==2){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                    }
                }
            }
            
            if(i>Nx-4){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-3){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                    }
                }
            }
        }
        
        
        /* compute indicators of smoothness of rou, qx, qy, E */
        for (n=0;n<4;n++){
            IS0_px[n]=(13.0/12.0)*(yp[1+n*7]-2.0*yp[2+n*7]+yp[3+n*7])*(yp[1+n*7]-2.0*yp[2+n*7]+yp[3+n*7])+(1.0/4.0)*(yp[1+n*7]-4.0*yp[2+n*7]+3.0*yp[3+n*7])*(yp[1+n*7]-4.0*yp[2+n*7]+3.0*yp[3+n*7]);
            
            IS1_px[n]=(13.0/12.0)*(yp[2+n*7]-2.0*yp[3+n*7]+yp[4+n*7])*(yp[2+n*7]-2.0*yp[3+n*7]+yp[4+n*7])+(1.0/4.0)*(yp[2+n*7]-yp[4+n*7])*(yp[2+n*7]-yp[4+n*7]);
            
            IS2_px[n]=(13.0/12.0)*(yp[3+n*7]-2.0*yp[4+n*7]+yp[5+n*7])*(yp[3+n*7]-2.0*yp[4+n*7]+yp[5+n*7])+(1.0/4.0)*(3.0*yp[3+n*7]-4.0*yp[4+n*7]+yp[5+n*7])*(3.0*yp[3+n*7]-4.0*yp[4+n*7]+yp[5+n*7]);
            
            IS0_nx[n]=(13.0/12.0)*(yn[4+n*7]-2.0*yn[5+n*7]+yn[6+n*7])*(yn[4+n*7]-2.0*yn[5+n*7]+yn[6+n*7])+(1.0/4.0)*(3.0*yn[4+n*7]-4.0*yn[5+n*7]+yn[6+n*7])*(3.0*yn[4+n*7]-4.0*yn[5+n*7]+yn[6+n*7]);
            
            IS1_nx[n]=(13.0/12.0)*(yn[3+n*7]-2.0*yn[4+n*7]+yn[5+n*7])*(yn[3+n*7]-2.0*yn[4+n*7]+yn[5+n*7])+(1.0/4.0)*(yn[3+n*7]-yn[5+n*7])*(yn[3+n*7]-yn[5+n*7]);
            
            IS2_nx[n]=(13.0/12.0)*(yn[2+n*7]-2.0*yn[3+n*7]+yn[4+n*7])*(yn[2+n*7]-2.0*yn[3+n*7]+yn[4+n*7])+(1.0/4.0)*(yn[2+n*7]-4.0*yn[3+n*7]+3.0*yn[4+n*7])*(yn[2+n*7]-4.0*yn[3+n*7]+3.0*yn[4+n*7]);
        }
        
        delete []yp;
        delete []yn;
        
        return 0;
    }
}

/* Fill in the indicators of smoothness on y direction */
static int SetISY(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny,realtype *yypdata, realtype *yyndata, long int i, long int j, long int Nx, long int Ny, int flag)
{
    /* consider Ny */
    if (Ny<4){
        cerr << "\nNy should be more than 3!\n";
        return 1;
    }
    else {
        /* declaration */
        long int k, n;
        realtype *yp = new realtype [28];
        realtype *yn = new realtype [28];
        
        if (j>2&&j<Ny-3){
            for (n=0;n<4;n++){
                for (k=0;k<7;k++){
                    yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                    yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                }
            }
        }
        
        if (flag == 0){
            if (j<3){
                for (n=0;n<4;n++){
                    for (k=3;k<7;k++){
                        yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                    }
                }
                
                if (j==0){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,Ny-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-3,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,Ny-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-3,Nx,Ny,n)];
                    }
                }
                
                if (j==1){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-2,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-2,Nx,Ny,n)];
                    }
                }
                
                if (j==2){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                    }
                }
            }
            
            if(j>Ny-4){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-3){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,0,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,0,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,0,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,0,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,1,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,1,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,0,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,0,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,2,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,2,Nx,Ny,n)];
                    }
                }
            }
        }
        
        
        if (flag == 1){
            /* consider the situations of different i */
            if (j<3){
                for (n=0;n<4;n++){
                    for (k=3;k<7;k++){
                        yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                    }
                }
                if (j==0){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yypdata[idx(i,j,Nx,Ny,3)];
                    //yp[1+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yp[0+3*7] = -yypdata[idx(i,j+2,Nx,Ny,3)];
                    //yn[2+3*7] = -yyndata[idx(i,j,Nx,Ny,3)];
                    //yn[1+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yn[0+3*7] = -yyndata[idx(i,j+2,Nx,Ny,3)];
                }
                if (j==1){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yypdata[idx(i,j-1,Nx,Ny,3)];
                    //yp[1+3*7] = -yypdata[idx(i,j-1,Nx,Ny,3)];
                    //yp[0+3*7] = -yypdata[idx(i,j,Nx,Ny,3)];
                    //yn[2+3*7] = -yyndata[idx(i,j-1,Nx,Ny,3)];
                    //yn[1+3*7] = -yyndata[idx(i,j-1,Nx,Ny,3)];
                    //yn[0+3*7] = -yyndata[idx(i,j,Nx,Ny,3)];
                }
                if (j==2){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yypdata[idx(i,j-1,Nx,Ny,3)];
                    //yp[1+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yp[0+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yn[2+3*7] = -yyndata[idx(i,j-1,Nx,Ny,3)];
                    //yn[1+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                    //yn[0+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                }
            }
            
            if(j>Ny-4){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-3){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yp[5+3*7] = -yypdata[idx(i,j+2,Nx,Ny,3)];
                    //yp[6+3*7] = -yypdata[idx(i,j+2,Nx,Ny,3)];
                    //yn[4+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yn[5+3*7] = -yyndata[idx(i,j+2,Nx,Ny,3)];
                    //yn[6+3*7] = -yyndata[idx(i,j+2,Nx,Ny,3)];
                }
                if (j==Ny-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yn[4+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yp[5+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yn[5+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yp[6+3*7] = -yypdata[idx(i,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yyndata[idx(i,j,Nx,Ny,3)];
                }
                if (j==Ny-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yypdata[idx(i,j,Nx,Ny,3)];
                    //yn[4+3*7] = -yyndata[idx(i,j,Nx,Ny,3)];
                    //yp[5+3*7] = -yypdata[idx(i,j-1,Nx,Ny,3)];
                    //yn[5+3*7] = -yyndata[idx(i,j-1,Nx,Ny,3)];
                    //yp[6+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yn[6+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                }
            }
        }
        
        
        if (flag == 2){
            if (j<3){
                for (n=0;n<4;n++){
                    for (k=3;k<7;k++){
                        yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                    }
                }
                
                if (j==0){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                    }
                }
                
                if (j==1){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                    }
                }
                
                if (j==2){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                    }
                }
            }
            
            if(j>Ny-4){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-3){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                    }
                }
            }
        }
        
        /* compute indicators of smoothness of rou, qx, qy, E */
        for (n=0;n<4;n++){
            IS0_py[n]=(13.0/12.0)*(yp[1+n*7]-2.0*yp[2+n*7]+yp[3+n*7])*(yp[1+n*7]-2.0*yp[2+n*7]+yp[3+n*7])+(1.0/4.0)*(yp[1+n*7]-4.0*yp[2+n*7]+3.0*yp[3+n*7])*(yp[1+n*7]-4.0*yp[2+n*7]+3.0*yp[3+n*7]);
            
            IS1_py[n]=(13.0/12.0)*(yp[2+n*7]-2.0*yp[3+n*7]+yp[4+n*7])*(yp[2+n*7]-2.0*yp[3+n*7]+yp[4+n*7])+(1.0/4.0)*(yp[2+n*7]-yp[4+n*7])*(yp[2+n*7]-yp[4+n*7]);
            
            IS2_py[n]=(13.0/12.0)*(yp[3+n*7]-2.0*yp[4+n*7]+yp[5+n*7])*(yp[3+n*7]-2.0*yp[4+n*7]+yp[5+n*7])+(1.0/4.0)*(3.0*yp[3+n*7]-4.0*yp[4+n*7]+yp[5+n*7])*(3.0*yp[3+n*7]-4.0*yp[4+n*7]+yp[5+n*7]);
            
            IS0_ny[n]=(13.0/12.0)*(yn[4+n*7]-2.0*yn[5+n*7]+yn[6+n*7])*(yn[4+n*7]-2.0*yn[5+n*7]+yn[6+n*7])+(1.0/4.0)*(3.0*yn[4+n*7]-4.0*yn[5+n*7]+yn[6+n*7])*(3.0*yn[4+n*7]-4.0*yn[5+n*7]+yn[6+n*7]);
            
            IS1_ny[n]=(13.0/12.0)*(yn[3+n*7]-2.0*yn[4+n*7]+yn[5+n*7])*(yn[3+n*7]-2.0*yn[4+n*7]+yn[5+n*7])+(1.0/4.0)*(yn[3+n*7]-yn[5+n*7])*(yn[3+n*7]-yn[5+n*7]);
            
            IS2_ny[n]=(13.0/12.0)*(yn[2+n*7]-2.0*yn[3+n*7]+yn[4+n*7])*(yn[2+n*7]-2.0*yn[3+n*7]+yn[4+n*7])+(1.0/4.0)*(yn[2+n*7]-4.0*yn[3+n*7]+3.0*yn[4+n*7])*(yn[2+n*7]-4.0*yn[3+n*7]+3.0*yn[4+n*7]);
        }
        
        delete []yp;
        delete []yn;
        
        return 0;
    }
}

/* Fill in the stencil weights on x direction */
static int Setalphawx(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx, realtype *alpha_0px, realtype *alpha_1px, realtype *alpha_2px, realtype *alpha_0nx, realtype *alpha_1nx, realtype *alpha_2nx, realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype Epsilon)
{
    /* compute the weights for rou, qx, qy, E */
    long int k;
    for(k=0;k<4;k++)
    {
        alpha_0px[k]=(1.0/10.0)*(1.0/(Epsilon+IS0_px[k]))*(1.0/(Epsilon+IS0_px[k]));
        alpha_1px[k]=(6.0/10.0)*(1.0/(Epsilon+IS1_px[k]))*(1.0/(Epsilon+IS1_px[k]));
        alpha_2px[k]=(3.0/10.0)*(1.0/(Epsilon+IS2_px[k]))*(1.0/(Epsilon+IS2_px[k]));
        alpha_0nx[k]=(1.0/10.0)*(1.0/(Epsilon+IS0_nx[k]))*(1.0/(Epsilon+IS0_nx[k]));
        alpha_1nx[k]=(6.0/10.0)*(1.0/(Epsilon+IS1_nx[k]))*(1.0/(Epsilon+IS1_nx[k]));
        alpha_2nx[k]=(3.0/10.0)*(1.0/(Epsilon+IS2_nx[k]))*(1.0/(Epsilon+IS2_nx[k]));
        
        w0_px[k]=alpha_0px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
        w1_px[k]=alpha_1px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
        w2_px[k]=alpha_2px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
        w0_nx[k]=alpha_0nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
        w1_nx[k]=alpha_1nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
        w2_nx[k]=alpha_2nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
    }
    return 0;
}

/* Fill in the stencil weights on y direction */
static int Setalphawy(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny, realtype *alpha_0py, realtype *alpha_1py, realtype *alpha_2py, realtype *alpha_0ny, realtype *alpha_1ny, realtype *alpha_2ny, realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype Epsilon)
{
    /* compute the weights for rou, qx, qy, E */
    long int k;
    for(k=0;k<4;k++)
    {
        alpha_0py[k]=(1.0/10.0)*(1.0/(Epsilon+IS0_py[k]))*(1.0/(Epsilon+IS0_py[k]));
        alpha_1py[k]=(6.0/10.0)*(1.0/(Epsilon+IS1_py[k]))*(1.0/(Epsilon+IS1_py[k]));
        alpha_2py[k]=(3.0/10.0)*(1.0/(Epsilon+IS2_py[k]))*(1.0/(Epsilon+IS2_py[k]));
        alpha_0ny[k]=(1.0/10.0)*(1.0/(Epsilon+IS0_ny[k]))*(1.0/(Epsilon+IS0_ny[k]));
        alpha_1ny[k]=(6.0/10.0)*(1.0/(Epsilon+IS1_ny[k]))*(1.0/(Epsilon+IS1_ny[k]));
        alpha_2ny[k]=(3.0/10.0)*(1.0/(Epsilon+IS2_ny[k]))*(1.0/(Epsilon+IS2_ny[k]));
        
        w0_py[k]=alpha_0py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
        w1_py[k]=alpha_1py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
        w2_py[k]=alpha_2py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
        w0_ny[k]=alpha_0ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
        w1_ny[k]=alpha_1ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
        w2_ny[k]=alpha_2ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
    }
    return 0;
}

/* Get the derivative on x direction */
static int SetUx(realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype *yxpdata, realtype *yxndata, realtype *u_tpphx, realtype *u_tnphx, realtype *u_tpnhx, realtype *u_tnnhx, long int i, long int j, long int Nx, long int Ny, int flag)
{
    /* consider Nx */
    if (Nx<4){
        cerr << "\nNx should be more than 3!\n";
        return 1;
    }
    else {
        /* declaration */
        long int k, n;
        realtype *yp = new realtype [28];
        realtype *yn = new realtype [28];
        
        if (i>2&&i<Nx-3){
            for (n=0;n<4;n++){
                for (k=0;k<7;k++){
                    yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                    yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                }
            }
        }
        
        if (flag == 0){
            if (i<3){
                for (n=0;n<4;n++){
                    for (k=3;k<7;k++){
                        yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==0){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(Nx-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-3,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(Nx-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-3,j,Nx,Ny,n)];
                    }
                }
                
                if (i==1){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-2,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-2,j,Nx,Ny,n)];
                    }
                }
                
                if (i==2){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                    }
                }
            }
            
            if(i>Nx-4){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-3){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(0,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(0,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(0,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(0,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(1,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(1,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(0,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(0,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(2,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(2,j,Nx,Ny,n)];
                    }
                }
            }
        }
        
        
        if (flag == 1){
            /* consider the situations of different i */
            if (i<3){
                for (n=0;n<4;n++){
                    for (k=3;k<7;k++){
                        yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                    }
                }
                if (i==0){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yxpdata[idx(i,j,Nx,Ny,3)];
                    //yp[1+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yp[0+3*7] = -yxpdata[idx(i+2,j,Nx,Ny,3)];
                    //yn[2+3*7] = -yxndata[idx(i,j,Nx,Ny,3)];
                    //yn[1+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yn[0+3*7] = -yxndata[idx(i+2,j,Nx,Ny,3)];
                }
                if (i==1){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yxpdata[idx(i-1,j,Nx,Ny,3)];
                    //yp[1+3*7] = -yxpdata[idx(i-1,j,Nx,Ny,3)];
                    //yp[0+3*7] = -yxpdata[idx(i,j,Nx,Ny,3)];
                    //yn[2+3*7] = -yxndata[idx(i-1,j,Nx,Ny,3)];
                    //yn[1+3*7] = -yxndata[idx(i-1,j,Nx,Ny,3)];
                    //yn[0+3*7] = -yxndata[idx(i,j,Nx,Ny,3)];
                }
                if (i==2){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yxpdata[idx(i-1,j,Nx,Ny,3)];
                    //yp[1+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yp[0+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yn[2+3*7] = -yxndata[idx(i-1,j,Nx,Ny,3)];
                    //yn[1+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                    //yn[0+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                }
            }
            
            if(i>Nx-4){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-3){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yp[5+3*7] = -yxpdata[idx(i+2,j,Nx,Ny,3)];
                    //yp[6+3*7] = -yxpdata[idx(i+2,j,Nx,Ny,3)];
                    //yn[4+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yn[5+3*7] = -yxndata[idx(i+2,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yxndata[idx(i+2,j,Nx,Ny,3)];
                }
                if (i==Nx-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yn[4+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yp[5+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yn[5+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yp[6+3*7] = -yxpdata[idx(i,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yxndata[idx(i,j,Nx,Ny,3)];
                }
                if (i==Nx-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yxpdata[idx(i,j,Nx,Ny,3)];
                    //yn[4+3*7] = -yxndata[idx(i,j,Nx,Ny,3)];
                    //yp[5+3*7] = -yxpdata[idx(i-1,j,Nx,Ny,3)];
                    //yn[5+3*7] = -yxndata[idx(i-1,j,Nx,Ny,3)];
                    //yp[6+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                }
            }
        }
        
        
        if (flag == 2){
            if (i<3){
                for (n=0;n<4;n++){
                    for (k=3;k<7;k++){
                        yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==0){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                    }
                }
                
                if (i==1){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                    }
                }
                
                if (i==2){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                    }
                }
            }
            
            if(i>Nx-4){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-3){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                    }
                }
            }
        }
        
        /* compute positive and negative solutions on the interface */
        for (n=0;n<4;n++){
            u_tpphx[n]=w0_px[n]*((2.0/6.0)*yp[1+n*7]-(7.0/6.0)*yp[2+n*7]+(11.0/6.0)*yp[3+n*7])+w1_px[n]*((-1.0/6.0)*yp[2+n*7]+(5.0/6.0)*yp[3+n*7]+(2.0/6.0)*yp[4+n*7])+w2_px[n]*((2.0/6.0)*yp[3+n*7]+(5.0/6.0)*yp[4+n*7]-(1.0/6.0)*yp[5+n*7]);
            
            u_tnphx[n]=w2_nx[n]*((-1.0/6.0)*yn[2+n*7]+(5.0/6.0)*yn[3+n*7]+(2.0/6.0)*yn[4+n*7])+w1_nx[n]*((2.0/6.0)*yn[3+n*7]+(5.0/6.0)*yn[4+n*7]-(1.0/6.0)*yn[5+n*7])+w0_nx[n]*((11.0/6.0)*yn[4+n*7]-(7.0/6.0)*yn[5+n*7]+(2.0/6.0)*yn[6+n*7]);
            
            u_tpnhx[n]=w0_px[n]*((2.0/6.0)*yp[0+n*7]-(7.0/6.0)*yp[1+n*7]+(11.0/6.0)*yp[2+n*7])+w1_px[n]*((-1.0/6.0)*yp[1+n*7]+(5.0/6.0)*yp[2+n*7]+(2.0/6.0)*yp[3+n*7])+w2_px[n]*((2.0/6.0)*yp[2+n*7]+(5.0/6.0)*yp[3+n*7]-(1.0/6.0)*yp[4+n*7]);
            
            u_tnnhx[n]=w2_nx[n]*((-1.0/6.0)*yn[1+n*7]+(5.0/6.0)*yn[2+n*7]+(2.0/6.0)*yn[3+n*7])+w1_nx[n]*((2.0/6.0)*yn[2+n*7]+(5.0/6.0)*yn[3+n*7]-(1.0/6.0)*yn[4+n*7])+w0_nx[n]*((11.0/6.0)*yn[3+n*7]-(7.0/6.0)*yn[4+n*7]+(2.0/6.0)*yn[5+n*7]);
        }
        
        delete []yp;
        delete []yn;
        
        return 0;
    }
}

/* Get the derivative on y direction */
static int SetUy(realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype *yypdata, realtype *yyndata, realtype *u_tpphy, realtype *u_tnphy, realtype *u_tpnhy, realtype *u_tnnhy, long int i, long int j, long int Nx, long int Ny, int flag)
{
    /* consider Ny */
    if (Ny<4){
        cerr << "\nNy should be more than 3!\n";
        return 1;
    }
    else {
        /* declaration */
        long int k, n;
        realtype *yp = new realtype [28];
        realtype *yn = new realtype [28];
        
        if (j>2&&j<Ny-3){
            for (n=0;n<4;n++){
                for (k=0;k<7;k++){
                    yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                    yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                }
            }
        }
        
        if (flag == 0){
            if (j<3){
                for (n=0;n<4;n++){
                    for (k=3;k<7;k++){
                        yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                    }
                }
                
                if (j==0){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,Ny-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-3,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,Ny-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-3,Nx,Ny,n)];
                    }
                }
                
                if (j==1){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-2,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-2,Nx,Ny,n)];
                    }
                }
                
                if (j==2){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                    }
                }
            }
            
            if(j>Ny-4){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-3){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,0,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,0,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,0,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,0,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,1,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,1,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,0,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,0,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,2,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,2,Nx,Ny,n)];
                    }
                }
            }
        }
        
        
        if (flag == 1){
            /* consider the situations of different i */
            if (j<3){
                for (n=0;n<4;n++){
                    for (k=3;k<7;k++){
                        yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                    }
                }
                if (j==0){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yypdata[idx(i,j,Nx,Ny,3)];
                    //yp[1+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yp[0+3*7] = -yypdata[idx(i,j+2,Nx,Ny,3)];
                    //yn[2+3*7] = -yyndata[idx(i,j,Nx,Ny,3)];
		    // yn[1+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yn[0+3*7] = -yyndata[idx(i,j+2,Nx,Ny,3)];
                }
                if (j==1){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yypdata[idx(i,j-1,Nx,Ny,3)];
                    //yp[1+3*7] = -yypdata[idx(i,j-1,Nx,Ny,3)];
                    //yp[0+3*7] = -yypdata[idx(i,j,Nx,Ny,3)];
                    //yn[2+3*7] = -yyndata[idx(i,j-1,Nx,Ny,3)];
                    //yn[1+3*7] = -yyndata[idx(i,j-1,Nx,Ny,3)];
                    //yn[0+3*7] = -yyndata[idx(i,j,Nx,Ny,3)];
                }
                if (j==2){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yypdata[idx(i,j-1,Nx,Ny,3)];
                    //yp[1+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yp[0+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yn[2+3*7] = -yyndata[idx(i,j-1,Nx,Ny,3)];
                    //yn[1+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                    //yn[0+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                }
            }
            
            if(j>Ny-4){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-3){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yp[5+3*7] = -yypdata[idx(i,j+2,Nx,Ny,3)];
                    //yp[6+3*7] = -yypdata[idx(i,j+2,Nx,Ny,3)];
                    //yn[4+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yn[5+3*7] = -yyndata[idx(i,j+2,Nx,Ny,3)];
                    //yn[6+3*7] = -yyndata[idx(i,j+2,Nx,Ny,3)];
                }
                if (j==Ny-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yn[4+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yp[5+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yn[5+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yp[6+3*7] = -yypdata[idx(i,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yyndata[idx(i,j,Nx,Ny,3)];
                }
                if (j==Ny-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yypdata[idx(i,j,Nx,Ny,3)];
                    //yn[4+3*7] = -yyndata[idx(i,j,Nx,Ny,3)];
                    //yp[5+3*7] = -yypdata[idx(i,j-1,Nx,Ny,3)];
                    //yn[5+3*7] = -yyndata[idx(i,j-1,Nx,Ny,3)];
                    //yp[6+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yn[6+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                }
            }
        }
        
        
        if (flag == 2){
            if (j<3){
                for (n=0;n<4;n++){
                    for (k=3;k<7;k++){
                        yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                    }
                }
                
                if (j==0){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                    }
                }
                
                if (j==1){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                    }
                }
                
                if (j==2){
                    for (n=0;n<4;n++){
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                    }
                }
            }
            
            if(j>Ny-4){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-3){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                    }
                }
            }
        }
        
        /* compute positive and negative solutions on the interface */
        for (n=0;n<4;n++){
            u_tpphy[n]=w0_py[n]*((2.0/6.0)*yp[1+n*7]-(7.0/6.0)*yp[2+n*7]+(11.0/6.0)*yp[3+n*7])+w1_py[n]*((-1.0/6.0)*yp[2+n*7]+(5.0/6.0)*yp[3+n*7]+(2.0/6.0)*yp[4+n*7])+w2_py[n]*((2.0/6.0)*yp[3+n*7]+(5.0/6.0)*yp[4+n*7]-(1.0/6.0)*yp[5+n*7]);
            
            u_tnphy[n]=w2_ny[n]*((-1.0/6.0)*yn[2+n*7]+(5.0/6.0)*yn[3+n*7]+(2.0/6.0)*yn[4+n*7])+w1_ny[n]*((2.0/6.0)*yn[3+n*7]+(5.0/6.0)*yn[4+n*7]-(1.0/6.0)*yn[5+n*7])+w0_ny[n]*((11.0/6.0)*yn[4+n*7]-(7.0/6.0)*yn[5+n*7]+(2.0/6.0)*yn[6+n*7]);
            
            u_tpnhy[n]=w0_py[n]*((2.0/6.0)*yp[0+n*7]-(7.0/6.0)*yp[1+n*7]+(11.0/6.0)*yp[2+n*7])+w1_py[n]*((-1.0/6.0)*yp[1+n*7]+(5.0/6.0)*yp[2+n*7]+(2.0/6.0)*yp[3+n*7])+w2_py[n]*((2.0/6.0)*yp[2+n*7]+(5.0/6.0)*yp[3+n*7]-(1.0/6.0)*yp[4+n*7]);
            
            u_tnnhy[n]=w2_ny[n]*((-1.0/6.0)*yn[1+n*7]+(5.0/6.0)*yn[2+n*7]+(2.0/6.0)*yn[3+n*7])+w1_ny[n]*((2.0/6.0)*yn[2+n*7]+(5.0/6.0)*yn[3+n*7]-(1.0/6.0)*yn[4+n*7])+w0_ny[n]*((11.0/6.0)*yn[3+n*7]-(7.0/6.0)*yn[4+n*7]+(2.0/6.0)*yn[5+n*7]);
        }
        
        delete []yp;
        delete []yn;
        
        return 0;
    }
}

/* Fill in values of tao in the whole domain */
//static int Gettao(N_Vector yxup, N_Vector yxdown, N_Vector yyup, N_Vector yydown, N_Vector taoup, N_Vector taodown, long int Nx, long int Ny, realtype gama)
static int Gettao(realtype *yxupdata, realtype *yxdowndata, realtype *yyupdata, realtype *yydowndata, N_Vector taoup, N_Vector taodown, long int Nx, long int Ny, realtype gama)
{
    /* declare parameters */
    long int i, j, NEQ;
    //realtype *vxupdata, *vyupdata, *vxdowndata, *vydowndata, *yxupdata, *yxdowndata, *yyupdata, *yydowndata, *taoupdata, *taodowndata;
    realtype *vxupdata, *vyupdata, *vxdowndata, *vydowndata, *taoupdata, *taodowndata;
    realtype pxup, pyup, pxdown, pydown;
    NEQ = Nx*Ny;
    
    /* create vectors */
    N_Vector vxup = NULL;
    N_Vector vyup = NULL;
    N_Vector vxdown = NULL;
    N_Vector vydown = NULL;
    
    /* Create serial vector of length N */
    vxup = N_VNew_Serial(NEQ);
    if (check_flag((void *) vxup, "N_VNew_Serial", 0)) return 1;
    vyup = N_VNew_Serial(NEQ);
    if (check_flag((void *) vyup, "N_VNew_Serial", 0)) return 1;
    vxdown = N_VNew_Serial(NEQ);
    if (check_flag((void *) vxdown, "N_VNew_Serial", 0)) return 1;
    vydown = N_VNew_Serial(NEQ);
    if (check_flag((void *) vydown, "N_VNew_Serial", 0)) return 1;
    
    /* Access data array for new NVector y, tao, vx, vy */
    vxupdata = N_VGetArrayPointer(vxup);
    if (check_flag((void *) vxupdata, "N_VGetArrayPointer", 0)) return 1;
    vyupdata = N_VGetArrayPointer(vyup);
    if (check_flag((void *) vyupdata, "N_VGetArrayPointer", 0)) return 1;
    vxdowndata = N_VGetArrayPointer(vxdown);
    if (check_flag((void *) vxdowndata, "N_VGetArrayPointer", 0)) return 1;
    vydowndata = N_VGetArrayPointer(vydown);
    if (check_flag((void *) vydowndata, "N_VGetArrayPointer", 0)) return 1;
    // yxupdata = N_VGetArrayPointer(yxup);
    //if (check_flag((void *) yxupdata, "N_VGetArrayPointer", 0)) return 1;
    // yyupdata = N_VGetArrayPointer(yyup);
    //if (check_flag((void *) yyupdata, "N_VGetArrayPointer", 0)) return 1;
    //yxdowndata = N_VGetArrayPointer(yxdown);
    //if (check_flag((void *) yxdowndata, "N_VGetArrayPointer", 0)) return 1;
    //yydowndata = N_VGetArrayPointer(yydown);
    //if (check_flag((void *) yydowndata, "N_VGetArrayPointer", 0)) return 1;
    taoupdata = N_VGetArrayPointer(taoup);
    if (check_flag((void *) taoupdata, "N_VGetArrayPointer", 0)) return 1;
    taodowndata = N_VGetArrayPointer(taodown);
    if (check_flag((void *) taodowndata, "N_VGetArrayPointer", 0)) return 1;
    
    /* Set values into vx and vy */
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            vxupdata[idx_v(i,j,Nx)]=yxupdata[idx(i, j, Nx, Ny, 1)]/yxupdata[idx(i, j, Nx, Ny, 0)];
	    //printf("i = %li, j=%li, yxupdata[idx(i, j, Nx, Ny, 1)]=%f,yxupdata[idx(i, j, Nx, Ny, 0)]=%f\n",i,j,yxupdata[idx(i, j, Nx, Ny, 1)],yxupdata[idx(i, j, Nx, Ny, 0)]);
            vyupdata[idx_v(i,j,Nx)]=yyupdata[idx(i, j, Nx, Ny, 2)]/yyupdata[idx(i, j, Nx, Ny, 0)];
            vxdowndata[idx_v(i,j,Nx)]=yxdowndata[idx(i, j, Nx, Ny, 1)]/yxdowndata[idx(i, j, Nx, Ny, 0)];
            vydowndata[idx_v(i,j,Nx)]=yydowndata[idx(i, j, Nx, Ny, 2)]/yydowndata[idx(i, j, Nx, Ny, 0)];
	    // printf("yydowndata[idx(i, j, Nx, Ny, 2)]=%f,yydowndata[idx(i, j, Nx, Ny, 0)]=%f\n",yydowndata[idx(i, j, Nx, Ny, 2)],yydowndata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    
    /* Compute the values of tao in the whole domain */
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            pxup = yxupdata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(yxupdata[idx(i, j, Nx, Ny, 3)]/yxupdata[idx(i, j, Nx, Ny, 0)]-0.5*(vxupdata[idx_v(i,j,Nx)]*vxupdata[idx_v(i,j,Nx)]+(yxupdata[idx(i, j, Nx, Ny, 2)]/yxupdata[idx(i, j, Nx, Ny, 0)])*(yxupdata[idx(i, j, Nx, Ny, 2)]/yxupdata[idx(i, j, Nx, Ny, 0)])));
            pyup = yyupdata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(yyupdata[idx(i, j, Nx, Ny, 3)]/yyupdata[idx(i, j, Nx, Ny, 0)]-0.5*((yyupdata[idx(i, j, Nx, Ny, 1)]/yyupdata[idx(i, j, Nx, Ny, 0)])*(yyupdata[idx(i, j, Nx, Ny, 1)]/yyupdata[idx(i, j, Nx, Ny, 0)])+vyupdata[idx_v(i,j,Nx)]*vyupdata[idx_v(i,j,Nx)]));
            
	    //printf("i = %li, j=%li, yxupdata[idx(i, j, Nx, Ny, 1)]=%f,vxupdata[idx_v(i,j,Nx)]=%f,pxup=%f\n",i,j,yxupdata[idx(i, j, Nx, Ny, 1)],vxupdata[idx_v(i,j,Nx)],pxup);
            taoupdata[idx(i, j, Nx, Ny, 0)] = yxupdata[idx(i, j, Nx, Ny, 1)]*vxupdata[idx_v(i,j,Nx)]+pxup;
            taoupdata[idx(i, j, Nx, Ny, 1)] = yyupdata[idx(i, j, Nx, Ny, 1)]*vyupdata[idx_v(i,j,Nx)];
            taoupdata[idx(i, j, Nx, Ny, 2)] = yxupdata[idx(i, j, Nx, Ny, 2)]*vxupdata[idx_v(i,j,Nx)];
            taoupdata[idx(i, j, Nx, Ny, 3)] = yyupdata[idx(i, j, Nx, Ny, 2)]*vyupdata[idx_v(i,j,Nx)]+pyup;
            
            pxdown = yxdowndata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(yxdowndata[idx(i, j, Nx, Ny, 3)]/yxdowndata[idx(i, j, Nx, Ny, 0)]-0.5*(vxdowndata[idx_v(i,j,Nx)]*vxdowndata[idx_v(i,j,Nx)]+(yxdowndata[idx(i, j, Nx, Ny, 2)]/yxdowndata[idx(i, j, Nx, Ny, 0)])*(yxdowndata[idx(i, j, Nx, Ny, 2)]/yxdowndata[idx(i, j, Nx, Ny, 0)])));
            pydown = yydowndata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(yydowndata[idx(i, j, Nx, Ny, 3)]/yydowndata[idx(i, j, Nx, Ny, 0)]-0.5*((yydowndata[idx(i, j, Nx, Ny, 1)]/yydowndata[idx(i, j, Nx, Ny, 0)])*(yydowndata[idx(i, j, Nx, Ny, 1)]/yydowndata[idx(i, j, Nx, Ny, 0)])+vydowndata[idx_v(i,j,Nx)]*vydowndata[idx_v(i,j,Nx)]));
            
            taodowndata[idx(i, j, Nx, Ny, 0)] = yxdowndata[idx(i, j, Nx, Ny, 1)]*vxdowndata[idx_v(i,j,Nx)]+pxdown;
	    //printf("yydowndata[idx(i, j, Nx, Ny, 1)]=%f,vydowndata[idx_v(i,j,Nx)]=%f\n",yydowndata[idx(i, j, Nx, Ny, 1)],vydowndata[idx_v(i,j,Nx)]); 
            taodowndata[idx(i, j, Nx, Ny, 1)] = yydowndata[idx(i, j, Nx, Ny, 1)]*vydowndata[idx_v(i,j,Nx)];
	    //printf("yxdowndata[idx(i, j, Nx, Ny, 2)]=%f,vxdowndata[idx_v(i,j,Nx)]=%f\n",yxdowndata[idx(i, j, Nx, Ny, 2)],vxdowndata[idx_v(i,j,Nx)]); 
            taodowndata[idx(i, j, Nx, Ny, 2)] = yxdowndata[idx(i, j, Nx, Ny, 2)]*vxdowndata[idx_v(i,j,Nx)];
            taodowndata[idx(i, j, Nx, Ny, 3)] = yydowndata[idx(i, j, Nx, Ny, 2)]*vydowndata[idx_v(i,j,Nx)]+pydown;
            
        }
    }
    
    /* Free vectors */
    N_VDestroy_Serial(vxup);
    N_VDestroy_Serial(vyup);
    N_VDestroy_Serial(vxdown);
    N_VDestroy_Serial(vydown);
    
    return 0;
}

/* Fill in the value of J in the whole domain */
//static int GetCj(N_Vector yxup, N_Vector yxdown, N_Vector yyup, N_Vector yydown, N_Vector Cjup, N_Vector Cjdown, long int Nx, long int Ny, realtype gama)
static int GetCj(realtype *yxupdata, realtype *yxdowndata, realtype *yyupdata, realtype *yydowndata, N_Vector Cjup, N_Vector Cjdown, long int Nx, long int Ny, realtype gama)
{
    /* declare parameters */
    long int i, j, NEQ;
    //realtype *vxupdata, *vyupdata, *vxdowndata, *vydowndata, *yxupdata, *yxdowndata, *yyupdata, *yydowndata, *Cjupdata, *Cjdowndata;
    realtype *vxupdata, *vyupdata, *vxdowndata, *vydowndata, *Cjupdata, *Cjdowndata;
    NEQ = Nx*Ny;    
    realtype pxup, pyup, pxdown, pydown;
        
    /* create vectors */
    N_Vector vxup = NULL;
    N_Vector vyup = NULL;
    N_Vector vxdown = NULL;
    N_Vector vydown = NULL;
    
    /* Create serial vector of length N */
    vxup = N_VNew_Serial(NEQ);
    if (check_flag((void *) vxup, "N_VNew_Serial", 0)) return 1;
    vyup = N_VNew_Serial(NEQ);
    if (check_flag((void *) vyup, "N_VNew_Serial", 0)) return 1;
    vxdown = N_VNew_Serial(NEQ);
    if (check_flag((void *) vxdown, "N_VNew_Serial", 0)) return 1;
    vydown = N_VNew_Serial(NEQ);
    if (check_flag((void *) vydown, "N_VNew_Serial", 0)) return 1;
    
    /* Access data array for new NVector y, tao, vx, vy */
    vxupdata = N_VGetArrayPointer(vxup);
    if (check_flag((void *) vxupdata, "N_VGetArrayPointer", 0)) return 1;
    vyupdata = N_VGetArrayPointer(vyup);
    if (check_flag((void *) vyupdata, "N_VGetArrayPointer", 0)) return 1;
    vxdowndata = N_VGetArrayPointer(vxdown);
    if (check_flag((void *) vxdowndata, "N_VGetArrayPointer", 0)) return 1;
    vydowndata = N_VGetArrayPointer(vydown);
    if (check_flag((void *) vydowndata, "N_VGetArrayPointer", 0)) return 1;
    //yxupdata = N_VGetArrayPointer(yxup);
    //if (check_flag((void *) yxupdata, "N_VGetArrayPointer", 0)) return 1;
    //yyupdata = N_VGetArrayPointer(yyup);
    //if (check_flag((void *) yyupdata, "N_VGetArrayPointer", 0)) return 1;
    // yxdowndata = N_VGetArrayPointer(yxdown);
    //if (check_flag((void *) yxdowndata, "N_VGetArrayPointer", 0)) return 1;
    //yydowndata = N_VGetArrayPointer(yydown);
    //if (check_flag((void *) yydowndata, "N_VGetArrayPointer", 0)) return 1;
    Cjupdata = N_VGetArrayPointer(Cjup);
    if (check_flag((void *) Cjupdata, "N_VGetArrayPointer", 0)) return 1;
    Cjdowndata = N_VGetArrayPointer(Cjdown);
    if (check_flag((void *) Cjdowndata, "N_VGetArrayPointer", 0)) return 1;
    
    /* Set values into vx and vy */
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            vxupdata[idx_v(i,j,Nx)]=yxupdata[idx(i, j, Nx, Ny, 1)]/yxupdata[idx(i, j, Nx, Ny, 0)];
            vyupdata[idx_v(i,j,Nx)]=yyupdata[idx(i, j, Nx, Ny, 2)]/yyupdata[idx(i, j, Nx, Ny, 0)];
            vxdowndata[idx_v(i,j,Nx)]=yxdowndata[idx(i, j, Nx, Ny, 1)]/yxdowndata[idx(i, j, Nx, Ny, 0)];
            vydowndata[idx_v(i,j,Nx)]=yydowndata[idx(i, j, Nx, Ny, 2)]/yydowndata[idx(i, j, Nx, Ny, 0)];
        }
    }
    
    /* Compute the values of Cj in the whole domain */
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            pxup = yxupdata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(yxupdata[idx(i, j, Nx, Ny, 3)]/yxupdata[idx(i, j, Nx, Ny, 0)]-0.5*(vxupdata[idx_v(i,j,Nx)]*vxupdata[idx_v(i,j,Nx)]+(yxupdata[idx(i, j, Nx, Ny, 2)]/yxupdata[idx(i, j, Nx, Ny, 0)])*(yxupdata[idx(i, j, Nx, Ny, 2)]/yxupdata[idx(i, j, Nx, Ny, 0)])));
            pyup = yyupdata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(yyupdata[idx(i, j, Nx, Ny, 3)]/yyupdata[idx(i, j, Nx, Ny, 0)]-0.5*((yyupdata[idx(i, j, Nx, Ny, 1)]/yyupdata[idx(i, j, Nx, Ny, 0)])*(yyupdata[idx(i, j, Nx, Ny, 1)]/yyupdata[idx(i, j, Nx, Ny, 0)])+vyupdata[idx_v(i,j,Nx)]*vyupdata[idx_v(i,j,Nx)]));
            Cjupdata[idx(i, j, Nx, Ny, 0)] = vxupdata[idx_v(i,j,Nx)]*(yxupdata[idx(i, j, Nx, Ny, 3)]+pxup);
            Cjupdata[idx(i, j, Nx, Ny, 1)] = vyupdata[idx_v(i,j,Nx)]*(yyupdata[idx(i, j, Nx, Ny, 3)]+pyup);
            
            pxdown = yxdowndata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(yxdowndata[idx(i, j, Nx, Ny, 3)]/yxdowndata[idx(i, j, Nx, Ny, 0)]-0.5*(vxdowndata[idx_v(i,j,Nx)]*vxdowndata[idx_v(i,j,Nx)]+(yxdowndata[idx(i, j, Nx, Ny, 2)]/yxdowndata[idx(i, j, Nx, Ny, 0)])*(yxdowndata[idx(i, j, Nx, Ny, 2)]/yxdowndata[idx(i, j, Nx, Ny, 0)])));
            pydown = yydowndata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(yydowndata[idx(i, j, Nx, Ny, 3)]/yydowndata[idx(i, j, Nx, Ny, 0)]-0.5*((yydowndata[idx(i, j, Nx, Ny, 1)]/yydowndata[idx(i, j, Nx, Ny, 0)])*(yydowndata[idx(i, j, Nx, Ny, 1)]/yydowndata[idx(i, j, Nx, Ny, 0)])+vydowndata[idx_v(i,j,Nx)]*vydowndata[idx_v(i,j,Nx)]));
            Cjdowndata[idx(i, j, Nx, Ny, 0)] = vxdowndata[idx_v(i,j,Nx)]*(yxdowndata[idx(i, j, Nx, Ny, 3)]+pxdown);
            Cjdowndata[idx(i, j, Nx, Ny, 1)] = vydowndata[idx_v(i,j,Nx)]*(yydowndata[idx(i, j, Nx, Ny, 3)]+pydown);
        }
    }
    
    /* Free vectors */
    N_VDestroy_Serial(vxup);
    N_VDestroy_Serial(vyup);
    N_VDestroy_Serial(vxdown);
    N_VDestroy_Serial(vydown);
    
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
