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
    int bcflag;    /* boundary conditions      */
} *UserData;

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/* Private function to check function return values */
static int check_flag(void *flagvalue, const string funcname, int opt);

/* Set value of tao in whole domain*/
static int Gettao(realtype *yxbackupdata, realtype *yxbackdowndata, realtype *yybackupdata, realtype *yybackdowndata, realtype *taoup, realtype *taodown, long int Nx, long int Ny, realtype gama);

/* Set value of J in whole domain*/
static int GetCj(realtype *yxbackupdata, realtype *yxbackdowndata, realtype *yybackupdata, realtype *yybackdowndata, realtype *Cjup, realtype *Cjdown, long int Nx, long int Ny, realtype gama);

static int Getmaxegvx(realtype *egvx, realtype *egxmax, long int Nx, long int Ny);

static int Getmaxegvy(realtype *egvy, realtype *egymax, long int Nx, long int Ny);

/* Split f to f_positve and f_negative parts for x component */
static int Splitfluxesx(realtype *Ydata, realtype *fp, realtype *fn, long int i, long int j, long int Nx, long int Ny, realtype gama);

/* Split g to g_positve and g_negative parts for y component */
static int Splitfluxesy(realtype *Ydata, realtype *gp, realtype *gn, long int i, long int j, long int Nx, long int Ny, realtype gama);

static int Getegvx(realtype *Ydata, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype *egvx, int flag);

static int Getegvy(realtype *Ydata, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype *egvy, int flag);

/* Set left eigenvectors in x direction and eigenvalues*/
static int Setlfxegm(realtype *Ydata, realtype **lfxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag);

/* Set right eigenvectors in x direction */
static int Setrhxegm(realtype *Ydata, realtype **rhxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag);

/* Set left eigenvectors in y direction and eigenvalues*/
static int Setlfyegm(realtype *Ydata, realtype **lfyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag);

/* Set right eigenvectors in y direction */
static int Setrhyegm(realtype *Ydata, realtype **rhyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag);

/* Fill in the indicators of smoothness on x direction */
static int SetISX(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx,realtype *yxpdata, realtype *yxndata, long int i, long int j, long int Nx, long int Ny);

/* Fill in the indicators of smoothness on y direction */
static int SetISY(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny,realtype *yypdata, realtype *yyndata, long int i, long int j, long int Nx, long int Ny);

/* Fill in the stencil weights on x direction */
static int Setalphawx(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx, realtype *alpha_0px, realtype *alpha_1px, realtype *alpha_2px, realtype *alpha_0nx, realtype *alpha_1nx, realtype *alpha_2nx, realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype Epsilon);

/* Fill in the stencil weights on y direction */
static int Setalphawy(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny, realtype *alpha_0py, realtype *alpha_1py, realtype *alpha_2py, realtype *alpha_0ny, realtype *alpha_1ny, realtype *alpha_2ny, realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype Epsilon);

/* Get the interface value on x direction */
static int SetUx(realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype *yxpdata, realtype *yxndata, realtype *u_tpphx, realtype *u_tnphx, realtype *u_tpnhx, realtype *u_tnnhx, long int i, long int j, long int Nx, long int Ny);

/* Get the interface value on y direction */
static int SetUy(realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype *yypdata, realtype *yyndata, realtype *u_tpphy, realtype *u_tnphy, realtype *u_tpnhy, realtype *u_tnnhy, long int i, long int j, long int Nx, long int Ny);

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
    int bcflag;
    
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
    flag = fscanf(FID," bcflag = %i\n", &bcflag);
    fclose(FID);
    
    /* store the inputs in the UserData structure */
    udata->Nx = Nx;
    udata->Ny = Ny;
    udata->Lx = Lx;
    udata->Ly = Ly;
    udata->gama = gama;
    udata->bcflag = bcflag;
    
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
    
    /* Set initial data */
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
    long int i, j, k;
    int flag;
    realtype p, Epsilon;
    Epsilon = 10.0e-6;
    
    /* create relative 2D WENO arrays */
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
    int bcflag = udata->bcflag;
    
    /* fill in relative data arrays */
    realtype *yxupdata = new realtype [4*(Nx)*Ny];
    realtype *yyupdata = new realtype [4*Nx*(Ny)];
    realtype *yxdowndata = new realtype [4*(Nx)*Ny];
    realtype *yydowndata = new realtype [4*Nx*(Ny)];
    realtype *yxpdata = new realtype [4*7];
    realtype *yxndata = new realtype [4*7];
    realtype *yypdata = new realtype [4*7];
    realtype *yyndata = new realtype [4*7];
    realtype *yxnewdata = new realtype [4*7];
    realtype *yynewdata = new realtype [4*7];
    realtype *yxnew = new realtype [4*7];
    realtype *yynew = new realtype [4*7];
    realtype *yxbackupdata = new realtype [4*(Nx)*Ny];
    realtype *yybackupdata = new realtype [4*Nx*(Ny)];
    realtype *yxbackdowndata = new realtype [4*(Nx)*Ny];
    realtype *yybackdowndata = new realtype [4*Nx*(Ny)];
    realtype *taoup = new realtype [4*(Nx)*Ny];
    realtype *Cjup = new realtype [4*(Nx)*Ny];
    realtype *taodown = new realtype [4*(Nx)*Ny];
    realtype *Cjdown = new realtype [4*(Nx)*Ny];
    realtype *egvx = new realtype [4*(Nx)*Ny];
    realtype *egvy = new realtype [4*Nx*(Ny)];
    //realtype *egvxmax = new realtype [Nx*Ny];
    //realtype *egvymax = new realtype [Nx*Ny];
    realtype *egxmax = new realtype [4];
    realtype *egymax = new realtype [4];
    
    /*
     realtype *yxptdata = new realtype [4*Nx*Ny];
     realtype *yxntdata = new realtype [4*Nx*Ny];
     realtype *yyptdata = new realtype [4*Nx*Ny];
     realtype *yyntdata = new realtype [4*Nx*Ny];
     */
    
    /* access data arrays */
    realtype *Ydata = N_VGetArrayPointer(y);
    if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;
    realtype *dYdata = N_VGetArrayPointer(ydot);
    if (check_flag((void *) dYdata, "N_VGetArrayPointer", 0)) return 1;
    /*
     for(j=0;j<Ny;j++){
     for(i=0;i<Nx;i++){
     flag = Splitfluxesx(Ydata, yxptdata, yxntdata, i, j, Nx, Ny, gama);
     if (flag!=0) printf("error in Splitfluxesx function \n");
     }
     }
     
     for(i=0;i<Nx;i++){
     for(j=0;j<Ny;j++){
     flag = Splitfluxesy(Ydata, yyptdata, yyntdata, i, j, Nx, Ny, gama);
     if (flag!=0) printf("error in Splitfluxesy function \n");
     }
     }
     */
    /* compute eigenvalues and fill in positive and negative characteristic variables yxpdata, yxndata, x direction */
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
            flag = Getegvx(Ydata, i, j, Nx, Ny, gama, egvx, bcflag);
            if (flag!=0) printf("error in Getegvx function \n");
        }
    }
    
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
            flag = Getmaxegvx(egvx, egxmax, Nx, Ny);
            if (flag!=0) printf("error in Getmaxegx function \n");
        }
    }
    
    //printf("egvmax: egxmax[0]=%lf,egxmax[1]=%lf,egxmax[2]=%lf,egxmax[3]=%lf\n",egxmax[0],egxmax[1],egxmax[2],egxmax[3]);
    /*
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
            printf(" Ydata : i = %li, j = %li, Ydata[idx(i, j, Nx, Ny, 0)] = %f, Ydata[idx(i, j, Nx, Ny, 1)] = %f, Ydata[idx(i, j, Nx, Ny, 2)] = %f, Ydata[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, Ydata[idx(i, j, Nx, Ny, 0)], Ydata[idx(i, j, Nx, Ny, 1)],Ydata[idx(i, j, Nx, Ny, 2)],Ydata[idx(i, j, Nx, Ny, 3)]);
        }
    }
    */

    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
            //printf(" Ydata : i = %li, j = %li, Ydata[idx(i, j, Nx, Ny, 0)] = %f, Ydata[idx(i, j, Nx, Ny, 1)] = %f, Ydata[idx(i, j, Nx, Ny, 2)] = %f, Ydata[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, Ydata[idx(i, j, Nx, Ny, 0)], Ydata[idx(i, j, Nx, Ny, 1)],Ydata[idx(i, j, Nx, Ny, 2)],Ydata[idx(i, j, Nx, Ny, 3)]);
            
            if(i>2&&i<Nx-3)
            {
                for(k=0;k<4;k++)
                {
                    yxnew[0+k*7]=Ydata[idx(i-3,j,Nx,Ny,k)];
                    yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                    yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                    yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                    yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                    yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                    yxnew[6+k*7]=Ydata[idx(i+3,j,Nx,Ny,k)];
                }
            }
            
            if (bcflag==0)
            {
                if(i==0){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(Nx-3,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(Nx-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(Nx-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i+3,j,Nx,Ny,k)];
                    }
                }
                
                if(i==1){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(Nx-2,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(Nx-1,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i+3,j,Nx,Ny,k)];
                    }
                }
            
                if(i==2){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(Nx-1,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i+3,j,Nx,Ny,k)];
                    }
                }
                
                if(i==Nx-1){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i-3,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(0,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(1,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(2,j,Nx,Ny,k)];
                    }
                }
                
                if(i==Nx-2){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i-3,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(0,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(1,j,Nx,Ny,k)];
                    }
                }
                
                if(i==Nx-3){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i-3,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(0,j,Nx,Ny,k)];
                    }
                }
            }
            
            if (bcflag==1)
            {
                if(i==0){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i+3,j,Nx,Ny,k)];
                    }
                    yxnew[0+1*7]=-Ydata[idx(i+2,j,Nx,Ny,1)];
                    yxnew[1+1*7]=-Ydata[idx(i+1,j,Nx,Ny,1)];
                    yxnew[2+1*7]=-Ydata[idx(i,j,Nx,Ny,1)];
                }
                
                if(i==1){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i+3,j,Nx,Ny,k)];
                    }
                    yxnew[0+1*7]=-Ydata[idx(i,j,Nx,Ny,1)];
                    yxnew[1+1*7]=-Ydata[idx(i-1,j,Nx,Ny,1)];
                }
                
                if(i==2){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i+3,j,Nx,Ny,k)];
                    }
                    yxnew[0+1*7]=-Ydata[idx(i-2,j,Nx,Ny,1)];
                }
                
                if(i==Nx-1){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i-3,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                    }
                    yxnew[4+1*7]=-Ydata[idx(i,j,Nx,Ny,1)];
                    yxnew[5+1*7]=-Ydata[idx(i-1,j,Nx,Ny,1)];
                    yxnew[6+1*7]=-Ydata[idx(i-2,j,Nx,Ny,1)];
                }
                
                if(i==Nx-2){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i-3,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                    }
                    yxnew[5+1*7]=-Ydata[idx(i+1,j,Nx,Ny,1)];
                    yxnew[6+1*7]=-Ydata[idx(i,j,Nx,Ny,1)];
                }
                
                if(i==Nx-3){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i-3,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                    }
                    yxnew[6+1*7]=-Ydata[idx(i+2,j,Nx,Ny,1)];
                }
            }
            
            if (bcflag==2)
            {
                if(i==0){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i+3,j,Nx,Ny,k)];
                    }
                   
                }
                
                if(i==1){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i+3,j,Nx,Ny,k)];
                    }
                   
                }
                
                if(i==2){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i+3,j,Nx,Ny,k)];
                    }
                  
                }
                
                if(i==Nx-1){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i-3,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                    }
                  
                }
                
                if(i==Nx-2){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i-3,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                    }
                  
                }
                
                if(i==Nx-3){
                    for(k=0;k<4;k++)
                    {
                        yxnew[0+k*7]=Ydata[idx(i-3,j,Nx,Ny,k)];
                        yxnew[1+k*7]=Ydata[idx(i-2,j,Nx,Ny,k)];
                        yxnew[2+k*7]=Ydata[idx(i-1,j,Nx,Ny,k)];
                        yxnew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yxnew[4+k*7]=Ydata[idx(i+1,j,Nx,Ny,k)];
                        yxnew[5+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                        yxnew[6+k*7]=Ydata[idx(i+2,j,Nx,Ny,k)];
                    }
               
                }
            }
            
            flag = Setlfxegm(Ydata, lfxegm, i, j, Nx, Ny, gama, bcflag);
            if (flag!=0) printf("error in Setlfxegm function \n");
            
            /* compute characteristic variables */
            for (k=0;k<4;k++){
                //printf(" lfxegm : i = %li, j = %li, lfxegm[k][0] = %f, lfxegm[k][1] = %f, lfxegm[k][2] = %f, lfxegm[k][3] = %f\n", i, j, lfxegm[k][0], lfxegm[k][1],lfxegm[k][2],lfxegm[k][3]);
                yxnewdata[0+k*7] = lfxegm[k][0]*yxnew[0+0*7]+lfxegm[k][1]*yxnew[0+1*7]+lfxegm[k][2]*yxnew[0+2*7]+lfxegm[k][3]*yxnew[0+3*7];
                yxnewdata[1+k*7] = lfxegm[k][0]*yxnew[1+0*7]+lfxegm[k][1]*yxnew[1+1*7]+lfxegm[k][2]*yxnew[1+2*7]+lfxegm[k][3]*yxnew[1+3*7];
                yxnewdata[2+k*7] = lfxegm[k][0]*yxnew[2+0*7]+lfxegm[k][1]*yxnew[2+1*7]+lfxegm[k][2]*yxnew[2+2*7]+lfxegm[k][3]*yxnew[2+3*7];
                yxnewdata[3+k*7] = lfxegm[k][0]*yxnew[3+0*7]+lfxegm[k][1]*yxnew[3+1*7]+lfxegm[k][2]*yxnew[3+2*7]+lfxegm[k][3]*yxnew[3+3*7];
                yxnewdata[4+k*7] = lfxegm[k][0]*yxnew[4+0*7]+lfxegm[k][1]*yxnew[4+1*7]+lfxegm[k][2]*yxnew[4+2*7]+lfxegm[k][3]*yxnew[4+3*7];
                yxnewdata[5+k*7] = lfxegm[k][0]*yxnew[5+0*7]+lfxegm[k][1]*yxnew[5+1*7]+lfxegm[k][2]*yxnew[5+2*7]+lfxegm[k][3]*yxnew[5+3*7];
                yxnewdata[6+k*7] = lfxegm[k][0]*yxnew[6+0*7]+lfxegm[k][1]*yxnew[6+1*7]+lfxegm[k][2]*yxnew[6+2*7]+lfxegm[k][3]*yxnew[6+3*7];
                //yxpdata[idx(i, j, Nx, Ny, k)] = lfxegm[k][0]*yxptdata[idx(i, j, Nx, Ny, 0)]+lfxegm[k][1]*yxptdata[idx(i, j, Nx, Ny, 1)]+lfxegm[k][2]*yxptdata[idx(i, j, Nx, Ny, 2)]+lfxegm[k][3]*yxptdata[idx(i, j, Nx, Ny, 3)];
                //yxndata[idx(i, j, Nx, Ny, k)] = lfxegm[k][0]*yxntdata[idx(i, j, Nx, Ny, 0)]+lfxegm[k][1]*yxntdata[idx(i, j, Nx, Ny, 1)]+lfxegm[k][2]*yxntdata[idx(i, j, Nx, Ny, 2)]+lfxegm[k][3]*yxntdata[idx(i, j, Nx, Ny, 3)];
            }
            
            for (k=0;k<4;k++){
                
                yxpdata[0+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[0+k*7]+egxmax[k]*yxnewdata[0+k*7]);
                yxndata[0+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[0+k*7]-egxmax[k]*yxnewdata[0+k*7]);
                yxpdata[1+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[1+k*7]+egxmax[k]*yxnewdata[1+k*7]);
                yxndata[1+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[1+k*7]-egxmax[k]*yxnewdata[1+k*7]);
                yxpdata[2+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[2+k*7]+egxmax[k]*yxnewdata[2+k*7]);
                yxndata[2+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[2+k*7]-egxmax[k]*yxnewdata[2+k*7]);
                yxpdata[3+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[3+k*7]+egxmax[k]*yxnewdata[3+k*7]);
                yxndata[3+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[3+k*7]-egxmax[k]*yxnewdata[3+k*7]);
                yxpdata[4+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[4+k*7]+egxmax[k]*yxnewdata[4+k*7]);
                yxndata[4+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[4+k*7]-egxmax[k]*yxnewdata[4+k*7]);
                yxpdata[5+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[5+k*7]+egxmax[k]*yxnewdata[5+k*7]);
                yxndata[5+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[5+k*7]-egxmax[k]*yxnewdata[5+k*7]);
                yxpdata[6+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[6+k*7]+egxmax[k]*yxnewdata[6+k*7]);
                yxndata[6+k*7]=0.5*(egvx[idx(i,j,Nx,Ny,k)]*yxnewdata[6+k*7]-egxmax[k]*yxnewdata[6+k*7]);
                
                /*
                 yxpdata[0+k*7]=0.5*(egvx[idx(i,j,Nx+1,Ny,k)]*yxnewdata[0+k*7]+fabs(egvx[idx(i,j,Nx+1,Ny,k)])*yxnewdata[0+k*7]);
                 yxndata[0+k*7]=0.5*(egvx[idx(i,j,Nx+1,Ny,k)]*yxnewdata[0+k*7]-fabs(egvx[idx(i,j,Nx+1,Ny,k)])*yxnewdata[0+k*7]);
                 yxpdata[1+k*7]=0.5*(egvx[idx(i,j,Nx+1,Ny,k)]*yxnewdata[1+k*7]+fabs(egvx[idx(i,j,Nx+1,Ny,k)])*yxnewdata[1+k*7]);
                 yxndata[1+k*7]=0.5*(egvx[idx(i,j,Nx+1,Ny,k)]*yxnewdata[1+k*7]-fabs(egvx[idx(i,j,Nx+1,Ny,k)])*yxnewdata[1+k*7]);
                 yxpdata[2+k*7]=0.5*(egvx[idx(i,j,Nx+1,Ny,k)]*yxnewdata[2+k*7]+fabs(egvx[idx(i,j,Nx+1,Ny,k)])*yxnewdata[2+k*7]);
                 yxndata[2+k*7]=0.5*(egvx[idx(i,j,Nx+1,Ny,k)]*yxnewdata[2+k*7]-fabs(egvx[idx(i,j,Nx+1,Ny,k)])*yxnewdata[2+k*7]);
                 yxpdata[3+k*7]=0.5*(egvx[idx(i,j,Nx+1,Ny,k)]*yxnewdata[3+k*7]+fabs(egvx[idx(i,j,Nx+1,Ny,k)])*yxnewdata[3+k*7]);
                 yxndata[3+k*7]=0.5*(egvx[idx(i,j,Nx+1,Ny,k)]*yxnewdata[3+k*7]-fabs(egvx[idx(i,j,Nx+1,Ny,k)])*yxnewdata[3+k*7]);
                 yxpdata[4+k*7]=0.5*(egvx[idx(i,j,Nx+1,Ny,k)]*yxnewdata[4+k*7]+fabs(egvx[idx(i,j,Nx+1,Ny,k)])*yxnewdata[4+k*7]);
                 yxndata[4+k*7]=0.5*(egvx[idx(i,j,Nx+1,Ny,k)]*yxnewdata[4+k*7]-fabs(egvx[idx(i,j,Nx+1,Ny,k)])*yxnewdata[4+k*7]);
                 yxpdata[5+k*7]=0.5*(egvx[idx(i,j,Nx+1,Ny,k)]*yxnewdata[5+k*7]+fabs(egvx[idx(i,j,Nx+1,Ny,k)])*yxnewdata[5+k*7]);
                 yxndata[5+k*7]=0.5*(egvx[idx(i,j,Nx+1,Ny,k)]*yxnewdata[5+k*7]-fabs(egvx[idx(i,j,Nx+1,Ny,k)])*yxnewdata[5+k*7]);
                 */
            }
            
            /* iterate over domain, using WENO to compute all equations on x direction */
            
            /* get indicator of smoothness on x direction */
            flag = SetISX(IS0_px, IS1_px, IS2_px, IS0_nx, IS1_nx, IS2_nx, yxpdata, yxndata, i, j, Nx, Ny);
            if (flag!=0) printf("error in SetISX function \n");
            //for (n=0;n<4;n++){
            // printf("i=%li,j=%li,n=%li,IS0_px[n]=%f,IS1_px[n]=%f,IS2_px[n]=%f,IS0_nx[n]=%f,IS1_nx[n]=%f,IS2_nx[n]=%f\n",i,j,n,IS0_px[n],IS1_px[n],IS2_px[n],IS0_nx[n],IS1_nx[n],IS2_nx[n]);
            // }
            /* get weight on x direction for four variables */
            flag = Setalphawx(IS0_px, IS1_px, IS2_px, IS0_nx, IS1_nx, IS2_nx, alpha_0px, alpha_1px, alpha_2px, alpha_0nx, alpha_1nx, alpha_2nx, w0_px, w1_px, w2_px, w0_nx, w1_nx, w2_nx, Epsilon);
            if (flag!=0) printf("error in Setalphawx function \n");
            
            /* compute the positive and negative parts of right interface on x direction */
            flag = SetUx(w0_px, w1_px, w2_px, w0_nx, w1_nx, w2_nx, yxpdata, yxndata, u_tpphx, u_tnphx, u_tpnhx, u_tnnhx, i, j, Nx, Ny);
            if (flag!=0) printf("error in SetUx function \n");
            
            /* get the interface values */
            for(k=0;k<4;k++){
                //printf("w: i=%li, j=%li, k=%li, w0_px[k]=%g,w1_px[k]=%g,w2_px[k]=%g,w0_nx[k]=%g,w1_nx[k]=%g,w2_nx[k]=%g\n",i,j,k,w0_px[k],w1_px[k],w2_px[k],w0_nx[k],w1_nx[k],w2_nx[k]);
                yxupdata[idx(i, j, Nx, Ny, k)]=(u_tpphx[k]+u_tnphx[k]);
                yxdowndata[idx(i, j, Nx, Ny, k)]=(u_tpnhx[k]+u_tnnhx[k]);
                //printf("yxdata: i=%li, j=%li, k=%li, yxdata[idx(i, j, Nx+1, Ny, k)]=%g\n",i,j,k,yxdata[idx(i, j, Nx+1, Ny, k)]);
                //printf("yx: i=%li, j=%li, k=%li, u_tpphx[k]=%f, u_tnphx[k]=%f\n",i,j, k,u_tpphx[k], u_tnphx[k]);
            }
            
            /* transform the interface values back to the physical ones on x direction*/
            flag = Setrhxegm(Ydata, rhxegm, i+1, j, Nx, Ny, gama, bcflag);
            if (flag!=0) printf("error in Setrhxegm function \n");
            for (k=0;k<4;k++){
                
                //printf(" rhxegm : i = %li, j = %li, k = %li, rhxegm[k][0] = %f, rhxegm[k][1] = %f, rhxegm[k][2] = %f, rhxegm[k][3] = %f\n", i, j, k, rhxegm[k][0], rhxegm[k][1],rhxegm[k][2],rhxegm[k][3]);
                
                yxbackupdata[idx(i, j, Nx, Ny, k)] = rhxegm[k][0]*yxupdata[idx(i, j, Nx, Ny, 0)]+rhxegm[k][1]*yxupdata[idx(i, j, Nx, Ny, 1)]+rhxegm[k][2]*yxupdata[idx(i, j, Nx, Ny, 2)]+rhxegm[k][3]*yxupdata[idx(i, j, Nx, Ny, 3)];
            
            }
            flag = Setrhxegm(Ydata, rhxegm, i, j, Nx, Ny, gama, bcflag);
            if (flag!=0) printf("error in Setrhxegm function \n");
            for (k=0;k<4;k++){
                
                //printf(" rhxegm : i = %li, j = %li, k = %li, rhxegm[k][0] = %f, rhxegm[k][1] = %f, rhxegm[k][2] = %f, rhxegm[k][3] = %f\n", i, j, k, rhxegm[k][0], rhxegm[k][1],rhxegm[k][2],rhxegm[k][3]);
                
                yxbackdowndata[idx(i, j, Nx, Ny, k)] = rhxegm[k][0]*yxdowndata[idx(i, j, Nx, Ny, 0)]+rhxegm[k][1]*yxdowndata[idx(i, j, Nx, Ny, 1)]+rhxegm[k][2]*yxdowndata[idx(i, j, Nx, Ny, 2)]+rhxegm[k][3]*yxdowndata[idx(i, j, Nx, Ny, 3)];
            }

            //printf(" yxdata : i = %li, j = %li, yxdata[idx(i, j, Nx+1, Ny, 0)]=%f, yxdata[idx(i, j, Nx+1, Ny, 1)]=%f, yxdata[idx(i, j, Nx+1, Ny, 2)]=%f, yxdata[idx(i, j, Nx+1, Ny, 3)]=%f\n",i,j,yxdata[idx(i, j, Nx+1, Ny, 0)],yxdata[idx(i, j, Nx+1, Ny, 1)],yxdata[idx(i, j, Nx+1, Ny, 2)],yxdata[idx(i, j, Nx+1, Ny, 3)]);
            //printf("yxbackdata : i = %li, j = %li, yxbackdata[idx(i, j, Nx+1, Ny, 0)]=%f, yxbackdata[idx(i, j, Nx+1, Ny, 1)]=%f, yxbackdata[idx(i, j, Nx+1, Ny, 2)]=%f, yxbackdata[idx(i, j, Nx+1, Ny, 3)]=%e\n",i,j,yxbackdata[idx(i, j, Nx+1, Ny, 0)],yxbackdata[idx(i, j, Nx+1, Ny, 1)],yxbackdata[idx(i, j, Nx+1, Ny, 2)],yxbackdata[idx(i, j, Nx+1, Ny, 3)]);
        }
    }
    
    /* compute eigenvalues and fill in positive and negative characteristic variables yxpdata, yxndata, x direction */
    for(i=0; i<Nx; i++){
        for (j=0; j<Ny; j++){
            flag =  Getegvy(Ydata, i, j, Nx, Ny, gama, egvy, bcflag);
            if (flag!=0) printf("error in Getegvy function \n");
        }
    }
    
    for(i=0; i<Nx; i++){
        for (j=0; j<Ny; j++){
            flag = Getmaxegvy(egvy, egymax, Nx, Ny);
            if (flag!=0) printf("error in Getmaxegy function \n");
        }
    }
    
    for(i=0; i<Nx; i++){
        for (j=0; j<Ny; j++){
            //printf(" Ydata : i = %li, j = %li, Ydata[idx(i, j, Nx, Ny, 0)] = %f, Ydata[idx(i, j, Nx, Ny, 1)] = %f, Ydata[idx(i, j, Nx, Ny, 2)] = %f, Ydata[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, Ydata[idx(i, j, Nx, Ny, 0)], Ydata[idx(i, j, Nx, Ny, 1)],Ydata[idx(i, j, Nx, Ny, 2)],Ydata[idx(i, j, Nx, Ny, 3)]);
            
            if(j>2&&j<Ny-3)
            {
                for(k=0;k<4;k++)
                {
                    yynew[0+k*7]=Ydata[idx(i,j-3,Nx,Ny,k)];
                    yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                    yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                    yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                    yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                    yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                    yynew[6+k*7]=Ydata[idx(i,j+3,Nx,Ny,k)];
                }
            }
            
            if (bcflag==0)
            {
                if(j==0){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,Ny-3,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,Ny-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,Ny-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j+3,Nx,Ny,k)];
                    }
                }
                
                if(j==1){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,Ny-2,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,Ny-1,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j+3,Nx,Ny,k)];
                    }
                }
                
                if(j==2){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,Ny-1,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j+3,Nx,Ny,k)];
                    }
                }
                
                if(j==Ny-1){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j-3,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,0,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,1,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,2,Nx,Ny,k)];
                    }
                }
                
                if(j==Ny-2){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j-3,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,0,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,1,Nx,Ny,k)];
                    }
                }
                
                if(j==Ny-3){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j-3,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,0,Nx,Ny,k)];
                    }
                }
            }
            
            if (bcflag==1)
            {
                if(j==0){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j+3,Nx,Ny,k)];
                    }
                    yynew[0+1*7]=-Ydata[idx(i,j+2,Nx,Ny,2)];
                    yynew[1+1*7]=-Ydata[idx(i,j+1,Nx,Ny,2)];
                    yynew[2+1*7]=-Ydata[idx(i,j,Nx,Ny,2)];
                }
                
                if(j==1){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j+3,Nx,Ny,k)];
                    }
                    yynew[0+1*7]=-Ydata[idx(i,j,Nx,Ny,2)];
                    yynew[1+1*7]=-Ydata[idx(i,j-1,Nx,Ny,2)];
                }
                
                if(j==2){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j+3,Nx,Ny,k)];
                    }
                    yynew[0+1*7]=-Ydata[idx(i,j-2,Nx,Ny,2)];
                }
                
                if(j==Ny-1){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j-3,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                    }
                    yynew[4+1*7]=-Ydata[idx(i,j,Nx,Ny,2)];
                    yynew[5+1*7]=-Ydata[idx(i,j-1,Nx,Ny,2)];
                    yynew[6+1*7]=-Ydata[idx(i,j-2,Nx,Ny,2)];
                }
                
                if(j==Ny-2){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j-3,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                    }
                    yynew[5+1*7]=-Ydata[idx(i,j+1,Nx,Ny,2)];
                    yynew[6+1*7]=-Ydata[idx(i,j,Nx,Ny,2)];
                }
                
                if(j==Ny-3){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j-3,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                    }
                    yynew[6+1*7]=-Ydata[idx(i,j+2,Nx,Ny,2)];
                }
            }
            
            if (bcflag==2)
            {
                if(j==0){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j+3,Nx,Ny,k)];
                    }
            }
                
                if(j==1){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j+3,Nx,Ny,k)];
                    }
            }
                
                if(j==2){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j+3,Nx,Ny,k)];
                    }
                }
                
                if(j==Ny-1){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j-3,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                    }
                }
                
                if(j==Ny-2){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j-3,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                    }
                }
                
                if(j==Ny-3){
                    for(k=0;k<4;k++)
                    {
                        yynew[0+k*7]=Ydata[idx(i,j-3,Nx,Ny,k)];
                        yynew[1+k*7]=Ydata[idx(i,j-2,Nx,Ny,k)];
                        yynew[2+k*7]=Ydata[idx(i,j-1,Nx,Ny,k)];
                        yynew[3+k*7]=Ydata[idx(i,j,Nx,Ny,k)];
                        yynew[4+k*7]=Ydata[idx(i,j+1,Nx,Ny,k)];
                        yynew[5+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                        yynew[6+k*7]=Ydata[idx(i,j+2,Nx,Ny,k)];
                    }
                }
            }

            
            flag = Setlfyegm(Ydata, lfyegm, i, j, Nx, Ny, gama, bcflag);
            if (flag!=0) printf("error in Setlfyegm function \n");
            
            /* compute characteristic variables */
            for (k=0;k<4;k++){
                //printf(" lfxegm : i = %li, j = %li, lfxegm[k][0] = %f, lfxegm[k][1] = %f, lfxegm[k][2] = %f, lfxegm[k][3] = %f\n", i, j, lfxegm[k][0], lfxegm[k][1],lfxegm[k][2],lfxegm[k][3]);
                yynewdata[0+k*7] = lfyegm[k][0]*yynew[0+0*7]+lfyegm[k][1]*yynew[0+1*7]+lfyegm[k][2]*yynew[0+2*7]+lfyegm[k][3]*yynew[0+3*7];
                yynewdata[1+k*7] = lfyegm[k][0]*yynew[1+0*7]+lfyegm[k][1]*yynew[1+1*7]+lfyegm[k][2]*yynew[1+2*7]+lfyegm[k][3]*yynew[1+3*7];
                yynewdata[2+k*7] = lfyegm[k][0]*yynew[2+0*7]+lfyegm[k][1]*yynew[2+1*7]+lfyegm[k][2]*yynew[2+2*7]+lfyegm[k][3]*yynew[2+3*7];
                yynewdata[3+k*7] = lfyegm[k][0]*yynew[3+0*7]+lfyegm[k][1]*yynew[3+1*7]+lfyegm[k][2]*yynew[3+2*7]+lfyegm[k][3]*yynew[3+3*7];
                yynewdata[4+k*7] = lfyegm[k][0]*yynew[4+0*7]+lfyegm[k][1]*yynew[4+1*7]+lfyegm[k][2]*yynew[4+2*7]+lfyegm[k][3]*yynew[4+3*7];
                yynewdata[5+k*7] = lfyegm[k][0]*yynew[5+0*7]+lfyegm[k][1]*yynew[5+1*7]+lfyegm[k][2]*yynew[5+2*7]+lfyegm[k][3]*yynew[5+3*7];
                yynewdata[6+k*7] = lfyegm[k][0]*yynew[6+0*7]+lfyegm[k][1]*yynew[6+1*7]+lfyegm[k][2]*yynew[6+2*7]+lfyegm[k][3]*yynew[6+3*7];
                //yxpdata[idx(i, j, Nx, Ny, k)] = lfxegm[k][0]*yxptdata[idx(i, j, Nx, Ny, 0)]+lfxegm[k][1]*yxptdata[idx(i, j, Nx, Ny, 1)]+lfxegm[k][2]*yxptdata[idx(i, j, Nx, Ny, 2)]+lfxegm[k][3]*yxptdata[idx(i, j, Nx, Ny, 3)];
                //yxndata[idx(i, j, Nx, Ny, k)] = lfxegm[k][0]*yxntdata[idx(i, j, Nx, Ny, 0)]+lfxegm[k][1]*yxntdata[idx(i, j, Nx, Ny, 1)]+lfxegm[k][2]*yxntdata[idx(i, j, Nx, Ny, 2)]+lfxegm[k][3]*yxntdata[idx(i, j, Nx, Ny, 3)];
            }
            
            for (k=0;k<4;k++){
                
                yypdata[0+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[0+k*7]+egymax[k]*yynewdata[0+k*7]);
                yyndata[0+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[0+k*7]-egymax[k]*yynewdata[0+k*7]);
                yypdata[1+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[1+k*7]+egymax[k]*yynewdata[1+k*7]);
                yyndata[1+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[1+k*7]-egymax[k]*yynewdata[1+k*7]);
                yypdata[2+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[2+k*7]+egymax[k]*yynewdata[2+k*7]);
                yyndata[2+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[2+k*7]-egymax[k]*yynewdata[2+k*7]);
                yypdata[3+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[3+k*7]+egymax[k]*yynewdata[3+k*7]);
                yyndata[3+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[3+k*7]-egymax[k]*yynewdata[3+k*7]);
                yypdata[4+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[4+k*7]+egymax[k]*yynewdata[4+k*7]);
                yyndata[4+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[4+k*7]-egymax[k]*yynewdata[4+k*7]);
                yypdata[5+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[5+k*7]+egymax[k]*yynewdata[5+k*7]);
                yyndata[5+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[5+k*7]-egymax[k]*yynewdata[5+k*7]);
                yypdata[6+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[6+k*7]+egymax[k]*yynewdata[6+k*7]);
                yyndata[6+k*7]=0.5*(egvy[idx(i,j,Nx,Ny,k)]*yynewdata[6+k*7]-egymax[k]*yynewdata[6+k*7]);
                
                /*
                 yypdata[0+k*7]=0.5*(egvy[idx(i,j,Nx,Ny+1,k)]*yynewdata[0+k*7]+fabs(egvy[idx(i,j,Nx,Ny+1,k)])*yynewdata[0+k*7]);
                 yyndata[0+k*7]=0.5*(egvy[idx(i,j,Nx,Ny+1,k)]*yynewdata[0+k*7]-fabs(egvy[idx(i,j,Nx,Ny+1,k)])*yynewdata[0+k*7]);
                 yypdata[1+k*7]=0.5*(egvy[idx(i,j,Nx,Ny+1,k)]*yynewdata[1+k*7]+fabs(egvy[idx(i,j,Nx,Ny+1,k)])*yynewdata[1+k*7]);
                 yyndata[1+k*7]=0.5*(egvy[idx(i,j,Nx,Ny+1,k)]*yynewdata[1+k*7]-fabs(egvy[idx(i,j,Nx,Ny+1,k)])*yynewdata[1+k*7]);
                 yypdata[2+k*7]=0.5*(egvy[idx(i,j,Nx,Ny+1,k)]*yynewdata[2+k*7]+fabs(egvy[idx(i,j,Nx,Ny+1,k)])*yynewdata[2+k*7]);
                 yyndata[2+k*7]=0.5*(egvy[idx(i,j,Nx,Ny+1,k)]*yynewdata[2+k*7]-fabs(egvy[idx(i,j,Nx,Ny+1,k)])*yynewdata[2+k*7]);
                 yypdata[3+k*7]=0.5*(egvy[idx(i,j,Nx,Ny+1,k)]*yynewdata[3+k*7]+fabs(egvy[idx(i,j,Nx,Ny+1,k)])*yynewdata[3+k*7]);
                 yyndata[3+k*7]=0.5*(egvy[idx(i,j,Nx,Ny+1,k)]*yynewdata[3+k*7]-fabs(egvy[idx(i,j,Nx,Ny+1,k)])*yynewdata[3+k*7]);
                 yypdata[4+k*7]=0.5*(egvy[idx(i,j,Nx,Ny+1,k)]*yynewdata[4+k*7]+fabs(egvy[idx(i,j,Nx,Ny+1,k)])*yynewdata[4+k*7]);
                 yyndata[4+k*7]=0.5*(egvy[idx(i,j,Nx,Ny+1,k)]*yynewdata[4+k*7]-fabs(egvy[idx(i,j,Nx,Ny+1,k)])*yynewdata[4+k*7]);
                 yypdata[5+k*7]=0.5*(egvy[idx(i,j,Nx,Ny+1,k)]*yynewdata[5+k*7]+fabs(egvy[idx(i,j,Nx,Ny+1,k)])*yynewdata[5+k*7]);
                 yyndata[5+k*7]=0.5*(egvy[idx(i,j,Nx,Ny+1,k)]*yynewdata[5+k*7]-fabs(egvy[idx(i,j,Nx,Ny+1,k)])*yynewdata[5+k*7]);
                 */
            }
            
            /* get indicator of smoothness on y direction */
            flag = SetISY(IS0_py, IS1_py, IS2_py, IS0_ny, IS1_ny, IS2_ny, yypdata, yyndata, i, j, Nx, Ny);
            if (flag!=0) printf("error in SetISY function \n");
            
            /* get weight on y direction for four variables */
            flag = Setalphawy(IS0_py, IS1_py, IS2_py, IS0_ny, IS1_ny, IS2_ny, alpha_0py, alpha_1py, alpha_2py, alpha_0ny, alpha_1ny, alpha_2ny, w0_py, w1_py, w2_py, w0_ny, w1_ny, w2_ny, Epsilon);
            if (flag!=0) printf("error in Setalphawy function \n");
            
            /* compute the positive and negative parts of right interface on y direction */
            flag = SetUy(w0_py, w1_py, w2_py, w0_ny, w1_ny, w2_ny, yypdata, yyndata, u_tpphy, u_tnphy, u_tpnhy, u_tnnhy, i, j, Nx, Ny);
            if (flag!=0) printf("error in SetUy function \n");
            
            /* get the interface values */
            for(k=0;k<4;k++){
                yyupdata[idx(i, j, Nx, Ny, k)]=u_tpphy[k]+u_tnphy[k];
                yydowndata[idx(i, j, Nx, Ny, k)]=u_tpnhy[k]+u_tnnhy[k];
                //printf("yy: i=%li, j=%li, k=%li, u_tpphy[k]=%f, u_tnphy[k]=%f,u_tpnhy[k]=%f,u_tnnhy[k]=%f\n",i,j, k,u_tpphy[k], u_tnphy[k], u_tpnhy[k],u_tnnhy[k]);
            }
            
            flag = Setrhyegm(Ydata, rhyegm, i, j+1, Nx, Ny, gama, bcflag);
            if (flag!=0) printf("error in Setrhyegm function \n");
            for (k=0;k<4;k++){
                
                //printf(" rhyegm : i = %li, j = %li, rhyegm[k][0] = %f, rhyegm[k][1] = %f, rhyegm[k][2] = %f, rhyegm[k][3] = %f\n", i, j, rhyegm[k][0], rhyegm[k][1],rhyegm[k][2],rhyegm[k][3]);
                
                yybackupdata[idx(i, j, Nx, Ny, k)] = rhyegm[k][0]*yyupdata[idx(i, j, Nx, Ny, 0)]+rhyegm[k][1]*yyupdata[idx(i, j, Nx, Ny, 1)]+rhyegm[k][2]*yyupdata[idx(i, j, Nx, Ny, 2)]+rhyegm[k][3]*yyupdata[idx(i, j, Nx, Ny, 3)];
            }
            
            flag = Setrhyegm(Ydata, rhyegm, i, j, Nx, Ny, gama, bcflag);
            if (flag!=0) printf("error in Setrhyegm function \n");
            for (k=0;k<4;k++){
                
                //printf(" rhyegm : i = %li, j = %li, rhyegm[k][0] = %f, rhyegm[k][1] = %f, rhyegm[k][2] = %f, rhyegm[k][3] = %f\n", i, j, rhyegm[k][0], rhyegm[k][1],rhyegm[k][2],rhyegm[k][3]);
                
                yybackdowndata[idx(i, j, Nx, Ny, k)] = rhyegm[k][0]*yydowndata[idx(i, j, Nx, Ny, 0)]+rhyegm[k][1]*yydowndata[idx(i, j, Nx, Ny, 1)]+rhyegm[k][2]*yydowndata[idx(i, j, Nx, Ny, 2)]+rhyegm[k][3]*yydowndata[idx(i, j, Nx, Ny, 3)];
            }

        }
    }
    
    /* fill in the value of tao and J in the whole domain */
    flag = Gettao(yxbackupdata, yxbackdowndata, yybackupdata, yybackdowndata, taoup, taodown, Nx, Ny, gama);
    if (flag!=0) printf("error in Gettao function \n");
    flag = GetCj(yxbackupdata, yxbackdowndata, yybackupdata, yybackdowndata, Cjup, Cjdown, Nx, Ny, gama);
    if (flag!=0) printf("error in GetCj function \n");
    
    /* get derivative both on x and y direction for the whole problem*/
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
            
            //for (k=0; k<4; k++){
            //  if(k==2)
            //        dYdata[idx(i, j, Nx, Ny, k)]=yxdata[idx(i, j, Nx, Ny, k)]+yydata[idx(i, j, Nx, Ny, k)]-0.1;
            // else
	  /*
            printf(" dYdata0 : i = %li, j = %li, yxbackupdata[idx(i, j, Nx, Ny, 1)] = %f, yxbackdowndata[idx(i, j, Nx, Ny, 1)] = %f, yybackupdata[idx(i, j, Nx, Ny, 2)]= %f, yybackdowndata[idx(i, j, Nx, Ny, 2)] = %f\n", i, j, yxbackupdata[idx(i, j, Nx, Ny, 1)], yxbackdowndata[idx(i, j, Nx, Ny, 1)], yybackupdata[idx(i, j, Nx, Ny, 2)], yybackdowndata[idx(i, j, Nx, Ny, 2)]);
                printf(" dYdata1 : i = %li, j = %li, taoup[idx(i, j, Nx, Ny, 0)] = %f, taodown[idx(i, j, Nx, Ny, 0)] = %f, taoup[idx(i, j, Nx, Ny, 2)] = %f, taodown[idx(i, j, Nx, Ny, 2)]] = %f\n", i, j, taoup[idx(i, j, Nx, Ny, 0)], taodown[idx(i, j, Nx, Ny, 0)], taoup[idx(i, j, Nx, Ny, 2)], taodown[idx(i, j, Nx, Ny, 2)]);
		printf(" dYdata2 : i = %li, j = %li, taoup[idx(i, j, Nx, Ny, 1)] = %f, taodown[idx(i, j, Nx, Ny, 1)] = %f, taoup[idx(i, j, Nx, Ny, 3)] = %f, taodown[idx(i, j, Nx, Ny, 3)]] = %f\n", i, j, taoup[idx(i, j, Nx, Ny, 1)], taodown[idx(i, j, Nx, Ny, 1)], taoup[idx(i, j, Nx, Ny, 3)], taodown[idx(i, j, Nx, Ny, 3)]);
                printf(" dYdata3 : i = %li, j = %li, Cjup[idx(i, j, Nx, Ny, 0)] = %f, Cjdown[idx(i, j, Nx, Ny, 0)] = %f, Cjup[idx(i, j, Nx, Ny, 1)] = %f, Cjdown[idx(i, j, Nx, Ny, 1)] = %f\n", i, j, Cjup[idx(i, j, Nx, Ny, 0)], Cjdown[idx(i, j, Nx, Ny, 0)], Cjup[idx(i, j, Nx, Ny, 1)], Cjdown[idx(i, j, Nx, Ny, 1)]);
	  */
            
            dYdata[idx(i, j, Nx, Ny, 0)]=(-1.0/dx)*(yxbackupdata[idx(i, j, Nx, Ny, 1)]-yxbackdowndata[idx(i, j, Nx, Ny, 1)])+(-1.0/dy)*(yybackupdata[idx(i, j, Nx, Ny, 2)]-yybackdowndata[idx(i, j, Nx, Ny, 2)]);
            dYdata[idx(i, j, Nx, Ny, 1)]=(-1.0/dx)*(taoup[idx(i, j, Nx, Ny, 0)]-taodown[idx(i, j, Nx, Ny, 0)])+(-1.0/dy)*(taoup[idx(i, j, Nx, Ny, 2)]-taodown[idx(i, j, Nx, Ny, 2)]);
            dYdata[idx(i, j, Nx, Ny, 2)]=(-1.0/dx)*(taoup[idx(i, j, Nx, Ny, 1)]-taodown[idx(i, j, Nx, Ny, 1)])+(-1.0/dy)*(taoup[idx(i, j, Nx, Ny, 3)]-taodown[idx(i, j, Nx, Ny, 3)]);
            dYdata[idx(i, j, Nx, Ny, 3)]=(-1.0/dx)*(Cjup[idx(i, j, Nx, Ny, 0)]-Cjdown[idx(i, j, Nx, Ny, 0)])+(-1.0/dy)*(Cjup[idx(i, j, Nx, Ny, 1)]-Cjdown[idx(i, j, Nx, Ny, 1)]);
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
    delete []yxbackupdata;
    delete []yybackupdata;
    delete []yxbackdowndata;
    delete []yybackdowndata;
    delete []taoup;
    delete []Cjup;
    delete []taodown;
    delete []Cjdown;
    delete []egvx;
    delete []egvy;
    //delete []egvxmax;
    //delete []egvymax;
    delete []egxmax;
    delete []egymax;
    /*
     delete []yxptdata;
     delete []yxntdata;
     delete []yyptdata;
     delete []yyntdata;
     */
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
    
    return 0;
}

static int Getmaxegvx(realtype *egvx, realtype *egxmax, long int Nx, long int Ny)
{
    long int i, j, k;
    
    for(k=0;k<4;k++){
        egxmax[k]=0.0;
    }
    
    for (j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            for(k=0;k<4;k++){
                if(fabs(egvx[idx(i,j,Nx,Ny,k)])>egxmax[k])
                {
                    egxmax[k]=fabs(egvx[idx(i,j,Nx,Ny,k)]);
                }
            }
        }
    }
    return 0;
}

static int Getmaxegvy(realtype *egvy, realtype *egymax, long int Nx, long int Ny)
{
    long int i, j, k;
    
    for(k=0;k<4;k++){
        egymax[k]=0.0;
    }
    
    for (i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            for(k=0;k<4;k++){
                if(fabs(egvy[idx(i,j,Nx,Ny,k)])>egymax[k])
                {
                    egymax[k]=fabs(egvy[idx(i,j,Nx,Ny,k)]);
                }
            }
        }
    }
    return 0;
}

/* split f into f_positive part and f_negative part for x component*/
static int Splitfluxesx(realtype *Ydata, realtype *fp, realtype *fn, long int i, long int j, long int Nx, long int Ny, realtype gama)
{
    /* declaration */
    realtype rou, vx, vy, p, a, h;
    int m, n;
    realtype eps = 1.0e-3;
    
    realtype **lfxegm = new realtype *[4];
    for (m=0;m<4;m++){
        lfxegm[m] = new realtype [4];
    }
    realtype **rhxegm = new realtype *[4];
    for (m=0;m<4;m++){
        rhxegm[m] = new realtype [4];
    }
    realtype **tempx1 = new realtype *[4];
    for (m=0;m<4;m++){
        tempx1[m] = new realtype [4];
    }
    realtype **tempx2 = new realtype *[4];
    for (m=0;m<4;m++){
        tempx2[m] = new realtype [4];
    }
    realtype **egvxp = new realtype *[4];
    for (m=0;m<4;m++){
        egvxp[m] = new realtype [4];
    }
    realtype **egvxn = new realtype *[4];
    for (m=0;m<4;m++){
        egvxn[m] = new realtype [4];
    }
    
    vx = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
    vy = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
    p = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vx*vx+vy*vy);
    a = sqrt(gama*p/Ydata[idx(i, j, Nx, Ny, 0)]);
    h = (Ydata[idx(i, j, Nx, Ny, 3)]+p)/Ydata[idx(i, j, Nx, Ny, 0)];
    
    for (m=0;m<4;m++)
    {
        for (n=0;n<4;n++)
        {
            egvxp[m][n]=0.0;
        }
    }
    
    for (m=0;m<4;m++)
    {
        for (n=0;n<4;n++)
        {
            egvxn[m][n]=0.0;
        }
    }
    /* set eigenvalues */
    egvxp[0][0]=(vx-a+sqrt((vx-a)*(vx-a)+eps))/2.0;
    egvxp[1][1]=(vx+sqrt(vx*vx+eps))/2.0;
    egvxp[2][2]=(vx+sqrt(vx*vx+eps))/2.0;
    egvxp[3][3]=(vx+a+sqrt((vx+a)*(vx+a)+eps))/2.0;
    egvxn[0][0]=(vx-a-sqrt((vx-a)*(vx-a)+eps))/2.0;
    egvxn[1][1]=(vx-sqrt(vx*vx+eps))/2.0;
    egvxn[2][2]=(vx-sqrt(vx*vx+eps))/2.0;
    egvxn[3][3]=(vx+a-sqrt((vx+a)*(vx+a)+eps))/2.0;
    
    /* left eigenvectors derivated from the book */
    lfxegm[0][0] = ((gama-1.0)/(2.0*a*a))*(h+(a/(gama-1.0))*(vx-a));
    lfxegm[1][0] = ((gama-1.0)/(2.0*a*a))*(-2.0*h+(4.0/(gama-1.0))*a*a);
    lfxegm[2][0] = ((gama-1.0)/(2.0*a*a))*(-2.0*vy*a*a/(gama-1.0));
    lfxegm[3][0] = ((gama-1.0)/(2.0*a*a))*(h-(a/(gama-1.0))*(vx+a));
    lfxegm[0][1] = ((gama-1.0)/(2.0*a*a))*(-(vx+a/(gama-1.0)));
    lfxegm[1][1] = ((gama-1.0)/(2.0*a*a))*(2.0*vx);
    lfxegm[2][1] = 0.0;
    lfxegm[3][1] = ((gama-1.0)/(2.0*a*a))*(-vx+a/(gama-1.0));
    lfxegm[0][2] = ((gama-1.0)/(2.0*a*a))*(-vy);
    lfxegm[1][2] = ((gama-1.0)/(2.0*a*a))*(2.0*vy);
    lfxegm[2][2] = ((gama-1.0)/(2.0*a*a))*(2.0*a*a/(gama-1.0));
    lfxegm[3][2] = ((gama-1.0)/(2.0*a*a))*(-vy);
    lfxegm[0][3] = ((gama-1.0)/(2.0*a*a))*1.0;
    lfxegm[1][3] = ((gama-1.0)/(2.0*a*a))*(-2.0);
    lfxegm[2][3] = 0.0;
    lfxegm[3][3] = ((gama-1.0)/(2.0*a*a))*1.0;
    /* left eigenvectors derivated by myself */
    /*
     lfxegm[0][0] = ((gama-1.0)*h)/(2.0*a*a)+vx/(2.0*a)-0.5;
     lfxegm[1][0] = ((1.0-gama)*h)/(a*a)+2.0;
     lfxegm[3][0] = ((gama-1.0)*h)/(2.0*a*a)-vx/(2.0*a)-0.5;
     lfxegm[2][0] = ((1.0-gama)*h*vy)/(a*a)+vy;
     lfxegm[0][1] = ((1.0-gama)*vx)/(2*a*a)-1.0/(2.0*a);
     lfxegm[1][1] = (gama-1.0)*vx/(a*a);
     lfxegm[3][1] = ((1.0-gama)*vx)/(2*a*a)+1.0/(2.0*a);
     lfxegm[2][1] = (gama-1.0)*vx*vy/(a*a);
     lfxegm[0][2] = ((1.0-gama)*vy)/(2*a*a);
     lfxegm[1][2] = (gama-1.0)*vy/(a*a);
     lfxegm[3][2] = ((1.0-gama)*vy)/(2*a*a);
     lfxegm[2][2] = (gama-1.0)*vy*vy/(a*a)+1.0;
     lfxegm[0][3] = (gama-1.0)/(2*a*a);
     lfxegm[1][3] = (1.0-gama)/(a*a);
     lfxegm[3][3] = (gama-1.0)/(2*a*a);
     lfxegm[2][3] = (1.0-gama)*vy/(a*a);
     */
    /* right eigenvectors derivated from the book */
    rhxegm[0][0] = 1.0;
    rhxegm[1][0] = vx-a;
    rhxegm[2][0] = vy;
    rhxegm[3][0] = h-a*vx;
    rhxegm[0][1] = 1.0;
    rhxegm[1][1] = vx;
    rhxegm[2][1] = vy;
    rhxegm[3][1] = 0.5*(vx*vx+vy*vy);
    rhxegm[0][2] = 0.0;
    rhxegm[1][2] = 0.0;
    rhxegm[2][2] = 1.0;
    rhxegm[3][2] = vy;
    rhxegm[0][3] = 1.0;
    rhxegm[1][3] = vx+a;
    rhxegm[2][3] = vy;
    rhxegm[3][3] = h+a*vx;
    /* right eigenvectors derivated by myself */
    /*
     rhxegm[0][0] = 1.0;
     rhxegm[1][0] = vx-a;
     rhxegm[2][0] = vy;
     rhxegm[3][0] = vx*vx+vy*vy-h-a*vx+(2*a*a)/(gama-1.0);
     rhxegm[0][1] = 1.0;
     rhxegm[1][1] = vx;
     rhxegm[2][1] = 0.0;
     rhxegm[3][1] = vx*vx+a*a/(gama-1.0)-h;
     rhxegm[0][2] = 0.0;
     rhxegm[1][2] = 0.0;
     rhxegm[2][2] = 1.0;
     rhxegm[3][2] = vy;
     rhxegm[0][3] = 1.0;
     rhxegm[1][3] = vx+a;
     rhxegm[2][3] = vy;
     rhxegm[3][3] = vx*vx+vy*vy-h+a*vx+(2*a*a)/(gama-1.0);
     */
    /* get fp */
    for (m=0;m<4;m++)
    {
        for (n=0;n<4;n++)
        {
            tempx1[m][n]=rhxegm[m][0]*egvxp[0][n]+rhxegm[m][1]*egvxp[1][n]+rhxegm[m][2]*egvxp[2][n]+rhxegm[m][3]*egvxp[3][n];
        }
    }
    
    for (m=0;m<4;m++)
    {
        for (n=0;n<4;n++)
        {
            tempx2[m][n]=tempx1[m][0]*lfxegm[0][n]+tempx1[m][1]*lfxegm[1][n]+tempx1[m][2]*lfxegm[2][n]+tempx1[m][3]*lfxegm[3][n];
        }
    }
    
    for (m=0;m<4;m++)
    {
        fp[idx(i, j, Nx, Ny, m)]=tempx2[m][0]*Ydata[idx(i, j, Nx, Ny, 0)]+tempx2[m][1]*Ydata[idx(i, j, Nx, Ny, 1)]+tempx2[m][2]*Ydata[idx(i, j, Nx, Ny, 2)]+tempx2[m][3]*Ydata[idx(i, j, Nx, Ny, 3)];
    }
    
    /* get fn */
    for (m=0;m<4;m++)
    {
        for (n=0;n<4;n++)
        {
            tempx1[m][n]=rhxegm[m][0]*egvxn[0][n]+rhxegm[m][1]*egvxn[1][n]+rhxegm[m][2]*egvxn[2][n]+rhxegm[m][3]*egvxn[3][n];
        }
    }
    
    for (m=0;m<4;m++)
    {
        for (n=0;n<4;n++)
        {
            tempx2[m][n]=tempx1[m][0]*lfxegm[0][n]+tempx1[m][1]*lfxegm[1][n]+tempx1[m][2]*lfxegm[2][n]+tempx1[m][3]*lfxegm[3][n];
        }
    }
    
    for (m=0;m<4;m++)
    {
        fn[idx(i, j, Nx, Ny, m)]=tempx2[m][0]*Ydata[idx(i, j, Nx, Ny, 0)]+tempx2[m][1]*Ydata[idx(i, j, Nx, Ny, 1)]+tempx2[m][2]*Ydata[idx(i, j, Nx, Ny, 2)]+tempx2[m][3]*Ydata[idx(i, j, Nx, Ny, 3)];
    }
    
    /* delete arrays */
    for (m=0;m<4;m++){
        delete[] lfxegm[m];
    }
    delete []lfxegm;
    
    for (m=0;m<4;m++){
        delete[] rhxegm[m];
    }
    delete []rhxegm;
    
    for (m=0;m<4;m++){
        delete[] tempx1[m];
    }
    delete []tempx1;
    
    for (m=0;m<4;m++){
        delete[] tempx2[m];
    }
    delete []tempx2;
    
    for (m=0;m<4;m++){
        delete[] egvxp[m];
    }
    delete []egvxp;
    
    for (m=0;m<4;m++){
        delete[] egvxn[m];
    }
    delete []egvxn;
    
    return 0;
}

/* split g into g_positive part and g_negative part for y component*/
static int Splitfluxesy(realtype *Ydata, realtype *gp, realtype *gn, long int i, long int j, long int Nx, long int Ny, realtype gama)
{
    /* declaration */
    realtype rou, vx, vy, p, a, h;
    int m, n;
    realtype eps = 1.0e-3;
    
    realtype **lfyegm = new realtype *[4];
    for (m=0;m<4;m++){
        lfyegm[m] = new realtype [4];
    }
    realtype **rhyegm = new realtype *[4];
    for (m=0;m<4;m++){
        rhyegm[m] = new realtype [4];
    }
    realtype **tempy1 = new realtype *[4];
    for (m=0;m<4;m++){
        tempy1[m] = new realtype [4];
    }
    realtype **tempy2 = new realtype *[4];
    for (m=0;m<4;m++){
        tempy2[m] = new realtype [4];
    }
    realtype **egvyp = new realtype *[4];
    for (m=0;m<4;m++){
        egvyp[m] = new realtype [4];
    }
    realtype **egvyn = new realtype *[4];
    for (m=0;m<4;m++){
        egvyn[m] = new realtype [4];
    }
    
    vx = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
    vy = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
    p = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vx*vx+vy*vy);
    a = sqrt(gama*p/Ydata[idx(i, j, Nx, Ny, 0)]);
    h = (Ydata[idx(i, j, Nx, Ny, 3)]+p)/Ydata[idx(i, j, Nx, Ny, 0)];
    
    for (m=0;m<4;m++)
    {
        for (n=0;n<4;n++)
        {
            egvyp[m][n]=0.0;
        }
    }
    
    for (m=0;m<4;m++)
    {
        for (n=0;n<4;n++)
        {
            egvyn[m][n]=0.0;
        }
    }
    /* eigenvalues */
    egvyp[0][0]=(vy-a+sqrt((vy-a)*(vy-a)+eps))/2.0;
    egvyp[1][1]=(vy+sqrt(vy*vy+eps))/2.0;
    egvyp[2][2]=(vy+sqrt(vy*vy+eps))/2.0;
    egvyp[3][3]=(vy+a+sqrt((vy+a)*(vy+a)+eps))/2.0;
    egvyn[0][0]=(vy-a-sqrt((vy-a)*(vy-a)+eps))/2.0;
    egvyn[1][1]=(vy-sqrt(vy*vy+eps))/2.0;
    egvyn[2][2]=(vy-sqrt(vy*vy+eps))/2.0;
    egvyn[3][3]=(vy+a-sqrt((vy+a)*(vy+a)+eps))/2.0;
    
    /* left eigenvectors derivated by myself */
    lfyegm[0][0] = ((gama-1.0)*h)/(2.0*a*a)+vy/(2.0*a)-0.5;
    lfyegm[1][0] = ((1.0-gama)*h)/(a*a)+2.0;
    lfyegm[2][0] = ((1.0-gama)*h*vx)/(a*a)+vx;
    lfyegm[3][0] = ((gama-1.0)*h)/(2.0*a*a)-vy/(2.0*a)-0.5;
    lfyegm[0][1] = ((1.0-gama)*vx)/(2*a*a);
    lfyegm[1][1] = (gama-1.0)*vx/(a*a);
    lfyegm[2][1] = (gama-1.0)*vx*vx/(a*a)+1.0;
    lfyegm[3][1] = ((1.0-gama)*vx)/(2*a*a);
    lfyegm[0][2] = ((1.0-gama)*vy)/(2*a*a)-1.0/(2.0*a);
    lfyegm[1][2] = (gama-1.0)*vy/(a*a);
    lfyegm[2][2] = (gama-1.0)*vx*vy/(a*a);
    lfyegm[3][2] = ((1.0-gama)*vy)/(2*a*a)+1.0/(2.0*a);
    lfyegm[0][3] = (gama-1.0)/(2*a*a);
    lfyegm[1][3] = (1.0-gama)/(a*a);
    lfyegm[2][3] = (1.0-gama)*vx/(a*a);
    lfyegm[3][3] = (gama-1.0)/(2*a*a);
    
    /* right eigenvectors derivated by myself */
    rhyegm[0][0] = 1.0;
    rhyegm[1][0] = vx;
    rhyegm[2][0] = vy-a;
    rhyegm[3][0] = vx*vx+vy*vy-h-vy*a+(2.0*a*a)/(gama-1.0);
    rhyegm[0][1] = 1.0;
    rhyegm[1][1] = 0.0;
    rhyegm[2][1] = vy;
    rhyegm[3][1] = vy*vy+a*a/(gama-1.0)-h;
    rhyegm[0][2] = 0.0;
    rhyegm[1][2] = 1.0;
    rhyegm[2][2] = 0.0;
    rhyegm[3][2] = vx;
    rhyegm[0][3] = 1.0;
    rhyegm[1][3] = vx;
    rhyegm[2][3] = vy+a;
    rhyegm[3][3] = vx*vx+vy*vy-h+vy*a+(2.0*a*a)/(gama-1.0);
    
    /* get gp */
    for (m=0;m<4;m++)
    {
        for (n=0;n<4;n++)
        {
            tempy1[m][n]=rhyegm[m][0]*egvyp[0][n]+rhyegm[m][1]*egvyp[1][n]+rhyegm[m][2]*egvyp[2][n]+rhyegm[m][3]*egvyp[3][n];
        }
    }
    
    for (m=0;m<4;m++)
    {
        for (n=0;n<4;n++)
        {
            tempy2[m][n]=tempy1[m][0]*lfyegm[0][n]+tempy1[m][1]*lfyegm[1][n]+tempy1[m][2]*lfyegm[2][n]+tempy1[m][3]*lfyegm[3][n];
        }
    }
    
    for (m=0;m<4;m++)
    {
        gp[idx(i, j, Nx, Ny, m)]=tempy2[m][0]*Ydata[idx(i, j, Nx, Ny, 0)]+tempy2[m][1]*Ydata[idx(i, j, Nx, Ny, 1)]+tempy2[m][2]*Ydata[idx(i, j, Nx, Ny, 2)]+tempy2[m][3]*Ydata[idx(i, j, Nx, Ny, 3)];
    }
    
    /* get gn */
    for (m=0;m<4;m++)
    {
        for (n=0;n<4;n++)
        {
            tempy1[m][n]=rhyegm[m][0]*egvyn[0][n]+rhyegm[m][1]*egvyn[1][n]+rhyegm[m][2]*egvyn[2][n]+rhyegm[m][3]*egvyn[3][n];
        }
    }
    
    for (m=0;m<4;m++)
    {
        for (n=0;n<4;n++)
        {
            tempy2[m][n]=tempy1[m][0]*lfyegm[0][n]+tempy1[m][1]*lfyegm[1][n]+tempy1[m][2]*lfyegm[2][n]+tempy1[m][3]*lfyegm[3][n];
        }
    }
    
    for (m=0;m<4;m++)
    {
        gn[idx(i, j, Nx, Ny, m)]=tempy2[m][0]*Ydata[idx(i, j, Nx, Ny, 0)]+tempy2[m][1]*Ydata[idx(i, j, Nx, Ny, 1)]+tempy2[m][2]*Ydata[idx(i, j, Nx, Ny, 2)]+tempy2[m][3]*Ydata[idx(i, j, Nx, Ny, 3)];
    }
    
    /* delete arrays */
    for (m=0;m<4;m++){
        delete[] lfyegm[m];
    }
    delete []lfyegm;
    
    for (m=0;m<4;m++){
        delete[] rhyegm[m];
    }
    delete []rhyegm;
    
    for (m=0;m<4;m++){
        delete[] tempy1[m];
    }
    delete []tempy1;
    
    for (m=0;m<4;m++){
        delete[] tempy2[m];
    }
    delete []tempy2;
    
    for (m=0;m<4;m++){
        delete[] egvyp[m];
    }
    delete []egvyp;
    
    for (m=0;m<4;m++){
        delete[] egvyn[m];
    }
    delete []egvyn;
    
    return 0;
}

/* fill in left eigenvectors for x component*/
static int Getegvx(realtype *Ydata, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype *egvx, int flag)
{
    /* declaration */
    realtype rou, vx, vy, p, a, h;
    
        vx = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
        vy = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
        p = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vx*vx+vy*vy);
        a = sqrt(gama*p/Ydata[idx(i, j, Nx, Ny, 0)]);
        h = (Ydata[idx(i, j, Nx, Ny, 3)]+p)/Ydata[idx(i, j, Nx, Ny, 0)];
    
    /* eigenvalues */
    egvx[idx(i,j,Nx,Ny,0)]=vx-a;
    egvx[idx(i,j,Nx,Ny,1)]=vx;
    egvx[idx(i,j,Nx,Ny,2)]=vx;
    egvx[idx(i,j,Nx,Ny,3)]=vx+a;
    
    return 0;
}

/* fill in left eigenvectors for x component*/
static int Getegvy(realtype *Ydata, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype *egvy, int flag)
{
    /* declaration */
    realtype rou, vx, vy, p, a, h;
    
        vx = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
        vy = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
        p = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vx*vx+vy*vy);
        a = sqrt(gama*p/Ydata[idx(i, j, Nx, Ny, 0)]);
        h = (Ydata[idx(i, j, Nx, Ny, 3)]+p)/Ydata[idx(i, j, Nx, Ny, 0)];
    
    /* eigenvalues */
    /* get eigenvalues */
    egvy[idx(i,j,Nx,Ny,0)]=vy-a;
    egvy[idx(i,j,Nx,Ny,1)]=vy;
    egvy[idx(i,j,Nx,Ny,2)]=vy;
    egvy[idx(i,j,Nx,Ny,3)]=vy+a;
    
    return 0;
}

/* fill in left eigenvectors for x component*/
static int Setlfxegm(realtype *Ydata, realtype **lfxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag)
{
    /* declaration */
    realtype rou, vx, vy, p, a, h;
    vx = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
    vy = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
    p = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vx*vx+vy*vy);
    a = sqrt(gama*p/Ydata[idx(i, j, Nx, Ny, 0)]);
    h = (Ydata[idx(i, j, Nx, Ny, 3)]+p)/Ydata[idx(i, j, Nx, Ny, 0)];
    
    /* left eigenvectors derivated by myself */
    /*
     lfxegm[0][0] = ((gama-1.0)*h)/(2.0*a*a)+vx/(2.0*a)-0.5;
     lfxegm[1][0] = ((1.0-gama)*h)/(a*a)+2.0;
     lfxegm[2][0] = ((1.0-gama)*h*vy)/(a*a)+vy;
     lfxegm[3][0] = ((gama-1.0)*h)/(2.0*a*a)-vx/(2.0*a)-0.5;
     lfxegm[0][1] = ((1.0-gama)*vx)/(2*a*a)-1.0/(2.0*a);
     lfxegm[1][1] = (gama-1.0)*vx/(a*a);
     lfxegm[2][1] = (gama-1.0)*vx*vy/(a*a);
     lfxegm[3][1] = ((1.0-gama)*vx)/(2*a*a)+1.0/(2.0*a);
     lfxegm[0][2] = ((1.0-gama)*vy)/(2*a*a);
     lfxegm[1][2] = (gama-1.0)*vy/(a*a);
     lfxegm[2][2] = (gama-1.0)*vy*vy/(a*a)+1.0;
     lfxegm[3][2] = ((1.0-gama)*vy)/(2*a*a);
     lfxegm[0][3] = (gama-1.0)/(2*a*a);
     lfxegm[1][3] = (1.0-gama)/(a*a);
     lfxegm[2][3] = (1.0-gama)*vy/(a*a);
     lfxegm[3][3] = (gama-1.0)/(2*a*a);
     */
    /* left eigenvectors derivated from the book */
    lfxegm[0][0] = ((gama-1.0)/(2.0*a*a))*(h+(a/(gama-1.0))*(vx-a));
    lfxegm[1][0] = ((gama-1.0)/(2.0*a*a))*(-2.0*h+(4.0/(gama-1.0))*a*a);
    lfxegm[2][0] = ((gama-1.0)/(2.0*a*a))*(-2.0*vy*a*a/(gama-1.0));
    lfxegm[3][0] = ((gama-1.0)/(2.0*a*a))*(h-(a/(gama-1.0))*(vx+a));
    lfxegm[0][1] = ((gama-1.0)/(2.0*a*a))*(-(vx+a/(gama-1.0)));
    lfxegm[1][1] = ((gama-1.0)/(2.0*a*a))*(2.0*vx);
    lfxegm[2][1] = 0.0;
    lfxegm[3][1] = ((gama-1.0)/(2.0*a*a))*(-vx+a/(gama-1.0));
    lfxegm[0][2] = ((gama-1.0)/(2.0*a*a))*(-vy);
    lfxegm[1][2] = ((gama-1.0)/(2.0*a*a))*(2.0*vy);
    lfxegm[2][2] = ((gama-1.0)/(2.0*a*a))*(2.0*a*a/(gama-1.0));
    lfxegm[3][2] = ((gama-1.0)/(2.0*a*a))*(-vy);
    lfxegm[0][3] = ((gama-1.0)/(2.0*a*a))*1.0;
    lfxegm[1][3] = ((gama-1.0)/(2.0*a*a))*(-2.0);
    lfxegm[2][3] = 0.0;
    lfxegm[3][3] = ((gama-1.0)/(2.0*a*a))*1.0;
    
    //printf("lfxegm : i=%li, j=%li, lfxegm[0][0]=%f, lfxegm[1][0]=%f, lfxegm[2][0]=%f, lfxegm[3][0]=%f, lfxegm[0][1]=%f, lfxegm[1][1]=%f, lfxegm[2][1]=%f, lfxegm[3][1]=%f, lfxegm[0][2]=%f, lfxegm[1][2]=%f, lfxegm[2][2]=%f, lfxegm[3][2]=%f, lfxegm[0][3]=%f, lfxegm[1][3]=%f, lfxegm[2][3]=%f, lfxegm[3][3]=%f\n", i, j, lfxegm[0][0], lfxegm[1][0], lfxegm[2][0], lfxegm[3][0], lfxegm[0][1], lfxegm[1][1], lfxegm[2][1], lfxegm[3][1], lfxegm[0][2], lfxegm[1][2], lfxegm[2][2], lfxegm[3][2], lfxegm[0][3], lfxegm[1][3], lfxegm[2][3], lfxegm[3][3]);
    
    return 0;
}

/* fill in right eigenvectors for x component*/
static int Setrhxegm(realtype *Ydata, realtype **rhxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag)
{
    /* declaration */
    realtype rou, vx, vy, p, a, h;
    /* without considering the boundary conditions */
    if (i!=Nx&&i!=0){
        vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i-1, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]);
        vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i-1, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]);
        p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i-1, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
        h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i-1, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]));
        rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]);
    }
    /* periodic boundary condition */
    if (i==Nx){
        if (flag == 0)
        {
            vx = (Ydata[idx(i-1, j, Nx, Ny, 1)]+Ydata[idx(0, j, Nx, Ny, 1)])/(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)]);
            vy = (Ydata[idx(i-1, j, Nx, Ny, 2)]+Ydata[idx(0, j, Nx, Ny, 2)])/(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)]);
            p = (gama-1.0)*0.5*(Ydata[idx(i-1, j, Nx, Ny, 3)]+Ydata[idx(0, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
            h = (0.5*(Ydata[idx(i-1, j, Nx, Ny, 3)]+Ydata[idx(0, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)]));
            rou = 0.5*(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(0, j, Nx, Ny, 0)]);
        }
        /* reflecting boundary condition */
        if (flag == 1)
        {
            vx = (-Ydata[idx(i-1, j, Nx, Ny, 1)]+Ydata[idx(i-1, j, Nx, Ny, 1)])/(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]);
            vy = (Ydata[idx(i-1, j, Nx, Ny, 2)]+Ydata[idx(i-1, j, Nx, Ny, 2)])/(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]);
            p = (gama-1.0)*0.5*(Ydata[idx(i-1, j, Nx, Ny, 3)]+Ydata[idx(i-1, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
            h = (0.5*(Ydata[idx(i-1, j, Nx, Ny, 3)]+Ydata[idx(i-1, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]));
            rou = 0.5*(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]);
        }
        /* nature or transmissive boundary condition */
        if (flag == 2)
        {
            vx = (Ydata[idx(i-1, j, Nx, Ny, 1)]+Ydata[idx(i-1, j, Nx, Ny, 1)])/(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]);
            vy = (Ydata[idx(i-1, j, Nx, Ny, 2)]+Ydata[idx(i-1, j, Nx, Ny, 2)])/(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]);
            p = (gama-1.0)*0.5*(Ydata[idx(i-1, j, Nx, Ny, 3)]+Ydata[idx(i-1, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
            h = (0.5*(Ydata[idx(i-1, j, Nx, Ny, 3)]+Ydata[idx(i-1, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]));
            rou = 0.5*(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]);
        }
    }
    
    if (i==0){
        if (flag == 0)
        {
            vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(Nx-1, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(Nx-1, j, Nx, Ny, 0)]);
            vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(Nx-1, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(Nx-1, j, Nx, Ny, 0)]);
            p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(Nx-1, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(Nx-1, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
            h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(Nx-1, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(Nx-1, j, Nx, Ny, 0)]));
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(Nx-1, j, Nx, Ny, 0)]);
        }
        /* reflecting boundary condition */
        if (flag == 1)
        {
            vx = (-Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
            vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
            p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
            h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]));
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
        }
        /* nature or transmissive boundary condition */
        if (flag == 2)
        {
            vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
            vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
            p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
            h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]));
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    a = sqrt(gama*p/rou);
    /* right eigenvectors derivated by myself */
    /*
     rhxegm[0][0] = 1.0;
     rhxegm[1][0] = vx-a;
     rhxegm[2][0] = vy;
     rhxegm[3][0] = vx*vx+vy*vy-h-a*vx+(2*a*a)/(gama-1.0);
     rhxegm[0][1] = 1.0;
     rhxegm[1][1] = vx;
     rhxegm[2][1] = 0.0;
     rhxegm[3][1] = vx*vx+a*a/(gama-1.0)-h;
     rhxegm[0][2] = 0.0;
     rhxegm[1][2] = 0.0;
     rhxegm[2][2] = 1.0;
     rhxegm[3][2] = vy;
     rhxegm[0][3] = 1.0;
     rhxegm[1][3] = vx+a;
     rhxegm[2][3] = vy;
     rhxegm[3][3] = vx*vx+vy*vy-h+a*vx+(2*a*a)/(gama-1.0);
     
     */
    /* right eigenvectors derivated from the book */
    rhxegm[0][0] = 1.0;
    rhxegm[1][0] = vx-a;
    rhxegm[2][0] = vy;
    rhxegm[3][0] = h-a*vx;
    rhxegm[0][1] = 1.0;
    rhxegm[1][1] = vx;
    rhxegm[2][1] = vy;
    rhxegm[3][1] = 0.5*(vx*vx+vy*vy);
    rhxegm[0][2] = 0.0;
    rhxegm[1][2] = 0.0;
    rhxegm[2][2] = 1.0;
    rhxegm[3][2] = vy;
    rhxegm[0][3] = 1.0;
    rhxegm[1][3] = vx+a;
    rhxegm[2][3] = vy;
    rhxegm[3][3] = h+a*vx;
    
    return 0;
}

/* fill in left eigenvectors for y component*/
static int Setlfyegm(realtype *Ydata, realtype **lfyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag)
{
    /* declaration */
    realtype rou, vx, vy, p, a, h;

     vx = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
     vy = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
     p = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vx*vx+vy*vy);
     a = sqrt(gama*p/Ydata[idx(i, j, Nx, Ny, 0)]);
     h = (Ydata[idx(i, j, Nx, Ny, 3)]+p)/Ydata[idx(i, j, Nx, Ny, 0)];
    
    /* left eigenvectors derivated by myself */
    /*
     lfyegm[0][0] = ((gama-1.0)*h)/(2.0*a*a)+vy/(2.0*a)-0.5;
     lfyegm[1][0] = ((1.0-gama)*h)/(a*a)+2.0;
     lfyegm[2][0] = ((1.0-gama)*h*vx)/(a*a)+vx;
     lfyegm[3][0] = ((gama-1.0)*h)/(2.0*a*a)-vy/(2.0*a)-0.5;
     lfyegm[0][1] = ((1.0-gama)*vx)/(2*a*a);
     lfyegm[1][1] = (gama-1.0)*vx/(a*a);
     lfyegm[2][1] = (gama-1.0)*vx*vx/(a*a)+1.0;
     lfyegm[3][1] = ((1.0-gama)*vx)/(2*a*a);
     lfyegm[0][2] = ((1.0-gama)*vy)/(2*a*a)-1.0/(2.0*a);
     lfyegm[1][2] = (gama-1.0)*vy/(a*a);
     lfyegm[2][2] = (gama-1.0)*vx*vy/(a*a);
     lfyegm[3][2] = ((1.0-gama)*vy)/(2*a*a)+1.0/(2.0*a);
     lfyegm[0][3] = (gama-1.0)/(2*a*a);
     lfyegm[1][3] = (1.0-gama)/(a*a);
     lfyegm[2][3] = (1.0-gama)*vx/(a*a);
     lfyegm[3][3] = (gama-1.0)/(2*a*a);
     */
    
    /* left eigenvectors derivated from the book */
    lfyegm[0][0] = ((gama-1.0)/(2.0*a*a))*(h+(a/(gama-1.0))*(vy-a));
    lfyegm[1][0] = ((gama-1.0)/(2.0*a*a))*(-2.0*h+(4.0/(gama-1.0))*a*a);
    lfyegm[2][0] = ((gama-1.0)/(2.0*a*a))*(-2.0*vx*a*a/(gama-1.0));
    lfyegm[3][0] = ((gama-1.0)/(2.0*a*a))*(h-(a/(gama-1.0))*(vy+a));
    lfyegm[0][1] = ((gama-1.0)/(2.0*a*a))*(-vx);
    lfyegm[1][1] = ((gama-1.0)/(2.0*a*a))*(2.0*vx);
    lfyegm[2][1] = ((gama-1.0)/(2.0*a*a))*(2.0*a*a/(gama-1.0));
    lfyegm[3][1] = ((gama-1.0)/(2.0*a*a))*(-vx);
    lfyegm[0][2] = ((gama-1.0)/(2.0*a*a))*(-(vy+a/(gama-1.0)));
    lfyegm[1][2] = ((gama-1.0)/(2.0*a*a))*(2.0*vy);
    lfyegm[2][2] = 0.0;
    lfyegm[3][2] = ((gama-1.0)/(2.0*a*a))*(-vy+a/(gama-1.0));
    lfyegm[0][3] = ((gama-1.0)/(2.0*a*a))*1.0;
    lfyegm[1][3] = ((gama-1.0)/(2.0*a*a))*(-2.0);
    lfyegm[2][3] = 0.0;
    lfyegm[3][3] = ((gama-1.0)/(2.0*a*a))*1.0;
    
    return 0;
}

/* fill in right eigenvectors for y component*/
static int Setrhyegm(realtype *Ydata, realtype **rhyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag)
{
    /* declaration */
    realtype rou, vx, vy, p, a, h;
    /* without considering boundary conditions */
    if (j!=Ny&&j!=0){
        vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i, j-1, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]);
        vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i, j-1, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]);
        p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j-1, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
        h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j-1, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]));
        rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]);
    }
    /* periodic boundary condition */
    if (j==Ny){
        if (flag == 0)
        {
            vx = (Ydata[idx(i, j-1, Nx, Ny, 1)]+Ydata[idx(i, 0, Nx, Ny, 1)])/(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, 0, Nx, Ny, 0)]);
            vy = (Ydata[idx(i, j-1, Nx, Ny, 2)]+Ydata[idx(i, 0, Nx, Ny, 2)])/(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, 0, Nx, Ny, 0)]);
            p = (gama-1.0)*0.5*(Ydata[idx(i, j-1, Nx, Ny, 3)]+Ydata[idx(i, 0, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, 0, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
            h = (0.5*(Ydata[idx(i, j-1, Nx, Ny, 3)]+Ydata[idx(i, 0, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, 0, Nx, Ny, 0)]));
            rou = 0.5*(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, 0, Nx, Ny, 0)]);
        }
        /* reflecting boundary condition */
        if (flag == 1)
        {
            vx = (Ydata[idx(i, j-1, Nx, Ny, 1)]+Ydata[idx(i, j-1, Nx, Ny, 1)])/(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]);
            vy = (-Ydata[idx(i, j-1, Nx, Ny, 2)]+Ydata[idx(i, j-1, Nx, Ny, 2)])/(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]);
            p = (gama-1.0)*0.5*(Ydata[idx(i, j-1, Nx, Ny, 3)]+Ydata[idx(i, j-1, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
            h = (0.5*(Ydata[idx(i, j-1, Nx, Ny, 3)]+Ydata[idx(i, j-1, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]));
            rou = 0.5*(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]);
        }
        /* nature or transmissive boundary condition */
        if (flag == 2)
        {
            vx = (Ydata[idx(i, j-1, Nx, Ny, 1)]+Ydata[idx(i, j-1, Nx, Ny, 1)])/(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]);
            vy = (Ydata[idx(i, j-1, Nx, Ny, 2)]+Ydata[idx(i, j-1, Nx, Ny, 2)])/(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]);
            p = (gama-1.0)*0.5*(Ydata[idx(i, j-1, Nx, Ny, 3)]+Ydata[idx(i, j-1, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
            h = (0.5*(Ydata[idx(i, j-1, Nx, Ny, 3)]+Ydata[idx(i, j-1, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]));
            rou = 0.5*(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]);
        }
    }
    
    if (j==0){
        if (flag == 0)
        {
            vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i, Ny-1, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, Ny-1, Nx, Ny, 0)]);
            vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i, Ny-1, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, Ny-1, Nx, Ny, 0)]);
            p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, Ny-1, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, Ny-1, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
            h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, Ny-1, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, Ny-1, Nx, Ny, 0)]));
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, Ny-1, Nx, Ny, 0)]);
        }
        /* reflecting boundary condition */
        if (flag == 1)
        {
            vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
            vy = (-Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
            p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
            h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]));
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
        }
        /* nature or transmissive boundary condition */
        if (flag == 2)
        {
            vx = (Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
            vy = (Ydata[idx(i, j, Nx, Ny, 2)]+Ydata[idx(i, j, Nx, Ny, 2)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
            p = (gama-1.0)*0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])-0.5*0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)])*(gama-1.0)*(vx*vx+vy*vy);
            h = (0.5*(Ydata[idx(i, j, Nx, Ny, 3)]+Ydata[idx(i, j, Nx, Ny, 3)])+p)/(0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]));
            rou = 0.5*(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    a = sqrt(gama*p/rou);
    
    /* right eigenvectors derivated by myself */
    /*
     rhyegm[0][0] = 1.0;
     rhyegm[1][0] = vx;
     rhyegm[2][0] = vy-a;
     rhyegm[3][0] = vx*vx+vy*vy-h-vy*a+(2.0*a*a)/(gama-1.0);
     rhyegm[0][1] = 1.0;
     rhyegm[1][1] = 0.0;
     rhyegm[2][1] = vy;
     rhyegm[3][1] = vy*vy+a*a/(gama-1.0)-h;
     rhyegm[0][2] = 0.0;
     rhyegm[1][2] = 1.0;
     rhyegm[2][2] = 0.0;
     rhyegm[3][2] = vx;
     rhyegm[0][3] = 1.0;
     rhyegm[1][3] = vx;
     rhyegm[2][3] = vy+a;
     rhyegm[3][3] = vx*vx+vy*vy-h+vy*a+(2.0*a*a)/(gama-1.0);
     */
    /* right eigenvectors derivated from the book */
    rhyegm[0][0] = 1.0;
    rhyegm[1][0] = vx;
    rhyegm[2][0] = vy-a;
    rhyegm[3][0] = h-vy*a;
    rhyegm[0][1] = 1.0;
    rhyegm[1][1] = vx;
    rhyegm[2][1] = vy;
    rhyegm[3][1] = 0.5*(vx*vx+vy*vy);
    rhyegm[0][2] = 0.0;
    rhyegm[1][2] = 1.0;
    rhyegm[2][2] = 0.0;
    rhyegm[3][2] = vx;
    rhyegm[0][3] = 1.0;
    rhyegm[1][3] = vx;
    rhyegm[2][3] = vy+a;
    rhyegm[3][3] = h+vy*a;
    
    return 0;
}

/* Fill in the indicators of smoothness on x direction */
static int SetISX(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx,realtype *yxpdata, realtype *yxndata, long int i, long int j, long int Nx, long int Ny)
{
    long int n;
    
    /* compute indicators of smoothness of rou, qx, qy, E */
    for (n=0;n<4;n++){
        //printf("i=%li,j=%li,n=%li,yxpdata[1+n*7]=%g,yxpdata[2+n*7]=%g,yxpdata[3+n*7]=%g,yxpdata[4+n*7]=%g,yxpdata[5+n*7]=%g,yxndata[2+n*7]=%g,yxndata[3+n*7]=%g,yxndata[4+n*7]=%g,yxndata[5+n*7]=%g,yxndata[6+n*7]=%g\n",i,j,n,yxpdata[1+n*7],yxpdata[2+n*7],yxpdata[3+n*7],yxpdata[4+n*7],yxpdata[5+n*7],yxndata[2+n*7],yxndata[3+n*7],yxndata[4+n*7],yxndata[5+n*7],yxndata[6+n*7]);
        IS0_px[n]=(13.0/12.0)*(yxpdata[1+n*7]-2.0*yxpdata[2+n*7]+yxpdata[3+n*7])*(yxpdata[1+n*7]-2.0*yxpdata[2+n*7]+yxpdata[3+n*7])+(1.0/4.0)*(yxpdata[1+n*7]-4.0*yxpdata[2+n*7]+3.0*yxpdata[3+n*7])*(yxpdata[1+n*7]-4.0*yxpdata[2+n*7]+3.0*yxpdata[3+n*7]);
        
        IS1_px[n]=(13.0/12.0)*(yxpdata[2+n*7]-2.0*yxpdata[3+n*7]+yxpdata[4+n*7])*(yxpdata[2+n*7]-2.0*yxpdata[3+n*7]+yxpdata[4+n*7])+(1.0/4.0)*(yxpdata[2+n*7]-yxpdata[4+n*7])*(yxpdata[2+n*7]-yxpdata[4+n*7]);
        
        IS2_px[n]=(13.0/12.0)*(yxpdata[3+n*7]-2.0*yxpdata[4+n*7]+yxpdata[5+n*7])*(yxpdata[3+n*7]-2.0*yxpdata[4+n*7]+yxpdata[5+n*7])+(1.0/4.0)*(3.0*yxpdata[3+n*7]-4.0*yxpdata[4+n*7]+yxpdata[5+n*7])*(3.0*yxpdata[3+n*7]-4.0*yxpdata[4+n*7]+yxpdata[5+n*7]);
        
        IS0_nx[n]=(13.0/12.0)*(yxndata[4+n*7]-2.0*yxndata[5+n*7]+yxndata[6+n*7])*(yxndata[4+n*7]-2.0*yxndata[5+n*7]+yxndata[6+n*7])+(1.0/4.0)*(3.0*yxndata[4+n*7]-4.0*yxndata[5+n*7]+yxndata[6+n*7])*(3.0*yxndata[4+n*7]-4.0*yxndata[5+n*7]+yxndata[6+n*7]);
        
        IS1_nx[n]=(13.0/12.0)*(yxndata[3+n*7]-2.0*yxndata[4+n*7]+yxndata[5+n*7])*(yxndata[3+n*7]-2.0*yxndata[4+n*7]+yxndata[5+n*7])+(1.0/4.0)*(yxndata[3+n*7]-yxndata[5+n*7])*(yxndata[3+n*7]-yxndata[5+n*7]);
        
        IS2_nx[n]=(13.0/12.0)*(yxndata[2+n*7]-2.0*yxndata[3+n*7]+yxndata[4+n*7])*(yxndata[2+n*7]-2.0*yxndata[3+n*7]+yxndata[4+n*7])+(1.0/4.0)*(yxndata[2+n*7]-4.0*yxndata[3+n*7]+3.0*yxndata[4+n*7])*(yxndata[2+n*7]-4.0*yxndata[3+n*7]+3.0*yxndata[4+n*7]);
        //printf("i=%li,j=%li,n=%li,IS0_px[n]=%g,IS1_px[n]=%g,IS2_px[n]=%g,IS0_nx[n]=%g,IS1_nx[n]=%g,IS2_nx[n]=%g\n",i,j,n,IS0_px[n],IS1_px[n],IS2_px[n],IS0_nx[n],IS1_nx[n],IS2_nx[n]);
    }
    
    return 0;
}

/* Fill in the indicators of smoothness on y direction */
static int SetISY(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny,realtype *yypdata, realtype *yyndata, long int i, long int j, long int Nx, long int Ny)
{
    long int n;
    
    /* compute indicators of smoothness of rou, qx, qy, E */
    for (n=0;n<4;n++){
        //printf("i=%li,j=%li,n=%li,yypdata[1+n*7]=%g,yypdata[2+n*7]=%g,yypdata[3+n*7]=%g,yypdata[4+n*7]=%g,yypdata[5+n*7]=%g,yyndata[2+n*7]=%g,yyndata[3+n*7]=%g,yyndata[4+n*7]=%g,yyndata[5+n*7]=%g,yyndata[6+n*7]=%g\n",i,j,n,yypdata[1+n*7],yypdata[2+n*7],yypdata[3+n*7],yypdata[4+n*7],yypdata[5+n*7],yyndata[2+n*7],yyndata[3+n*7],yyndata[4+n*7],yyndata[5+n*7],yyndata[6+n*7]);
        IS0_py[n]=(13.0/12.0)*(yypdata[1+n*7]-2.0*yypdata[2+n*7]+yypdata[3+n*7])*(yypdata[1+n*7]-2.0*yypdata[2+n*7]+yypdata[3+n*7])+(1.0/4.0)*(yypdata[1+n*7]-4.0*yypdata[2+n*7]+3.0*yypdata[3+n*7])*(yypdata[1+n*7]-4.0*yypdata[2+n*7]+3.0*yypdata[3+n*7]);
        
        IS1_py[n]=(13.0/12.0)*(yypdata[2+n*7]-2.0*yypdata[3+n*7]+yypdata[4+n*7])*(yypdata[2+n*7]-2.0*yypdata[3+n*7]+yypdata[4+n*7])+(1.0/4.0)*(yypdata[2+n*7]-yypdata[4+n*7])*(yypdata[2+n*7]-yypdata[4+n*7]);
        
        IS2_py[n]=(13.0/12.0)*(yypdata[3+n*7]-2.0*yypdata[4+n*7]+yypdata[5+n*7])*(yypdata[3+n*7]-2.0*yypdata[4+n*7]+yypdata[5+n*7])+(1.0/4.0)*(3.0*yypdata[3+n*7]-4.0*yypdata[4+n*7]+yypdata[5+n*7])*(3.0*yypdata[3+n*7]-4.0*yypdata[4+n*7]+yypdata[5+n*7]);
        
        IS0_ny[n]=(13.0/12.0)*(yyndata[4+n*7]-2.0*yyndata[5+n*7]+yyndata[6+n*7])*(yyndata[4+n*7]-2.0*yyndata[5+n*7]+yyndata[6+n*7])+(1.0/4.0)*(3.0*yyndata[4+n*7]-4.0*yyndata[5+n*7]+yyndata[6+n*7])*(3.0*yyndata[4+n*7]-4.0*yyndata[5+n*7]+yyndata[6+n*7]);
        
        IS1_ny[n]=(13.0/12.0)*(yyndata[3+n*7]-2.0*yyndata[4+n*7]+yyndata[5+n*7])*(yyndata[3+n*7]-2.0*yyndata[4+n*7]+yyndata[5+n*7])+(1.0/4.0)*(yyndata[3+n*7]-yyndata[5+n*7])*(yyndata[3+n*7]-yyndata[5+n*7]);
        
        IS2_ny[n]=(13.0/12.0)*(yyndata[2+n*7]-2.0*yyndata[3+n*7]+yyndata[4+n*7])*(yyndata[2+n*7]-2.0*yyndata[3+n*7]+yyndata[4+n*7])+(1.0/4.0)*(yyndata[2+n*7]-4.0*yyndata[3+n*7]+3.0*yyndata[4+n*7])*(yyndata[2+n*7]-4.0*yyndata[3+n*7]+3.0*yyndata[4+n*7]);
        //printf("i=%li,j=%li,n=%li,IS0_py[n]=%g,IS1_py[n]=%g,IS2_py[n]=%g,IS0_ny[n]=%g,IS1_ny[n]=%g,IS2_ny[n]=%g\n",i,j,n,IS0_py[n],IS1_py[n],IS2_py[n],IS0_ny[n],IS1_ny[n],IS2_ny[n]);
    }
    
    return 0;
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

/* Get the positive and negative parts of right interface on x direction */
static int SetUx(realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype *yxpdata, realtype *yxndata, realtype *u_tpphx, realtype *u_tnphx, realtype *u_tpnhx, realtype *u_tnnhx, long int i, long int j, long int Nx, long int Ny)
{       
    long int n;
    /* compute positive and negative solutions on the interface */
    for (n=0;n<4;n++){
        u_tpphx[n]=w0_px[n]*((2.0/6.0)*yxpdata[1+n*7]-(7.0/6.0)*yxpdata[2+n*7]+(11.0/6.0)*yxpdata[3+n*7])+w1_px[n]*((-1.0/6.0)*yxpdata[2+n*7]+(5.0/6.0)*yxpdata[3+n*7]+(2.0/6.0)*yxpdata[4+n*7])+w2_px[n]*((2.0/6.0)*yxpdata[3+n*7]+(5.0/6.0)*yxpdata[4+n*7]-(1.0/6.0)*yxpdata[5+n*7]);
        
        u_tnphx[n]=w2_nx[n]*((-1.0/6.0)*yxndata[2+n*7]+(5.0/6.0)*yxndata[3+n*7]+(2.0/6.0)*yxndata[4+n*7])+w1_nx[n]*((2.0/6.0)*yxndata[3+n*7]+(5.0/6.0)*yxndata[4+n*7]-(1.0/6.0)*yxndata[5+n*7])+w0_nx[n]*((11.0/6.0)*yxndata[4+n*7]-(7.0/6.0)*yxndata[5+n*7]+(2.0/6.0)*yxndata[6+n*7]);
        
        u_tpnhx[n]=w0_px[n]*((2.0/6.0)*yxpdata[0+n*7]-(7.0/6.0)*yxpdata[1+n*7]+(11.0/6.0)*yxpdata[2+n*7])+w1_px[n]*((-1.0/6.0)*yxpdata[1+n*7]+(5.0/6.0)*yxpdata[2+n*7]+(2.0/6.0)*yxpdata[3+n*7])+w2_px[n]*((2.0/6.0)*yxpdata[2+n*7]+(5.0/6.0)*yxpdata[3+n*7]-(1.0/6.0)*yxpdata[4+n*7]);
        
        u_tnnhx[n]=w2_nx[n]*((-1.0/6.0)*yxndata[1+n*7]+(5.0/6.0)*yxndata[2+n*7]+(2.0/6.0)*yxndata[3+n*7])+w1_nx[n]*((2.0/6.0)*yxndata[2+n*7]+(5.0/6.0)*yxndata[3+n*7]-(1.0/6.0)*yxndata[4+n*7])+w0_nx[n]*((11.0/6.0)*yxndata[3+n*7]-(7.0/6.0)*yxndata[4+n*7]+(2.0/6.0)*yxndata[5+n*7]);
    }
    return 0;
}

/* Get the positive and negative parts of right interface  on y direction */
static int SetUy(realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype *yypdata, realtype *yyndata, realtype *u_tpphy, realtype *u_tnphy, realtype *u_tpnhy, realtype *u_tnnhy, long int i, long int j, long int Nx, long int Ny)
{
    long int n;
    /* compute positive and negative solutions on the interface */
    for (n=0;n<4;n++){
        u_tpphy[n]=w0_py[n]*((2.0/6.0)*yypdata[1+n*7]-(7.0/6.0)*yypdata[2+n*7]+(11.0/6.0)*yypdata[3+n*7])+w1_py[n]*((-1.0/6.0)*yypdata[2+n*7]+(5.0/6.0)*yypdata[3+n*7]+(2.0/6.0)*yypdata[4+n*7])+w2_py[n]*((2.0/6.0)*yypdata[3+n*7]+(5.0/6.0)*yypdata[4+n*7]-(1.0/6.0)*yypdata[5+n*7]);
        
        u_tnphy[n]=w2_ny[n]*((-1.0/6.0)*yyndata[2+n*7]+(5.0/6.0)*yyndata[3+n*7]+(2.0/6.0)*yyndata[4+n*7])+w1_ny[n]*((2.0/6.0)*yyndata[3+n*7]+(5.0/6.0)*yyndata[4+n*7]-(1.0/6.0)*yyndata[5+n*7])+w0_ny[n]*((11.0/6.0)*yyndata[4+n*7]-(7.0/6.0)*yyndata[5+n*7]+(2.0/6.0)*yyndata[6+n*7]);
        
        u_tpnhy[n]=w0_py[n]*((2.0/6.0)*yypdata[0+n*7]-(7.0/6.0)*yypdata[1+n*7]+(11.0/6.0)*yypdata[2+n*7])+w1_py[n]*((-1.0/6.0)*yypdata[1+n*7]+(5.0/6.0)*yypdata[2+n*7]+(2.0/6.0)*yypdata[3+n*7])+w2_py[n]*((2.0/6.0)*yypdata[2+n*7]+(5.0/6.0)*yypdata[3+n*7]-(1.0/6.0)*yypdata[4+n*7]);
        
        u_tnnhy[n]=w2_ny[n]*((-1.0/6.0)*yyndata[1+n*7]+(5.0/6.0)*yyndata[2+n*7]+(2.0/6.0)*yyndata[3+n*7])+w1_ny[n]*((2.0/6.0)*yyndata[2+n*7]+(5.0/6.0)*yyndata[3+n*7]-(1.0/6.0)*yyndata[4+n*7])+w0_ny[n]*((11.0/6.0)*yyndata[3+n*7]-(7.0/6.0)*yyndata[4+n*7]+(2.0/6.0)*yyndata[5+n*7]);
    }
    
    return 0;
}

/* Fill in values of tao in the whole domain */
static int Gettao(realtype *yxbackupdata, realtype *yxbackdowndata, realtype *yybackupdata, realtype *yybackdowndata, realtype *taoup, realtype *taodown, long int Nx, long int Ny, realtype gama)
{
    /* declare parameters */
    long int i, j;
    
    realtype *vxxup = new realtype [(Nx)*Ny];
    realtype *vyxup = new realtype [(Nx)*Ny];
    realtype *vxyup = new realtype [Nx*(Ny)];
    realtype *vyyup = new realtype [Nx*(Ny)];
    realtype *pxup = new realtype [(Nx)*Ny];
    realtype *pyup = new realtype [Nx*(Ny)];
    
    realtype *vxxdown = new realtype [(Nx)*Ny];
    realtype *vyxdown = new realtype [(Nx)*Ny];
    realtype *vxydown = new realtype [Nx*(Ny)];
    realtype *vyydown = new realtype [Nx*(Ny)];
    realtype *pxdown = new realtype [(Nx)*Ny];
    realtype *pydown = new realtype [Nx*(Ny)];
    
    /* Set values into vx and vy on x direction*/
    for(j=0;j<Ny;j++){
        for(i=0;i<(Nx);i++){
            vxxup[idx_v(i,j,(Nx))]=yxbackupdata[idx(i, j, (Nx), Ny, 1)]/yxbackupdata[idx(i, j, (Nx), Ny, 0)];
            //printf("i = %li, j=%li, yxbackdata[idx(i, j, (Nx+1), Ny, 1)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",i,j,yxbackdata[idx(i, j, (Nx+1), Ny, 1)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
            vyxup[idx_v(i,j,(Nx))]=yxbackupdata[idx(i, j, (Nx), Ny, 2)]/yxbackupdata[idx(i, j, (Nx), Ny, 0)];
            //printf("yxbackdata[idx(i, j, (Nx+1), Ny, 2)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",yxbackdata[idx(i, j, (Nx+1), Ny, 2)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
            //printf("i = %li, j=%li, vxx[idx_v(i,j,(Nx+1))]=%f,vyx[idx_v(i,j,(Nx+1))]=%f\n",i,j,vxx[idx_v(i,j,(Nx+1))],vyx[idx_v(i,j,(Nx+1))]);
        }
    }
    
    /* Set values into vx and vy on x direction*/
    for(j=0;j<Ny;j++){
        for(i=0;i<(Nx);i++){
            vxxdown[idx_v(i,j,(Nx))]=yxbackdowndata[idx(i, j, (Nx), Ny, 1)]/yxbackdowndata[idx(i, j, (Nx), Ny, 0)];
            //printf("i = %li, j=%li, yxbackdata[idx(i, j, (Nx+1), Ny, 1)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",i,j,yxbackdata[idx(i, j, (Nx+1), Ny, 1)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
            vyxdown[idx_v(i,j,(Nx))]=yxbackdowndata[idx(i, j, (Nx), Ny, 2)]/yxbackdowndata[idx(i, j, (Nx), Ny, 0)];
            //printf("yxbackdata[idx(i, j, (Nx+1), Ny, 2)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",yxbackdata[idx(i, j, (Nx+1), Ny, 2)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
            //printf("i = %li, j=%li, vxx[idx_v(i,j,(Nx+1))]=%f,vyx[idx_v(i,j,(Nx+1))]=%f\n",i,j,vxx[idx_v(i,j,(Nx+1))],vyx[idx_v(i,j,(Nx+1))]);
        }
    }
    
    /* Set values into vxy and vyy on y direction*/
    for(i=0;i<Nx;i++){
        for(j=0;j<(Ny);j++){
            vxyup[idx_v(i,j,Nx)]=yybackupdata[idx(i, j, Nx, (Ny), 1)]/yybackupdata[idx(i, j, Nx, (Ny), 0)];
            //printf("i = %li, j=%li, yybackdata[idx(i, j, Nx, (Ny+1), 1)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",i,j,yybackdata[idx(i, j, Nx, (Ny+1), 1)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
            vyyup[idx_v(i,j,Nx)]=yybackupdata[idx(i, j, Nx, (Ny), 2)]/yybackupdata[idx(i, j, Nx, (Ny), 0)];
            //printf("yybackdata[idx(i, j, Nx, (Ny+1), 2)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",yybackdata[idx(i, j, Nx, (Ny+1), 2)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
        }
    }
    
    /* Set values into vxy and vyy on y direction*/
    for(i=0;i<Nx;i++){
        for(j=0;j<(Ny);j++){
            vxydown[idx_v(i,j,Nx)]=yybackdowndata[idx(i, j, Nx, (Ny), 1)]/yybackdowndata[idx(i, j, Nx, (Ny), 0)];
            //printf("i = %li, j=%li, yybackdata[idx(i, j, Nx, (Ny+1), 1)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",i,j,yybackdata[idx(i, j, Nx, (Ny+1), 1)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
            vyydown[idx_v(i,j,Nx)]=yybackdowndata[idx(i, j, Nx, (Ny), 2)]/yybackdowndata[idx(i, j, Nx, (Ny), 0)];
            //printf("yybackdata[idx(i, j, Nx, (Ny+1), 2)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",yybackdata[idx(i, j, Nx, (Ny+1), 2)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
        }
    }
    
    /* Compute the values of tao in the whole domain on x direction*/
    for(j=0;j<Ny;j++){
        for(i=0;i<(Nx);i++){
            pxup[idx_v(i,j,(Nx))] =  yxbackupdata[idx(i, j, (Nx), Ny, 0)]*(gama-1.0)*(yxbackupdata[idx(i, j, (Nx), Ny, 3)]/yxbackupdata[idx(i, j, (Nx), Ny, 0)]-0.5*(vxxup[idx_v(i,j,(Nx))]*vxxup[idx_v(i,j,(Nx))]+vyxup[idx_v(i,j,(Nx))]*vyxup[idx_v(i,j,(Nx))]));
            taoup[idx(i, j, Nx, Ny, 0)] = yxbackupdata[idx(i, j, (Nx), Ny, 1)]*vxxup[idx_v(i,j,(Nx))]+pxup[idx_v(i,j,(Nx))];
            taoup[idx(i, j, Nx, Ny, 1)] = yxbackupdata[idx(i, j, (Nx), Ny, 1)]*vyxup[idx_v(i,j,(Nx))];
        } 
	}
    
    /* Compute the values of tao in the whole domain on x direction*/
    for(j=0;j<Ny;j++){
        for(i=0;i<(Nx);i++){
            pxdown[idx_v(i,j,(Nx))] =  yxbackdowndata[idx(i, j, (Nx), Ny, 0)]*(gama-1.0)*(yxbackdowndata[idx(i, j, (Nx), Ny, 3)]/yxbackdowndata[idx(i, j, (Nx), Ny, 0)]-0.5*(vxxdown[idx_v(i,j,(Nx))]*vxxdown[idx_v(i,j,(Nx))]+vyxdown[idx_v(i,j,(Nx))]*vyxdown[idx_v(i,j,(Nx))]));
            taodown[idx(i, j, Nx, Ny, 0)] = yxbackdowndata[idx(i, j, (Nx), Ny, 1)]*vxxdown[idx_v(i,j,(Nx))]+pxdown[idx_v(i,j,(Nx))];
            taodown[idx(i, j, Nx, Ny, 1)] = yxbackdowndata[idx(i, j, (Nx), Ny, 1)]*vyxdown[idx_v(i,j,(Nx))];
        }
	}
    
    /* Compute the values of tao in the whole domain on y direction*/
    for (i=0;i<Nx;i++){
        for (j=0;j<(Ny);j++){
            pyup[idx_v(i,j,Nx)] =  yybackupdata[idx(i, j, Nx, (Ny), 0)]*(gama-1.0)*(yybackupdata[idx(i, j, Nx, (Ny), 3)]/yybackupdata[idx(i, j, Nx, (Ny), 0)]-0.5*(vyxup[idx_v(i,j,Nx)]*vyxup[idx_v(i,j,Nx)]+vyyup[idx_v(i,j,Nx)]*vyyup[idx_v(i,j,Nx)]));
            taoup[idx(i, j, Nx, Ny, 2)] = yybackupdata[idx(i, j, Nx, (Ny), 2)]*vxyup[idx_v(i,j,Nx)];
            taoup[idx(i, j, Nx, Ny, 3)] = yybackupdata[idx(i, j, Nx, (Ny), 2)]*vyyup[idx_v(i,j,Nx)]+pyup[idx_v(i,j,Nx)];
        }
    }
    
    for (i=0;i<Nx;i++){
        for (j=0;j<(Ny);j++){
            pydown[idx_v(i,j,Nx)] =  yybackdowndata[idx(i, j, Nx, (Ny), 0)]*(gama-1.0)*(yybackdowndata[idx(i, j, Nx, (Ny), 3)]/yybackdowndata[idx(i, j, Nx, (Ny), 0)]-0.5*(vyxdown[idx_v(i,j,Nx)]*vyxdown[idx_v(i,j,Nx)]+vyydown[idx_v(i,j,Nx)]*vyydown[idx_v(i,j,Nx)]));
            taodown[idx(i, j, Nx, Ny, 2)] = yybackdowndata[idx(i, j, Nx, (Ny), 2)]*vxydown[idx_v(i,j,Nx)];
            taodown[idx(i, j, Nx, Ny, 3)] = yybackdowndata[idx(i, j, Nx, (Ny), 2)]*vyydown[idx_v(i,j,Nx)]+pydown[idx_v(i,j,Nx)];
        }
    }
    
    delete []vxxup;
    delete []vxyup;
    delete []vyxup;
    delete []vyyup;
    delete []pxup;
    delete []pyup;
    delete []vxxdown;
    delete []vxydown;
    delete []vyxdown;
    delete []vyydown;
    delete []pxdown;
    delete []pydown;
    
    return 0;
}

/* Fill in the value of J in the whole domain */
static int GetCj(realtype *yxbackupdata, realtype *yxbackdowndata, realtype *yybackupdata, realtype *yybackdowndata, realtype *Cjup, realtype *Cjdown, long int Nx, long int Ny, realtype gama)
{
    /* declare parameters */
    long int i, j;
    
    realtype *vxxup = new realtype [(Nx)*Ny];
    realtype *vyxup = new realtype [(Nx)*Ny];
    realtype *vxyup = new realtype [Nx*(Ny)];
    realtype *vyyup = new realtype [Nx*(Ny)];
    realtype *pxup = new realtype [(Nx)*Ny];
    realtype *pyup = new realtype [Nx*(Ny)];
    
    realtype *vxxdown = new realtype [(Nx)*Ny];
    realtype *vyxdown = new realtype [(Nx)*Ny];
    realtype *vxydown = new realtype [Nx*(Ny)];
    realtype *vyydown = new realtype [Nx*(Ny)];
    realtype *pxdown = new realtype [(Nx)*Ny];
    realtype *pydown = new realtype [Nx*(Ny)];
    
    /* Set values into vx and vy on x direction*/
    for(j=0;j<Ny;j++){
        for(i=0;i<(Nx);i++){
            vxxup[idx_v(i,j,(Nx))]=yxbackupdata[idx(i, j, (Nx), Ny, 1)]/yxbackupdata[idx(i, j, (Nx), Ny, 0)];
            //printf("i = %li, j=%li, yxbackdata[idx(i, j, (Nx+1), Ny, 1)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",i,j,yxbackdata[idx(i, j, (Nx+1), Ny, 1)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
            vyxup[idx_v(i,j,(Nx))]=yxbackupdata[idx(i, j, (Nx), Ny, 2)]/yxbackupdata[idx(i, j, (Nx), Ny, 0)];
            //printf("yxbackdata[idx(i, j, (Nx+1), Ny, 2)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",yxbackdata[idx(i, j, (Nx+1), Ny, 2)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
            //printf("i = %li, j=%li, vxx[idx_v(i,j,(Nx+1))]=%f,vyx[idx_v(i,j,(Nx+1))]=%f\n",i,j,vxx[idx_v(i,j,(Nx+1))],vyx[idx_v(i,j,(Nx+1))]);
        }
    }
    
    /* Set values into vx and vy on x direction*/
    for(j=0;j<Ny;j++){
        for(i=0;i<(Nx);i++){
            vxxdown[idx_v(i,j,(Nx))]=yxbackdowndata[idx(i, j, (Nx), Ny, 1)]/yxbackdowndata[idx(i, j, (Nx), Ny, 0)];
            //printf("i = %li, j=%li, yxbackdata[idx(i, j, (Nx+1), Ny, 1)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",i,j,yxbackdata[idx(i, j, (Nx+1), Ny, 1)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
            vyxdown[idx_v(i,j,(Nx))]=yxbackdowndata[idx(i, j, (Nx), Ny, 2)]/yxbackdowndata[idx(i, j, (Nx), Ny, 0)];
            //printf("yxbackdata[idx(i, j, (Nx+1), Ny, 2)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",yxbackdata[idx(i, j, (Nx+1), Ny, 2)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
            //printf("i = %li, j=%li, vxx[idx_v(i,j,(Nx+1))]=%f,vyx[idx_v(i,j,(Nx+1))]=%f\n",i,j,vxx[idx_v(i,j,(Nx+1))],vyx[idx_v(i,j,(Nx+1))]);
        }
    }
    
    /* Set values into vxy and vyy on y direction*/
    for(i=0;i<Nx;i++){
        for(j=0;j<(Ny);j++){
            vxyup[idx_v(i,j,Nx)]=yybackupdata[idx(i, j, Nx, (Ny), 1)]/yybackupdata[idx(i, j, Nx, (Ny), 0)];
            //printf("i = %li, j=%li, yybackdata[idx(i, j, Nx, (Ny+1), 1)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",i,j,yybackdata[idx(i, j, Nx, (Ny+1), 1)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
            vyyup[idx_v(i,j,Nx)]=yybackupdata[idx(i, j, Nx, (Ny), 2)]/yybackupdata[idx(i, j, Nx, (Ny), 0)];
            //printf("yybackdata[idx(i, j, Nx, (Ny+1), 2)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",yybackdata[idx(i, j, Nx, (Ny+1), 2)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
        }
    }
    
    /* Set values into vxy and vyy on y direction*/
    for(i=0;i<Nx;i++){
        for(j=0;j<(Ny);j++){
            vxydown[idx_v(i,j,Nx)]=yybackdowndata[idx(i, j, Nx, (Ny), 1)]/yybackdowndata[idx(i, j, Nx, (Ny), 0)];
            //printf("i = %li, j=%li, yybackdata[idx(i, j, Nx, (Ny+1), 1)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",i,j,yybackdata[idx(i, j, Nx, (Ny+1), 1)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
            vyydown[idx_v(i,j,Nx)]=yybackdowndata[idx(i, j, Nx, (Ny), 2)]/yybackdowndata[idx(i, j, Nx, (Ny), 0)];
            //printf("yybackdata[idx(i, j, Nx, (Ny+1), 2)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",yybackdata[idx(i, j, Nx, (Ny+1), 2)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
        }
    }
    
    /* Compute the values of tao in the whole domain on x direction*/
    for(j=0;j<Ny;j++){
        for(i=0;i<(Nx);i++){
            pxup[idx_v(i,j,(Nx))] =  yxbackupdata[idx(i, j, (Nx), Ny, 0)]*(gama-1.0)*(yxbackupdata[idx(i, j, (Nx), Ny, 3)]/yxbackupdata[idx(i, j, (Nx), Ny, 0)]-0.5*(vxxup[idx_v(i,j,(Nx))]*vxxup[idx_v(i,j,(Nx))]+vyxup[idx_v(i,j,(Nx))]*vyxup[idx_v(i,j,(Nx))]));  
            Cjup[idx(i, j, Nx, Ny, 0)] = (yxbackupdata[idx(i, j, (Nx), Ny, 3)]+pxup[idx_v(i,j,(Nx))])*vxxup[idx_v(i,j,(Nx))];
        } 
    }
    
    /* Compute the values of tao in the whole domain on x direction*/
    for(j=0;j<Ny;j++){
        for(i=0;i<(Nx);i++){
            pxdown[idx_v(i,j,(Nx))] =  yxbackdowndata[idx(i, j, (Nx), Ny, 0)]*(gama-1.0)*(yxbackdowndata[idx(i, j, (Nx), Ny, 3)]/yxbackdowndata[idx(i, j, (Nx), Ny, 0)]-0.5*(vxxdown[idx_v(i,j,(Nx))]*vxxdown[idx_v(i,j,(Nx))]+vyxdown[idx_v(i,j,(Nx))]*vyxdown[idx_v(i,j,(Nx))]));
            Cjdown[idx(i, j, Nx, Ny, 0)] = (yxbackdowndata[idx(i, j, (Nx), Ny, 3)]+pxdown[idx_v(i,j,(Nx))])*vxxdown[idx_v(i,j,(Nx))];
        }
    }

    /* Compute the values of tao in the whole domain on y direction*/
    for (i=0;i<Nx;i++){
        for (j=0;j<(Ny);j++){
            pyup[idx_v(i,j,Nx)] =  yybackupdata[idx(i, j, Nx, (Ny), 0)]*(gama-1.0)*(yybackupdata[idx(i, j, Nx, (Ny), 3)]/yybackupdata[idx(i, j, Nx, (Ny), 0)]-0.5*(vyxup[idx_v(i,j,Nx)]*vyxup[idx_v(i,j,Nx)]+vyyup[idx_v(i,j,Nx)]*vyyup[idx_v(i,j,Nx)]));
            Cjup[idx(i, j, Nx, Ny, 1)] = (yybackupdata[idx(i, j, Nx, (Ny), 3)]+pyup[idx_v(i,j,Nx)])*vyyup[idx_v(i,j,Nx)];
        }
    }
    
    /* Compute the values of tao in the whole domain on y direction*/
    for (i=0;i<Nx;i++){
        for (j=0;j<(Ny);j++){
            pydown[idx_v(i,j,Nx)] =  yybackdowndata[idx(i, j, Nx, (Ny), 0)]*(gama-1.0)*(yybackdowndata[idx(i, j, Nx, (Ny), 3)]/yybackdowndata[idx(i, j, Nx, (Ny), 0)]-0.5*(vyxdown[idx_v(i,j,Nx)]*vyxdown[idx_v(i,j,Nx)]+vyydown[idx_v(i,j,Nx)]*vyydown[idx_v(i,j,Nx)]));
            Cjdown[idx(i, j, Nx, Ny, 1)] = (yybackdowndata[idx(i, j, Nx, (Ny), 3)]+pydown[idx_v(i,j,Nx)])*vyydown[idx_v(i,j,Nx)];
        }
    }
    
    delete []vxxup;
    delete []vxyup;
    delete []vyxup;
    delete []vyyup;
    delete []pxup;
    delete []pyup;
    delete []vxxdown;
    delete []vxydown;
    delete []vyxdown;
    delete []vyydown;
    delete []pxdown;
    delete []pydown;
    
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
