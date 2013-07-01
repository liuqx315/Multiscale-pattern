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
static int Gettao(N_Vector y, N_Vector tao, long int Nx, long int Ny, realtype gama);

/* Set value of J in whole domain*/
static int GetCj(N_Vector y, N_Vector Cj, long int Nx, long int Ny, realtype gama);

/* Fill in the indicators of smoothness on x direction */
static int SetISX(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx,realtype *ypdata, realtype *yndata, realtype *taopdata, realtype *taondata, realtype *Cjpdata, realtype *Cjndata, long int i, long int j, long int Nx, long int Ny);

/* Fill in the indicators of smoothness on y direction */
static int SetISY(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny,realtype *ypdata, realtype *yndata, realtype *taopdata, realtype *taondata, realtype *Cjpdata, realtype *Cjndata, long int i, long int j, long int Nx, long int Ny);

/* Fill in the stencil weights on x direction */
static int Setalphawx(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx, realtype *alpha_0px, realtype *alpha_1px, realtype *alpha_2px, realtype *alpha_0nx, realtype *alpha_1nx, realtype *alpha_2nx, realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype Epsilon);

/* Fill in the stencil weights on y direction */
static int Setalphawy(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny, realtype *alpha_0py, realtype *alpha_1py, realtype *alpha_2py, realtype *alpha_0ny, realtype *alpha_1ny, realtype *alpha_2ny, realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype Epsilon);

/* Get the derivative on x direction */
static int SetUx(realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype *ypdata, realtype *yndata, realtype *taopdata, realtype *taondata, realtype *Cjpdata, realtype *Cjndata, realtype *u_tpphx, realtype *u_tnphx, realtype *u_tpnhx, realtype *u_tnnhx, long int i, long int j, long int Nx, long int Ny);

/* Get the derivative on y direction */
static int SetUy(realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype *ypdata, realtype *yndata, realtype *taopdata, realtype *taondata, realtype *Cjpdata, realtype *Cjndata, realtype *u_tpphy, realtype *u_tnphy, realtype *u_tpnhy, realtype *u_tnnhy, long int i, long int j, long int Nx, long int Ny);

int main(int argc, const char * argv[])
{
    /* general problem parameters */
    realtype T0 = RCONST(0.0);
    realtype Tf = RCONST(0.3);
    int Nt = 6;
    int Nvar = 4;
    UserData udata = NULL;
    realtype *data;
    long int Nx, Ny, NEQ, i, j;
    realtype Lx, Ly, gama;
    
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
            data[idx(i,j,Nx,Ny,3)] =  0.3/((udata->gama-1.0)*0.5323);  /* E */
	  }
	  if (udata->dx*(i+0.5)>0.5&&udata->dy*(j+0.5)>0.5){
            data[idx(i,j,Nx,Ny,0)] =  1.5;  /* rou */
            data[idx(i,j,Nx,Ny,1)] =  1.5*0.0;  /* qx */
            data[idx(i,j,Nx,Ny,2)] =  1.5*0.0;  /* qy */
            data[idx(i,j,Nx,Ny,3)] =  1.5/((udata->gama-1.0)*1.5);  /* E */
	  }
	  if (udata->dx*(i+0.5)<0.5&&udata->dy*(j+0.5)<0.5){
            data[idx(i,j,Nx,Ny,0)] =  0.138;  /* rou */
            data[idx(i,j,Nx,Ny,1)] =  1.206*0.138;  /* qx */
            data[idx(i,j,Nx,Ny,2)] =  1.206*0.138;  /* qy */
            data[idx(i,j,Nx,Ny,3)] =  0.029/((udata->gama-1.0)*0.138);  /* E */
	  }
	  if (udata->dx*(i+0.5)>0.5&&udata->dy*(j+0.5)<0.5){
            data[idx(i,j,Nx,Ny,0)] =  0.5323;  /* rou */
            data[idx(i,j,Nx,Ny,1)] =  1.206*0.0;  /* qx */
            data[idx(i,j,Nx,Ny,2)] =  1.206*0.5323;  /* qy */
            data[idx(i,j,Nx,Ny,3)] =  0.3/((udata->gama-1.0)*0.5323);  /* E */
	  }
        }   
    }
    
    /* Call ARKodeCreate to create the solver memory */
    arkode_mem = ARKodeCreate();
    if (check_flag((void *) arkode_mem, "ARKodeCreate", 0)) return 1;
    
    /* Set solver parameters */
    realtype reltol  = 1.e-3;
    realtype abstol  = 1.e-6;
    
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
    flag = ARKodeSetInitStep(arkode_mem, 0.05);
    if (check_flag(&flag, "ARKodeSetInitStep", 1)) return 1;

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
    realtype egv1, egv2, egv3, egvmax;
    int flag;
    realtype Epsilon = 1.e-6;

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
    N_Vector yp = NULL;
    N_Vector yn = NULL;
    N_Vector tao = NULL;
    N_Vector Cj = NULL;
    N_Vector yx = NULL;
    N_Vector yy = NULL;
    N_Vector taop = NULL;
    N_Vector taon = NULL;
    N_Vector Cjp = NULL;
    N_Vector Cjn = NULL;
    
    /* Create serial vector of length NEQ and NEQS */
    yp = N_VNew_Serial(NEQ);
    if (check_flag((void *) yp, "N_VNew_Serial", 0)) return 1;
    yn = N_VNew_Serial(NEQ);
    if (check_flag((void *) yn, "N_VNew_Serial", 0)) return 1;
    tao = N_VNew_Serial(NEQ);
    if (check_flag((void *) tao, "N_VNew_Serial", 0)) return 1;
    Cj = N_VNew_Serial(NEQS);
    if (check_flag((void *) Cj, "N_VNew_Serial", 0)) return 1;
    yx = N_VNew_Serial(NEQ);
    if (check_flag((void *) yx, "N_VNew_Serial", 0)) return 1;
    yy = N_VNew_Serial(NEQ);
    if (check_flag((void *) yy, "N_VNew_Serial", 0)) return 1;
    taop = N_VNew_Serial(NEQ);
    if (check_flag((void *) taop, "N_VNew_Serial", 0)) return 1;
    taon = N_VNew_Serial(NEQ);
    if (check_flag((void *) taon, "N_VNew_Serial", 0)) return 1;
    Cjp = N_VNew_Serial(NEQS);
    if (check_flag((void *) Cjp, "N_VNew_Serial", 0)) return 1;
    Cjn = N_VNew_Serial(NEQS);
    if (check_flag((void *) Cjn, "N_VNew_Serial", 0)) return 1;
    
    /* fill in the value of tao and J in the whole domain */
    flag=Gettao(y, tao, Nx, Ny, gama);
    if (flag!=0) printf("error in Gettao function \n");
    flag=GetCj(y, Cj, Nx, Ny, gama);
    if (flag!=0) printf("error in GetCj function \n");
    
    /* access data arrays */
    realtype *Ydata = N_VGetArrayPointer(y);
    if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;
    realtype *dYdata = N_VGetArrayPointer(ydot);
    if (check_flag((void *) dYdata, "N_VGetArrayPointer", 0)) return 1;
    realtype *ypdata = N_VGetArrayPointer(yp);
    if (check_flag((void *) ypdata, "N_VGetArrayPointer", 0)) return 1;
    realtype *yndata = N_VGetArrayPointer(yn);
    if (check_flag((void *) yndata, "N_VGetArrayPointer", 0)) return 1;
    realtype *taodata = N_VGetArrayPointer(tao);
    if (check_flag((void *) taodata, "N_VGetArrayPointer", 0)) return 1;
    realtype *Cjdata = N_VGetArrayPointer(Cj);
    if (check_flag((void *) Cjdata, "N_VGetArrayPointer", 0)) return 1;
    realtype *yxdata = N_VGetArrayPointer(yx);
    if (check_flag((void *) yxdata, "N_VGetArrayPointer", 0)) return 1;
    realtype *yydata = N_VGetArrayPointer(yy);
    if (check_flag((void *) yydata, "N_VGetArrayPointer", 0)) return 1;
    realtype *taopdata = N_VGetArrayPointer(taop);
    if (check_flag((void *) taopdata, "N_VGetArrayPointer", 0)) return 1;
    realtype *taondata = N_VGetArrayPointer(taon);
    if (check_flag((void *) taondata, "N_VGetArrayPointer", 0)) return 1;
    realtype *Cjpdata = N_VGetArrayPointer(Cjp);
    if (check_flag((void *) Cjpdata, "N_VGetArrayPointer", 0)) return 1;
    realtype *Cjndata = N_VGetArrayPointer(Cjn);
    if (check_flag((void *) Cjndata, "N_VGetArrayPointer", 0)) return 1;
    
    /* compute max absolue eigenvalue and fill in ypdata, yndata, taopdata, taondata, Cjpdata, Cjpdata */
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
            /* get different eigenvalues */
            egv1 = Lx*(Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)])+Ly*(Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)]);
            egv2 = Lx*(Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)])+Ly*(Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)])-(gama*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]/Ydata[idx(i, j, Nx, Ny, 0)])*sqrt(Lx*Lx+Ly*Ly);
            egv3 = Lx*(Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)])+Ly*(Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)])+(gama*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]/Ydata[idx(i, j, Nx, Ny, 0)])*sqrt(Lx*Lx+Ly*Ly);
           /* get max absolute eigenvalue */
           egvmax = (fabs(egv1)>fabs(egv2))? fabs(egv1) : fabs(egv2);
           egvmax = (egvmax>fabs(egv3))? egvmax : fabs(egv3);
            /* fill in ypdata, yndata, taopdata, taondata, Cjpdata, Cjpdata */
            ypdata[idx(i, j, Nx, Ny, 1)]=0.5*(Ydata[idx(i, j, Nx, Ny, 1)]+egvmax*Ydata[idx(i, j, Nx, Ny, 0)]);
            yndata[idx(i, j, Nx, Ny, 1)]=0.5*(Ydata[idx(i, j, Nx, Ny, 1)]-egvmax*Ydata[idx(i, j, Nx, Ny, 0)]);
            ypdata[idx(i, j, Nx, Ny, 2)]=0.5*(Ydata[idx(i, j, Nx, Ny, 2)]+egvmax*Ydata[idx(i, j, Nx, Ny, 0)]);
            yndata[idx(i, j, Nx, Ny, 2)]=0.5*(Ydata[idx(i, j, Nx, Ny, 2)]-egvmax*Ydata[idx(i, j, Nx, Ny, 0)]);
            taopdata[idx(i, j, Nx, Ny, 0)]=0.5*(taodata[idx(i, j, Nx, Ny, 0)]+egvmax*Ydata[idx(i, j, Nx, Ny, 1)]);
            taondata[idx(i, j, Nx, Ny, 0)]=0.5*(taodata[idx(i, j, Nx, Ny, 0)]-egvmax*Ydata[idx(i, j, Nx, Ny, 1)]);
            taopdata[idx(i, j, Nx, Ny, 1)]=0.5*(taodata[idx(i, j, Nx, Ny, 1)]+egvmax*Ydata[idx(i, j, Nx, Ny, 1)]);
            taondata[idx(i, j, Nx, Ny, 1)]=0.5*(taodata[idx(i, j, Nx, Ny, 1)]-egvmax*Ydata[idx(i, j, Nx, Ny, 1)]);
            taopdata[idx(i, j, Nx, Ny, 2)]=0.5*(taodata[idx(i, j, Nx, Ny, 2)]+egvmax*Ydata[idx(i, j, Nx, Ny, 2)]);
            taondata[idx(i, j, Nx, Ny, 2)]=0.5*(taodata[idx(i, j, Nx, Ny, 2)]-egvmax*Ydata[idx(i, j, Nx, Ny, 2)]);
            taopdata[idx(i, j, Nx, Ny, 3)]=0.5*(taodata[idx(i, j, Nx, Ny, 3)]+egvmax*Ydata[idx(i, j, Nx, Ny, 2)]);
            taondata[idx(i, j, Nx, Ny, 3)]=0.5*(taodata[idx(i, j, Nx, Ny, 3)]-egvmax*Ydata[idx(i, j, Nx, Ny, 2)]);
            Cjpdata[idx(i, j, Nx, Ny, 0)]=0.5*(Cjdata[idx(i, j, Nx, Ny, 0)]+egvmax*Ydata[idx(i, j, Nx, Ny, 3)]);
            Cjndata[idx(i, j, Nx, Ny, 0)]=0.5*(Cjdata[idx(i, j, Nx, Ny, 0)]-egvmax*Ydata[idx(i, j, Nx, Ny, 3)]);
            Cjpdata[idx(i, j, Nx, Ny, 1)]=0.5*(Cjdata[idx(i, j, Nx, Ny, 1)]+egvmax*Ydata[idx(i, j, Nx, Ny, 3)]);
            Cjndata[idx(i, j, Nx, Ny, 1)]=0.5*(Cjdata[idx(i, j, Nx, Ny, 1)]-egvmax*Ydata[idx(i, j, Nx, Ny, 3)]);
        }
    }

    /* iterate over domain, computing all equations */
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
        
            /* get derivative on x direction */
            flag = SetISX(IS0_px, IS1_px, IS2_px, IS0_nx, IS1_nx, IS2_nx, ypdata, yndata, taopdata, taondata, Cjpdata, Cjndata, i, j, Nx, Ny);
            if (flag!=0) printf("error in SetISX function \n");
            
            flag = Setalphawx(IS0_px, IS1_px, IS2_px, IS0_nx, IS1_nx, IS2_nx, alpha_0px, alpha_1px, alpha_2px, alpha_0nx, alpha_1nx, alpha_2nx, w0_px, w1_px, w2_px, w0_nx, w1_nx, w2_nx, Epsilon);
            if (flag!=0) printf("error in Setalphawx function \n");
            
            flag = SetUx(w0_px, w1_px, w2_px, w0_nx, w1_nx, w2_nx, ypdata, yndata, taopdata, taondata, Cjpdata, Cjndata, u_tpphx, u_tnphx, u_tpnhx, u_tnnhx, i, j, Nx, Ny);
            if (flag!=0) printf("error in SetUx function \n");
            
            for(k=0;k<4;k++){
                yxdata[idx(i, j, Nx, Ny, k)]=-(1/dx)*((u_tpphx[k]-u_tpnhx[k])+(u_tnphx[k]-u_tnnhx[k]));
            }
            
            
            /* get derivative on y direction */
            flag = SetISY(IS0_py, IS1_py, IS2_py, IS0_ny, IS1_ny, IS2_ny, ypdata, yndata, taopdata, taondata, Cjpdata, Cjndata, i, j, Nx, Ny);
            if (flag!=0) printf("error in SetISY function \n");
            
            flag = Setalphawy(IS0_py, IS1_py, IS2_py, IS0_ny, IS1_ny, IS2_ny, alpha_0py, alpha_1py, alpha_2py, alpha_0ny, alpha_1ny, alpha_2ny, w0_py, w1_py, w2_py, w0_ny, w1_ny, w2_ny, Epsilon);
            if (flag!=0) printf("error in Setalphawy function \n");
            
            flag = SetUy(w0_py, w1_py, w2_py, w0_ny, w1_ny, w2_ny, ypdata, yndata, taopdata, taondata, Cjpdata, Cjndata, u_tpphy, u_tnphy, u_tpnhy, u_tnnhy, i, j, Nx, Ny);
            if (flag!=0) printf("error in SetUy function \n");
            
            for(k=0;k<4;k++){
                yydata[idx(i, j, Nx, Ny, k)]=-(1/dy)*((u_tpphy[k]-u_tpnhy[k])+(u_tnphy[k]-u_tnnhy[k]));
            }
    
            /* get derivative both on x and y direction */
            for (k=0; k<4; k++){
                dYdata[idx(i, j, Nx, Ny, k)]=yxdata[idx(i, j, Nx, Ny, k)]+yydata[idx(i, j, Nx, Ny, k)];
            }
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
    
    /* Free vectors */
    N_VDestroy_Serial(yp);
    N_VDestroy_Serial(yn);
    N_VDestroy_Serial(tao);
    N_VDestroy_Serial(Cj);
    N_VDestroy_Serial(yx);
    N_VDestroy_Serial(yy);
    N_VDestroy_Serial(taop);
    N_VDestroy_Serial(taon);
    N_VDestroy_Serial(Cjp);
    N_VDestroy_Serial(Cjn);
    
    return 0;
}

/* Fill in the indicators of smoothness on x direction */
static int SetISX(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx,realtype *ypdata, realtype *yndata, realtype *taopdata, realtype *taondata, realtype *Cjpdata, realtype *Cjndata, long int i, long int j, long int Nx, long int Ny)
{
  /* declaration */
  long int idown1, idown2, idown3, iup1, iup2, iup3;

    /* consider Nx */
    if (Nx<4){
        cerr << "\nNx should be more than 3!\n";
        return 1;
    }
    else {
      /* consider the situations of different i */  
      idown1 = i;
      idown2 = i;
      idown3 = i;
      iup1 = i;
      iup2 = i;
      iup3 = i;
      if (i<3||i>Nx-4){      
	if (i==0){
	  idown1 = Nx;
	  idown2 = Nx;
	  idown3 = Nx;
	}
    
	if (i==1){ 
	  idown2 = Nx+1;
	  idown3 = Nx+1;
	}
    
	if (i==2){
	  idown3 = Nx+2;
	}
    
	if (i==Nx-3){
	  iup3 = -3;
	}
    
	if (i==Nx-2){
	  iup2 = -2;
	  iup3 = -2;
	}
    
	if (i==Nx-1){
	  iup1 = -1;
	  iup2 = -1;
	  iup3 = -1;
	}
      }
      else {
    /* compute indicators of smoothness of rou, qx, qy, E */ 
    IS0_px[0]=(13.0/12.0)*(ypdata[idx(idown2-2, j, Nx, Ny, 1)]-2.0*ypdata[idx(idown1-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(idown2-2, j, Nx, Ny, 1)]-2.0*ypdata[idx(idown1-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])+(1.0/4.0)*(ypdata[idx(idown2-2, j, Nx, Ny, 1)]-4.0*ypdata[idx(idown1-1, j, Nx, Ny, 1)]+3.0*ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(idown2-2, j, Nx, Ny, 1)]-4.0*ypdata[idx(idown1-1, j, Nx, Ny, 1)]+3.0*ypdata[idx(i, j, Nx, Ny, 1)]);
    
    IS1_px[0]=(13.0/12.0)*(ypdata[idx(idown1-1, j, Nx, Ny, 1)]-2.0*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(iup1+1, j, Nx, Ny, 1)])*(ypdata[idx(idown1-1, j, Nx, Ny, 1)]-2.0*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(iup1+1, j, Nx, Ny, 1)])+(1.0/4.0)*(ypdata[idx(idown1-1, j, Nx, Ny, 1)]-ypdata[idx(iup1+1, j, Nx, Ny, 1)])*(ypdata[idx(idown1-1, j, Nx, Ny, 1)]-ypdata[idx(iup1+1, j, Nx, Ny, 1)]);
    
    IS2_px[0]=(13.0/12.0)*(ypdata[idx(i, j, Nx, Ny, 1)]-2.0*ypdata[idx(iup1+1, j, Nx, Ny, 1)]+ypdata[idx(iup2+2, j, Nx, Ny, 1)])*(ypdata[idx(i, j, Nx, Ny, 1)]-2.0*ypdata[idx(iup1+1, j, Nx, Ny, 1)]+ypdata[idx(iup2+2, j, Nx, Ny, 1)])+(1.0/4.0)*(3.0*ypdata[idx(i, j, Nx, Ny, 1)]-4.0*ypdata[idx(iup1+1, j, Nx, Ny, 1)]+ypdata[idx(iup2+2, j, Nx, Ny, 1)])*(3.0*ypdata[idx(i, j, Nx, Ny, 1)]-4.0*ypdata[idx(iup1+1, j, Nx, Ny, 1)]+ypdata[idx(iup2+2, j, Nx, Ny, 1)]);
    
    IS0_nx[0]=(13.0/12.0)*(yndata[idx(iup1+1, j, Nx, Ny, 1)]-2.0*yndata[idx(iup2+2, j, Nx, Ny, 1)]+yndata[idx(iup3+3, j, Nx, Ny, 1)])*(yndata[idx(iup1+1, j, Nx, Ny, 1)]-2.0*yndata[idx(iup2+2, j, Nx, Ny, 1)]+yndata[idx(iup3+3, j, Nx, Ny, 1)])+(1.0/4.0)*(3.0*yndata[idx(iup1+1, j, Nx, Ny, 1)]-4.0*yndata[idx(iup2+2, j, Nx, Ny, 1)]+yndata[idx(iup3+3, j, Nx, Ny, 1)])*(3.0*yndata[idx(iup1+1, j, Nx, Ny, 1)]-4.0*yndata[idx(iup2+2, j, Nx, Ny, 1)]+yndata[idx(iup3+3, j, Nx, Ny, 1)]);
    
    IS1_nx[0]=(13.0/12.0)*(yndata[idx(i, j, Nx, Ny, 1)]-2.0*yndata[idx(iup1+1, j, Nx, Ny, 1)]+yndata[idx(iup2+2, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-2.0*yndata[idx(iup1+1, j, Nx, Ny, 1)]+yndata[idx(iup2+2, j, Nx, Ny, 1)])+(1.0/4.0)*(yndata[idx(i, j, Nx, Ny, 1)]-yndata[idx(iup2+2, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-yndata[idx(iup2+2, j, Nx, Ny, 1)]);
    
    IS2_nx[0]=(13.0/12.0)*(yndata[idx(idown1-1, j, Nx, Ny, 1)]-2.0*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(iup1+1, j, Nx, Ny, 1)])*(yndata[idx(idown1-1, j, Nx, Ny, 1)]-2.0*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(iup1+1, j, Nx, Ny, 1)])+(1.0/4.0)*(yndata[idx(idown1-1, j, Nx, Ny, 1)]-4.0*yndata[idx(i, j, Nx, Ny, 1)]+3.0*yndata[idx(iup1+1, j, Nx, Ny, 1)])*(yndata[idx(idown1-1, j, Nx, Ny, 1)]-4.0*yndata[idx(i, j, Nx, Ny, 1)]+3.0*yndata[idx(iup1+1, j, Nx, Ny, 1)]);
    
    
    IS0_px[1]=(13.0/12.0)*(taopdata[idx(idown2-2, j, Nx, Ny, 0)]-2.0*taopdata[idx(idown1-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(idown2-2, j, Nx, Ny, 0)]-2.0*taopdata[idx(idown1-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])+(1.0/4.0)*(taopdata[idx(idown2-2, j, Nx, Ny, 0)]-4.0*taopdata[idx(idown1-1, j, Nx, Ny, 0)]+3.0*taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(idown2-2, j, Nx, Ny, 0)]-4.0*taopdata[idx(idown1-1, j, Nx, Ny, 0)]+3.0*taopdata[idx(i, j, Nx, Ny, 0)]);
    
    IS1_px[1]=(13.0/12.0)*(taopdata[idx(idown1-1, j, Nx, Ny, 0)]-2.0*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(iup1+1, j, Nx, Ny, 0)])*(taopdata[idx(idown1-1, j, Nx, Ny, 0)]-2.0*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(iup1+1, j, Nx, Ny, 0)])+(1.0/4.0)*(taopdata[idx(idown1-1, j, Nx, Ny, 0)]-taopdata[idx(iup1+1, j, Nx, Ny, 0)])*(taopdata[idx(idown1-1, j, Nx, Ny, 0)]-taopdata[idx(iup1+1, j, Nx, Ny, 0)]);
    
    IS2_px[1]=(13.0/12.0)*(taopdata[idx(i, j, Nx, Ny, 0)]-2.0*taopdata[idx(iup1+1, j, Nx, Ny, 0)]+taopdata[idx(iup2+2, j, Nx, Ny, 0)])*(taopdata[idx(i, j, Nx, Ny, 0)]-2.0*taopdata[idx(iup1+1, j, Nx, Ny, 0)]+taopdata[idx(iup2+2, j, Nx, Ny, 0)])+(1.0/4.0)*(3.0*taopdata[idx(i, j, Nx, Ny, 0)]-4.0*taopdata[idx(iup1+1, j, Nx, Ny, 0)]+taopdata[idx(iup2+2, j, Nx, Ny, 0)])*(3.0*taopdata[idx(i, j, Nx, Ny, 0)]-4.0*taopdata[idx(iup1+1, j, Nx, Ny, 0)]+taopdata[idx(iup2+2, j, Nx, Ny, 0)]);
    
    IS0_nx[1]=(13.0/12.0)*(taondata[idx(iup1+1, j, Nx, Ny, 0)]-2.0*taondata[idx(iup2+2, j, Nx, Ny, 0)]+taondata[idx(iup3+3, j, Nx, Ny, 0)])*(taondata[idx(iup1+1, j, Nx, Ny, 0)]-2.0*taondata[idx(iup2+2, j, Nx, Ny, 0)]+taondata[idx(iup3+3, j, Nx, Ny, 0)])+(1.0/4.0)*(3.0*taondata[idx(iup1+1, j, Nx, Ny, 0)]-4.0*taondata[idx(iup2+2, j, Nx, Ny, 0)]+taondata[idx(iup3+3, j, Nx, Ny, 0)])*(3.0*taondata[idx(iup1+1, j, Nx, Ny, 0)]-4.0*taondata[idx(iup2+2, j, Nx, Ny, 0)]+taondata[idx(iup3+3, j, Nx, Ny, 0)]);
    
    IS1_nx[1]=(13.0/12.0)*(taondata[idx(i, j, Nx, Ny, 0)]-2.0*taondata[idx(iup1+1, j, Nx, Ny, 0)]+taondata[idx(iup2+2, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-2.0*taondata[idx(iup1+1, j, Nx, Ny, 0)]+taondata[idx(iup2+2, j, Nx, Ny, 0)])+(1.0/4.0)*(taondata[idx(i, j, Nx, Ny, 0)]-taondata[idx(iup2+2, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-taondata[idx(iup2+2, j, Nx, Ny, 0)]);
    
    IS2_nx[1]=(13.0/12.0)*(taondata[idx(idown1-1, j, Nx, Ny, 0)]-2.0*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(iup1+1, j, Nx, Ny, 0)])*(taondata[idx(idown1-1, j, Nx, Ny, 0)]-2.0*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(iup1+1, j, Nx, Ny, 0)])+(1.0/4.0)*(taondata[idx(idown1-1, j, Nx, Ny, 0)]-4.0*taondata[idx(i, j, Nx, Ny, 0)]+3.0*taondata[idx(iup1+1, j, Nx, Ny, 0)])*(taondata[idx(idown1-1, j, Nx, Ny, 0)]-4.0*taondata[idx(i, j, Nx, Ny, 0)]+3.0*taondata[idx(iup1+1, j, Nx, Ny, 0)]);
    
    
    IS0_px[2]=(13.0/12.0)*(taopdata[idx(idown2-2, j, Nx, Ny, 2)]-2.0*taopdata[idx(idown1-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(idown2-2, j, Nx, Ny, 2)]-2.0*taopdata[idx(idown1-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])+(1.0/4.0)*(taopdata[idx(idown2-2, j, Nx, Ny, 2)]-4.0*taopdata[idx(idown1-1, j, Nx, Ny, 2)]+3.0*taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(idown2-2, j, Nx, Ny, 2)]-4.0*taopdata[idx(idown1-1, j, Nx, Ny, 2)]+3.0*taopdata[idx(i, j, Nx, Ny, 2)]);
    
    IS1_px[2]=(13.0/12.0)*(taopdata[idx(idown1-1, j, Nx, Ny, 2)]-2.0*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(iup1+1, j, Nx, Ny, 2)])*(taopdata[idx(idown1-1, j, Nx, Ny, 2)]-2.0*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(iup1+1, j, Nx, Ny, 2)])+(1.0/4.0)*(taopdata[idx(idown1-1, j, Nx, Ny, 2)]-taopdata[idx(iup1+1, j, Nx, Ny, 2)])*(taopdata[idx(idown1-1, j, Nx, Ny, 2)]-taopdata[idx(iup1+1, j, Nx, Ny, 2)]);
    
    IS2_px[2]=(13.0/12.0)*(taopdata[idx(i, j, Nx, Ny, 2)]-2.0*taopdata[idx(iup1+1, j, Nx, Ny, 2)]+taopdata[idx(iup2+2, j, Nx, Ny, 2)])*(taopdata[idx(i, j, Nx, Ny, 2)]-2.0*taopdata[idx(iup1+1, j, Nx, Ny, 2)]+taopdata[idx(iup2+2, j, Nx, Ny, 2)])+(1.0/4.0)*(3.0*taopdata[idx(i, j, Nx, Ny, 2)]-4.0*taopdata[idx(iup1+1, j, Nx, Ny, 2)]+taopdata[idx(iup2+2, j, Nx, Ny, 2)])*(3.0*taopdata[idx(i, j, Nx, Ny, 2)]-4.0*taopdata[idx(iup1+1, j, Nx, Ny, 2)]+taopdata[idx(iup2+2, j, Nx, Ny, 2)]);
    
    IS0_nx[2]=(13.0/12.0)*(taondata[idx(iup1+1, j, Nx, Ny, 2)]-2.0*taondata[idx(iup2+2, j, Nx, Ny, 2)]+taondata[idx(iup3+3, j, Nx, Ny, 2)])*(taondata[idx(iup1+1, j, Nx, Ny, 2)]-2.0*taondata[idx(iup2+2, j, Nx, Ny, 2)]+taondata[idx(iup3+3, j, Nx, Ny, 2)])+(1.0/4.0)*(3.0*taondata[idx(iup1+1, j, Nx, Ny, 2)]-4.0*taondata[idx(iup2+2, j, Nx, Ny, 2)]+taondata[idx(iup3+3, j, Nx, Ny, 2)])*(3.0*taondata[idx(iup1+1, j, Nx, Ny, 2)]-4.0*taondata[idx(iup2+2, j, Nx, Ny, 2)]+taondata[idx(iup3+3, j, Nx, Ny, 2)]);
    
    IS1_nx[2]=(13.0/12.0)*(taondata[idx(i, j, Nx, Ny, 2)]-2.0*taondata[idx(iup1+1, j, Nx, Ny, 2)]+taondata[idx(iup2+2, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-2.0*taondata[idx(iup1+1, j, Nx, Ny, 2)]+taondata[idx(iup2+2, j, Nx, Ny, 2)])+(1.0/4.0)*(taondata[idx(i, j, Nx, Ny, 2)]-taondata[idx(iup2+2, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-taondata[idx(iup2+2, j, Nx, Ny, 2)]);
    
    IS2_nx[2]=(13.0/12.0)*(taondata[idx(idown1-1, j, Nx, Ny, 2)]-2.0*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(iup1+1, j, Nx, Ny, 2)])*(taondata[idx(idown1-1, j, Nx, Ny, 2)]-2.0*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(iup1+1, j, Nx, Ny, 2)])+(1.0/4.0)*(taondata[idx(idown1-1, j, Nx, Ny, 2)]-4.0*taondata[idx(i, j, Nx, Ny, 2)]+3.0*taondata[idx(iup1+1, j, Nx, Ny, 2)])*(taondata[idx(idown1-1, j, Nx, Ny, 2)]-4.0*taondata[idx(i, j, Nx, Ny, 2)]+3.0*taondata[idx(iup1+1, j, Nx, Ny, 2)]);
    
    
    IS0_px[3]=(13.0/12.0)*(Cjpdata[idx(idown2-2, j, Nx, Ny, 0)]-2.0*Cjpdata[idx(idown1-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(idown2-2, j, Nx, Ny, 0)]-2.0*Cjpdata[idx(idown1-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])+(1.0/4.0)*(Cjpdata[idx(idown2-2, j, Nx, Ny, 0)]-4.0*Cjpdata[idx(idown1-1, j, Nx, Ny, 0)]+3.0*Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(idown2-2, j, Nx, Ny, 0)]-4.0*Cjpdata[idx(idown1-1, j, Nx, Ny, 0)]+3.0*Cjpdata[idx(i, j, Nx, Ny, 0)]);
    
    IS1_px[3]=(13.0/12.0)*(Cjpdata[idx(idown1-1, j, Nx, Ny, 0)]-2.0*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(iup1+1, j, Nx, Ny, 0)])*(Cjpdata[idx(idown1-1, j, Nx, Ny, 0)]-2.0*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(iup1+1, j, Nx, Ny, 0)])+(1.0/4.0)*(Cjpdata[idx(idown1-1, j, Nx, Ny, 0)]-Cjpdata[idx(iup1+1, j, Nx, Ny, 0)])*(Cjpdata[idx(idown1-1, j, Nx, Ny, 0)]-Cjpdata[idx(iup1+1, j, Nx, Ny, 0)]);
    
    IS2_px[3]=(13.0/12.0)*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2.0*Cjpdata[idx(iup1+1, j, Nx, Ny, 0)]+Cjpdata[idx(iup2+2, j, Nx, Ny, 0)])*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2.0*Cjpdata[idx(iup1+1, j, Nx, Ny, 0)]+Cjpdata[idx(iup2+2, j, Nx, Ny, 0)])+(1.0/4.0)*(3.0*Cjpdata[idx(i, j, Nx, Ny, 0)]-4.0*Cjpdata[idx(iup1+1, j, Nx, Ny, 0)]+Cjpdata[idx(iup2+2, j, Nx, Ny, 0)])*(3.0*Cjpdata[idx(i, j, Nx, Ny, 0)]-4.0*Cjpdata[idx(iup1+1, j, Nx, Ny, 0)]+Cjpdata[idx(iup2+2, j, Nx, Ny, 0)]);
    
    IS0_nx[3]=(13.0/12.0)*(Cjndata[idx(iup1+1, j, Nx, Ny, 0)]-2.0*Cjndata[idx(iup2+2, j, Nx, Ny, 0)]+Cjndata[idx(iup3+3, j, Nx, Ny, 0)])*(Cjndata[idx(iup1+1, j, Nx, Ny, 0)]-2.0*Cjndata[idx(iup2+2, j, Nx, Ny, 0)]+Cjndata[idx(iup3+3, j, Nx, Ny, 0)])+(1.0/4.0)*(3.0*Cjndata[idx(iup1+1, j, Nx, Ny, 0)]-4.0*Cjndata[idx(iup2+2, j, Nx, Ny, 0)]+Cjndata[idx(iup3+3, j, Nx, Ny, 0)])*(3.0*Cjndata[idx(iup1+1, j, Nx, Ny, 0)]-4.0*Cjndata[idx(iup2+2, j, Nx, Ny, 0)]+Cjndata[idx(iup3+3, j, Nx, Ny, 0)]);
    
    IS1_nx[3]=(13.0/12.0)*(Cjndata[idx(i, j, Nx, Ny, 0)]-2.0*Cjndata[idx(iup1+1, j, Nx, Ny, 0)]+Cjndata[idx(iup2+2, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-2.0*Cjndata[idx(iup1+1, j, Nx, Ny, 0)]+Cjndata[idx(iup2+2, j, Nx, Ny, 0)])+(1.0/4.0)*(Cjndata[idx(i, j, Nx, Ny, 0)]-Cjndata[idx(iup2+2, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-Cjndata[idx(iup2+2, j, Nx, Ny, 0)]);
    
    IS2_nx[3]=(13.0/12.0)*(Cjndata[idx(idown1-1, j, Nx, Ny, 0)]-2.0*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(iup1+1, j, Nx, Ny, 0)])*(Cjndata[idx(idown1-1, j, Nx, Ny, 0)]-2.0*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(iup1+1, j, Nx, Ny, 0)])+(1.0/4.0)*(Cjndata[idx(idown1-1, j, Nx, Ny, 0)]-4.0*Cjndata[idx(i, j, Nx, Ny, 0)]+3.0*Cjndata[idx(iup1+1, j, Nx, Ny, 0)])*(Cjndata[idx(idown1-1, j, Nx, Ny, 0)]-4.0*Cjndata[idx(i, j, Nx, Ny, 0)]+3.0*Cjndata[idx(iup1+1, j, Nx, Ny, 0)]);
   
      }
    return 0;
    }
}

/* Fill in the indicators of smoothness on y direction */
static int SetISY(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny,realtype *ypdata, realtype *yndata, realtype *taopdata, realtype *taondata, realtype *Cjpdata, realtype *Cjndata, long int i, long int j, long int Nx, long int Ny)
{
  /* declaration */
  long int jdown1, jdown2, jdown3, jup1, jup2, jup3;

    /* consider Ny */
    if (Ny<4){
        cerr << "\nNy should be more than 3!\n";
        return 1;
    }
    else {
      /* consider situations of different j */
      jdown1 = j;
      jdown2 = j;
      jdown3 = j;
      jup1 = j;
      jup2 = j;
      jup3 = j;
      if (j<3||j>Ny-4){
        if (j==0){
            jdown1 = Ny;
            jdown2 = Ny;
            jdown3 = Ny;
        }
        
        if (j==1){
            jdown2 = Ny+1;
            jdown3 = Ny+1;
        }
        
        if (j==2){
            jdown3 = Ny+2;
        }
        
        if (j==Ny-3){
            jup3 = -3;
        }
        
        if (j==Ny-2){
            jup2 = -2;
            jup3 = -2;
        }
        
        if (j==Ny-1){
            jup1 = -1;
            jup2 = -1;
            jup3 = -1;
        }
      }
      else
	{
	/*compute indicators of smoothness of rou, qx, qy, E */
    IS0_py[0]=(13.0/12.0)*(ypdata[idx(i, jdown2-2, Nx, Ny, 2)]-2.0*ypdata[idx(i, jdown1-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, jdown2-2, Nx, Ny, 2)]-2.0*ypdata[idx(i, jdown1-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])+(1.0/4.0)*(ypdata[idx(i, jdown2-2, Nx, Ny, 2)]-4.0*ypdata[idx(i, jdown1-1, Nx, Ny, 2)]+3.0*ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, jdown2-2, Nx, Ny, 2)]-4.0*ypdata[idx(i, jdown1-1, Nx, Ny, 2)]+3.0*ypdata[idx(i, j, Nx, Ny, 2)]);
    
    IS1_py[0]=(13.0/12.0)*(ypdata[idx(i, jdown1-1, Nx, Ny, 2)]-2.0*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, jup1+1, Nx, Ny, 2)])*(ypdata[idx(i, jdown1-1, Nx, Ny, 2)]-2.0*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, jup1+1, Nx, Ny, 2)])+(1.0/4.0)*(ypdata[idx(i, jdown1-1, Nx, Ny, 2)]-ypdata[idx(i, jup1+1, Nx, Ny, 2)])*(ypdata[idx(i, jdown1-1, Nx, Ny, 2)]-ypdata[idx(i, jup1+1, Nx, Ny, 2)]);
    
    IS2_py[0]=(13.0/12.0)*(ypdata[idx(i, j, Nx, Ny, 2)]-2.0*ypdata[idx(i, jup1+1, Nx, Ny, 2)]+ypdata[idx(i, jup2+2, Nx, Ny, 2)])*(ypdata[idx(i, j, Nx, Ny, 2)]-2.0*ypdata[idx(i, jup1+1, Nx, Ny, 2)]+ypdata[idx(i, jup2+2, Nx, Ny, 2)])+(1.0/4.0)*(3.0*ypdata[idx(i, j, Nx, Ny, 2)]-4.0*ypdata[idx(i, jup1+1, Nx, Ny, 2)]+ypdata[idx(i, jup2+2, Nx, Ny, 2)])*(3.0*ypdata[idx(i, j, Nx, Ny, 2)]-4.0*ypdata[idx(i, jup1+1, Nx, Ny, 2)]+ypdata[idx(i, jup2+2, Nx, Ny, 2)]);
    
    IS0_ny[0]=(13.0/12.0)*(yndata[idx(i, jup1+1, Nx, Ny, 2)]-2.0*yndata[idx(i, jup2+2, Nx, Ny, 2)]+yndata[idx(i, jup3+3, Nx, Ny, 2)])*(yndata[idx(i, jup1+1, Nx, Ny, 2)]-2.0*yndata[idx(i, jup2+2, Nx, Ny, 2)]+yndata[idx(i, jup3+3, Nx, Ny, 2)])+(1.0/4.0)*(3.0*yndata[idx(i, jup1+1, Nx, Ny, 2)]-4.0*yndata[idx(i, jup2+2, Nx, Ny, 2)]+yndata[idx(i, jup3+3, Nx, Ny, 2)])*(3.0*yndata[idx(i, jup1+1, Nx, Ny, 2)]-4.0*yndata[idx(i, jup2+2, Nx, Ny, 2)]+yndata[idx(i, jup3+3, Nx, Ny, 2)]);
    
    IS1_ny[0]=(13.0/12.0)*(yndata[idx(i, j, Nx, Ny, 2)]-2.0*yndata[idx(i, jup1+1, Nx, Ny, 2)]+yndata[idx(i, jup2+2, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-2.0*yndata[idx(i, jup1+1, Nx, Ny, 2)]+yndata[idx(i, jup2+2, Nx, Ny, 2)])+(1.0/4.0)*(yndata[idx(i, j, Nx, Ny, 2)]-yndata[idx(i, jup2+2, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-yndata[idx(i, jup2+2, Nx, Ny, 2)]);
    
    IS2_ny[0]=(13.0/12.0)*(yndata[idx(i, jdown1-1, Nx, Ny, 2)]-2.0*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, jup1+1, Nx, Ny, 2)])*(yndata[idx(i, jdown1-1, Nx, Ny, 2)]-2.0*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, jup1+1, Nx, Ny, 2)])+(1.0/4.0)*(yndata[idx(i, jdown1-1, Nx, Ny, 2)]-4.0*yndata[idx(i, j, Nx, Ny, 2)]+3.0*yndata[idx(i, jup1+1, Nx, Ny, 2)])*(yndata[idx(i, jdown1-1, Nx, Ny, 2)]-4.0*yndata[idx(i, j, Nx, Ny, 2)]+3.0*yndata[idx(i, jup1+1, Nx, Ny, 2)]);
    
    
    IS0_py[1]=(13.0/12.0)*(taopdata[idx(i, jdown2-2, Nx, Ny, 1)]-2.0*taopdata[idx(i, jdown1-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, jdown2-2, Nx, Ny, 1)]-2.0*taopdata[idx(i, jdown1-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])+(1.0/4.0)*(taopdata[idx(i, jdown2-2, Nx, Ny, 1)]-4.0*taopdata[idx(i, jdown1-1, Nx, Ny, 1)]+3.0*taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, jdown2-2, Nx, Ny, 1)]-4.0*taopdata[idx(i, jdown1-1, Nx, Ny, 1)]+3.0*taopdata[idx(i, j, Nx, Ny, 1)]);
    
    IS1_py[1]=(13.0/12.0)*(taopdata[idx(i, jdown1-1, Nx, Ny, 1)]-2.0*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, jup1+1, Nx, Ny, 1)])*(taopdata[idx(i, jdown1-1, Nx, Ny, 1)]-2.0*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, jup1+1, Nx, Ny, 1)])+(1.0/4.0)*(taopdata[idx(i, jdown1-1, Nx, Ny, 1)]-taopdata[idx(i, jup1+1, Nx, Ny, 1)])*(taopdata[idx(i, jdown1-1, Nx, Ny, 1)]-taopdata[idx(i, jup1+1, Nx, Ny, 1)]);
    
    IS2_py[1]=(13.0/12.0)*(taopdata[idx(i, j, Nx, Ny, 1)]-2.0*taopdata[idx(i, jup1+1, Nx, Ny, 1)]+taopdata[idx(i, jup2+2, Nx, Ny, 1)])*(taopdata[idx(i, j, Nx, Ny, 1)]-2.0*taopdata[idx(i, jup1+1, Nx, Ny, 1)]+taopdata[idx(i, jup2+2, Nx, Ny, 1)])+(1.0/4.0)*(3.0*taopdata[idx(i, j, Nx, Ny, 1)]-4.0*taopdata[idx(i, jup1+1, Nx, Ny, 1)]+taopdata[idx(i, jup2+2, Nx, Ny, 1)])*(3.0*taopdata[idx(i, j, Nx, Ny, 1)]-4.0*taopdata[idx(i, jup1+1, Nx, Ny, 1)]+taopdata[idx(i, jup2+2, Nx, Ny, 1)]);
    
    IS0_ny[1]=(13.0/12.0)*(taondata[idx(i, jup1+1, Nx, Ny, 1)]-2.0*taondata[idx(i, jup2+2, Nx, Ny, 1)]+taondata[idx(i, jup3+3, Nx, Ny, 1)])*(taondata[idx(i, jup1+1, Nx, Ny, 1)]-2.0*taondata[idx(i, jup2+2, Nx, Ny, 1)]+taondata[idx(i, jup3+3, Nx, Ny, 1)])+(1.0/4.0)*(3.0*taondata[idx(i, jup1+1, Nx, Ny, 1)]-4.0*taondata[idx(i, jup2+2, Nx, Ny, 1)]+taondata[idx(i, jup3+3, Nx, Ny, 1)])*(3.0*taondata[idx(i, jup1+1, Nx, Ny, 1)]-4.0*taondata[idx(i, jup2+2, Nx, Ny, 1)]+taondata[idx(i, jup3+3, Nx, Ny, 1)]);
    
    IS1_ny[1]=(13.0/12.0)*(taondata[idx(i, j, Nx, Ny, 1)]-2.0*taondata[idx(i, jup1+1, Nx, Ny, 1)]+taondata[idx(i, jup2+2, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-2.0*taondata[idx(i, jup1+1, Nx, Ny, 1)]+taondata[idx(i, jup2+2, Nx, Ny, 1)])+(1.0/4.0)*(taondata[idx(i, j, Nx, Ny, 1)]-taondata[idx(i, jup2+2, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-taondata[idx(i, jup2+2, Nx, Ny, 1)]);
    
    IS2_ny[1]=(13.0/12.0)*(taondata[idx(i, jdown1-1, Nx, Ny, 1)]-2.0*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, jup1+1, Nx, Ny, 1)])*(taondata[idx(i, jdown1-1, Nx, Ny, 1)]-2.0*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, jup1+1, Nx, Ny, 1)])+(1.0/4.0)*(taondata[idx(i, jdown1-1, Nx, Ny, 1)]-4.0*taondata[idx(i, j, Nx, Ny, 1)]+3.0*taondata[idx(i, jup1+1, Nx, Ny, 1)])*(taondata[idx(i, jdown1-1, Nx, Ny, 1)]-4.0*taondata[idx(i, j, Nx, Ny, 1)]+3.0*taondata[idx(i, jup1+1, Nx, Ny, 1)]);
    
    
    IS0_py[2]=(13.0/12.0)*(taopdata[idx(i, jdown2-2, Nx, Ny, 3)]-2.0*taopdata[idx(i, jdown1-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, jdown2-2, Nx, Ny, 3)]-2.0*taopdata[idx(i, jdown1-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])+(1.0/4.0)*(taopdata[idx(i, jdown2-2, Nx, Ny, 3)]-4.0*taopdata[idx(i, jdown1-1, Nx, Ny, 3)]+3.0*taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, jdown2-2, Nx, Ny, 3)]-4.0*taopdata[idx(i, jdown1-1, Nx, Ny, 3)]+3.0*taopdata[idx(i, j, Nx, Ny, 3)]);
    
    IS1_py[2]=(13.0/12.0)*(taopdata[idx(i, jdown1-1, Nx, Ny, 3)]-2.0*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, jup1+1, Nx, Ny, 3)])*(taopdata[idx(i, jdown1-1, Nx, Ny, 3)]-2.0*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, jup1+1, Nx, Ny, 3)])+(1.0/4.0)*(taopdata[idx(i, jdown1-1, Nx, Ny, 3)]-taopdata[idx(i, jup1+1, Nx, Ny, 3)])*(taopdata[idx(i, jdown1-1, Nx, Ny, 3)]-taopdata[idx(i, jup1+1, Nx, Ny, 3)]);
    
    IS2_py[2]=(13.0/12.0)*(taopdata[idx(i, j, Nx, Ny, 3)]-2.0*taopdata[idx(i, jup1+1, Nx, Ny, 3)]+taopdata[idx(i, jup2+2, Nx, Ny, 3)])*(taopdata[idx(i, j, Nx, Ny, 3)]-2.0*taopdata[idx(i, jup1+1, Nx, Ny, 3)]+taopdata[idx(i, jup2+2, Nx, Ny, 3)])+(1.0/4.0)*(3.0*taopdata[idx(i, j, Nx, Ny, 3)]-4.0*taopdata[idx(i, jup1+1, Nx, Ny, 3)]+taopdata[idx(i, jup2+2, Nx, Ny, 3)])*(3.0*taopdata[idx(i, j, Nx, Ny, 3)]-4.0*taopdata[idx(i, jup1+1, Nx, Ny, 3)]+taopdata[idx(i, jup2+2, Nx, Ny, 3)]);
    
    IS0_ny[2]=(13.0/12.0)*(taondata[idx(i, jup1+1, Nx, Ny, 3)]-2.0*taondata[idx(i, jup2+2, Nx, Ny, 3)]+taondata[idx(i, jup3+3, Nx, Ny, 3)])*(taondata[idx(i, jup1+1, Nx, Ny, 3)]-2.0*taondata[idx(i, jup2+2, Nx, Ny, 3)]+taondata[idx(i, jup3+3, Nx, Ny, 3)])+(1.0/4.0)*(3.0*taondata[idx(i, jup1+1, Nx, Ny, 3)]-4.0*taondata[idx(i, jup2+2, Nx, Ny, 3)]+taondata[idx(i, jup3+3, Nx, Ny, 3)])*(3.0*taondata[idx(i, jup1+1, Nx, Ny, 3)]-4.0*taondata[idx(i, jup2+2, Nx, Ny, 3)]+taondata[idx(i, jup3+3, Nx, Ny, 3)]);
    
    IS1_ny[2]=(13.0/12.0)*(taondata[idx(i, j, Nx, Ny, 3)]-2.0*taondata[idx(i, jup1+1, Nx, Ny, 3)]+taondata[idx(i, jup2+2, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-2.0*taondata[idx(i, jup1+1, Nx, Ny, 3)]+taondata[idx(i, jup2+2, Nx, Ny, 3)])+(1.0/4.0)*(taondata[idx(i, j, Nx, Ny, 3)]-taondata[idx(i, jup2+2, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-taondata[idx(i, jup2+2, Nx, Ny, 3)]);
    
    IS2_ny[2]=(13.0/12.0)*(taondata[idx(i, jdown1-1, Nx, Ny, 3)]-2.0*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, jup1+1, Nx, Ny, 3)])*(taondata[idx(i, jdown1-1, Nx, Ny, 3)]-2.0*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, jup1+1, Nx, Ny, 3)])+(1.0/4.0)*(taondata[idx(i, jdown1-1, Nx, Ny, 3)]-4.0*taondata[idx(i, j, Nx, Ny, 3)]+3.0*taondata[idx(i, jup1+1, Nx, Ny, 3)])*(taondata[idx(i, jdown1-1, Nx, Ny, 3)]-4.0*taondata[idx(i, j, Nx, Ny, 3)]+3.0*taondata[idx(i, jup1+1, Nx, Ny, 3)]);
    
    
    IS0_py[3]=(13.0/12.0)*(Cjpdata[idx(i, jdown2-2, Nx, Ny, 1)]-2.0*Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, jdown2-2, Nx, Ny, 1)]-2.0*Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])+(1.0/4.0)*(Cjpdata[idx(i, jdown2-2, Nx, Ny, 1)]-4.0*Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)]+3.0*Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, jdown2-2, Nx, Ny, 1)]-4.0*Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)]+3.0*Cjpdata[idx(i, j, Nx, Ny, 1)]);
    
    IS1_py[3]=(13.0/12.0)*(Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)]-2.0*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, jup1+1, Nx, Ny, 1)])*(Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)]-2.0*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, jup1+1, Nx, Ny, 1)])+(1.0/4.0)*(Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)]-Cjpdata[idx(i, jup1+1, Nx, Ny, 1)])*(Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)]-Cjpdata[idx(i, jup1+1, Nx, Ny, 1)]);
    
    IS2_py[3]=(13.0/12.0)*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2.0*Cjpdata[idx(i, jup1+1, Nx, Ny, 1)]+Cjpdata[idx(i, jup2+2, Nx, Ny, 1)])*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2.0*Cjpdata[idx(i, jup1+1, Nx, Ny, 1)]+Cjpdata[idx(i, jup2+2, Nx, Ny, 1)])+(1.0/4.0)*(3.0*Cjpdata[idx(i, j, Nx, Ny, 1)]-4.0*Cjpdata[idx(i, jup1+1, Nx, Ny, 1)]+Cjpdata[idx(i, jup2+2, Nx, Ny, 1)])*(3.0*Cjpdata[idx(i, j, Nx, Ny, 1)]-4.0*Cjpdata[idx(i, jup1+1, Nx, Ny, 1)]+Cjpdata[idx(i, jup2+2, Nx, Ny, 1)]);
    
    IS0_ny[3]=(13.0/12.0)*(Cjndata[idx(i, jup1+1, Nx, Ny, 1)]-2.0*Cjndata[idx(i, jup2+2, Nx, Ny, 1)]+Cjndata[idx(i, jup3+3, Nx, Ny, 1)])*(Cjndata[idx(i, jup1+1, Nx, Ny, 1)]-2.0*Cjndata[idx(i, jup2+2, Nx, Ny, 1)]+Cjndata[idx(i, jup3+3, Nx, Ny, 1)])+(1.0/4.0)*(3.0*Cjndata[idx(i, jup1+1, Nx, Ny, 1)]-4.0*Cjndata[idx(i, jup2+2, Nx, Ny, 1)]+Cjndata[idx(i, jup3+3, Nx, Ny, 1)])*(3.0*Cjndata[idx(i, jup1+1, Nx, Ny, 1)]-4.0*Cjndata[idx(i, jup2+2, Nx, Ny, 1)]+Cjndata[idx(i, jup3+3, Nx, Ny, 1)]);
    
    IS1_ny[3]=(13.0/12.0)*(Cjndata[idx(i, j, Nx, Ny, 1)]-2.0*Cjndata[idx(i, jup1+1, Nx, Ny, 1)]+Cjndata[idx(i, jup2+2, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-2.0*Cjndata[idx(i, jup1+1, Nx, Ny, 1)]+Cjndata[idx(i, jup2+2, Nx, Ny, 1)])+(1.0/4.0)*(Cjndata[idx(i, j, Nx, Ny, 1)]-Cjndata[idx(i, jup2+2, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-Cjndata[idx(i, jup2+2, Nx, Ny, 1)]);
    
    IS2_ny[3]=(13.0/12.0)*(Cjndata[idx(i, jdown1-1, Nx, Ny, 1)]-2.0*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, jup1+1, Nx, Ny, 1)])*(Cjndata[idx(i, jdown1-1, Nx, Ny, 1)]-2.0*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, jup1+1, Nx, Ny, 1)])+(1.0/4.0)*(Cjndata[idx(i, jdown1-1, Nx, Ny, 1)]-4.0*Cjndata[idx(i, j, Nx, Ny, 1)]+3.0*Cjndata[idx(i, jup1+1, Nx, Ny, 1)])*(Cjndata[idx(i, jdown1-1, Nx, Ny, 1)]-4.0*Cjndata[idx(i, j, Nx, Ny, 1)]+3.0*Cjndata[idx(i, jup1+1, Nx, Ny, 1)]);
	}
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
static int SetUx(realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype *ypdata, realtype *yndata, realtype *taopdata, realtype *taondata, realtype *Cjpdata, realtype *Cjndata, realtype *u_tpphx, realtype *u_tnphx, realtype *u_tpnhx, realtype *u_tnnhx, long int i, long int j, long int Nx, long int Ny)
{
  /* declaration */
  long int idown1, idown2, idown3, iup1, iup2, iup3;

  /* consider Nx */
    if (Nx<4){
        cerr << "\nNx should be more than 3!\n";
        return 1;
    }
    else {
      /* consider situations of different i */
      idown1 = i;
      idown2 = i;
      idown3 = i;
      iup1 = i;
      iup2 = i;
      iup3 = i;

      if(i<3||i>Nx-4){
        if (i==0){
            idown1 = Nx;
            idown2 = Nx;
            idown3 = Nx;
        }
        
        if (i==1){
            idown2 = Nx+1;
            idown3 = Nx+1;
        }
        
        if (i==2){
            idown3 = Nx+2;
        }
        
        if (i==Nx-3){
            iup3 = -3;
        }
        
        if (i==Nx-2){
            iup2 = -2;
            iup3 = -2;
        }
        
        if (i==Nx-1){
            iup1 = -1;
            iup2 = -1;
            iup3 = -1;
        }
      }
      else {
	/* compute positive and negative solutions on the interface */
    u_tpphx[0]=w0_px[0]*((2.0/6.0)*ypdata[idx(idown2-2, j, Nx, Ny, 1)]-(7.0/6.0)*ypdata[idx(idown1-1, j, Nx, Ny, 1)]+(11.0/6.0)*ypdata[idx(i, j, Nx, Ny, 1)])+w1_px[0]*((-1.0/6.0)*ypdata[idx(idown1-1, j, Nx, Ny, 1)]+(5.0/6.0)*ypdata[idx(i, j, Nx, Ny, 1)]+(2.0/6.0)*ypdata[idx(iup1+1, j, Nx, Ny, 1)])+w2_px[0]*((2.0/6.0)*ypdata[idx(i, j, Nx, Ny, 1)]+(5.0/6.0)*ypdata[idx(iup1+1, j, Nx, Ny, 1)]-(1.0/6.0)*ypdata[idx(iup2+2, j, Nx, Ny, 1)]);
    
    u_tnphx[0]=w2_nx[0]*((-1.0/6.0)*yndata[idx(idown1-1, j, Nx, Ny, 1)]+(5.0/6.0)*yndata[idx(i, j, Nx, Ny, 1)]+(2.0/6.0)*yndata[idx(iup1+1, j, Nx, Ny, 1)])+w1_nx[0]*((2.0/6.0)*yndata[idx(i, j, Nx, Ny, 1)]+(5.0/6.0)*yndata[idx(iup1+1, j, Nx, Ny, 1)]-(1.0/6.0)*yndata[idx(iup2+2, j, Nx, Ny, 1)])+w0_nx[0]*((11.0/6.0)*yndata[idx(iup1+1, j, Nx, Ny, 1)]-(7.0/6.0)*yndata[idx(iup2+2, j, Nx, Ny, 1)]+(2.0/6.0)*yndata[idx(iup3+3, j, Nx, Ny, 1)]);
    
    u_tpnhx[0]=w0_px[0]*((2.0/6.0)*ypdata[idx(idown3-3, j, Nx, Ny, 1)]-(7.0/6.0)*ypdata[idx(idown2-2, j, Nx, Ny, 1)]+(11.0/6.0)*ypdata[idx(idown1-1, j, Nx, Ny, 1)])+w1_px[0]*((-1.0/6.0)*ypdata[idx(idown2-2, j, Nx, Ny, 1)]+(5.0/6.0)*ypdata[idx(idown1-1, j, Nx, Ny, 1)]+(2.0/6.0)*ypdata[idx(i, j, Nx, Ny, 1)])+w2_px[0]*((2.0/6.0)*ypdata[idx(idown1-1, j, Nx, Ny, 1)]+(5.0/6.0)*ypdata[idx(i, j, Nx, Ny, 1)]-(1.0/6.0)*ypdata[idx(iup1+1, j, Nx, Ny, 1)]);
    
    u_tnnhx[0]=w2_nx[0]*((-1.0/6.0)*yndata[idx(idown2-2, j, Nx, Ny, 1)]+(5.0/6.0)*yndata[idx(idown1-1, j, Nx, Ny, 1)]+(2.0/6.0)*yndata[idx(i, j, Nx, Ny, 1)])+w1_nx[0]*((2.0/6.0)*yndata[idx(idown1-1, j, Nx, Ny, 1)]+(5.0/6.0)*yndata[idx(i, j, Nx, Ny, 1)]-(1.0/6.0)*yndata[idx(iup1+1, j, Nx, Ny, 1)])+w0_nx[0]*((11.0/6.0)*yndata[idx(i, j, Nx, Ny, 1)]-(7.0/6.0)*yndata[idx(iup1+1, j, Nx, Ny, 1)]+(2.0/6.0)*yndata[idx(iup2+2, j, Nx, Ny, 1)]);
    
    
    u_tpphx[1]=w0_px[1]*((2.0/6.0)*taopdata[idx(idown2-2, j, Nx, Ny, 0)]-(7.0/6.0)*taopdata[idx(idown1-1, j, Nx, Ny, 0)]+(11.0/6.0)*taopdata[idx(i, j, Nx, Ny, 0)])+w1_px[1]*((-1.0/6.0)*taopdata[idx(idown1-1, j, Nx, Ny, 0)]+(5.0/6.0)*taopdata[idx(i, j, Nx, Ny, 0)]+(2.0/6.0)*taopdata[idx(iup1+1, j, Nx, Ny, 0)])+w2_px[1]*((2.0/6.0)*taopdata[idx(i, j, Nx, Ny, 0)]+(5.0/6.0)*taopdata[idx(iup1+1, j, Nx, Ny, 0)]-(1.0/6.0)*taopdata[idx(iup2+2, j, Nx, Ny, 0)]);
    
    u_tnphx[1]=w2_nx[1]*((-1.0/6.0)*taondata[idx(idown1-1, j, Nx, Ny, 0)]+(5.0/6.0)*taondata[idx(i, j, Nx, Ny, 0)]+(2.0/6.0)*taondata[idx(iup1+1, j, Nx, Ny, 0)])+w1_nx[1]*((2.0/6.0)*taondata[idx(i, j, Nx, Ny, 0)]+(5.0/6.0)*taondata[idx(iup1+1, j, Nx, Ny, 0)]-(1.0/6.0)*taondata[idx(iup2+2, j, Nx, Ny, 0)])+w0_nx[1]*((11.0/6.0)*taondata[idx(iup1+1, j, Nx, Ny, 0)]-(7.0/6.0)*taondata[idx(iup2+2, j, Nx, Ny, 0)]+(2.0/6.0)*taondata[idx(iup3+3, j, Nx, Ny, 0)]);
    
    u_tpnhx[1]=w0_px[1]*((2.0/6.0)*taopdata[idx(idown3-3, j, Nx, Ny, 0)]-(7.0/6.0)*taopdata[idx(idown2-2, j, Nx, Ny, 0)]+(11.0/6.0)*taopdata[idx(idown1-1, j, Nx, Ny, 0)])+w1_px[1]*((-1.0/6.0)*taopdata[idx(idown2-2, j, Nx, Ny, 0)]+(5.0/6.0)*taopdata[idx(idown1-1, j, Nx, Ny, 0)]+(2.0/6.0)*taopdata[idx(i, j, Nx, Ny, 0)])+w2_px[1]*((2.0/6.0)*taopdata[idx(idown1-1, j, Nx, Ny, 0)]+(5.0/6.0)*taopdata[idx(i, j, Nx, Ny, 0)]-(1.0/6.0)*taopdata[idx(iup1+1, j, Nx, Ny, 0)]);
    
    u_tnnhx[1]=w2_nx[1]*((-1.0/6.0)*taondata[idx(idown2-2, j, Nx, Ny, 0)]+(5.0/6.0)*taondata[idx(idown1-1, j, Nx, Ny, 0)]+(2.0/6.0)*taondata[idx(i, j, Nx, Ny, 0)])+w1_nx[1]*((2.0/6.0)*taondata[idx(idown1-1, j, Nx, Ny, 0)]+(5.0/6.0)*taondata[idx(i, j, Nx, Ny, 0)]-(1.0/6.0)*taondata[idx(iup1+1, j, Nx, Ny, 0)])+w0_nx[1]*((11.0/6.0)*taondata[idx(i, j, Nx, Ny, 0)]-(7.0/6.0)*taondata[idx(iup1+1, j, Nx, Ny, 0)]+(2.0/6.0)*taondata[idx(iup2+2, j, Nx, Ny, 0)]);
    
    
    u_tpphx[2]=w0_px[2]*((2.0/6.0)*taopdata[idx(idown2-2, j, Nx, Ny, 2)]-(7.0/6.0)*taopdata[idx(idown1-1, j, Nx, Ny, 2)]+(11.0/6.0)*taopdata[idx(i, j, Nx, Ny, 2)])+w1_px[2]*((-1.0/6.0)*taopdata[idx(idown1-1, j, Nx, Ny, 2)]+(5.0/6.0)*taopdata[idx(i, j, Nx, Ny, 2)]+(2.0/6.0)*taopdata[idx(iup1+1, j, Nx, Ny, 2)])+w2_px[2]*((2.0/6.0)*taopdata[idx(i, j, Nx, Ny, 2)]+(5.0/6.0)*taopdata[idx(iup1+1, j, Nx, Ny, 2)]-(1.0/6.0)*taopdata[idx(iup2+2, j, Nx, Ny, 2)]);
    
    u_tnphx[2]=w2_nx[2]*((-1.0/6.0)*taondata[idx(idown1-1, j, Nx, Ny, 2)]+(5.0/6.0)*taondata[idx(i, j, Nx, Ny, 2)]+(2.0/6.0)*taondata[idx(iup1+1, j, Nx, Ny, 2)])+w1_nx[2]*((2.0/6.0)*taondata[idx(i, j, Nx, Ny, 2)]+(5.0/6.0)*taondata[idx(iup1+1, j, Nx, Ny, 2)]-(1.0/6.0)*taondata[idx(iup2+2, j, Nx, Ny, 2)])+w0_nx[2]*((11.0/6.0)*taondata[idx(iup1+1, j, Nx, Ny, 2)]-(7.0/6.0)*taondata[idx(iup2+2, j, Nx, Ny, 2)]+(2.0/6.0)*taondata[idx(iup3+3, j, Nx, Ny, 2)]);
    
    u_tpnhx[2]=w0_px[2]*((2.0/6.0)*taopdata[idx(idown3-3, j, Nx, Ny, 2)]-(7.0/6.0)*taopdata[idx(idown2-2, j, Nx, Ny, 2)]+(11.0/6.0)*taopdata[idx(idown1-1, j, Nx, Ny, 2)])+w1_px[2]*((-1.0/6.0)*taopdata[idx(idown2-2, j, Nx, Ny, 2)]+(5.0/6.0)*taopdata[idx(idown1-1, j, Nx, Ny, 2)]+(2.0/6.0)*taopdata[idx(i, j, Nx, Ny, 2)])+w2_px[2]*((2.0/6.0)*taopdata[idx(idown1-1, j, Nx, Ny, 2)]+(5.0/6.0)*taopdata[idx(i, j, Nx, Ny, 2)]-(1.0/6.0)*taopdata[idx(iup1+1, j, Nx, Ny, 2)]);
    
    u_tnnhx[2]=w2_nx[2]*((-1.0/6.0)*taondata[idx(idown2-2, j, Nx, Ny, 2)]+(5.0/6.0)*taondata[idx(idown1-1, j, Nx, Ny, 2)]+(2.0/6.0)*taondata[idx(i, j, Nx, Ny, 2)])+w1_nx[2]*((2.0/6.0)*taondata[idx(idown1-1, j, Nx, Ny, 2)]+(5.0/6.0)*taondata[idx(i, j, Nx, Ny, 2)]-(1.0/6.0)*taondata[idx(iup1+1, j, Nx, Ny, 2)])+w0_nx[2]*((11.0/6.0)*taondata[idx(i, j, Nx, Ny, 2)]-(7.0/6.0)*taondata[idx(iup1+1, j, Nx, Ny, 2)]+(2.0/6.0)*taondata[idx(iup2+2, j, Nx, Ny, 2)]);
    
    
    u_tpphx[3]=w0_px[3]*((2.0/6.0)*Cjpdata[idx(idown2-2, j, Nx, Ny, 0)]-(7.0/6.0)*Cjpdata[idx(idown1-1, j, Nx, Ny, 0)]+(11.0/6.0)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w1_px[3]*((-1.0/6.0)*Cjpdata[idx(idown1-1, j, Nx, Ny, 0)]+(5.0/6.0)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(2.0/6.0)*Cjpdata[idx(iup1+1, j, Nx, Ny, 0)])+w2_px[3]*((2.0/6.0)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(5.0/6.0)*Cjpdata[idx(iup1+1, j, Nx, Ny, 0)]-(1.0/6.0)*Cjpdata[idx(iup2+2, j, Nx, Ny, 0)]);
    
    u_tnphx[3]=w2_nx[3]*((-1.0/6.0)*Cjndata[idx(idown1-1, j, Nx, Ny, 0)]+(5.0/6.0)*Cjndata[idx(i, j, Nx, Ny, 0)]+(2.0/6.0)*Cjndata[idx(iup1+1, j, Nx, Ny, 0)])+w1_nx[3]*((2.0/6.0)*Cjndata[idx(i, j, Nx, Ny, 0)]+(5.0/6.0)*Cjndata[idx(iup1+1, j, Nx, Ny, 0)]-(1.0/6.0)*Cjndata[idx(iup2+2, j, Nx, Ny, 0)])+w0_nx[3]*((11.0/6.0)*Cjndata[idx(iup1+1, j, Nx, Ny, 0)]-(7.0/6.0)*Cjndata[idx(iup2+2, j, Nx, Ny, 0)]+(2.0/6.0)*Cjndata[idx(iup3+3, j, Nx, Ny, 0)]);
    
    u_tpnhx[3]=w0_px[3]*((2.0/6.0)*Cjpdata[idx(idown3-3, j, Nx, Ny, 0)]-(7.0/6.0)*Cjpdata[idx(idown2-2, j, Nx, Ny, 0)]+(11.0/6.0)*Cjpdata[idx(idown1-1, j, Nx, Ny, 0)])+w1_px[3]*((-1.0/6.0)*Cjpdata[idx(idown2-2, j, Nx, Ny, 0)]+(5.0/6.0)*Cjpdata[idx(idown1-1, j, Nx, Ny, 0)]+(2.0/6.0)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w2_px[3]*((2.0/6.0)*Cjpdata[idx(idown1-1, j, Nx, Ny, 0)]+(5.0/6.0)*Cjpdata[idx(i, j, Nx, Ny, 0)]-(1.0/6.0)*Cjpdata[idx(iup1+1, j, Nx, Ny, 0)]);
    
    u_tnnhx[3]=w2_nx[3]*((-1.0/6.0)*Cjndata[idx(idown2-2, j, Nx, Ny, 0)]+(5.0/6.0)*Cjndata[idx(idown1-1, j, Nx, Ny, 0)]+(2.0/6.0)*Cjndata[idx(i, j, Nx, Ny, 0)])+w1_nx[3]*((2.0/6.0)*Cjndata[idx(idown1-1, j, Nx, Ny, 0)]+(5.0/6.0)*Cjndata[idx(i, j, Nx, Ny, 0)]-(1.0/6.0)*Cjndata[idx(iup1+1, j, Nx, Ny, 0)])+w0_nx[3]*((11.0/6.0)*Cjndata[idx(i, j, Nx, Ny, 0)]-(7.0/6.0)*Cjndata[idx(iup1+1, j, Nx, Ny, 0)]+(2.0/6.0)*Cjndata[idx(iup2+2, j, Nx, Ny, 0)]);
      }
    return 0;
    }
}

/* Get the derivative on y direction */
static int SetUy(realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype *ypdata, realtype *yndata, realtype *taopdata, realtype *taondata, realtype *Cjpdata, realtype *Cjndata, realtype *u_tpphy, realtype *u_tnphy, realtype *u_tpnhy, realtype *u_tnnhy, long int i, long int j, long int Nx, long int Ny)
{
  /* declaration */
  long int jdown1, jdown2, jdown3, jup1, jup2, jup3;

  /* consider Ny */
    if (Ny<4){
        cerr << "\nNy should be more than 3!\n";
        return 1;
    }
    else {
      jdown1 = j;
      jdown2 = j;
      jdown3 = j;
      jup1 = j;
      jup2 = j;
      jup3 = j;
      
      if (j<3||j>Ny-4){
      /* consider situations of different j */
        if (j==0){
            jdown1 = Ny;
            jdown2 = Ny;
            jdown3 = Ny;
        }
        
        if (j==1){
            jdown2 = Ny+1;
            jdown3 = Ny+1;
        }
        
        if (j==2){
            jdown3 = Ny+2;
        }
        
        if (j==Ny-3){
            jup3 = -3;
        }
        
        if (j==Ny-2){
            jup2 = -2;
            jup3 = -2;
        }
        
        if (j==Ny-1){
            jup1 = -1;
            jup2 = -1;
            jup3 = -1;
        }
      }
      else {
	/* compute positive and negative solutions on the interface */
    u_tpphy[0]=w0_py[0]*((2.0/6.0)*ypdata[idx(i, jdown2-2, Nx, Ny, 2)]-(7.0/6.0)*ypdata[idx(i, jdown1-1, Nx, Ny, 2)]+(11.0/6.0)*ypdata[idx(i, j, Nx, Ny, 2)])+w1_py[0]*((-1.0/6.0)*ypdata[idx(i, jdown1-1, Nx, Ny, 2)]+(5.0/6.0)*ypdata[idx(i, j, Nx, Ny, 2)]+(2.0/6.0)*ypdata[idx(i, jup1+1, Nx, Ny, 2)])+w2_py[0]*((2.0/6.0)*ypdata[idx(i, j, Nx, Ny, 2)]+(5.0/6.0)*ypdata[idx(i, jup1+1, Nx, Ny, 2)]-(1.0/6.0)*ypdata[idx(i, jup2+2, Nx, Ny, 2)]);
    
    u_tnphy[0]=w2_ny[0]*((-1.0/6.0)*yndata[idx(i, jdown1-1, Nx, Ny, 2)]+(5.0/6.0)*yndata[idx(i, j, Nx, Ny, 2)]+(2.0/6.0)*yndata[idx(i, jup1+1, Nx, Ny, 2)])+w1_ny[0]*((2.0/6.0)*yndata[idx(i, j, Nx, Ny, 2)]+(5.0/6.0)*yndata[idx(i, jup1+1, Nx, Ny, 2)]-(1.0/6.0)*yndata[idx(i, jup2+2, Nx, Ny, 2)])+w0_ny[0]*((11.0/6.0)*yndata[idx(i, jup1+1, Nx, Ny, 2)]-(7.0/6.0)*yndata[idx(i, jup2+2, Nx, Ny, 2)]+(2.0/6.0)*yndata[idx(i, jup3+3, Nx, Ny, 2)]);
    
    u_tpnhy[0]=w0_py[0]*((2.0/6.0)*ypdata[idx(i, jdown3-3, Nx, Ny, 2)]-(7.0/6.0)*ypdata[idx(i, jdown2-2, Nx, Ny, 2)]+(11.0/6.0)*ypdata[idx(i, jdown1-1, Nx, Ny, 2)])+w1_py[0]*((-1.0/6.0)*ypdata[idx(i, jdown2-2, Nx, Ny, 2)]+(5.0/6.0)*ypdata[idx(i, jdown1-1, Nx, Ny, 2)]+(2.0/6.0)*ypdata[idx(i, j, Nx, Ny, 2)])+w2_py[0]*((2.0/6.0)*ypdata[idx(i, jdown1-1, Nx, Ny, 2)]+(5.0/6.0)*ypdata[idx(i, j, Nx, Ny, 2)]-(1.0/6.0)*ypdata[idx(i, jup1+1, Nx, Ny, 2)]);
    
    u_tnnhy[0]=w2_ny[0]*((-1.0/6.0)*yndata[idx(i, jdown2-2, Nx, Ny, 2)]+(5.0/6.0)*yndata[idx(i, jdown1-1, Nx, Ny, 2)]+(2.0/6.0)*yndata[idx(i, j, Nx, Ny, 2)])+w1_ny[0]*((2.0/6.0)*yndata[idx(i, jdown1-1, Nx, Ny, 2)]+(5.0/6.0)*yndata[idx(i, j, Nx, Ny, 2)]-(1.0/6.0)*yndata[idx(i, jup1+1, Nx, Ny, 2)])+w0_ny[0]*((11.0/6.0)*yndata[idx(i, j, Nx, Ny, 2)]-(7.0/6.0)*yndata[idx(i, jup1+1, Nx, Ny, 2)]+(2.0/6.0)*yndata[idx(i, jup2+2, Nx, Ny, 2)]);
    
    
    u_tpphy[1]=w0_py[1]*((2.0/6.0)*taopdata[idx(i, jdown2-2, Nx, Ny, 1)]-(7.0/6.0)*taopdata[idx(i, jdown1-1, Nx, Ny, 1)]+(11.0/6.0)*taopdata[idx(i, j, Nx, Ny, 1)])+w1_py[1]*((-1.0/6.0)*taopdata[idx(i, jdown1-1, Nx, Ny, 1)]+(5.0/6.0)*taopdata[idx(i, j, Nx, Ny, 1)]+(2.0/6.0)*taopdata[idx(i, jup1+1, Nx, Ny, 1)])+w2_py[1]*((2.0/6.0)*taopdata[idx(i, j, Nx, Ny, 1)]+(5.0/6.0)*taopdata[idx(i, jup1+1, Nx, Ny, 1)]-(1.0/6.0)*taopdata[idx(i, jup2+2, Nx, Ny, 1)]);
    
    u_tnphy[1]=w2_ny[1]*((-1.0/6.0)*taondata[idx(i, jdown1-1, Nx, Ny, 1)]+(5.0/6.0)*taondata[idx(i, j, Nx, Ny, 1)]+(2.0/6.0)*taondata[idx(i, jup1+1, Nx, Ny, 1)])+w1_ny[1]*((2.0/6.0)*taondata[idx(i, j, Nx, Ny, 1)]+(5.0/6.0)*taondata[idx(i, jup1+1, Nx, Ny, 1)]-(1.0/6.0)*taondata[idx(i, jup2+2, Nx, Ny, 1)])+w0_ny[1]*((11.0/6.0)*taondata[idx(i, jup1+1, Nx, Ny, 1)]-(7.0/6.0)*taondata[idx(i, jup2+2, Nx, Ny, 1)]+(2.0/6.0)*taondata[idx(i, jup3+3, Nx, Ny, 1)]);
    
    u_tpnhy[1]=w0_py[1]*((2.0/6.0)*taopdata[idx(i, jdown3-3, Nx, Ny, 1)]-(7.0/6.0)*taopdata[idx(i, jdown2-2, Nx, Ny, 1)]+(11.0/6.0)*taopdata[idx(i, jdown1-1, Nx, Ny, 1)])+w1_py[1]*((-1.0/6.0)*taopdata[idx(i, jdown2-2, Nx, Ny, 1)]+(5.0/6.0)*taopdata[idx(i, jdown1-1, Nx, Ny, 1)]+(2.0/6.0)*taopdata[idx(i, j, Nx, Ny, 1)])+w2_py[1]*((2.0/6.0)*taopdata[idx(i, jdown1-1, Nx, Ny, 1)]+(5.0/6.0)*taopdata[idx(i, j, Nx, Ny, 1)]-(1.0/6.0)*taopdata[idx(i, jup1+1, Nx, Ny, 1)]);
    
    u_tnnhy[1]=w2_ny[1]*((-1.0/6.0)*taondata[idx(i, jdown2-2, Nx, Ny, 1)]+(5.0/6.0)*taondata[idx(i, jdown1-1, Nx, Ny, 1)]+(2.0/6.0)*taondata[idx(i, j, Nx, Ny, 1)])+w1_ny[1]*((2.0/6.0)*taondata[idx(i, jdown1-1, Nx, Ny, 1)]+(5.0/6.0)*taondata[idx(i, j, Nx, Ny, 1)]-(1.0/6.0)*taondata[idx(i, jup1+1, Nx, Ny, 1)])+w0_ny[1]*((11.0/6.0)*taondata[idx(i, j, Nx, Ny, 1)]-(7.0/6.0)*taondata[idx(i, jup1+1, Nx, Ny, 1)]+(2.0/6.0)*taondata[idx(i, jup2+2, Nx, Ny, 1)]);
    
    
    u_tpphy[2]=w0_py[2]*((2.0/6.0)*taopdata[idx(i, jdown2-2, Nx, Ny, 3)]-(7.0/6.0)*taopdata[idx(i, jdown1-1, Nx, Ny, 3)]+(11.0/6.0)*taopdata[idx(i, j, Nx, Ny, 3)])+w1_py[2]*((-1.0/6.0)*taopdata[idx(i, jdown1-1, Nx, Ny, 3)]+(5.0/6.0)*taopdata[idx(i, j, Nx, Ny, 3)]+(2.0/6.0)*taopdata[idx(i, jup1+1, Nx, Ny, 3)])+w2_py[2]*((2.0/6.0)*taopdata[idx(i, j, Nx, Ny, 3)]+(5.0/6.0)*taopdata[idx(i, jup1+1, Nx, Ny, 3)]-(1.0/6.0)*taopdata[idx(i, jup2+2, Nx, Ny, 3)]);
    
    u_tnphy[2]=w2_ny[2]*((-1.0/6.0)*taondata[idx(i, jdown1-1, Nx, Ny, 3)]+(5.0/6.0)*taondata[idx(i, j, Nx, Ny, 3)]+(2.0/6.0)*taondata[idx(i, jup1+1, Nx, Ny, 3)])+w1_ny[2]*((2.0/6.0)*taondata[idx(i, j, Nx, Ny, 3)]+(5.0/6.0)*taondata[idx(i, jup1+1, Nx, Ny, 3)]-(1.0/6.0)*taondata[idx(i, jup2+2, Nx, Ny, 3)])+w0_ny[2]*((11.0/6.0)*taondata[idx(i, jup1+1, Nx, Ny, 3)]-(7.0/6.0)*taondata[idx(i, jup2+2, Nx, Ny, 3)]+(2.0/6.0)*taondata[idx(i, jup3+3, Nx, Ny, 3)]);
    
    u_tpnhy[2]=w0_py[2]*((2.0/6.0)*taopdata[idx(i, jdown3-3, Nx, Ny, 3)]-(7.0/6.0)*taopdata[idx(i, jdown2-2, Nx, Ny, 3)]+(11.0/6.0)*taopdata[idx(i, jdown1-1, Nx, Ny, 3)])+w1_py[2]*((-1.0/6.0)*taopdata[idx(i, jdown2-2, Nx, Ny, 3)]+(5.0/6.0)*taopdata[idx(i, jdown1-1, Nx, Ny, 3)]+(2.0/6.0)*taopdata[idx(i, j, Nx, Ny, 3)])+w2_py[2]*((2.0/6.0)*taopdata[idx(i, jdown1-1, Nx, Ny, 3)]+(5.0/6.0)*taopdata[idx(i, j, Nx, Ny, 3)]-(1.0/6.0)*taopdata[idx(i, jup1+1, Nx, Ny, 3)]);
    
    u_tnnhy[2]=w2_ny[2]*((-1.0/6.0)*taondata[idx(i, jdown2-2, Nx, Ny, 3)]+(5.0/6.0)*taondata[idx(i, jdown1-1, Nx, Ny, 3)]+(2.0/6.0)*taondata[idx(i, j, Nx, Ny, 3)])+w1_ny[2]*((2.0/6.0)*taondata[idx(i, jdown1-1, Nx, Ny, 3)]+(5.0/6.0)*taondata[idx(i, j, Nx, Ny, 3)]-(1.0/6.0)*taondata[idx(i, jup1+1, Nx, Ny, 3)])+w0_ny[2]*((11.0/6.0)*taondata[idx(i, j, Nx, Ny, 3)]-(7.0/6.0)*taondata[idx(i, jup1+1, Nx, Ny, 3)]+(2.0/6.0)*taondata[idx(i, jup2+2, Nx, Ny, 3)]);
    
    
    u_tpphy[3]=w0_py[3]*((2.0/6.0)*Cjpdata[idx(i, jdown2-2, Nx, Ny, 1)]-(7.0/6.0)*Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)]+(11.0/6.0)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w1_py[3]*((-1.0/6.0)*Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)]+(5.0/6.0)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(2.0/6.0)*Cjpdata[idx(i, jup1+1, Nx, Ny, 1)])+w2_py[3]*((2.0/6.0)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(5.0/6.0)*Cjpdata[idx(i, jup1+1, Nx, Ny, 1)]-(1.0/6.0)*Cjpdata[idx(i, jup2+2, Nx, Ny, 1)]);
    
    u_tnphy[3]=w2_ny[3]*((-1.0/6.0)*Cjndata[idx(i, jdown1-1, Nx, Ny, 1)]+(5.0/6.0)*Cjndata[idx(i, j, Nx, Ny, 1)]+(2.0/6.0)*Cjndata[idx(i, jup1+1, Nx, Ny, 1)])+w1_ny[3]*((2.0/6.0)*Cjndata[idx(i, j, Nx, Ny, 1)]+(5.0/6.0)*Cjndata[idx(i, jup1+1, Nx, Ny, 1)]-(1.0/6.0)*Cjndata[idx(i, jup2+2, Nx, Ny, 1)])+w0_ny[3]*((11.0/6.0)*Cjndata[idx(i, jup1+1, Nx, Ny, 1)]-(7.0/6.0)*Cjndata[idx(i, jup2+2, Nx, Ny, 1)]+(2.0/6.0)*Cjndata[idx(i, jup3+3, Nx, Ny, 1)]);
    
    u_tpnhy[3]=w0_py[3]*((2.0/6.0)*Cjpdata[idx(i, jdown3-3, Nx, Ny, 1)]-(7.0/6.0)*Cjpdata[idx(i, jdown2-2, Nx, Ny, 1)]+(11.0/6.0)*Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)])+w1_py[3]*((-1.0/6.0)*Cjpdata[idx(i, jdown2-2, Nx, Ny, 1)]+(5.0/6.0)*Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)]+(2.0/6.0)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w2_py[3]*((2.0/6.0)*Cjpdata[idx(i, jdown1-1, Nx, Ny, 1)]+(5.0/6.0)*Cjpdata[idx(i, j, Nx, Ny, 1)]-(1.0/6.0)*Cjpdata[idx(i, jup1+1, Nx, Ny, 1)]);
    
    u_tnnhy[3]=w2_ny[3]*((-1.0/6.0)*Cjndata[idx(i, jdown2-2, Nx, Ny, 1)]+(5.0/6.0)*Cjndata[idx(i, jdown1-1, Nx, Ny, 1)]+(2.0/6.0)*Cjndata[idx(i, j, Nx, Ny, 1)])+w1_ny[3]*((2.0/6.0)*Cjndata[idx(i, jdown1-1, Nx, Ny, 1)]+(5.0/6.0)*Cjndata[idx(i, j, Nx, Ny, 1)]-(1.0/6.0)*Cjndata[idx(i, jup1+1, Nx, Ny, 1)])+w0_ny[3]*((11.0/6.0)*Cjndata[idx(i, j, Nx, Ny, 1)]-(7.0/6.0)*Cjndata[idx(i, jup1+1, Nx, Ny, 1)]+(2.0/6.0)*Cjndata[idx(i, jup2+2, Nx, Ny, 1)]);
      }
    return 0;
    }
}

/* Fill in values of tao in the whole domain */
static int Gettao(N_Vector y, N_Vector tao, long int Nx, long int Ny, realtype gama)
{
    /* declare parameters */
    long int i, j, NEQ;
    realtype *vx_data, *vy_data, *data, *tao_data;
    NEQ = Nx*Ny;

    /* create vectors */
    N_Vector vx = NULL;
    N_Vector vy = NULL;
    
    /* Create serial vector of length N */
    vx = N_VNew_Serial(NEQ);
    if (check_flag((void *) vx, "N_VNew_Serial", 0)) return 1;
    vy = N_VNew_Serial(NEQ);
    if (check_flag((void *) vy, "N_VNew_Serial", 0)) return 1;
    
    /* Access data array for new NVector y, tao, vx, vy */
    vx_data = N_VGetArrayPointer(vx);
    if (check_flag((void *) vx_data, "N_VGetArrayPointer", 0)) return 1;
    vy_data = N_VGetArrayPointer(vy);
    if (check_flag((void *) vy_data, "N_VGetArrayPointer", 0)) return 1;
    data = N_VGetArrayPointer(y);
    if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
    tao_data = N_VGetArrayPointer(tao);
    if (check_flag((void *) tao_data, "N_VGetArrayPointer", 0)) return 1;
    
    /* Set values into vx and vy */
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            vx_data[idx_v(i,j,Nx)]=data[idx(i, j, Nx, Ny, 1)]/data[idx(i, j, Nx, Ny, 0)];
            vy_data[idx_v(i,j,Nx)]=data[idx(i, j, Nx, Ny, 2)]/data[idx(i, j, Nx, Ny, 0)];
        }
    }
    
    /* Compute the values of tao in the whole domain */
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            tao_data[idx(i, j, Nx, Ny, 0)] = data[idx(i, j, Nx, Ny, 1)]*vx_data[idx_v(i,j,Nx)]+data[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*data[idx(i, j, Nx, Ny, 3)];
            tao_data[idx(i, j, Nx, Ny, 1)] = data[idx(i, j, Nx, Ny, 1)]*vy_data[idx_v(i,j,Nx)];
            tao_data[idx(i, j, Nx, Ny, 2)] = data[idx(i, j, Nx, Ny, 2)]*vx_data[idx_v(i,j,Nx)];
            tao_data[idx(i, j, Nx, Ny, 3)] = data[idx(i, j, Nx, Ny, 2)]*vy_data[idx_v(i,j,Nx)]+data[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*data[idx(i, j, Nx, Ny, 3)];
        }
    }

    /* Free vectors */
    N_VDestroy_Serial(vx);
    N_VDestroy_Serial(vy);
    
    return 0;
}

/* Fill in the value of J in the whole domain */
static int GetCj(N_Vector y, N_Vector Cj, long int Nx, long int Ny, realtype gama)
{
    /* declare parameters */
    long int i, j, NEQ;
    realtype *vx_data, *vy_data, *data, *Cj_data;
    NEQ = Nx*Ny;

    /* create vectors */
    N_Vector vx = NULL;
    N_Vector vy = NULL;
    
    /* Create serial vector of length N */
    vx = N_VNew_Serial(NEQ);
    if (check_flag((void *) vx, "N_VNew_Serial", 0)) return 1;
    vy = N_VNew_Serial(NEQ);
    if (check_flag((void *) vy, "N_VNew_Serial", 0)) return 1;
    
    /* Access data array for new NVector y, Cj, vx, vy */
    vx_data = N_VGetArrayPointer(vx);
    if (check_flag((void *) vx_data, "N_VGetArrayPointer", 0)) return 1;
    vy_data = N_VGetArrayPointer(vy);
    if (check_flag((void *) vy_data, "N_VGetArrayPointer", 0)) return 1;
    data = N_VGetArrayPointer(y);
    if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
    Cj_data = N_VGetArrayPointer(Cj);
    if (check_flag((void *) Cj_data, "N_VGetArrayPointer", 0)) return 1;

    /* Set values into vx and vy */
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            vx_data[idx_v(i,j,Nx)]=data[idx(i, j, Nx, Ny, 1)]/data[idx(i, j, Nx, Ny, 0)];
            vy_data[idx_v(i,j,Nx)]=data[idx(i, j, Nx, Ny, 2)]/data[idx(i, j, Nx, Ny, 0)];
        }
    }    

    /* Compute the values of Cj in the whole domain */
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            Cj_data[idx(i, j, Nx, Ny, 0)] = data[idx(i, j, Nx, Ny, 1)]*data[idx(i, j, Nx, Ny, 1)]*vx_data[idx_v(i,j,Nx)]+data[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*data[idx(i, j, Nx, Ny, 3)]*data[idx(i, j, Nx, Ny, 1)]+data[idx(i, j, Nx, Ny, 1)]*data[idx(i, j, Nx, Ny, 2)]*vy_data[idx_v(i,j,Nx)];
            Cj_data[idx(i, j, Nx, Ny, 1)] = data[idx(i, j, Nx, Ny, 2)]*data[idx(i, j, Nx, Ny, 2)]*vy_data[idx_v(i,j,Nx)]+data[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*data[idx(i, j, Nx, Ny, 3)]*data[idx(i, j, Nx, Ny, 2)]+data[idx(i, j, Nx, Ny, 1)]*data[idx(i, j, Nx, Ny, 2)]*vx_data[idx_v(i,j,Nx)];
        }
    }

    /* Free vectors */
    N_VDestroy_Serial(vx);
    N_VDestroy_Serial(vy);
    
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
