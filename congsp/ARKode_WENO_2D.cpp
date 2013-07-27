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
static int Gettao(realtype *yxbackdata, realtype *yybackdata, realtype *tao, long int Nx, long int Ny, realtype gama);

/* Set value of J in whole domain*/
static int GetCj(realtype *yxbackdata, realtype *yybackdata, realtype *Cj, long int Nx, long int Ny, realtype gama);

/* Set left eigenvectors in x direction */
static int Setlfxegm(realtype *Ydata, realtype **lfxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype *egvx);

/* Set right eigenvectors in x direction */
static int Setrhxegm(realtype *Ydata, realtype **rhxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag);

/* Set left eigenvectors in y direction */
static int Setlfyegm(realtype *Ydata, realtype **lfyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype *egvy);

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
//static int SetUx(realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype *yxpdata, realtype *yxndata, realtype *u_tpphx, realtype *u_tnphx, realtype *u_tpnhx, realtype *u_tnnhx, long int i, long int j, long int Nx, long int Ny, int flag);
static int SetUx(realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype *yxpdata, realtype *yxndata, realtype *u_tpphx, realtype *u_tnphx, long int i, long int j, long int Nx, long int Ny, int flag);

/* Get the derivative on y direction */
//static int SetUy(realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype *yypdata, realtype *yyndata, realtype *u_tpphy, realtype *u_tnphy, realtype *u_tpnhy, realtype *u_tnnhy, long int i, long int j, long int Nx, long int Ny, int flag);
static int SetUy(realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype *yypdata, realtype *yyndata, realtype *u_tpphy, realtype *u_tnphy, long int i, long int j, long int Nx, long int Ny, int flag);

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
    int flag;
    realtype p, Epsilon;
    Epsilon = 10.0e-6;
    
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
    
    realtype *yxdata = new realtype [4*(Nx+1)*Ny];
    realtype *yydata = new realtype [4*Nx*(Ny+1)];
    realtype *yxpdata = new realtype [4*Nx*Ny];
    realtype *yxndata = new realtype [4*Nx*Ny];
    realtype *yypdata = new realtype [4*Nx*Ny];
    realtype *yyndata = new realtype [4*Nx*Ny];
    realtype *yxnewdata = new realtype [4*Nx*Ny];
    realtype *yynewdata = new realtype [4*Nx*Ny];
    realtype *yxbackdata = new realtype [4*(Nx+1)*Ny];
    realtype *yybackdata = new realtype [4*Nx*(Ny+1)];
    realtype *tao = new realtype [2*(Nx+1)*Ny+2*Nx*(Ny+1)];
    realtype *Cj = new realtype [(Nx+1)*Ny+Nx*(Ny+1)];
    realtype *egvx = new realtype [4*Nx*Ny];
    realtype *egvy = new realtype [4*Nx*Ny];
    realtype *egvxmax = new realtype [Nx*Ny];
    realtype *egvymax = new realtype [Nx*Ny];

    /* access data arrays */
    realtype *Ydata = N_VGetArrayPointer(y);
    if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;
    realtype *dYdata = N_VGetArrayPointer(ydot);
    if (check_flag((void *) dYdata, "N_VGetArrayPointer", 0)) return 1;
    
    /* compute max absolue eigenvalue and fill in ypdata, yndata, taopdata, taondata, Cjpdata, Cjpdata */
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
	    flag = Setlfxegm(Ydata, lfxegm, i, j, Nx, Ny, gama, egvx);
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
            //printf("   yxnew problem parameters: i = %li, j = %li, yxnew0 = %e,  yxnew1 = %e, yxnew2 = %e,  yxnew3 = %e\n", i, j, yxnewdata[idx(i, j, Nx, Ny, 0)], yxnewdata[idx(i, j, Nx, Ny, 1)], yxnewdata[idx(i, j, Nx, Ny, 2)], yxnewdata[idx(i, j, Nx, Ny, 3)]);
           
            egvxmax[idx_v(i,j,Nx)] = (fabs(egvx[idx(i,j,Nx,Ny,0)])>fabs(egvx[idx(i,j,Nx,Ny,2)])) ? fabs(egvx[idx(i,j,Nx,Ny,0)]) : fabs(egvx[idx(i,j,Nx,Ny,2)]);
	 
	    //if (i==0&&j==0){    
	    //printf(" maxeigen : i = %li, j = %li, egv1 = %f, egv2 = %f, egv3 = %f, egv4 = %f, egvmaxx = %f\n", i, j, egvx[idx(i,j,Nx,Ny,0)], egvx[idx(i,j,Nx,Ny,1)],egvx[idx(i,j,Nx,Ny,2)],egvx[idx(i,j,Nx,Ny,3)], egvxmax[idx_v(i,j,Nx)]);
	    //}
            yxpdata[idx(i, j, Nx, Ny, 0)]=0.5*(egvx[idx(i,j,Nx,Ny,0)]*yxnewdata[idx(i, j, Nx, Ny, 0)]+egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 0)]);
            yxndata[idx(i, j, Nx, Ny, 0)]=0.5*(egvx[idx(i,j,Nx,Ny,0)]*yxnewdata[idx(i, j, Nx, Ny, 0)]-egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 0)]);
	    //printf("i = %li, j = %li,yxpdata[idx(i, j, Nx, Ny, 0)] = %g, yxndata[idx(i, j, Nx, Ny, 0)] = %g\n", i,j,yxpdata[idx(i, j, Nx, Ny, 0)], yxndata[idx(i, j, Nx, Ny, 0)]);
            yxpdata[idx(i, j, Nx, Ny, 1)]=0.5*(egvx[idx(i,j,Nx,Ny,1)]*yxnewdata[idx(i, j, Nx, Ny, 1)]+egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 1)]);
            yxndata[idx(i, j, Nx, Ny, 1)]=0.5*(egvx[idx(i,j,Nx,Ny,1)]*yxnewdata[idx(i, j, Nx, Ny, 1)]-egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 1)]);
	    //printf("i = %li, j = %li,yxpdata[idx(i, j, Nx, Ny, 1)] = %g, yxndata[idx(i, j, Nx, Ny, 1)] = %g\n", i,j,yxpdata[idx(i, j, Nx, Ny, 1)], yxndata[idx(i, j, Nx, Ny, 1)]);
            yxpdata[idx(i, j, Nx, Ny, 2)]=0.5*(egvx[idx(i,j,Nx,Ny,2)]*yxnewdata[idx(i, j, Nx, Ny, 2)]+egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 2)]);
            yxndata[idx(i, j, Nx, Ny, 2)]=0.5*(egvx[idx(i,j,Nx,Ny,2)]*yxnewdata[idx(i, j, Nx, Ny, 2)]-egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 2)]);
	    //printf("i = %li, j = %li,yxpdata[idx(i, j, Nx, Ny, 2)] = %g, yxndata[idx(i, j, Nx, Ny, 2)] = %g\n", i,j,yxpdata[idx(i, j, Nx, Ny, 2)], yxndata[idx(i, j, Nx, Ny, 2)]);
            yxpdata[idx(i, j, Nx, Ny, 3)]=0.5*(egvx[idx(i,j,Nx,Ny,3)]*yxnewdata[idx(i, j, Nx, Ny, 3)]+egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 3)]);
            yxndata[idx(i, j, Nx, Ny, 3)]=0.5*(egvx[idx(i,j,Nx,Ny,3)]*yxnewdata[idx(i, j, Nx, Ny, 3)]-egvxmax[idx_v(i,j,Nx)]*yxnewdata[idx(i, j, Nx, Ny, 3)]);
	    //printf("i = %li, j = %li,yxpdata[idx(i, j, Nx, Ny, 3)] = %g, yxndata[idx(i, j, Nx, Ny, 3)] = %g\n", i,j,yxpdata[idx(i, j, Nx, Ny, 3)], yxndata[idx(i, j, Nx, Ny, 3)]);
        }
    }
    
    for(i=0; i<Nx; i++){
        for (j=0; j<Ny; j++){
	  //flag = Setlfyegm(Ydata, lfyegm, i, j, Nx, Ny, gama, egv1y, egv2y, egv3y, egv4y, 2);
	  flag = Setlfyegm(Ydata, lfyegm, i, j, Nx, Ny, gama, egvy);
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
      for (i=0; i<(Nx+1); i++){
	long int n;
            /* get derivative on x direction */
            flag = SetISX(IS0_px, IS1_px, IS2_px, IS0_nx, IS1_nx, IS2_nx, yxpdata, yxndata, i, j, Nx, Ny, 0);
            if (flag!=0) printf("error in SetISX function \n");
	    for (n=0;n<4;n++){
	      // printf("i=%li,j=%li,n=%li,IS0_px[n]=%f,IS1_px[n]=%f,IS2_px[n]=%f,IS0_nx[n]=%f,IS1_nx[n]=%f,IS2_nx[n]=%f\n",i,j,n,IS0_px[n],IS1_px[n],IS2_px[n],IS0_nx[n],IS1_nx[n],IS2_nx[n]);
	      }
            flag = Setalphawx(IS0_px, IS1_px, IS2_px, IS0_nx, IS1_nx, IS2_nx, alpha_0px, alpha_1px, alpha_2px, alpha_0nx, alpha_1nx, alpha_2nx, w0_px, w1_px, w2_px, w0_nx, w1_nx, w2_nx, Epsilon);
            if (flag!=0) printf("error in Setalphawx function \n");
            
            flag = SetUx(w0_px, w1_px, w2_px, w0_nx, w1_nx, w2_nx, yxpdata, yxndata, u_tpphx, u_tnphx, i, j, Nx, Ny, 0);
            if (flag!=0) printf("error in SetUx function \n");
            
	    //if (i==0&j==0){
            for(k=0;k<4;k++){
	      //printf("w: i=%li, j=%li, k=%li, w0_px[k]=%g,w1_px[k]=%g,w2_px[k]=%g,w0_nx[k]=%g,w1_nx[k]=%g,w2_nx[k]=%g\n",i,j,k,w0_px[k],w1_px[k],w2_px[k],w0_nx[k],w1_nx[k],w2_nx[k]);
	      yxdata[idx(i, j, Nx+1, Ny, k)]=(u_tpphx[k]+u_tnphx[k]);
	      //yxdata[idx(i, j, Nx+1, Ny, k)]=(u_tpnhx[k]+u_tnnhx[k]);
	      //printf("yxdata: i=%li, j=%li, k=%li, yxdata[idx(i, j, Nx+1, Ny, k)]=%g\n",i,j,k,yxdata[idx(i, j, Nx+1, Ny, k)]);
	      //printf("yx: i=%li, j=%li, k=%li, u_tpphx[k]=%f, u_tnphx[k]=%f\n",i,j, k,u_tpphx[k], u_tnphx[k]);
            }
	    
	    //	    }
        }
    }
    
    for (i=0;i<Nx;i++){
      for(j=0;j<(Ny+1);j++){
            
            /* get derivative on y direction */
            flag = SetISY(IS0_py, IS1_py, IS2_py, IS0_ny, IS1_ny, IS2_ny, yypdata, yyndata, i, j, Nx, Ny, 0);
            if (flag!=0) printf("error in SetISY function \n");
            
            flag = Setalphawy(IS0_py, IS1_py, IS2_py, IS0_ny, IS1_ny, IS2_ny, alpha_0py, alpha_1py, alpha_2py, alpha_0ny, alpha_1ny, alpha_2ny, w0_py, w1_py, w2_py, w0_ny, w1_ny, w2_ny, Epsilon);
            if (flag!=0) printf("error in Setalphawy function \n");
            
            flag = SetUy(w0_py, w1_py, w2_py, w0_ny, w1_ny, w2_ny, yypdata, yyndata, u_tpphy, u_tnphy, i, j, Nx, Ny, 0);
            if (flag!=0) printf("error in SetUy function \n");
            
            for(k=0;k<4;k++){
                yydata[idx(i, j, Nx, Ny+1, k)]=u_tpphy[k]+u_tnphy[k];
                //yydata[idx(i, j, Nx, Ny+1, k)]=u_tpnhy[k]+u_tnnhy[k];
		//printf("yy: i=%li, j=%li, k=%li, u_tpphy[k]=%f, u_tnphy[k]=%f,u_tpnhy[k]=%f,u_tnnhy[k]=%f\n",i,j, k,u_tpphy[k], u_tnphy[k], u_tpnhy[k],u_tnnhy[k]);
            }
        }
    }
    
    for(j=0; j<Ny; j++){
      for (i=0; i<(Nx+1); i++){
            flag = Setrhxegm(Ydata, rhxegm, i, j, Nx, Ny, gama, 0);
            if (flag!=0) printf("error in Setrhxegm function \n");
            for (k=0;k<4;k++){
	      //if (i==0&&j==0){
	      //printf(" rhxegm : i = %li, j = %li, k = %li, rhxegm[k][0] = %f, rhxegm[k][1] = %f, rhxegm[k][2] = %f, rhxegm[k][3] = %f\n", i, j, k, rhxegm[k][0], rhxegm[k][1],rhxegm[k][2],rhxegm[k][3]);
	      //}
                yxbackdata[idx(i, j, Nx+1, Ny, k)] = rhxegm[k][0]*yxdata[idx(i, j, Nx+1, Ny, 0)]+rhxegm[k][1]*yxdata[idx(i, j, Nx+1, Ny, 1)]+rhxegm[k][2]*yxdata[idx(i, j, Nx+1, Ny, 2)]+rhxegm[k][3]*yxdata[idx(i, j, Nx+1, Ny, 3)];
            }
	    //printf(" yxdata : i = %li, j = %li, yxdata[idx(i, j, Nx+1, Ny, 0)]=%f, yxdata[idx(i, j, Nx+1, Ny, 1)]=%f, yxdata[idx(i, j, Nx+1, Ny, 2)]=%f, yxdata[idx(i, j, Nx+1, Ny, 3)]=%f\n",i,j,yxdata[idx(i, j, Nx+1, Ny, 0)],yxdata[idx(i, j, Nx+1, Ny, 1)],yxdata[idx(i, j, Nx+1, Ny, 2)],yxdata[idx(i, j, Nx+1, Ny, 3)]);
	    //printf("yxbackdata : i = %li, j = %li, yxbackdata[idx(i, j, Nx+1, Ny, 0)]=%f, yxbackdata[idx(i, j, Nx+1, Ny, 1)]=%f, yxbackdata[idx(i, j, Nx+1, Ny, 2)]=%f, yxbackdata[idx(i, j, Nx+1, Ny, 3)]=%e\n",i,j,yxbackdata[idx(i, j, Nx+1, Ny, 0)],yxbackdata[idx(i, j, Nx+1, Ny, 1)],yxbackdata[idx(i, j, Nx+1, Ny, 2)],yxbackdata[idx(i, j, Nx+1, Ny, 3)]);
        }
    }
    
    for(i=0; i<Nx; i++){
      for (j=0; j<(Ny+1); j++){
            flag = Setrhyegm(Ydata, rhyegm, i, j, Nx, Ny, gama, 0);
            if (flag!=0) printf("error in Setrhyegm function \n");
            for (k=0;k<4;k++){
	      //  if (i==0&&j==0){
	      //printf(" rhyegm : i = %li, j = %li, rhyegm[k][0] = %f, rhyegm[k][1] = %f, rhyegm[k][2] = %f, rhyegm[k][3] = %f\n", i, j, rhyegm[k][0], rhyegm[k][1],rhyegm[k][2],rhyegm[k][3]);
	      //  }
                yybackdata[idx(i, j, Nx, Ny+1, k)] = rhyegm[k][0]*yydata[idx(i, j, Nx, Ny+1, 0)]+rhyegm[k][1]*yydata[idx(i, j, Nx, Ny+1, 1)]+rhyegm[k][2]*yydata[idx(i, j, Nx, Ny+1, 2)]+rhyegm[k][3]*yydata[idx(i, j, Nx, Ny+1, 3)];
            }
	    
	    //if(i==0&&j==0){
	    //	printf(" yydownbackdata : i = %li, j = %li, yydownbackdata[idx(i, j, Nx, Ny, 0)] = %f, yydownbackdata[idx(i, j, Nx, Ny, 1)] = %f, yydownbackdata[idx(i, j, Nx, Ny, 2)] = %f, yydownbackdata[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, yydownbackdata[idx(i, j, Nx, Ny, 0)], yydownbackdata[idx(i, j, Nx, Ny, 1)],yydownbackdata[idx(i, j, Nx, Ny, 2)],yydownbackdata[idx(i, j, Nx, Ny, 3)]);
	    //	}
        }
    }
    
     /* fill in the value of tao and J in the whole domain */
    flag = Gettao(yxbackdata, yybackdata, tao, Nx, Ny, gama);
    if (flag!=0) printf("error in Gettao function \n");
    flag = GetCj(yxbackdata, yybackdata, Cj, Nx, Ny, gama);
    if (flag!=0) printf("error in GetCj function \n");
    
    for(j=0; j<Ny; j++){
      for (i=0; i<Nx; i++){
            /* get derivative both on x and y direction */
            //for (k=0; k<4; k++){
                //  if(k==2)
                //        dYdata[idx(i, j, Nx, Ny, k)]=yxdata[idx(i, j, Nx, Ny, k)]+yydata[idx(i, j, Nx, Ny, k)]-0.1;
                // else
	      //if (i==0&&j==0){
	      //printf(" dYdata0 : i = %li, j = %li, yxbackdata[idx(i+1, j, (Nx+1), Ny, 1)] = %f, yxbackdata[idx(i, j, (Nx+1), Ny, 1)] = %f, yybackdata[idx(i, j+1, Nx, (Ny+1), 2)] = %f, yybackdata[idx(i, j, Nx, (Ny+1), 2)] = %f\n", i, j, yxbackdata[idx(i+1, j, (Nx+1), Ny, 1)], yxbackdata[idx(i, j, (Nx+1), Ny, 1)], yybackdata[idx(i, j+1, Nx, (Ny+1), 2)],yybackdata[idx(i, j, Nx, (Ny+1), 2)]);
	      //printf(" dYdata1 : i = %li, j = %li, tao[idx_v(i+1,j,Nx+1)] = %f, tao[idx_v(i,j,Nx+1)] = %f, tao[idx_v(i,j+1,Nx)+2*((Nx+1)*Ny)] = %f, tao[idx_v(i,j,Nx)+2*((Nx+1)*Ny)] = %f\n", i, j, tao[idx_v(i+1,j,Nx+1)], tao[idx_v(i,j,Nx+1)], tao[idx_v(i,j+1,Nx)+2*((Nx+1)*Ny)], tao[idx_v(i,j,Nx)+2*((Nx+1)*Ny)]);
	      //printf(" dYdata2 : i = %li, j = %li, tao[idx_v(i+1,j,Nx+1)+(Nx+1)*Ny] = %f, tao[idx_v(i,j,Nx+1)+(Nx+1)*Ny] = %f, tao[idx_v(i,j+1,Nx)+2*((Nx+1)*Ny)+Nx*(Ny+1)] = %f, tao[idx_v(i,j,Nx)+2*((Nx+1)*Ny)+Nx*(Ny+1)] = %f\n", i, j, tao[idx_v(i+1,j,Nx+1)+(Nx+1)*Ny], tao[idx_v(i,j,Nx+1)+(Nx+1)*Ny], tao[idx_v(i,j+1,Nx)+2*((Nx+1)*Ny)+Nx*(Ny+1)], tao[idx_v(i,j,Nx)+2*((Nx+1)*Ny)+Nx*(Ny+1)]);
	      //xprintf(" dYdata3 : i = %li, j = %li, Cj[idx_v(i+1,j,Nx+1)] = %f, Cj[idx_v(i,j,Nx+1)] = %f, Cj[idx_v(i,j+1,Nx)+(Nx+1)*Ny] = %f, Cj[idx_v(i,j,Nx)+(Nx+1)*Ny] = %f\n", i, j, Cj[idx_v(i+1,j,Nx+1)], Cj[idx_v(i,j,Nx+1)], Cj[idx_v(i,j+1,Nx)+(Nx+1)*Ny], Cj[idx_v(i,j,Nx)+(Nx+1)*Ny]);
	      //}
	dYdata[idx(i, j, Nx, Ny, 0)]=(-1.0/dx)*(yxbackdata[idx(i+1, j, (Nx+1), Ny, 1)]-yxbackdata[idx(i, j, (Nx+1), Ny, 1)])+(-1.0/dy)*(yybackdata[idx(i, j+1, Nx, (Ny+1), 2)]-yybackdata[idx(i, j, Nx, (Ny+1), 2)]);
	dYdata[idx(i, j, Nx, Ny, 1)]=(-1.0/dx)*(tao[idx_v(i+1,j,Nx+1)]-tao[idx_v(i,j,Nx+1)])+(-1.0/dy)*(tao[idx_v(i,j+1,Nx)+2*((Nx+1)*Ny)]-tao[idx_v(i,j,Nx)+2*((Nx+1)*Ny)]);
	dYdata[idx(i, j, Nx, Ny, 2)]=(-1.0/dx)*(tao[idx_v(i+1,j,Nx+1)+(Nx+1)*Ny]-tao[idx_v(i,j,Nx+1)+(Nx+1)*Ny])+(-1.0/dy)*(tao[idx_v(i,j+1,Nx)+2*((Nx+1)*Ny)+Nx*(Ny+1)]-tao[idx_v(i,j,Nx)+2*((Nx+1)*Ny)+Nx*(Ny+1)]);
        dYdata[idx(i, j, Nx, Ny, 3)]=(-1.0/dx)*(Cj[idx_v(i+1,j,Nx+1)]-Cj[idx_v(i,j,Nx+1)])+(-1.0/dy)*(Cj[idx_v(i,j+1,Nx)+(Nx+1)*Ny]-Cj[idx_v(i,j,Nx)+(Nx+1)*Ny]);
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
    
    delete []yxdata;
    delete []yydata;
    delete []yxpdata;
    delete []yxndata;
    delete []yypdata;
    delete []yyndata;
    delete []yxnewdata;
    delete []yynewdata;
    delete []yxbackdata;
    delete []yybackdata;
    delete []tao;
    delete []Cj;
    delete []egvx;
    delete []egvy;
    delete []egvxmax;
    delete []egvymax;

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

/* fill in left eigenvectors for x component*/
static int Setlfxegm(realtype *Ydata, realtype **lfxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype *egvx)
{
    realtype rou, vx, vy, p, a, vxnext, vynext, pnext, vxcur, vycur, pcur, h, hcur, hnext;
   
    vx = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
    vy = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
    p = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vx*vx+vy*vy);
    a = sqrt(gama*p/Ydata[idx(i, j, Nx, Ny, 0)]);
    h = (Ydata[idx(i, j, Nx, Ny, 3)]+p)/Ydata[idx(i, j, Nx, Ny, 0)];
   
    egvx[idx(i,j,Nx,Ny,0)]=vx-a;
    egvx[idx(i,j,Nx,Ny,1)]=vx;
    egvx[idx(i,j,Nx,Ny,2)]=vx;
    egvx[idx(i,j,Nx,Ny,3)]=vx+a;
    /*
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
    */
    /*
    lfxegm[0][0] = vx/(2.0*a)+(vx*vx+vy*vy)/(4.0*h-2.0*(vx*vx+vy*vy));
    lfxegm[1][0] = 1.0-(vx*vx+vy*vy)/(2.0*h-1.0*(vx*vx+vy*vy));
    lfxegm[2][0] = -vx/(2.0*a)+(vx*vx+vy*vy)/(4.0*h-2.0*(vx*vx+vy*vy));
    lfxegm[3][0] = -vy;
    lfxegm[0][1] = -1.0/(2.0*a)+vx/(2.0*h-1.0*(vx*vx+vy*vy));
    lfxegm[1][1] = 2.0*vx/(2.0*h-1.0*(vx*vx+vy*vy));
    lfxegm[2][1] = 1.0/(2.0*a)+vx/(2.0*h-1.0*(vx*vx+vy*vy));
    lfxegm[3][1] = 0.0;
    lfxegm[0][2] = -vy/(2.0*h-1.0*(vx*vx+vy*vy));
    lfxegm[1][2] = 2.0*vy/(2.0*h-1.0*(vx*vx+vy*vy));
    lfxegm[2][2] = -vy/(2.0*h-1.0*(vx*vx+vy*vy));
    lfxegm[3][2] = 1.0;
    lfxegm[0][3] = 1.0/(2.0*h-1.0*(vx*vx+vy*vy));
    lfxegm[1][3] = -2.0/(2.0*h-1.0*(vx*vx+vy*vy));
    lfxegm[2][3] = 1.0/(2.0*h-1.0*(vx*vx+vy*vy));
    lfxegm[3][3] = 0.0;
    */
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
    realtype rou, vx, vy, p, a, vxnext, vynext, pnext, vxcur, vycur, pcur, h, hcur, hnext;
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
      vx = -(Ydata[idx(i-1, j, Nx, Ny, 1)]+Ydata[idx(i-1, j, Nx, Ny, 1)])/(Ydata[idx(i-1, j, Nx, Ny, 0)]+Ydata[idx(i-1, j, Nx, Ny, 0)]);
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
      vx = -(Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
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
    /*
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
    */
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
static int Setlfyegm(realtype *Ydata, realtype **lfyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, realtype *egvy)
{
    realtype rou, vx, vy, p, a, vxnext, vynext, pnext, vxcur, vycur, pcur, h, hcur, hnext;
    
    vx = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
    vy = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
    p = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vx*vx+vy*vy);
    a = sqrt(gama*p/Ydata[idx(i, j, Nx, Ny, 0)]);
    h = (Ydata[idx(i, j, Nx, Ny, 3)]+p)/Ydata[idx(i, j, Nx, Ny, 0)];
    
    egvy[idx(i,j,Nx,Ny,0)]=vy-a;
    egvy[idx(i,j,Nx,Ny,1)]=vy;
    egvy[idx(i,j,Nx,Ny,2)]=vy;
    egvy[idx(i,j,Nx,Ny,3)]=vy+a;

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
    
    
    return 0;
}

/* fill in right eigenvectors for y component*/
static int Setrhyegm(realtype *Ydata, realtype **rhyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag)
{
    realtype rou, vx, vy, p, a, vxnext, vynext, pnext, vxcur, vycur, pcur, h, hcur, hnext;
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
      vx = -(Ydata[idx(i, j-1, Nx, Ny, 1)]+Ydata[idx(i, j-1, Nx, Ny, 1)])/(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]);
      vy = (Ydata[idx(i, j-1, Nx, Ny, 2)]+Ydata[idx(i, j-1, Nx, Ny, 2)])/(Ydata[idx(i, j-1, Nx, Ny, 0)]+Ydata[idx(i, j-1, Nx, Ny, 0)]);
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
      vx = -(Ydata[idx(i, j, Nx, Ny, 1)]+Ydata[idx(i, j, Nx, Ny, 1)])/(Ydata[idx(i, j, Nx, Ny, 0)]+Ydata[idx(i, j, Nx, Ny, 0)]);
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
        
        if (i>3&&i<Nx-2){
            for (n=0;n<4;n++){
                for (k=0;k<7;k++){
                    yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                    yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                }
            }
        }
        
        if (flag == 0){
            if (i<4){
                for (n=0;n<4;n++){
                    for (k=4;k<7;k++){
                        yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==0){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(Nx-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(Nx-3,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-4,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(Nx-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(Nx-3,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-4,j,Nx,Ny,n)];
                    }
                }
                
                if (i==1){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(Nx-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-3,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(Nx-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-3,j,Nx,Ny,n)];
                    }
                }
                
                if (i==2){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-2,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-2,j,Nx,Ny,n)];
                    }
                }

		if (i==3){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                    }
                }
            }
            
            if(i>Nx-3){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(0,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(0,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(0,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(0,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(1,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(1,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx){
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
            if (i<4){
                for (n=0;n<4;n++){
                    for (k=4;k<7;k++){
                        yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                    }
                }
                if (i==0){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i+3,j,Nx,Ny,n)];
                        yn[3+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
			yn[2+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i+3,j,Nx,Ny,n)];
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
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
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
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yxpdata[idx(i-1,j,Nx,Ny,3)];
                    //yp[1+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yp[0+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yn[2+3*7] = -yxndata[idx(i-1,j,Nx,Ny,3)];
                    //yn[1+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                    //yn[0+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                }
		if (i==3){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yxpdata[idx(i-1,j,Nx,Ny,3)];
                    //yp[1+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yp[0+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yn[2+3*7] = -yxndata[idx(i-1,j,Nx,Ny,3)];
                    //yn[1+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                    //yn[0+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                }
            }
            
            if(i>Nx-3){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yp[5+3*7] = -yxpdata[idx(i+2,j,Nx,Ny,3)];
                    //yp[6+3*7] = -yxpdata[idx(i+2,j,Nx,Ny,3)];
                    //yn[4+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yn[5+3*7] = -yxndata[idx(i+2,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yxndata[idx(i+2,j,Nx,Ny,3)];
                }
                if (i==Nx-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yn[4+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yp[5+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yn[5+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yp[6+3*7] = -yxpdata[idx(i,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yxndata[idx(i,j,Nx,Ny,3)];
                }
                if (i==Nx){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
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
            
            if (i<4){
                for (n=0;n<4;n++){
                    for (k=4;k<7;k++){
                        yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                    }
                }
                if (i==0){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i+3,j,Nx,Ny,n)];
                        yn[3+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
			yn[2+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i+3,j,Nx,Ny,n)];
                    }
                }
                if (i==1){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                    }
                }
                if (i==2){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                    }
                }
		if (i==3){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
                    }
                }
            }
            
            if(i>Nx-3){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                    }
                }
                if (i==Nx-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                    }
                }
                if (i==Nx){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
                    }
                }
            }
        }
        
        
        /* compute indicators of smoothness of rou, qx, qy, E */
        for (n=0;n<4;n++){
	  //printf("i=%li,j=%li,n=%li,yp[1+n*7]=%g,yp[2+n*7]=%g,yp[3+n*7]=%g,yp[4+n*7]=%g,yp[5+n*7]=%g,yn[2+n*7]=%g,yn[3+n*7]=%g,yn[4+n*7]=%g,yn[5+n*7]=%g,yn[6+n*7]=%g\n",i,j,n,yp[1+n*7],yp[2+n*7],yp[3+n*7],yp[4+n*7],yp[5+n*7],yn[2+n*7],yn[3+n*7],yn[4+n*7],yn[5+n*7],yn[6+n*7]);
            IS0_px[n]=(13.0/12.0)*(yp[1+n*7]-2.0*yp[2+n*7]+yp[3+n*7])*(yp[1+n*7]-2.0*yp[2+n*7]+yp[3+n*7])+(1.0/4.0)*(yp[1+n*7]-4.0*yp[2+n*7]+3.0*yp[3+n*7])*(yp[1+n*7]-4.0*yp[2+n*7]+3.0*yp[3+n*7]);
            
            IS1_px[n]=(13.0/12.0)*(yp[2+n*7]-2.0*yp[3+n*7]+yp[4+n*7])*(yp[2+n*7]-2.0*yp[3+n*7]+yp[4+n*7])+(1.0/4.0)*(yp[2+n*7]-yp[4+n*7])*(yp[2+n*7]-yp[4+n*7]);
            
            IS2_px[n]=(13.0/12.0)*(yp[3+n*7]-2.0*yp[4+n*7]+yp[5+n*7])*(yp[3+n*7]-2.0*yp[4+n*7]+yp[5+n*7])+(1.0/4.0)*(3.0*yp[3+n*7]-4.0*yp[4+n*7]+yp[5+n*7])*(3.0*yp[3+n*7]-4.0*yp[4+n*7]+yp[5+n*7]);
            
            IS0_nx[n]=(13.0/12.0)*(yn[4+n*7]-2.0*yn[5+n*7]+yn[6+n*7])*(yn[4+n*7]-2.0*yn[5+n*7]+yn[6+n*7])+(1.0/4.0)*(3.0*yn[4+n*7]-4.0*yn[5+n*7]+yn[6+n*7])*(3.0*yn[4+n*7]-4.0*yn[5+n*7]+yn[6+n*7]);
            
            IS1_nx[n]=(13.0/12.0)*(yn[3+n*7]-2.0*yn[4+n*7]+yn[5+n*7])*(yn[3+n*7]-2.0*yn[4+n*7]+yn[5+n*7])+(1.0/4.0)*(yn[3+n*7]-yn[5+n*7])*(yn[3+n*7]-yn[5+n*7]);
            
            IS2_nx[n]=(13.0/12.0)*(yn[2+n*7]-2.0*yn[3+n*7]+yn[4+n*7])*(yn[2+n*7]-2.0*yn[3+n*7]+yn[4+n*7])+(1.0/4.0)*(yn[2+n*7]-4.0*yn[3+n*7]+3.0*yn[4+n*7])*(yn[2+n*7]-4.0*yn[3+n*7]+3.0*yn[4+n*7]);
	    //printf("i=%li,j=%li,n=%li,IS0_px[n]=%g,IS1_px[n]=%g,IS2_px[n]=%g,IS0_nx[n]=%g,IS1_nx[n]=%g,IS2_nx[n]=%g\n",i,j,n,IS0_px[n],IS1_px[n],IS2_px[n],IS0_nx[n],IS1_nx[n],IS2_nx[n]);
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
        
        if (j>3&&j<Ny-2){
            for (n=0;n<4;n++){
                for (k=0;k<7;k++){
                    yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                    yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                }
            }
        }
        
        if (flag == 0){
            if (j<4){
                for (n=0;n<4;n++){
                    for (k=4;k<7;k++){
                        yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                    }
                }
                
                if (j==0){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,Ny-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,Ny-3,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-4,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,Ny-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,Ny-3,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-4,Nx,Ny,n)];
                    }
                }
                
                if (j==1){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,Ny-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-3,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,Ny-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-3,Nx,Ny,n)];
                    }
                }
                
                if (j==2){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-2,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-2,Nx,Ny,n)];
                    }
                }

		if (j==3){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                    }
                }
            }
            
            if(j>Ny-3){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,0,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,0,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,0,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,0,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,1,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,1,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny){
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
            if (j<4){
                for (n=0;n<4;n++){
                    for (k=4;k<7;k++){
                        yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                    }
                }
                if (j==0){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j+3,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j+3,Nx,Ny,n)];
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
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
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
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yypdata[idx(i,j-1,Nx,Ny,3)];
                    //yp[1+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yp[0+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yn[2+3*7] = -yyndata[idx(i,j-1,Nx,Ny,3)];
                    //yn[1+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                    //yn[0+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                }
		if (j==3){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yypdata[idx(i,j-1,Nx,Ny,3)];
                    //yp[1+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yp[0+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yn[2+3*7] = -yyndata[idx(i,j-1,Nx,Ny,3)];
                    //yn[1+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                    //yn[0+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                }
            }
            
            if(j>Ny-3){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yp[5+3*7] = -yypdata[idx(i,j+2,Nx,Ny,3)];
                    //yp[6+3*7] = -yypdata[idx(i,j+2,Nx,Ny,3)];
                    //yn[4+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yn[5+3*7] = -yyndata[idx(i,j+2,Nx,Ny,3)];
                    //yn[6+3*7] = -yyndata[idx(i,j+2,Nx,Ny,3)];
                }
                if (j==Ny-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yn[4+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yp[5+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yn[5+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yp[6+3*7] = -yypdata[idx(i,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yyndata[idx(i,j,Nx,Ny,3)];
                }
                if (j==Ny){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
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
            /* consider the situations of different i */
            if (j<4){
                for (n=0;n<4;n++){
                    for (k=4;k<7;k++){
                        yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                    }
                }
                if (j==0){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j+3,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j+3,Nx,Ny,n)];
                    }
                }
                if (j==1){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                    }
                }
                if (j==2){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                    }
                }
		if (j==3){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
                    }
                }
            }
            
            if(j>Ny-3){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                    }
                }
                if (j==Ny-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                    }
                }
                if (j==Ny){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
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
//static int SetUx(realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype *yxpdata, realtype *yxndata, realtype *u_tpphx, realtype *u_tnphx, realtype *u_tpnhx, realtype *u_tnnhx, long int i, long int j, long int Nx, long int Ny, int flag)
static int SetUx(realtype *w0_px, realtype *w1_px, realtype *w2_px, realtype *w0_nx, realtype *w1_nx, realtype *w2_nx, realtype *yxpdata, realtype *yxndata, realtype *u_tpphx, realtype *u_tnphx, long int i, long int j, long int Nx, long int Ny, int flag)
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
        
        if (i>3&&i<Nx-2){
            for (n=0;n<4;n++){
                for (k=0;k<7;k++){
                    yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                    yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                }
            }
        }
        
        if (flag == 0){
            if (i<4){
                for (n=0;n<4;n++){
                    for (k=4;k<7;k++){
                        yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==0){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(Nx-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(Nx-3,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-4,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(Nx-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(Nx-3,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-4,j,Nx,Ny,n)];
                    }
                }
                
                if (i==1){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(Nx-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-3,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(Nx-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-3,j,Nx,Ny,n)];
                    }
                }
                
                if (i==2){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-2,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-2,j,Nx,Ny,n)];
                    }
                }

		if (i==3){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(Nx-1,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(Nx-1,j,Nx,Ny,n)];
                    }
                }
            }
            
            if(i>Nx-3){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(0,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(0,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(0,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(0,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(1,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(1,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx){
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
            if (i<4){
                for (n=0;n<4;n++){
                    for (k=4;k<7;k++){
                        yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                    }
                }
                if (i==0){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i+3,j,Nx,Ny,n)];
                        yn[3+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
			yn[2+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i+3,j,Nx,Ny,n)];
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
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
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
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yxpdata[idx(i-1,j,Nx,Ny,3)];
                    //yp[1+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yp[0+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yn[2+3*7] = -yxndata[idx(i-1,j,Nx,Ny,3)];
                    //yn[1+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                    //yn[0+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                }
		if (i==3){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yxpdata[idx(i-1,j,Nx,Ny,3)];
                    //yp[1+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yp[0+3*7] = -yxpdata[idx(i-2,j,Nx,Ny,3)];
                    //yn[2+3*7] = -yxndata[idx(i-1,j,Nx,Ny,3)];
                    //yn[1+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                    //yn[0+3*7] = -yxndata[idx(i-2,j,Nx,Ny,3)];
                }
            }
            
            if(i>Nx-3){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yp[5+3*7] = -yxpdata[idx(i+2,j,Nx,Ny,3)];
                    //yp[6+3*7] = -yxpdata[idx(i+2,j,Nx,Ny,3)];
                    //yn[4+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yn[5+3*7] = -yxndata[idx(i+2,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yxndata[idx(i+2,j,Nx,Ny,3)];
                }
                if (i==Nx-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yn[4+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yp[5+3*7] = -yxpdata[idx(i+1,j,Nx,Ny,3)];
                    //yn[5+3*7] = -yxndata[idx(i+1,j,Nx,Ny,3)];
                    //yp[6+3*7] = -yxpdata[idx(i,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yxndata[idx(i,j,Nx,Ny,3)];
                }
                if (i==Nx){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
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
            if (i<4){
                for (n=0;n<4;n++){
                    for (k=4;k<7;k++){
                        yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                    }
                }
                if (i==0){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i+2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i+3,j,Nx,Ny,n)];
                        yn[3+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
			yn[2+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i+2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i+3,j,Nx,Ny,n)];
                    }
                }
                if (i==1){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                    }
                }
                if (i==2){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                    }
                }
		if (i==3){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yp[2+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yp[1+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
                        yp[0+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
			yn[3+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yn[2+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yn[1+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
                        yn[0+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
                    }
                }
            }
            
            if(i>Nx-3){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yxpdata[idx(i-4+k,j,Nx,Ny,n)];
                        yn[k+n*7] = yxndata[idx(i-4+k,j,Nx,Ny,n)];
                    }
                }
                
                if (i==Nx-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i+1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i+1,j,Nx,Ny,n)];
                    }
                }
                if (i==Nx-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                    }
                }
                if (i==Nx){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yxpdata[idx(i-1,j,Nx,Ny,n)];
                        yn[4+n*7] = yxndata[idx(i-1,j,Nx,Ny,n)];
                        yp[5+n*7] = yxpdata[idx(i-2,j,Nx,Ny,n)];
                        yn[5+n*7] = yxndata[idx(i-2,j,Nx,Ny,n)];
                        yp[6+n*7] = yxpdata[idx(i-3,j,Nx,Ny,n)];
                        yn[6+n*7] = yxndata[idx(i-3,j,Nx,Ny,n)];
                    }
                }
            }
        }
        
        /* compute positive and negative solutions on the interface */
        for (n=0;n<4;n++){
            u_tpphx[n]=w0_px[n]*((2.0/6.0)*yp[1+n*7]-(7.0/6.0)*yp[2+n*7]+(11.0/6.0)*yp[3+n*7])+w1_px[n]*((-1.0/6.0)*yp[2+n*7]+(5.0/6.0)*yp[3+n*7]+(2.0/6.0)*yp[4+n*7])+w2_px[n]*((2.0/6.0)*yp[3+n*7]+(5.0/6.0)*yp[4+n*7]-(1.0/6.0)*yp[5+n*7]);
            
            u_tnphx[n]=w2_nx[n]*((-1.0/6.0)*yn[2+n*7]+(5.0/6.0)*yn[3+n*7]+(2.0/6.0)*yn[4+n*7])+w1_nx[n]*((2.0/6.0)*yn[3+n*7]+(5.0/6.0)*yn[4+n*7]-(1.0/6.0)*yn[5+n*7])+w0_nx[n]*((11.0/6.0)*yn[4+n*7]-(7.0/6.0)*yn[5+n*7]+(2.0/6.0)*yn[6+n*7]);
            
            //u_tpnhx[n]=w0_px[n]*((2.0/6.0)*yp[0+n*7]-(7.0/6.0)*yp[1+n*7]+(11.0/6.0)*yp[2+n*7])+w1_px[n]*((-1.0/6.0)*yp[1+n*7]+(5.0/6.0)*yp[2+n*7]+(2.0/6.0)*yp[3+n*7])+w2_px[n]*((2.0/6.0)*yp[2+n*7]+(5.0/6.0)*yp[3+n*7]-(1.0/6.0)*yp[4+n*7]);
            
            //u_tnnhx[n]=w2_nx[n]*((-1.0/6.0)*yn[1+n*7]+(5.0/6.0)*yn[2+n*7]+(2.0/6.0)*yn[3+n*7])+w1_nx[n]*((2.0/6.0)*yn[2+n*7]+(5.0/6.0)*yn[3+n*7]-(1.0/6.0)*yn[4+n*7])+w0_nx[n]*((11.0/6.0)*yn[3+n*7]-(7.0/6.0)*yn[4+n*7]+(2.0/6.0)*yn[5+n*7]);
	    
        }
        
        delete []yp;
        delete []yn;
        
        return 0;
    }
}

/* Get the derivative on y direction */
//static int SetUy(realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype *yypdata, realtype *yyndata, realtype *u_tpphy, realtype *u_tnphy, realtype *u_tpnhy, realtype *u_tnnhy, long int i, long int j, long int Nx, long int Ny, int flag)
static int SetUy(realtype *w0_py, realtype *w1_py, realtype *w2_py, realtype *w0_ny, realtype *w1_ny, realtype *w2_ny, realtype *yypdata, realtype *yyndata, realtype *u_tpphy, realtype *u_tnphy, long int i, long int j, long int Nx, long int Ny, int flag)
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
        
        if (j>3&&j<Ny-2){
            for (n=0;n<4;n++){
                for (k=0;k<7;k++){
                    yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                    yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                }
            }
        }
        
        if (flag == 0){
            if (j<4){
                for (n=0;n<4;n++){
                    for (k=4;k<7;k++){
                        yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                    }
                }
                
                if (j==0){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,Ny-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,Ny-3,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-4,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,Ny-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,Ny-3,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-4,Nx,Ny,n)];
                    }
                }
                
                if (j==1){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,Ny-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-3,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,Ny-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-3,Nx,Ny,n)];
                    }
                }
                
                if (j==2){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-2,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-2,Nx,Ny,n)];
                    }
                }

		if (j==3){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,Ny-1,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,Ny-1,Nx,Ny,n)];
                    }
                }
            }
            
            if(j>Ny-3){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,0,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,0,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,0,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,0,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,1,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,1,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny){
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
            if (j<4){
                for (n=0;n<4;n++){
                    for (k=4;k<7;k++){
                        yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                    }
                }
                if (j==0){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j+3,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j+3,Nx,Ny,n)];
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
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
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
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yypdata[idx(i,j-1,Nx,Ny,3)];
                    //yp[1+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yp[0+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yn[2+3*7] = -yyndata[idx(i,j-1,Nx,Ny,3)];
                    //yn[1+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                    //yn[0+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                }
		if (j==3){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
                    }
                    //yp[2+3*7] = -yypdata[idx(i,j-1,Nx,Ny,3)];
                    //yp[1+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yp[0+3*7] = -yypdata[idx(i,j-2,Nx,Ny,3)];
                    //yn[2+3*7] = -yyndata[idx(i,j-1,Nx,Ny,3)];
                    //yn[1+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                    //yn[0+3*7] = -yyndata[idx(i,j-2,Nx,Ny,3)];
                }
            }
            
            if(j>Ny-3){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yp[5+3*7] = -yypdata[idx(i,j+2,Nx,Ny,3)];
                    //yp[6+3*7] = -yypdata[idx(i,j+2,Nx,Ny,3)];
                    //yn[4+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yn[5+3*7] = -yyndata[idx(i,j+2,Nx,Ny,3)];
                    //yn[6+3*7] = -yyndata[idx(i,j+2,Nx,Ny,3)];
                }
                if (j==Ny-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                    }
                    //yp[4+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yn[4+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yp[5+3*7] = -yypdata[idx(i,j+1,Nx,Ny,3)];
                    //yn[5+3*7] = -yyndata[idx(i,j+1,Nx,Ny,3)];
                    //yp[6+3*7] = -yypdata[idx(i,j,Nx,Ny,3)];
                    //yn[6+3*7] = -yyndata[idx(i,j,Nx,Ny,3)];
                }
                if (j==Ny){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
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
            /* consider the situations of different i */
            if (j<4){
                for (n=0;n<4;n++){
                    for (k=4;k<7;k++){
                        yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                    }
                }
                if (j==0){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j+2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j+3,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j+2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j+3,Nx,Ny,n)];
                    }
                }
                if (j==1){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                    }
                }
                if (j==2){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                    }
                }
		if (j==3){
                    for (n=0;n<4;n++){
		        yp[3+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yp[2+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yp[1+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
                        yp[0+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
			yn[3+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yn[2+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yn[1+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
                        yn[0+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
                    }
                }
            }
            
            if(j>Ny-3){
                for (n=0;n<4;n++){
                    for (k=0;k<4;k++){
                        yp[k+n*7] = yypdata[idx(i,j-4+k,Nx,Ny,n)];
                        yn[k+n*7] = yyndata[idx(i,j-4+k,Nx,Ny,n)];
                    }
                }
                
                if (j==Ny-2){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j+1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j+1,Nx,Ny,n)];
                    }
                }
                if (j==Ny-1){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                    }
                }
                if (j==Ny){
                    for (n=0;n<4;n++){
                        yp[4+n*7] = yypdata[idx(i,j-1,Nx,Ny,n)];
                        yn[4+n*7] = yyndata[idx(i,j-1,Nx,Ny,n)];
                        yp[5+n*7] = yypdata[idx(i,j-2,Nx,Ny,n)];
                        yn[5+n*7] = yyndata[idx(i,j-2,Nx,Ny,n)];
                        yp[6+n*7] = yypdata[idx(i,j-3,Nx,Ny,n)];
                        yn[6+n*7] = yyndata[idx(i,j-3,Nx,Ny,n)];
                    }
                }
            }
        }
        
        /* compute positive and negative solutions on the interface */
        for (n=0;n<4;n++){
            u_tpphy[n]=w0_py[n]*((2.0/6.0)*yp[1+n*7]-(7.0/6.0)*yp[2+n*7]+(11.0/6.0)*yp[3+n*7])+w1_py[n]*((-1.0/6.0)*yp[2+n*7]+(5.0/6.0)*yp[3+n*7]+(2.0/6.0)*yp[4+n*7])+w2_py[n]*((2.0/6.0)*yp[3+n*7]+(5.0/6.0)*yp[4+n*7]-(1.0/6.0)*yp[5+n*7]);
            
            u_tnphy[n]=w2_ny[n]*((-1.0/6.0)*yn[2+n*7]+(5.0/6.0)*yn[3+n*7]+(2.0/6.0)*yn[4+n*7])+w1_ny[n]*((2.0/6.0)*yn[3+n*7]+(5.0/6.0)*yn[4+n*7]-(1.0/6.0)*yn[5+n*7])+w0_ny[n]*((11.0/6.0)*yn[4+n*7]-(7.0/6.0)*yn[5+n*7]+(2.0/6.0)*yn[6+n*7]);
            
            //u_tpnhy[n]=w0_py[n]*((2.0/6.0)*yp[0+n*7]-(7.0/6.0)*yp[1+n*7]+(11.0/6.0)*yp[2+n*7])+w1_py[n]*((-1.0/6.0)*yp[1+n*7]+(5.0/6.0)*yp[2+n*7]+(2.0/6.0)*yp[3+n*7])+w2_py[n]*((2.0/6.0)*yp[2+n*7]+(5.0/6.0)*yp[3+n*7]-(1.0/6.0)*yp[4+n*7]);
            
            //u_tnnhy[n]=w2_ny[n]*((-1.0/6.0)*yn[1+n*7]+(5.0/6.0)*yn[2+n*7]+(2.0/6.0)*yn[3+n*7])+w1_ny[n]*((2.0/6.0)*yn[2+n*7]+(5.0/6.0)*yn[3+n*7]-(1.0/6.0)*yn[4+n*7])+w0_ny[n]*((11.0/6.0)*yn[3+n*7]-(7.0/6.0)*yn[4+n*7]+(2.0/6.0)*yn[5+n*7]);
        }
        
        delete []yp;
        delete []yn;
        
        return 0;
    }
}

/* Fill in values of tao in the whole domain */
static int Gettao(realtype *yxbackdata, realtype *yybackdata, realtype *tao, long int Nx, long int Ny, realtype gama)
{
    /* declare parameters */
    long int i, j;
    //realtype px, py;
   
    realtype *vxx = new realtype [(Nx+1)*Ny];
    realtype *vyx = new realtype [(Nx+1)*Ny];
    realtype *vxy = new realtype [Nx*(Ny+1)];
    realtype *vyy = new realtype [Nx*(Ny+1)];
    realtype *px = new realtype [(Nx+1)*Ny];
    realtype *py = new realtype [Nx*(Ny+1)];

    /* Set values into vx and vy */
    for(j=0;j<Ny;j++){
      for(i=0;i<(Nx+1);i++){
	vxx[idx_v(i,j,(Nx+1))]=yxbackdata[idx(i, j, (Nx+1), Ny, 1)]/yxbackdata[idx(i, j, (Nx+1), Ny, 0)];
	//printf("i = %li, j=%li, yxbackdata[idx(i, j, (Nx+1), Ny, 1)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",i,j,yxbackdata[idx(i, j, (Nx+1), Ny, 1)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
	vyx[idx_v(i,j,(Nx+1))]=yxbackdata[idx(i, j, (Nx+1), Ny, 2)]/yxbackdata[idx(i, j, (Nx+1), Ny, 0)];
	//printf("yxbackdata[idx(i, j, (Nx+1), Ny, 2)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",yxbackdata[idx(i, j, (Nx+1), Ny, 2)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
	//printf("i = %li, j=%li, vxx[idx_v(i,j,(Nx+1))]=%f,vyx[idx_v(i,j,(Nx+1))]=%f\n",i,j,vxx[idx_v(i,j,(Nx+1))],vyx[idx_v(i,j,(Nx+1))]);
        }
    }
    
    for(i=0;i<Nx;i++){
      for(j=0;j<(Ny+1);j++){
	vxy[idx_v(i,j,Nx)]=yybackdata[idx(i, j, Nx, (Ny+1), 1)]/yybackdata[idx(i, j, Nx, (Ny+1), 0)];
	//printf("i = %li, j=%li, yybackdata[idx(i, j, Nx, (Ny+1), 1)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",i,j,yybackdata[idx(i, j, Nx, (Ny+1), 1)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
        vyy[idx_v(i,j,Nx)]=yybackdata[idx(i, j, Nx, (Ny+1), 2)]/yybackdata[idx(i, j, Nx, (Ny+1), 0)];
	//printf("yybackdata[idx(i, j, Nx, (Ny+1), 2)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",yybackdata[idx(i, j, Nx, (Ny+1), 2)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
        }
    }

    /* Compute the values of tao in the whole domain */
    for(j=0;j<Ny;j++){
      for(i=0;i<(Nx+1);i++){
	px[idx_v(i,j,(Nx+1))] =  yxbackdata[idx(i, j, (Nx+1), Ny, 0)]*(gama-1.0)*(yxbackdata[idx(i, j, (Nx+1), Ny, 3)]/yxbackdata[idx(i, j, (Nx+1), Ny, 0)]-0.5*(vxx[idx_v(i,j,(Nx+1))]*vxx[idx_v(i,j,(Nx+1))]+vyx[idx_v(i,j,(Nx+1))]*vyx[idx_v(i,j,(Nx+1))]));	   
	tao[idx_v(i,j,(Nx+1))] = yxbackdata[idx(i, j, (Nx+1), Ny, 1)]*vxx[idx_v(i,j,(Nx+1))]+px[idx_v(i,j,(Nx+1))];
	tao[idx_v(i,j,(Nx+1))+(Nx+1)*Ny] = yxbackdata[idx(i, j, (Nx+1), Ny, 1)]*vyx[idx_v(i,j,(Nx+1))];
	     } 
	}
	     for (i=0;i<Nx;i++){
	       for (j=0;j<(Ny+1);j++){
		 py[idx_v(i,j,Nx)] =  yybackdata[idx(i, j, Nx, (Ny+1), 0)]*(gama-1.0)*(yybackdata[idx(i, j, Nx, (Ny+1), 3)]/yybackdata[idx(i, j, Nx, (Ny+1), 0)]-0.5*(vyx[idx_v(i,j,Nx)]*vyx[idx_v(i,j,Nx)]+vyy[idx_v(i,j,Nx)]*vyy[idx_v(i,j,Nx)]));	   
		 tao[idx_v(i,j,Nx)+2*((Nx+1)*Ny)] = yybackdata[idx(i, j, Nx, (Ny+1), 2)]*vxy[idx_v(i,j,Nx)];
		 tao[idx_v(i,j,Nx)+2*((Nx+1)*Ny)+Nx*(Ny+1)] = yybackdata[idx(i, j, Nx, (Ny+1), 2)]*vyy[idx_v(i,j,Nx)]+py[idx_v(i,j,Nx)];
	     }
	}
    
	       delete []vxx;
	       delete []vxy;
	       delete []vyx;
	       delete []vyy;
	       delete []px;
	       delete []py;
    
	       return 0;
}

/* Fill in the value of J in the whole domain */
static int GetCj(realtype *yxbackdata, realtype *yybackdata, realtype *Cj, long int Nx, long int Ny, realtype gama)
{
    /* declare parameters */
    long int i, j;
    
    realtype *vxx = new realtype [(Nx+1)*Ny];
    realtype *vyx = new realtype [(Nx+1)*Ny];
    realtype *vxy = new realtype [Nx*(Ny+1)];
    realtype *vyy = new realtype [Nx*(Ny+1)];
    realtype *px = new realtype [(Nx+1)*Ny];
    realtype *py = new realtype [Nx*(Ny+1)];

    /* Set values into vx and vy */
    for(j=0;j<Ny;j++){
      for(i=0;i<(Nx+1);i++){
	vxx[idx_v(i,j,(Nx+1))]=yxbackdata[idx(i, j, (Nx+1), Ny, 1)]/yxbackdata[idx(i, j, (Nx+1), Ny, 0)];
	    //printf("i = %li, j=%li, yxupdata[idx(i, j, Nx, Ny, 1)]=%f,yxupdata[idx(i, j, Nx, Ny, 0)]=%f\n",i,j,yxupdata[idx(i, j, Nx, Ny, 1)],yxupdata[idx(i, j, Nx, Ny, 0)]);
	vyx[idx_v(i,j,(Nx+1))]=yxbackdata[idx(i, j, (Nx+1), Ny, 2)]/yxbackdata[idx(i, j, (Nx+1), Ny, 0)];
	    // printf("yydowndata[idx(i, j, Nx, Ny, 2)]=%f,yydowndata[idx(i, j, Nx, Ny, 0)]=%f\n",yydowndata[idx(i, j, Nx, Ny, 2)],yydowndata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    
    for(i=0;i<Nx;i++){
      for(j=0;j<(Ny+1);j++){
	vxy[idx_v(i,j,Nx)]=yybackdata[idx(i, j, Nx, (Ny+1), 1)]/yybackdata[idx(i, j, Nx, (Ny+1), 0)];
	    //printf("i = %li, j=%li, yxupdata[idx(i, j, Nx, Ny, 1)]=%f,yxupdata[idx(i, j, Nx, Ny, 0)]=%f\n",i,j,yxupdata[idx(i, j, Nx, Ny, 1)],yxupdata[idx(i, j, Nx, Ny, 0)]);
        vyy[idx_v(i,j,Nx)]=yybackdata[idx(i, j, Nx, (Ny+1), 2)]/yybackdata[idx(i, j, Nx, (Ny+1), 0)];
	    // printf("yydowndata[idx(i, j, Nx, Ny, 2)]=%f,yydowndata[idx(i, j, Nx, Ny, 0)]=%f\n",yydowndata[idx(i, j, Nx, Ny, 2)],yydowndata[idx(i, j, Nx, Ny, 0)]);
        }
    }
    
     /* Compute the values of tao in the whole domain */
    for(j=0;j<Ny;j++){
      for(i=0;i<(Nx+1);i++){
	px[idx_v(i,j,(Nx+1))] =  yxbackdata[idx(i, j, (Nx+1), Ny, 0)]*(gama-1.0)*(yxbackdata[idx(i, j, (Nx+1), Ny, 3)]/yxbackdata[idx(i, j, (Nx+1), Ny, 0)]-0.5*(vxx[idx_v(i,j,(Nx+1))]*vxx[idx_v(i,j,(Nx+1))]+vyx[idx_v(i,j,(Nx+1))]*vyx[idx_v(i,j,(Nx+1))]));	   
	Cj[idx_v(i,j,(Nx+1))] = (yxbackdata[idx(i, j, (Nx+1), Ny, 3)]+px[idx_v(i,j,(Nx+1))])*vxx[idx_v(i,j,(Nx+1))];
	     } 
	}
	     for (i=0;i<Nx;i++){
	       for (j=0;j<(Ny+1);j++){
		 py[idx_v(i,j,Nx)] =  yybackdata[idx(i, j, Nx, (Ny+1), 0)]*(gama-1.0)*(yybackdata[idx(i, j, Nx, (Ny+1), 3)]/yybackdata[idx(i, j, Nx, (Ny+1), 0)]-0.5*(vyx[idx_v(i,j,Nx)]*vyx[idx_v(i,j,Nx)]+vyy[idx_v(i,j,Nx)]*vyy[idx_v(i,j,Nx)]));	   
		 Cj[idx_v(i,j,Nx)+(Nx+1)*Ny] = (yybackdata[idx(i, j, Nx, (Ny+1), 3)]+py[idx_v(i,j,Nx)])*vyy[idx_v(i,j,Nx)];
	     }
	}

	       delete []vxx;
	       delete []vxy;
	       delete []vyx;
	       delete []vyy;
	       delete []px;
	       delete []py;
    
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
