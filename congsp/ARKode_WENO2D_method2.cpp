/*
Programmer: Cong Zhang
Example:
u_t+f(u)_x+g(u)_y = 0，f(u) and g(u) come from the Euler equations (fluid dynamics) with 2D Riemann test initial conditions (different initial conditions in different quadrants) and natural boundary conditions.
This program solves the problem with WENO method in paper (Jun Luo, Lijun Xuan and Kun Xu, Comparison of Fifth-Order WENO Scheme and Finite Volume WENO-Gas-Kinetic Scheme for Inviscid and Viscous Flow Simulation), especially section 3.2.1 and section 2. 
All the parameters are provided in the input file input_WENO2D.txt.
*/
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
static int Gettao(realtype *fupback, realtype *fdownback, realtype *gupback, realtype *gdownback, realtype *taoup, realtype *taodown, long int Nx, long int Ny, realtype gama);

/* Set value of J in whole domain*/
static int GetCj(realtype *fupback, realtype *fdownback, realtype *gupback, realtype *gdownback, realtype *Cjup, realtype *Cjdown, long int Nx, long int Ny, realtype gama);

/* Split f to f_positve and f_negative parts for x component */
static int Splitfluxesx(realtype *Ydata, realtype *fp, realtype *fn, long int i, long int j, long int Nx, long int Ny, realtype gama);

/* Split g to g_positve and g_negative parts for y component */
static int Splitfluxesy(realtype *Ydata, realtype *gp, realtype *gn, long int i, long int j, long int Nx, long int Ny, realtype gama);

/* Set left eigenvectors in x direction */
static int Setlfxegm(realtype *Ydata, realtype **lfxegm, long int i, long int j, long int Nx, long int Ny, realtype gama);

/* Set right eigenvectors in x direction */
static int Setrhxegm(realtype *Ydata, realtype **rhxegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag);

/* Set left eigenvectors in y direction */
static int Setlfyegm(realtype *Ydata, realtype **lfyegm, long int i, long int j, long int Nx, long int Ny, realtype gama);

/* Set right eigenvectors in y direction */
static int Setrhyegm(realtype *Ydata, realtype **rhyegm, long int i, long int j, long int Nx, long int Ny, realtype gama, int flag);

/* Fill in the indicators of smoothness on x direction */
static int SetISX(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx, realtype *yxpdata, realtype *yxndata, long int i, long int j, long int Nx, long int Ny, int flag);

/* Fill in the indicators of smoothness on y direction */
static int SetISY(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny, realtype *yypdata, realtype *yyndata, long int i, long int j, long int Nx, long int Ny, int flag);

/* Fill in the stencil weights on x direction */
static int Setalphawx(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx, realtype *alpha_0prx, realtype *alpha_1prx, realtype *alpha_2prx, realtype *alpha_0plx, realtype *alpha_1plx, realtype *alpha_2plx, realtype *alpha_0nrx, realtype *alpha_1nrx, realtype *alpha_2nrx, realtype *alpha_0nlx, realtype *alpha_1nlx, realtype *alpha_2nlx, realtype *w0_prx, realtype *w1_prx, realtype *w2_prx, realtype *w0_plx, realtype *w1_plx, realtype *w2_plx, realtype *w0_nrx, realtype *w1_nrx, realtype *w2_nrx, realtype *w0_nlx, realtype *w1_nlx, realtype *w2_nlx, realtype Epsilon);

/* Fill in the stencil weights on y direction */
static int Setalphawy(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny, realtype *alpha_0pry, realtype *alpha_1pry, realtype *alpha_2pry, realtype *alpha_0ply, realtype *alpha_1ply, realtype *alpha_2ply, realtype *alpha_0nry, realtype *alpha_1nry, realtype *alpha_2nry, realtype *alpha_0nly, realtype *alpha_1nly, realtype *alpha_2nly, realtype *w0_pry, realtype *w1_pry, realtype *w2_pry, realtype *w0_ply, realtype *w1_ply, realtype *w2_ply, realtype *w0_nry, realtype *w1_nry, realtype *w2_nry, realtype *w0_nly, realtype *w1_nly, realtype *w2_nly, realtype Epsilon);

/* Get the derivative on x direction */
static int SetUx(realtype *w0_prx, realtype *w1_prx, realtype *w2_prx, realtype *w0_plx, realtype *w1_plx, realtype *w2_plx, realtype *w0_nrx, realtype *w1_nrx, realtype *w2_nrx, realtype *w0_nlx, realtype *w1_nlx, realtype *w2_nlx, realtype *yxpdata, realtype *yxndata, realtype *u_tprhx, realtype *u_tnrhx, realtype *u_tplhx, realtype *u_tnlhx, long int i, long int j, long int Nx, long int Ny, int flag);

/* Get the derivative on y direction */
static int SetUy(realtype *w0_pry, realtype *w1_pry, realtype *w2_pry, realtype *w0_ply, realtype *w1_ply, realtype *w2_ply, realtype *w0_nry, realtype *w1_nry, realtype *w2_nry, realtype *w0_nly, realtype *w1_nly, realtype *w2_nly, realtype *yypdata, realtype *yyndata, realtype *u_tprhy, realtype *u_tnrhy, realtype *u_tplhy, realtype *u_tnlhy, long int i, long int j, long int Nx, long int Ny, int flag);

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
    long int i, j, k;
    int flag;
    realtype  Epsilon;
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
    
    realtype *alpha_0prx = new realtype [4];
    realtype *alpha_1prx = new realtype [4];
    realtype *alpha_2prx = new realtype [4];
    realtype *alpha_0nrx = new realtype [4];
    realtype *alpha_1nrx = new realtype [4];
    realtype *alpha_2nrx = new realtype [4];
    realtype *alpha_0pry = new realtype [4];
    realtype *alpha_1pry = new realtype [4];
    realtype *alpha_2pry = new realtype [4];
    realtype *alpha_0nry = new realtype [4];
    realtype *alpha_1nry = new realtype [4];
    realtype *alpha_2nry = new realtype [4];
    realtype *alpha_0plx = new realtype [4];
    realtype *alpha_1plx = new realtype [4];
    realtype *alpha_2plx = new realtype [4];
    realtype *alpha_0nlx = new realtype [4];
    realtype *alpha_1nlx = new realtype [4];
    realtype *alpha_2nlx = new realtype [4];
    realtype *alpha_0ply = new realtype [4];
    realtype *alpha_1ply = new realtype [4];
    realtype *alpha_2ply = new realtype [4];
    realtype *alpha_0nly = new realtype [4];
    realtype *alpha_1nly = new realtype [4];
    realtype *alpha_2nly = new realtype [4];

    realtype *w0_prx = new realtype [4];
    realtype *w1_prx = new realtype [4];
    realtype *w2_prx = new realtype [4];
    realtype *w0_nrx = new realtype [4];
    realtype *w1_nrx = new realtype [4];
    realtype *w2_nrx = new realtype [4];
    realtype *w0_pry = new realtype [4];
    realtype *w1_pry = new realtype [4];
    realtype *w2_pry = new realtype [4];
    realtype *w0_nry = new realtype [4];
    realtype *w1_nry = new realtype [4];
    realtype *w2_nry = new realtype [4];
    realtype *w0_plx = new realtype [4];
    realtype *w1_plx = new realtype [4];
    realtype *w2_plx = new realtype [4];
    realtype *w0_nlx = new realtype [4];
    realtype *w1_nlx = new realtype [4];
    realtype *w2_nlx = new realtype [4];
    realtype *w0_ply = new realtype [4];
    realtype *w1_ply = new realtype [4];
    realtype *w2_ply = new realtype [4];
    realtype *w0_nly = new realtype [4];
    realtype *w1_nly = new realtype [4];
    realtype *w2_nly = new realtype [4];

    realtype *u_tprhx = new realtype [4];
    realtype *u_tplhx = new realtype [4];
    realtype *u_tnrhx = new realtype [4];
    realtype *u_tnlhx = new realtype [4];
    realtype *u_tprhy = new realtype [4];
    realtype *u_tplhy = new realtype [4];
    realtype *u_tnrhy = new realtype [4];
    realtype *u_tnlhy = new realtype [4];
    
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

    /* access data arrays */
    realtype *Ydata = N_VGetArrayPointer(y);
    if (check_flag((void *) Ydata, "N_VGetArrayPointer", 0)) return 1;
    realtype *dYdata = N_VGetArrayPointer(ydot);
    if (check_flag((void *) dYdata, "N_VGetArrayPointer", 0)) return 1;

    /* declare relative arrays */
    /* up means the right interface and down means the left interface for the same grid */
    realtype *fp = new realtype [4*Nx*Ny];
    realtype *fn = new realtype [4*Nx*Ny];
    realtype *gp = new realtype [4*Nx*Ny];
    realtype *gn = new realtype [4*Nx*Ny];
    realtype *fpnew = new realtype [4*Nx*Ny];
    realtype *fnnew = new realtype [4*Nx*Ny];
    realtype *gpnew = new realtype [4*Nx*Ny];
    realtype *gnnew = new realtype [4*Nx*Ny];
    realtype *fup = new realtype [4*Nx*Ny];
    realtype *fdown = new realtype [4*Nx*Ny];
    realtype *gup = new realtype [4*Nx*Ny];
    realtype *gdown = new realtype [4*Nx*Ny];
    realtype *fupback = new realtype [4*Nx*Ny];
    realtype *fdownback = new realtype [4*Nx*Ny];
    realtype *gupback = new realtype [4*Nx*Ny];
    realtype *gdownback = new realtype [4*Nx*Ny];
    realtype *taoup = new realtype [4*Nx*Ny];
    realtype *taodown = new realtype [4*Nx*Ny];
    realtype *Cjup = new realtype [2*Nx*Ny];
    realtype *Cjdown = new realtype [2*Nx*Ny];

    /* Split f to f_positve and f_negative parts for x component */
    for(j=0;j<Ny;j++){
      for(i=0;i<Nx;i++){
	flag = Splitfluxesx(Ydata, fp, fn, i, j, Nx, Ny, gama);
	if (flag!=0) printf("error in Splitfluxesx function \n");
	//printf(" fp : i = %li, j = %li, fp[idx(i, j, Nx, Ny, 0)] = %f, fp[idx(i, j, Nx, Ny, 1)] = %f, fp[idx(i, j, Nx, Ny, 2)] = %f, fp[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, fp[idx(i, j, Nx, Ny, 0)], fp[idx(i, j, Nx, Ny, 1)],fp[idx(i, j, Nx, Ny, 2)],fp[idx(i, j, Nx, Ny, 3)]);
	//printf(" fn : i = %li, j = %li, fn[idx(i, j, Nx, Ny, 0)] = %f, fn[idx(i, j, Nx, Ny, 1)] = %f, fn[idx(i, j, Nx, Ny, 2)] = %f, fn[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, fn[idx(i, j, Nx, Ny, 0)], fn[idx(i, j, Nx, Ny, 1)],fn[idx(i, j, Nx, Ny, 2)],fn[idx(i, j, Nx, Ny, 3)]);
      }
    }

    /* Split g to g_positve and g_negative parts for y component */
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	flag = Splitfluxesy(Ydata, gp, gn, i, j, Nx, Ny, gama);
	if (flag!=0) printf("error in Splitfluxesy function \n");
      }
    }

    /* compute inverse eigenvector matrix and fill in fpnew, fnnew */
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
	    flag = Setlfxegm(Ydata, lfxegm, i, j, Nx, Ny, gama);
            if (flag!=0) printf("error in Setlfxegm function \n");
	   
	    //printf(" Ydata : i = %li, j = %li, Ydata[idx(i, j, Nx, Ny, 0)] = %f, Ydata[idx(i, j, Nx, Ny, 1)] = %f, Ydata[idx(i, j, Nx, Ny, 2)] = %f, Ydata[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, Ydata[idx(i, j, Nx, Ny, 0)], Ydata[idx(i, j, Nx, Ny, 1)],Ydata[idx(i, j, Nx, Ny, 2)],Ydata[idx(i, j, Nx, Ny, 3)]);
	       
            for (k=0;k<4;k++){
	     
	      //printf(" lfxegm : i = %li, j = %li, lfxegm[k][0] = %f, lfxegm[k][1] = %f, lfxegm[k][2] = %f, lfxegm[k][3] = %f\n", i, j, lfxegm[k][0], lfxegm[k][1],lfxegm[k][2],lfxegm[k][3]);
	     
                fpnew[idx(i, j, Nx, Ny, k)] = lfxegm[k][0]*fp[idx(i, j, Nx, Ny, 0)]+lfxegm[k][1]*fp[idx(i, j, Nx, Ny, 1)]+lfxegm[k][2]*fp[idx(i, j, Nx, Ny, 2)]+lfxegm[k][3]*fp[idx(i, j, Nx, Ny, 3)];
		fnnew[idx(i, j, Nx, Ny, k)] = lfxegm[k][0]*fn[idx(i, j, Nx, Ny, 0)]+lfxegm[k][1]*fn[idx(i, j, Nx, Ny, 1)]+lfxegm[k][2]*fn[idx(i, j, Nx, Ny, 2)]+lfxegm[k][3]*fn[idx(i, j, Nx, Ny, 3)];
            }
            //printf("   yxnew problem parameters: i = %li, j = %li, yxnew0 = %e,  yxnew1 = %e, yxnew2 = %e,  yxnew3 = %e\n", i, j, yxnewdata[idx(i, j, Nx, Ny, 0)], yxnewdata[idx(i, j, Nx, Ny, 1)], yxnewdata[idx(i, j, Nx, Ny, 2)], yxnewdata[idx(i, j, Nx, Ny, 3)]);
           
        }
    }
    
    /* compute inverse eigenvector matrix and fill in gpnew, gnnew */
    for(i=0; i<Nx; i++){
        for (j=0; j<Ny; j++){
	  flag = Setlfyegm(Ydata, lfyegm, i, j, Nx, Ny, gama);
            if (flag!=0) printf("error in Setlfyegm function \n");

            for (k=0;k<4;k++){
                gpnew[idx(i, j, Nx, Ny, k)] = lfyegm[k][0]*gp[idx(i, j, Nx, Ny, 0)]+lfyegm[k][1]*gp[idx(i, j, Nx, Ny, 1)]+lfyegm[k][2]*gp[idx(i, j, Nx, Ny, 2)]+lfyegm[k][3]*gp[idx(i, j, Nx, Ny, 3)];
		gnnew[idx(i, j, Nx, Ny, k)] = lfyegm[k][0]*gn[idx(i, j, Nx, Ny, 0)]+lfyegm[k][1]*gn[idx(i, j, Nx, Ny, 1)]+lfyegm[k][2]*gn[idx(i, j, Nx, Ny, 2)]+lfyegm[k][3]*gn[idx(i, j, Nx, Ny, 3)];
            }
	    // printf("    problem parameters:  i = %li, j=%li, yypdata0 = %g,  yypdata1 = %g, yypdata2 = %g,  yypdata3 = %g\n", i, j, yypdata[idx(i, j, Nx, Ny, 0)], yypdata[idx(i, j, Nx, Ny, 1)], yypdata[idx(i, j, Nx, Ny, 2)], yypdata[idx(i, j, Nx, Ny, 3)]);
        }
    }
    
    /* iterate over domain, computing all equations */
    for(j=0; j<Ny; j++){
      for (i=0; i<Nx; i++){
            /* get indicator of smoothness on x direction */
            flag = SetISX(IS0_px, IS1_px, IS2_px, IS0_nx, IS1_nx, IS2_nx, fpnew, fnnew, i, j, Nx, Ny, 0);
            if (flag!=0) printf("error in SetISX function \n");
	    //for (n=0;n<4;n++){
	      // printf("i=%li,j=%li,n=%li,IS0_px[n]=%f,IS1_px[n]=%f,IS2_px[n]=%f,IS0_nx[n]=%f,IS1_nx[n]=%f,IS2_nx[n]=%f\n",i,j,n,IS0_px[n],IS1_px[n],IS2_px[n],IS0_nx[n],IS1_nx[n],IS2_nx[n]);
	    // }
	    /* get weight on x direction for each variable */
            flag = Setalphawx(IS0_px, IS1_px, IS2_px, IS0_nx, IS1_nx, IS2_nx, alpha_0prx, alpha_1prx, alpha_2prx, alpha_0plx, alpha_1plx, alpha_2plx, alpha_0nrx, alpha_1nrx, alpha_2nrx, alpha_0nlx, alpha_1nlx, alpha_2nlx, w0_prx, w1_prx, w2_prx, w0_plx, w1_plx, w2_plx, w0_nrx, w1_nrx, w2_nrx, w0_nlx, w1_nlx, w2_nlx, Epsilon);
            if (flag!=0) printf("error in Setalphawx function \n");
            
	    /* compute the positive and negative parts of right interface on x direction */
            flag = SetUx(w0_prx, w1_prx, w2_prx, w0_plx, w1_plx, w2_plx, w0_nrx, w1_nrx, w2_nrx, w0_nlx, w1_nlx, w2_nlx, fpnew, fnnew, u_tprhx, u_tnrhx, u_tplhx, u_tnlhx, i, j, Nx, Ny, 0);
            if (flag!=0) printf("error in SetUx function \n");
            
            for(k=0;k<4;k++){
	      //printf("w: i=%li, j=%li, k=%li, w0_px[k]=%g,w1_px[k]=%g,w2_px[k]=%g,w0_nx[k]=%g,w1_nx[k]=%g,w2_nx[k]=%g\n",i,j,k,w0_px[k],w1_px[k],w2_px[k],w0_nx[k],w1_nx[k],w2_nx[k]);
	      /* get right interface value for x component */
	      fup[idx(i, j, Nx, Ny, k)]=(u_tprhx[k]+u_tnrhx[k]);
	      /* get left interface value for x component */
	      fdown[idx(i, j, Nx, Ny, k)]=(u_tplhx[k]+u_tnlhx[k]);
	      //printf("f: i=%li, j=%li, k=%li, fup[idx(i, j, Nx, Ny, k)]=%g, fdown[idx(i, j, Nx, Ny, k)]=%g\n",i,j,k,fup[idx(i, j, Nx, Ny, k)], fdown[idx(i, j, Nx, Ny, k)]);
	      //printf("yx: i=%li, j=%li, k=%li, u_tpphx[k]=%f, u_tnphx[k]=%f\n",i,j, k,u_tpphx[k], u_tnphx[k]);
            }
        }
    }
    
    /* iterate over domain, computing all equations */
    for (i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){        
            /* get indicator of smoothness on y direction */
            flag = SetISY(IS0_py, IS1_py, IS2_py, IS0_ny, IS1_ny, IS2_ny, gpnew, gnnew, i, j, Nx, Ny, 0);
            if (flag!=0) printf("error in SetISY function \n");

            /* get weight on y direction for each variable */
            flag = Setalphawy(IS0_py, IS1_py, IS2_py, IS0_ny, IS1_ny, IS2_ny, alpha_0pry, alpha_1pry, alpha_2pry, alpha_0ply, alpha_1ply, alpha_2ply, alpha_0nry, alpha_1nry, alpha_2nry, alpha_0nly, alpha_1nly, alpha_2nly, w0_pry, w1_pry, w2_pry, w0_ply, w1_ply, w2_ply, w0_nry, w1_nry, w2_nry, w0_nly, w1_nly, w2_nly, Epsilon);
            if (flag!=0) printf("error in Setalphawy function \n");
            
	    /* compute the positive and negative parts of right interface on y direction */
            flag = SetUy(w0_pry, w1_pry, w2_pry, w0_ply, w1_ply, w2_ply, w0_nry, w1_nry, w2_nry, w0_nly, w1_nly, w2_nly, gpnew, gnnew, u_tprhy, u_tnrhy, u_tplhy, u_tnlhy, i, j, Nx, Ny, 0);
            if (flag!=0) printf("error in SetUy function \n");
            
            for(k=0;k<4;k++){
	        /* get right interface value for y component */
                gup[idx(i, j, Nx, Ny, k)]=u_tprhy[k]+u_tnrhy[k];
		/* get left interface value for y component */
                gdown[idx(i, j, Nx, Ny, k)]=u_tplhy[k]+u_tnlhy[k];
		//printf("yy: i=%li, j=%li, k=%li, u_tpphy[k]=%f, u_tnphy[k]=%f,u_tpnhy[k]=%f,u_tnnhy[k]=%f\n",i,j, k,u_tpphy[k], u_tnphy[k], u_tpnhy[k],u_tnnhy[k]);
            }
        }
    }
    
    /* transform the interface values back to the physical ones on x direction*/
    for(j=0; j<Ny; j++){
      for (i=0; i<Nx; i++){
	/* the left interface */
            flag = Setrhxegm(Ydata, rhxegm, i, j, Nx, Ny, gama, 0);
            if (flag!=0) printf("error in Setrhxegm function \n");
            for (k=0;k<4;k++){
	      //if (i==0&&j==0){
	      //printf(" rhxegm : i = %li, j = %li, k = %li, rhxegm[k][0] = %f, rhxegm[k][1] = %f, rhxegm[k][2] = %f, rhxegm[k][3] = %f\n", i, j, k, rhxegm[k][0], rhxegm[k][1],rhxegm[k][2],rhxegm[k][3]);
	      //}
                fdownback[idx(i, j, Nx, Ny, k)] = rhxegm[k][0]*fdown[idx(i, j, Nx, Ny, 0)]+rhxegm[k][1]*fdown[idx(i, j, Nx, Ny, 1)]+rhxegm[k][2]*fdown[idx(i, j, Nx, Ny, 2)]+rhxegm[k][3]*fdown[idx(i, j, Nx, Ny, 3)];
            }

	    /* the right interface */
	    flag = Setrhxegm(Ydata, rhxegm, i+1, j, Nx, Ny, gama, 0);
            if (flag!=0) printf("error in Setrhxegm function \n");
            for (k=0;k<4;k++){
	      //if (i==0&&j==0){
	      //printf(" rhxegm : i = %li, j = %li, k = %li, rhxegm[k][0] = %f, rhxegm[k][1] = %f, rhxegm[k][2] = %f, rhxegm[k][3] = %f\n", i, j, k, rhxegm[k][0], rhxegm[k][1],rhxegm[k][2],rhxegm[k][3]);
	      //}
                fupback[idx(i, j, Nx, Ny, k)] = rhxegm[k][0]*fup[idx(i, j, Nx, Ny, 0)]+rhxegm[k][1]*fup[idx(i, j, Nx, Ny, 1)]+rhxegm[k][2]*fup[idx(i, j, Nx, Ny, 2)]+rhxegm[k][3]*fup[idx(i, j, Nx, Ny, 3)];
            }
	    //printf(" fupback : i = %li, j = %li, fupback[idx(i, j, Nx, Ny, 0)]=%f, fupback[idx(i, j, Nx, Ny, 1)]=%f, fupback[idx(i, j, Nx, Ny, 2)]=%f, fupback[idx(i, j, Nx, Ny, 3)]=%f\n",i,j,fupback[idx(i, j, Nx, Ny, 0)],fupback[idx(i, j, Nx, Ny, 1)],fupback[idx(i, j, Nx, Ny, 2)],fupback[idx(i, j, Nx, Ny, 3)]);
	    //printf("fdownback : i = %li, j = %li, fdownback[idx(i, j, Nx, Ny, 0)]=%f, fdownback[idx(i, j, Nx, Ny, 1)]=%f, fdownback[idx(i, j, Nx, Ny, 2)]=%f, fdownback[idx(i, j, Nx, Ny, 3)]=%e\n",i,j,fdownback[idx(i, j, Nx, Ny, 0)],fdownback[idx(i, j, Nx, Ny, 1)],fdownback[idx(i, j, Nx, Ny, 2)],fdownback[idx(i, j, Nx, Ny, 3)]);
        }
    }
   
    /* transform the interface values back to the physical ones on y direction*/
    for(i=0; i<Nx; i++){
      for (j=0; j<Ny; j++){
	/* left interface */
            flag = Setrhyegm(Ydata, rhyegm, i, j, Nx, Ny, gama, 0);
            if (flag!=0) printf("error in Setrhyegm function \n");
            for (k=0;k<4;k++){
	      //  if (i==0&&j==0){
	      //printf(" rhyegm : i = %li, j = %li, rhyegm[k][0] = %f, rhyegm[k][1] = %f, rhyegm[k][2] = %f, rhyegm[k][3] = %f\n", i, j, rhyegm[k][0], rhyegm[k][1],rhyegm[k][2],rhyegm[k][3]);
	      //  }
                gdownback[idx(i, j, Nx, Ny, k)] = rhyegm[k][0]*gdown[idx(i, j, Nx, Ny, 0)]+rhyegm[k][1]*gdown[idx(i, j, Nx, Ny, 1)]+rhyegm[k][2]*gdown[idx(i, j, Nx, Ny, 2)]+rhyegm[k][3]*gdown[idx(i, j, Nx, Ny, 3)];
            }
	    
	    /* right interafce */
	    flag = Setrhyegm(Ydata, rhyegm, i, j+1, Nx, Ny, gama, 0);
            if (flag!=0) printf("error in Setrhyegm function \n");
            for (k=0;k<4;k++){
	      //  if (i==0&&j==0){
	      //printf(" rhyegm : i = %li, j = %li, rhyegm[k][0] = %f, rhyegm[k][1] = %f, rhyegm[k][2] = %f, rhyegm[k][3] = %f\n", i, j, rhyegm[k][0], rhyegm[k][1],rhyegm[k][2],rhyegm[k][3]);
	      //  }
                gupback[idx(i, j, Nx, Ny, k)] = rhyegm[k][0]*gup[idx(i, j, Nx, Ny, 0)]+rhyegm[k][1]*gup[idx(i, j, Nx, Ny, 1)]+rhyegm[k][2]*gup[idx(i, j, Nx, Ny, 2)]+rhyegm[k][3]*gup[idx(i, j, Nx, Ny, 3)];
            }
	    //if(i==0&&j==0){
	    //	printf(" yydownbackdata : i = %li, j = %li, yydownbackdata[idx(i, j, Nx, Ny, 0)] = %f, yydownbackdata[idx(i, j, Nx, Ny, 1)] = %f, yydownbackdata[idx(i, j, Nx, Ny, 2)] = %f, yydownbackdata[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, yydownbackdata[idx(i, j, Nx, Ny, 0)], yydownbackdata[idx(i, j, Nx, Ny, 1)],yydownbackdata[idx(i, j, Nx, Ny, 2)],yydownbackdata[idx(i, j, Nx, Ny, 3)]);
	    //	}
        }
    }
    
     /* fill in the value of tao and J in the whole domain */
    flag = Gettao(fupback, fdownback, gupback, gdownback, taoup, taodown, Nx, Ny, gama);
    if (flag!=0) printf("error in Gettao function \n");
    flag = GetCj(fupback, fdownback, gupback, gdownback, Cjup, Cjdown, Nx, Ny, gama);
    if (flag!=0) printf("error in GetCj function \n");
    
    for(j=0; j<Ny; j++){
      for (i=0; i<Nx; i++){
            /* get derivative both on x and y direction for the whole problem*/

            //for (k=0; k<4; k++){
                //  if(k==2)
                //        dYdata[idx(i, j, Nx, Ny, k)]=yxdata[idx(i, j, Nx, Ny, k)]+yydata[idx(i, j, Nx, Ny, k)]-0.1;
                // else
	      //if (i==0&&j==0){
	      //printf(" dYdata0 : i = %li, j = %li, fupback[idx(i, j, Nx, Ny, 1)] = %f, fdownback[idx(i, j, Nx, Ny, 1)] = %f, gupback[idx(i, j, Nx, Ny, 2)] = %f, gdownback[idx(i, j, Nx, Ny, 2)] = %f\n", i, j, fupback[idx(i, j, Nx, Ny, 1)], fdownback[idx(i, j, Nx, Ny, 1)], gupback[idx(i, j, Nx, Ny, 2)], gdownback[idx(i, j, Nx, Ny, 2)]);
	      //printf(" dYdata1 : i = %li, j = %li, taoup[idx(i, j, Nx, Ny, 0)] = %f, taodown[idx(i, j, Nx, Ny, 0)] = %f, taoup[idx(i, j, Nx, Ny, 2)] = %f, taodown[idx(i, j, Nx, Ny, 2)] = %f\n", i, j, taoup[idx(i, j, Nx, Ny, 0)], taodown[idx(i, j, Nx, Ny, 0)], taoup[idx(i, j, Nx, Ny, 2)], taodown[idx(i, j, Nx, Ny, 2)]);
	      //printf(" dYdata2 : i = %li, j = %li, taoup[idx(i, j, Nx, Ny, 1)] = %f, taodown[idx(i, j, Nx, Ny, 1)] = %f, taoup[idx(i, j, Nx, Ny, 2)] = %f, taodown[idx(i, j, Nx, Ny, 3)] = %f\n", i, j, taoup[idx(i, j, Nx, Ny, 1)], taodown[idx(i, j, Nx, Ny, 1)], taoup[idx(i, j, Nx, Ny, 2)], taodown[idx(i, j, Nx, Ny, 3)]);
	      //printf(" dYdata3 : i = %li, j = %li, Cjup[idx(i, j, Nx, Ny, 0)] = %f, Cjdown[idx(i, j, Nx, Ny, 0)] = %f, Cjup[idx(i, j, Nx, Ny, 1)] = %f, Cjdown[idx(i, j, Nx, Ny, 1)] = %f\n", i, j, Cjup[idx(i, j, Nx, Ny, 0)], Cjdown[idx(i, j, Nx, Ny, 0)], Cjup[idx(i, j, Nx, Ny, 1)], Cjdown[idx(i, j, Nx, Ny, 1)]);
	      //}
	dYdata[idx(i, j, Nx, Ny, 0)]=(-1.0/dx)*(fupback[idx(i, j, Nx, Ny, 1)]-fdownback[idx(i, j, Nx, Ny, 1)])+(-1.0/dy)*(gupback[idx(i, j, Nx, Ny, 2)]-gdownback[idx(i, j, Nx, Ny, 2)]);
	dYdata[idx(i, j, Nx, Ny, 1)]=(-1.0/dx)*(taoup[idx(i, j, Nx, Ny, 0)]-taodown[idx(i, j, Nx, Ny, 0)])+(-1.0/dy)*(taoup[idx(i, j, Nx, Ny, 2)]-taodown[idx(i, j, Nx, Ny, 2)]);
	dYdata[idx(i, j, Nx, Ny, 2)]=(-1.0/dx)*(taoup[idx(i, j, Nx, Ny, 1)]-taodown[idx(i, j, Nx, Ny, 1)])+(-1.0/dy)*(taoup[idx(i, j, Nx, Ny, 3)]-taodown[idx(i, j, Nx, Ny, 3)]);
        dYdata[idx(i, j, Nx, Ny, 3)]=(-1.0/dx)*(Cjup[idx(i, j, Nx, Ny, 0)]-Cjdown[idx(i, j, Nx, Ny, 0)])+(-1.0/dy)*(Cjup[idx(i, j, Nx, Ny, 1)]-Cjdown[idx(i, j, Nx, Ny, 1)]);
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
    
    delete []alpha_0prx;
    delete []alpha_1prx;
    delete []alpha_2prx;
    delete []alpha_0nrx;
    delete []alpha_1nrx;
    delete []alpha_2nrx;
    delete []alpha_0pry;
    delete []alpha_1pry;
    delete []alpha_2pry;
    delete []alpha_0nry;
    delete []alpha_1nry;
    delete []alpha_2nry;
    delete []alpha_0plx;
    delete []alpha_1plx;
    delete []alpha_2plx;
    delete []alpha_0nlx;
    delete []alpha_1nlx;
    delete []alpha_2nlx;
    delete []alpha_0ply;
    delete []alpha_1ply;
    delete []alpha_2ply;
    delete []alpha_0nly;
    delete []alpha_1nly;
    delete []alpha_2nly;

    delete []w0_prx;
    delete []w1_prx;
    delete []w2_prx;
    delete []w0_nrx;
    delete []w1_nrx;
    delete []w2_nrx;
    delete []w0_pry;
    delete []w1_pry;
    delete []w2_pry;
    delete []w0_nry;
    delete []w1_nry;
    delete []w2_nry;
    delete []w0_plx;
    delete []w1_plx;
    delete []w2_plx;
    delete []w0_nlx;
    delete []w1_nlx;
    delete []w2_nlx;
    delete []w0_ply;
    delete []w1_ply;
    delete []w2_ply;
    delete []w0_nly;
    delete []w1_nly;
    delete []w2_nly;

    delete []u_tprhx;
    delete []u_tplhx;
    delete []u_tnrhx;
    delete []u_tnlhx;
    delete []u_tprhy;
    delete []u_tplhy;
    delete []u_tnrhy;
    delete []u_tnlhy;
   
    delete []fp;
    delete []fn;
    delete []gp;
    delete []gn;
    delete []fpnew;
    delete []fnnew;
    delete []gpnew;
    delete []gnnew;
    delete []fup;
    delete []fdown;
    delete []gup;
    delete []gdown;
    delete []fupback;
    delete []fdownback;
    delete []gupback;
    delete []gdownback;
    delete []taoup;
    delete []taodown;
    delete []Cjup;
    delete []Cjdown;

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

/* split f into f_positive part and f_negative part x component*/ 
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
    realtype **temp = new realtype *[4];
    for (m=0;m<4;m++){
        temp[m] = new realtype [4];
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
    /* eigenvalues */
    egvxp[0][0]=(vx-a+sqrt((vx-a)*(vx-a)+eps))/2.0;
    egvxp[1][1]=(vx+sqrt(vx*vx+eps))/2.0;
    egvxp[2][2]=(vx+sqrt(vx*vx+eps))/2.0;
    egvxp[3][3]=(vx+a+sqrt((vx+a)*(vx+a)+eps))/2.0;
    egvxn[0][0]=(vx-a-sqrt((vx-a)*(vx-a)+eps))/2.0;
    egvxn[1][1]=(vx-sqrt(vx*vx+eps))/2.0;
    egvxn[2][2]=(vx-sqrt(vx*vx+eps))/2.0;
    egvxn[3][3]=(vx+a-sqrt((vx+a)*(vx+a)+eps))/2.0;

    /* left eigenvectors derivated from the book */
    /*
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
    */
    /* left eigenvectors derivated by myself */
    
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
    
    /* right eigenvectors derivated from the book */
    /*
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
    */
    /* right eigenvectors derivated by myself */
    
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
    
    //printf("lfxegm : i=%li, j=%li, lfxegm[0][0]=%f, lfxegm[1][0]=%f, lfxegm[2][0]=%f, lfxegm[3][0]=%f, lfxegm[0][1]=%f, lfxegm[1][1]=%f, lfxegm[2][1]=%f, lfxegm[3][1]=%f, lfxegm[0][2]=%f, lfxegm[1][2]=%f, lfxegm[2][2]=%f, lfxegm[3][2]=%f, lfxegm[0][3]=%f, lfxegm[1][3]=%f, lfxegm[2][3]=%f, lfxegm[3][3]=%f\n", i, j, lfxegm[0][0], lfxegm[1][0], lfxegm[2][0], lfxegm[3][0], lfxegm[0][1], lfxegm[1][1], lfxegm[2][1], lfxegm[3][1], lfxegm[0][2], lfxegm[1][2], lfxegm[2][2], lfxegm[3][2], lfxegm[0][3], lfxegm[1][3], lfxegm[2][3], lfxegm[3][3]);

    //printf("rhxegm : i=%li, j=%li, rhxegm[0][0]=%f, rhxegm[1][0]=%f, rhxegm[2][0]=%f, rhxegm[3][0]=%f, rhxegm[0][1]=%f, rhxegm[1][1]=%f, rhxegm[2][1]=%f, rhxegm[3][1]=%f, rhxegm[0][2]=%f, rhxegm[1][2]=%f, rhxegm[2][2]=%f, rhxegm[3][2]=%f, rhxegm[0][3]=%f, rhxegm[1][3]=%f, rhxegm[2][3]=%f, rhxegm[3][3]=%f\n", i, j, rhxegm[0][0], rhxegm[1][0], rhxegm[2][0], rhxegm[3][0], rhxegm[0][1], rhxegm[1][1], rhxegm[2][1], rhxegm[3][1], rhxegm[0][2], rhxegm[1][2], rhxegm[2][2], rhxegm[3][2], rhxegm[0][3], rhxegm[1][3], rhxegm[2][3], rhxegm[3][3]);

    for (m=0;m<4;m++)
      {
	for (n=0;n<4;n++)
	  {
	    temp[m][n]=lfxegm[m][0]*rhxegm[0][n]+lfxegm[m][1]*rhxegm[1][n]+lfxegm[m][2]*rhxegm[2][n]+lfxegm[m][3]*rhxegm[3][n];
	  }
      }

    //printf("temp : i=%li, j=%li, temp[0][0]=%f, temp[1][0]=%f, temp[2][0]=%f, temp[3][0]=%f, temp[0][1]=%f, temp[1][1]=%f, temp[2][1]=%f, temp[3][1]=%f, temp[0][2]=%f, temp[1][2]=%f, temp[2][2]=%f, temp[3][2]=%f, temp[0][3]=%f, temp[1][3]=%f, temp[2][3]=%f, temp[3][3]=%f\n", i, j, temp[0][0], temp[1][0], temp[2][0], temp[3][0], temp[0][1], temp[1][1], temp[2][1], temp[3][1], temp[0][2], temp[1][2], temp[2][2], temp[3][2], temp[0][3], temp[1][3], temp[2][3], temp[3][3]);

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

    for (m=0;m<4;m++){
        delete[] temp[m];
    }
    delete []temp;  
    
    return 0;
}

/* split g into g_positive part and g_negative part y component*/ 
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
static int Setlfxegm(realtype *Ydata, realtype **lfxegm, long int i, long int j, long int Nx, long int Ny, realtype gama)
{
    /* declaration */
    realtype rou, vx, vy, p, a, vxnext, vynext, pnext, vxcur, vycur, pcur, h, hcur, hnext;
  
    vx = Ydata[idx(i, j, Nx, Ny, 1)]/Ydata[idx(i, j, Nx, Ny, 0)];
    vy = Ydata[idx(i, j, Nx, Ny, 2)]/Ydata[idx(i, j, Nx, Ny, 0)];
    p = (gama-1.0)*Ydata[idx(i, j, Nx, Ny, 3)]-0.5*Ydata[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(vx*vx+vy*vy);
    a = sqrt(gama*p/Ydata[idx(i, j, Nx, Ny, 0)]);
    h = (Ydata[idx(i, j, Nx, Ny, 3)]+p)/Ydata[idx(i, j, Nx, Ny, 0)];
    
    /* left eigenvectors derivated by myself */
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
    realtype rou, vx, vy, p, a, vxnext, vynext, pnext, vxcur, vycur, pcur, h, hcur, hnext;
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
    rhxegm[0][2] = 1.0;
    rhxegm[1][2] = vx+a;
    rhxegm[2][2] = vy;
    rhxegm[3][2] = vx*vx+vy*vy-h+a*vx+(2*a*a)/(gama-1.0);
    rhxegm[0][3] = 0.0;
    rhxegm[1][3] = 0.0;
    rhxegm[2][3] = 1.0;
    rhxegm[3][3] = vy;
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
static int Setlfyegm(realtype *Ydata, realtype **lfyegm, long int i, long int j, long int Nx, long int Ny, realtype gama)
{
    /* declarations */
    realtype rou, vx, vy, p, a, vxnext, vynext, pnext, vxcur, vycur, pcur, h, hcur, hnext;
    
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
static int SetISX(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx, realtype *yxpdata, realtype *yxndata, long int i, long int j, long int Nx, long int Ny, int flag)
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
        
	/* without condisering BCs */
        if (i>2&&i<Nx-3){
            for (n=0;n<4;n++){
                for (k=0;k<7;k++){
                    yp[k+n*7] = yxpdata[idx(i-3+k,j,Nx,Ny,n)];
                    yn[k+n*7] = yxndata[idx(i-3+k,j,Nx,Ny,n)];
                }
            }
        }
        
	/* periodic BC */
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
        
        /* reflecting BC */
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
        
        /* nature or transmissive BC */
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
	  //printf("i=%li,j=%li,n=%li,yp[1+n*7]=%g,yp[2+n*7]=%g,yp[3+n*7]=%g,yp[4+n*7]=%g,yp[5+n*7]=%g,yn[2+n*7]=%g,yn[3+n*7]=%g,yn[4+n*7]=%g,yn[5+n*7]=%g,yn[6+n*7]=%g\n",i,j,n,yp[1+n*7],yp[2+n*7],yp[3+n*7],yp[4+n*7],yp[5+n*7],yn[2+n*7],yn[3+n*7],yn[4+n*7],yn[5+n*7],yn[6+n*7]);
            IS2_px[n]=(13.0/12.0)*(yp[1+n*7]-2.0*yp[2+n*7]+yp[3+n*7])*(yp[1+n*7]-2.0*yp[2+n*7]+yp[3+n*7])+(1.0/4.0)*(yp[1+n*7]-4.0*yp[2+n*7]+3.0*yp[3+n*7])*(yp[1+n*7]-4.0*yp[2+n*7]+3.0*yp[3+n*7]);
            
            IS1_px[n]=(13.0/12.0)*(yp[2+n*7]-2.0*yp[3+n*7]+yp[4+n*7])*(yp[2+n*7]-2.0*yp[3+n*7]+yp[4+n*7])+(1.0/4.0)*(yp[2+n*7]-yp[4+n*7])*(yp[2+n*7]-yp[4+n*7]);
            
            IS0_px[n]=(13.0/12.0)*(yp[3+n*7]-2.0*yp[4+n*7]+yp[5+n*7])*(yp[3+n*7]-2.0*yp[4+n*7]+yp[5+n*7])+(1.0/4.0)*(3.0*yp[3+n*7]-4.0*yp[4+n*7]+yp[5+n*7])*(3.0*yp[3+n*7]-4.0*yp[4+n*7]+yp[5+n*7]);

	    IS2_nx[n]=(13.0/12.0)*(yn[1+n*7]-2.0*yn[2+n*7]+yn[3+n*7])*(yn[1+n*7]-2.0*yn[2+n*7]+yn[3+n*7])+(1.0/4.0)*(yn[1+n*7]-4.0*yn[2+n*7]+3.0*yn[3+n*7])*(yn[1+n*7]-4.0*yn[2+n*7]+3.0*yn[3+n*7]);
            
            IS1_nx[n]=(13.0/12.0)*(yn[2+n*7]-2.0*yn[3+n*7]+yn[4+n*7])*(yn[2+n*7]-2.0*yn[3+n*7]+yn[4+n*7])+(1.0/4.0)*(yn[2+n*7]-yn[4+n*7])*(yn[2+n*7]-yn[4+n*7]);
            
            IS0_nx[n]=(13.0/12.0)*(yn[3+n*7]-2.0*yn[4+n*7]+yn[5+n*7])*(yn[3+n*7]-2.0*yn[4+n*7]+yn[5+n*7])+(1.0/4.0)*(3.0*yn[3+n*7]-4.0*yn[4+n*7]+yn[5+n*7])*(3.0*yn[3+n*7]-4.0*yn[4+n*7]+yn[5+n*7]);
            
	    //printf("i=%li,j=%li,n=%li,IS0_px[n]=%g,IS1_px[n]=%g,IS2_px[n]=%g,IS0_nx[n]=%g,IS1_nx[n]=%g,IS2_nx[n]=%g\n",i,j,n,IS0_px[n],IS1_px[n],IS2_px[n],IS0_nx[n],IS1_nx[n],IS2_nx[n]);
        }
       
        delete []yp;
        delete []yn;

        return 0;
    }
}

/* Fill in the indicators of smoothness on y direction */
static int SetISY(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny, realtype *yypdata, realtype *yyndata, long int i, long int j, long int Nx, long int Ny, int flag)
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
       
	/* without considering BCs */
        if (j>2&&j<Ny-3){
            for (n=0;n<4;n++){
                for (k=0;k<7;k++){
                    yp[k+n*7] = yypdata[idx(i,j-3+k,Nx,Ny,n)];
                    yn[k+n*7] = yyndata[idx(i,j-3+k,Nx,Ny,n)];
                }
            }
        }
        
	/* periodic BC */
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
        
        /* reflecting BC */
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
        
        /* nature or transmissive BC */
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
            IS2_py[n]=(13.0/12.0)*(yp[1+n*7]-2.0*yp[2+n*7]+yp[3+n*7])*(yp[1+n*7]-2.0*yp[2+n*7]+yp[3+n*7])+(1.0/4.0)*(yp[1+n*7]-4.0*yp[2+n*7]+3.0*yp[3+n*7])*(yp[1+n*7]-4.0*yp[2+n*7]+3.0*yp[3+n*7]);
            
            IS1_py[n]=(13.0/12.0)*(yp[2+n*7]-2.0*yp[3+n*7]+yp[4+n*7])*(yp[2+n*7]-2.0*yp[3+n*7]+yp[4+n*7])+(1.0/4.0)*(yp[2+n*7]-yp[4+n*7])*(yp[2+n*7]-yp[4+n*7]);
            
            IS0_py[n]=(13.0/12.0)*(yp[3+n*7]-2.0*yp[4+n*7]+yp[5+n*7])*(yp[3+n*7]-2.0*yp[4+n*7]+yp[5+n*7])+(1.0/4.0)*(3.0*yp[3+n*7]-4.0*yp[4+n*7]+yp[5+n*7])*(3.0*yp[3+n*7]-4.0*yp[4+n*7]+yp[5+n*7]);

	    IS2_ny[n]=(13.0/12.0)*(yn[1+n*7]-2.0*yn[2+n*7]+yn[3+n*7])*(yn[1+n*7]-2.0*yn[2+n*7]+yn[3+n*7])+(1.0/4.0)*(yn[1+n*7]-4.0*yn[2+n*7]+3.0*yn[3+n*7])*(yn[1+n*7]-4.0*yn[2+n*7]+3.0*yn[3+n*7]);
            
            IS1_ny[n]=(13.0/12.0)*(yn[2+n*7]-2.0*yn[3+n*7]+yn[4+n*7])*(yn[2+n*7]-2.0*yn[3+n*7]+yn[4+n*7])+(1.0/4.0)*(yn[2+n*7]-yn[4+n*7])*(yn[2+n*7]-yn[4+n*7]);
            
            IS0_ny[n]=(13.0/12.0)*(yn[3+n*7]-2.0*yn[4+n*7]+yn[5+n*7])*(yn[3+n*7]-2.0*yn[4+n*7]+yn[5+n*7])+(1.0/4.0)*(3.0*yn[3+n*7]-4.0*yn[4+n*7]+yn[5+n*7])*(3.0*yn[3+n*7]-4.0*yn[4+n*7]+yn[5+n*7]);
        }
        
        delete []yp;
        delete []yn;

        return 0;
    }
}

/* Fill in the stencil weights on x direction */
static int Setalphawx(realtype *IS0_px, realtype *IS1_px, realtype *IS2_px, realtype *IS0_nx, realtype *IS1_nx, realtype *IS2_nx, realtype *alpha_0prx, realtype *alpha_1prx, realtype *alpha_2prx, realtype *alpha_0plx, realtype *alpha_1plx, realtype *alpha_2plx, realtype *alpha_0nrx, realtype *alpha_1nrx, realtype *alpha_2nrx, realtype *alpha_0nlx, realtype *alpha_1nlx, realtype *alpha_2nlx, realtype *w0_prx, realtype *w1_prx, realtype *w2_prx, realtype *w0_plx, realtype *w1_plx, realtype *w2_plx, realtype *w0_nrx, realtype *w1_nrx, realtype *w2_nrx, realtype *w0_nlx, realtype *w1_nlx, realtype *w2_nlx, realtype Epsilon)
{
    /* compute the weights for rou, qx, qy, E */
    long int k;
    for(k=0;k<4;k++)
    {
        alpha_0prx[k]=(3.0/10.0)*(1.0/(Epsilon+IS0_px[k]))*(1.0/(Epsilon+IS0_px[k]));
        alpha_1prx[k]=(6.0/10.0)*(1.0/(Epsilon+IS1_px[k]))*(1.0/(Epsilon+IS1_px[k]));
        alpha_2prx[k]=(1.0/10.0)*(1.0/(Epsilon+IS2_px[k]))*(1.0/(Epsilon+IS2_px[k]));
	alpha_0plx[k]=(1.0/10.0)*(1.0/(Epsilon+IS0_px[k]))*(1.0/(Epsilon+IS0_px[k]));
        alpha_1plx[k]=(6.0/10.0)*(1.0/(Epsilon+IS1_px[k]))*(1.0/(Epsilon+IS1_px[k]));
        alpha_2plx[k]=(3.0/10.0)*(1.0/(Epsilon+IS2_px[k]))*(1.0/(Epsilon+IS2_px[k]));
        alpha_0nrx[k]=(3.0/10.0)*(1.0/(Epsilon+IS0_nx[k]))*(1.0/(Epsilon+IS0_nx[k]));
        alpha_1nrx[k]=(6.0/10.0)*(1.0/(Epsilon+IS1_nx[k]))*(1.0/(Epsilon+IS1_nx[k]));
        alpha_2nrx[k]=(1.0/10.0)*(1.0/(Epsilon+IS2_nx[k]))*(1.0/(Epsilon+IS2_nx[k]));
	alpha_0nlx[k]=(1.0/10.0)*(1.0/(Epsilon+IS0_nx[k]))*(1.0/(Epsilon+IS0_nx[k]));
        alpha_1nlx[k]=(6.0/10.0)*(1.0/(Epsilon+IS1_nx[k]))*(1.0/(Epsilon+IS1_nx[k]));
        alpha_2nlx[k]=(3.0/10.0)*(1.0/(Epsilon+IS2_nx[k]))*(1.0/(Epsilon+IS2_nx[k]));
        
        w0_prx[k]=alpha_0prx[k]/(alpha_0prx[k]+alpha_1prx[k]+alpha_2prx[k]);
        w1_prx[k]=alpha_1prx[k]/(alpha_0prx[k]+alpha_1prx[k]+alpha_2prx[k]);
        w2_prx[k]=alpha_2prx[k]/(alpha_0prx[k]+alpha_1prx[k]+alpha_2prx[k]);
	w0_plx[k]=alpha_0plx[k]/(alpha_0plx[k]+alpha_1plx[k]+alpha_2plx[k]);
        w1_plx[k]=alpha_1plx[k]/(alpha_0plx[k]+alpha_1plx[k]+alpha_2plx[k]);
        w2_plx[k]=alpha_2plx[k]/(alpha_0plx[k]+alpha_1plx[k]+alpha_2plx[k]);
        w0_nrx[k]=alpha_0nrx[k]/(alpha_0nrx[k]+alpha_1nrx[k]+alpha_2nrx[k]);
        w1_nrx[k]=alpha_1nrx[k]/(alpha_0nrx[k]+alpha_1nrx[k]+alpha_2nrx[k]);
        w2_nrx[k]=alpha_2nrx[k]/(alpha_0nrx[k]+alpha_1nrx[k]+alpha_2nrx[k]);
	w0_nlx[k]=alpha_0nlx[k]/(alpha_0nlx[k]+alpha_1nlx[k]+alpha_2nlx[k]);
        w1_nlx[k]=alpha_1nlx[k]/(alpha_0nlx[k]+alpha_1nlx[k]+alpha_2nlx[k]);
        w2_nlx[k]=alpha_2nlx[k]/(alpha_0nlx[k]+alpha_1nlx[k]+alpha_2nlx[k]);
    }
    return 0;
}

/* Fill in the stencil weights on y direction */
static int Setalphawy(realtype *IS0_py, realtype *IS1_py, realtype *IS2_py, realtype *IS0_ny, realtype *IS1_ny, realtype *IS2_ny, realtype *alpha_0pry, realtype *alpha_1pry, realtype *alpha_2pry, realtype *alpha_0ply, realtype *alpha_1ply, realtype *alpha_2ply, realtype *alpha_0nry, realtype *alpha_1nry, realtype *alpha_2nry, realtype *alpha_0nly, realtype *alpha_1nly, realtype *alpha_2nly, realtype *w0_pry, realtype *w1_pry, realtype *w2_pry, realtype *w0_ply, realtype *w1_ply, realtype *w2_ply, realtype *w0_nry, realtype *w1_nry, realtype *w2_nry, realtype *w0_nly, realtype *w1_nly, realtype *w2_nly, realtype Epsilon)
{
    /* compute the weights for rou, qx, qy, E */
    long int k;
    for(k=0;k<4;k++)
    {
        alpha_0pry[k]=(3.0/10.0)*(1.0/(Epsilon+IS0_py[k]))*(1.0/(Epsilon+IS0_py[k]));
        alpha_1pry[k]=(6.0/10.0)*(1.0/(Epsilon+IS1_py[k]))*(1.0/(Epsilon+IS1_py[k]));
        alpha_2pry[k]=(1.0/10.0)*(1.0/(Epsilon+IS2_py[k]))*(1.0/(Epsilon+IS2_py[k]));
	alpha_0ply[k]=(1.0/10.0)*(1.0/(Epsilon+IS0_py[k]))*(1.0/(Epsilon+IS0_py[k]));
        alpha_1ply[k]=(6.0/10.0)*(1.0/(Epsilon+IS1_py[k]))*(1.0/(Epsilon+IS1_py[k]));
        alpha_2ply[k]=(3.0/10.0)*(1.0/(Epsilon+IS2_py[k]))*(1.0/(Epsilon+IS2_py[k]));
        alpha_0nry[k]=(3.0/10.0)*(1.0/(Epsilon+IS0_ny[k]))*(1.0/(Epsilon+IS0_ny[k]));
        alpha_1nry[k]=(6.0/10.0)*(1.0/(Epsilon+IS1_ny[k]))*(1.0/(Epsilon+IS1_ny[k]));
        alpha_2nry[k]=(1.0/10.0)*(1.0/(Epsilon+IS2_ny[k]))*(1.0/(Epsilon+IS2_ny[k]));
	alpha_0nly[k]=(1.0/10.0)*(1.0/(Epsilon+IS0_ny[k]))*(1.0/(Epsilon+IS0_ny[k]));
        alpha_1nly[k]=(6.0/10.0)*(1.0/(Epsilon+IS1_ny[k]))*(1.0/(Epsilon+IS1_ny[k]));
        alpha_2nly[k]=(3.0/10.0)*(1.0/(Epsilon+IS2_ny[k]))*(1.0/(Epsilon+IS2_ny[k]));
        
        w0_pry[k]=alpha_0pry[k]/(alpha_0pry[k]+alpha_1pry[k]+alpha_2pry[k]);
        w1_pry[k]=alpha_1pry[k]/(alpha_0pry[k]+alpha_1pry[k]+alpha_2pry[k]);
        w2_pry[k]=alpha_2pry[k]/(alpha_0pry[k]+alpha_1pry[k]+alpha_2pry[k]);
	w0_ply[k]=alpha_0ply[k]/(alpha_0ply[k]+alpha_1ply[k]+alpha_2ply[k]);
        w1_ply[k]=alpha_1ply[k]/(alpha_0ply[k]+alpha_1ply[k]+alpha_2ply[k]);
        w2_ply[k]=alpha_2ply[k]/(alpha_0ply[k]+alpha_1ply[k]+alpha_2ply[k]);
        w0_nry[k]=alpha_0nry[k]/(alpha_0nry[k]+alpha_1nry[k]+alpha_2nry[k]);
        w1_nry[k]=alpha_1nry[k]/(alpha_0nry[k]+alpha_1nry[k]+alpha_2nry[k]);
        w2_nry[k]=alpha_2nry[k]/(alpha_0nry[k]+alpha_1nry[k]+alpha_2nry[k]);
	w0_nly[k]=alpha_0nly[k]/(alpha_0nly[k]+alpha_1nly[k]+alpha_2nly[k]);
        w1_nly[k]=alpha_1nly[k]/(alpha_0nly[k]+alpha_1nly[k]+alpha_2nly[k]);
        w2_nly[k]=alpha_2nly[k]/(alpha_0nly[k]+alpha_1nly[k]+alpha_2nly[k]);
    }
    return 0;
}

/* Get the positive and negative parts of right interface on x direction */
static int SetUx(realtype *w0_prx, realtype *w1_prx, realtype *w2_prx, realtype *w0_plx, realtype *w1_plx, realtype *w2_plx, realtype *w0_nrx, realtype *w1_nrx, realtype *w2_nrx, realtype *w0_nlx, realtype *w1_nlx, realtype *w2_nlx, realtype *yxpdata, realtype *yxndata, realtype *u_tprhx, realtype *u_tnrhx, realtype *u_tplhx, realtype *u_tnlhx, long int i, long int j, long int Nx, long int Ny, int flag)
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
        
        /* compute positive and negative solutions on the interface */
        for (n=0;n<4;n++){
            u_tprhx[n]=w2_prx[n]*((2.0/6.0)*yp[1+n*7]-(7.0/6.0)*yp[2+n*7]+(11.0/6.0)*yp[3+n*7])+w1_prx[n]*((-1.0/6.0)*yp[2+n*7]+(5.0/6.0)*yp[3+n*7]+(2.0/6.0)*yp[4+n*7])+w0_prx[n]*((2.0/6.0)*yp[3+n*7]+(5.0/6.0)*yp[4+n*7]-(1.0/6.0)*yp[5+n*7]);
            
            u_tplhx[n]=w2_plx[n]*((-1.0/6.0)*yp[1+n*7]+(5.0/6.0)*yp[2+n*7]+(2.0/6.0)*yp[3+n*7])+w1_plx[n]*((2.0/6.0)*yp[2+n*7]+(5.0/6.0)*yp[3+n*7]-(1.0/6.0)*yp[4+n*7])+w0_plx[n]*((11.0/6.0)*yp[3+n*7]-(7.0/6.0)*yp[4+n*7]+(2.0/6.0)*yp[5+n*7]);
            
            u_tnrhx[n]=w2_nrx[n]*((2.0/6.0)*yn[1+n*7]-(7.0/6.0)*yn[2+n*7]+(11.0/6.0)*yn[3+n*7])+w1_nrx[n]*((-1.0/6.0)*yn[2+n*7]+(5.0/6.0)*yn[3+n*7]+(2.0/6.0)*yn[4+n*7])+w0_nrx[n]*((2.0/6.0)*yn[3+n*7]+(5.0/6.0)*yn[4+n*7]-(1.0/6.0)*yn[5+n*7]);
            
            u_tnlhx[n]=w2_nlx[n]*((-1.0/6.0)*yn[1+n*7]+(5.0/6.0)*yn[2+n*7]+(2.0/6.0)*yn[3+n*7])+w1_nlx[n]*((2.0/6.0)*yn[2+n*7]+(5.0/6.0)*yn[3+n*7]-(1.0/6.0)*yn[4+n*7])+w0_nlx[n]*((11.0/6.0)*yn[3+n*7]-(7.0/6.0)*yn[4+n*7]+(2.0/6.0)*yn[5+n*7]);
	    
        }
        
        delete []yp;
        delete []yn;
        
        return 0;
    }
}

/* Get the positive and negative parts of right interface  on y direction */
static int SetUy(realtype *w0_pry, realtype *w1_pry, realtype *w2_pry, realtype *w0_ply, realtype *w1_ply, realtype *w2_ply, realtype *w0_nry, realtype *w1_nry, realtype *w2_nry, realtype *w0_nly, realtype *w1_nly, realtype *w2_nly, realtype *yypdata, realtype *yyndata, realtype *u_tprhy, realtype *u_tnrhy, realtype *u_tplhy, realtype *u_tnlhy, long int i, long int j, long int Nx, long int Ny, int flag)
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
        
        /* compute positive and negative solutions on the interface */
        for (n=0;n<4;n++){
            u_tprhy[n]=w2_pry[n]*((2.0/6.0)*yp[1+n*7]-(7.0/6.0)*yp[2+n*7]+(11.0/6.0)*yp[3+n*7])+w1_pry[n]*((-1.0/6.0)*yp[2+n*7]+(5.0/6.0)*yp[3+n*7]+(2.0/6.0)*yp[4+n*7])+w0_pry[n]*((2.0/6.0)*yp[3+n*7]+(5.0/6.0)*yp[4+n*7]-(1.0/6.0)*yp[5+n*7]);
            
            u_tplhy[n]=w2_ply[n]*((-1.0/6.0)*yp[1+n*7]+(5.0/6.0)*yp[2+n*7]+(2.0/6.0)*yp[3+n*7])+w1_ply[n]*((2.0/6.0)*yp[2+n*7]+(5.0/6.0)*yp[3+n*7]-(1.0/6.0)*yp[4+n*7])+w0_ply[n]*((11.0/6.0)*yp[3+n*7]-(7.0/6.0)*yp[4+n*7]+(2.0/6.0)*yp[5+n*7]);
            
            u_tnrhy[n]=w2_nry[n]*((2.0/6.0)*yn[1+n*7]-(7.0/6.0)*yn[2+n*7]+(11.0/6.0)*yn[3+n*7])+w1_nry[n]*((-1.0/6.0)*yn[2+n*7]+(5.0/6.0)*yn[3+n*7]+(2.0/6.0)*yn[4+n*7])+w0_nry[n]*((2.0/6.0)*yn[3+n*7]+(5.0/6.0)*yn[4+n*7]-(1.0/6.0)*yn[5+n*7]);
            
            u_tnlhy[n]=w2_nly[n]*((-1.0/6.0)*yn[1+n*7]+(5.0/6.0)*yn[2+n*7]+(2.0/6.0)*yn[3+n*7])+w1_nly[n]*((2.0/6.0)*yn[2+n*7]+(5.0/6.0)*yn[3+n*7]-(1.0/6.0)*yn[4+n*7])+w0_nly[n]*((11.0/6.0)*yn[3+n*7]-(7.0/6.0)*yn[4+n*7]+(2.0/6.0)*yn[5+n*7]);
        }
        
        delete []yp;
        delete []yn;
        
        return 0;
    }
}

/* Fill in values of tao in the whole domain */
static int Gettao(realtype *fupback, realtype *fdownback, realtype *gupback, realtype *gdownback, realtype *taoup, realtype *taodown, long int Nx, long int Ny, realtype gama)
{
    /* declare parameters */
    long int i, j;
   
    realtype *vxxup = new realtype [Nx*Ny];
    realtype *vyxup = new realtype [Nx*Ny];
    realtype *vxyup = new realtype [Nx*Ny];
    realtype *vyyup = new realtype [Nx*Ny];
    realtype *pxup = new realtype [Nx*Ny];
    realtype *pyup = new realtype [Nx*Ny];
    realtype *vxxdown = new realtype [Nx*Ny];
    realtype *vyxdown = new realtype [Nx*Ny];
    realtype *vxydown = new realtype [Nx*Ny];
    realtype *vyydown = new realtype [Nx*Ny];
    realtype *pxdown = new realtype [Nx*Ny];
    realtype *pydown = new realtype [Nx*Ny];

    /* Set values into vx and vy on x direction*/
    for(j=0;j<Ny;j++){
      for(i=0;i<Nx;i++){
	vxxup[idx_v(i,j,Nx)]=fupback[idx(i, j, Nx, Ny, 1)]/fupback[idx(i, j, Nx, Ny, 0)];
	vxxdown[idx_v(i,j,Nx)]=fdownback[idx(i, j, Nx, Ny, 1)]/fdownback[idx(i, j, Nx, Ny, 0)];
	//printf("i = %li, j=%li, yxbackdata[idx(i, j, (Nx+1), Ny, 1)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",i,j,yxbackdata[idx(i, j, (Nx+1), Ny, 1)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
	vyxup[idx_v(i,j,Nx)]=fupback[idx(i, j, Nx, Ny, 2)]/fupback[idx(i, j, Nx, Ny, 0)];
	vyxdown[idx_v(i,j,Nx)]=fdownback[idx(i, j, Nx, Ny, 2)]/fdownback[idx(i, j, Nx, Ny, 0)];
	//printf("yxbackdata[idx(i, j, (Nx+1), Ny, 2)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",yxbackdata[idx(i, j, (Nx+1), Ny, 2)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
	//printf("i = %li, j=%li, vxx[idx_v(i,j,(Nx+1))]=%f,vyx[idx_v(i,j,(Nx+1))]=%f\n",i,j,vxx[idx_v(i,j,(Nx+1))],vyx[idx_v(i,j,(Nx+1))]);
        }
    }
    
    /* Set values into vxy and vyy on y direction*/
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	vxyup[idx_v(i,j,Nx)]=gupback[idx(i, j, Nx, Ny, 1)]/gupback[idx(i, j, Nx, Ny, 0)];
	vxydown[idx_v(i,j,Nx)]=gdownback[idx(i, j, Nx, Ny, 1)]/gdownback[idx(i, j, Nx, Ny, 0)];
	//printf("i = %li, j=%li, yybackdata[idx(i, j, Nx, (Ny+1), 1)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",i,j,yybackdata[idx(i, j, Nx, (Ny+1), 1)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
        vyyup[idx_v(i,j,Nx)]=gupback[idx(i, j, Nx, Ny, 2)]/gupback[idx(i, j, Nx, Ny, 0)];
	vyydown[idx_v(i,j,Nx)]=gdownback[idx(i, j, Nx, Ny, 2)]/gdownback[idx(i, j, Nx, Ny, 0)];
	//printf("yybackdata[idx(i, j, Nx, (Ny+1), 2)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",yybackdata[idx(i, j, Nx, (Ny+1), 2)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
        }
    }

    /* Compute the left and right interface values of tao in the whole domain x direction */
    for(j=0;j<Ny;j++){
      for(i=0;i<Nx;i++){
	pxup[idx_v(i,j,Nx)] =  fupback[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(fupback[idx(i, j, Nx, Ny, 3)]/fupback[idx(i, j, Nx, Ny, 0)]-0.5*(vxxup[idx_v(i,j,Nx)]*vxxup[idx_v(i,j,Nx)]+vyxup[idx_v(i,j,Nx)]*vyxup[idx_v(i,j,Nx)]));	   
	taoup[idx(i, j, Nx, Ny, 0)] = fupback[idx(i, j, Nx, Ny, 1)]*vxxup[idx_v(i,j,Nx)]+pxup[idx_v(i,j,Nx)];
	taoup[idx(i, j, Nx, Ny, 1)] = fupback[idx(i, j, Nx, Ny, 1)]*vyxup[idx_v(i,j,Nx)];

	pxdown[idx_v(i,j,Nx)] =  fdownback[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(fdownback[idx(i, j, Nx, Ny, 3)]/fdownback[idx(i, j, Nx, Ny, 0)]-0.5*(vxxdown[idx_v(i,j,Nx)]*vxxdown[idx_v(i,j,Nx)]+vyxdown[idx_v(i,j,Nx)]*vyxdown[idx_v(i,j,Nx)]));	  
	taodown[idx(i, j, Nx, Ny, 0)] = fdownback[idx(i, j, Nx, Ny, 1)]*vxxdown[idx_v(i,j,Nx)]+pxdown[idx_v(i,j,Nx)];
	taodown[idx(i, j, Nx, Ny, 1)] = fdownback[idx(i, j, Nx, Ny, 1)]*vyxdown[idx_v(i,j,Nx)];
	     } 
	}

    /* Compute the left and right interaface values of tao in the whole domain on y direction*/
    for (i=0;i<Nx;i++){
      for (j=0;j<Ny;j++){
	pyup[idx_v(i,j,Nx)] =  gupback[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(gupback[idx(i, j, Nx, Ny, 3)]/gupback[idx(i, j, Nx, Ny, 0)]-0.5*(vyxup[idx_v(i,j,Nx)]*vyxup[idx_v(i,j,Nx)]+vyyup[idx_v(i,j,Nx)]*vyyup[idx_v(i,j,Nx)]));	   
	taoup[idx(i, j, Nx, Ny, 2)] = gupback[idx(i, j, Nx, Ny, 2)]*vxyup[idx_v(i,j,Nx)];
	taoup[idx(i, j, Nx, Ny, 3)] = gupback[idx(i, j, Nx, Ny, 2)]*vyyup[idx_v(i,j,Nx)]+pyup[idx_v(i,j,Nx)];

	pydown[idx_v(i,j,Nx)] =  gdownback[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(gdownback[idx(i, j, Nx, Ny, 3)]/gdownback[idx(i, j, Nx, Ny, 0)]-0.5*(vyxdown[idx_v(i,j,Nx)]*vyxdown[idx_v(i,j,Nx)]+vyydown[idx_v(i,j,Nx)]*vyydown[idx_v(i,j,Nx)]));	  
	taodown[idx(i, j, Nx, Ny, 2)] = gdownback[idx(i, j, Nx, Ny, 2)]*vxydown[idx_v(i,j,Nx)];
	taodown[idx(i, j, Nx, Ny, 3)] = gdownback[idx(i, j, Nx, Ny, 2)]*vyydown[idx_v(i,j,Nx)]+pydown[idx_v(i,j,Nx)];
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
static int GetCj(realtype *fupback, realtype *fdownback, realtype *gupback, realtype *gdownback, realtype *Cjup, realtype *Cjdown, long int Nx, long int Ny, realtype gama)
{
    /* declare parameters */
    long int i, j;
    
    realtype *vxxup = new realtype [Nx*Ny];
    realtype *vyxup = new realtype [Nx*Ny];
    realtype *vxyup = new realtype [Nx*Ny];
    realtype *vyyup = new realtype [Nx*Ny];
    realtype *pxup = new realtype [Nx*Ny];
    realtype *pyup = new realtype [Nx*Ny];
    realtype *vxxdown = new realtype [Nx*Ny];
    realtype *vyxdown = new realtype [Nx*Ny];
    realtype *vxydown = new realtype [Nx*Ny];
    realtype *vyydown = new realtype [Nx*Ny];
    realtype *pxdown = new realtype [Nx*Ny];
    realtype *pydown = new realtype [Nx*Ny];

     /* Set values into vx and vy on x direction */
    for(j=0;j<Ny;j++){
      for(i=0;i<Nx;i++){
	vxxup[idx_v(i,j,Nx)]=fupback[idx(i, j, Nx, Ny, 1)]/fupback[idx(i, j, Nx, Ny, 0)];
	vxxdown[idx_v(i,j,Nx)]=fdownback[idx(i, j, Nx, Ny, 1)]/fdownback[idx(i, j, Nx, Ny, 0)];
	//printf("i = %li, j=%li, yxbackdata[idx(i, j, (Nx+1), Ny, 1)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",i,j,yxbackdata[idx(i, j, (Nx+1), Ny, 1)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
	vyxup[idx_v(i,j,Nx)]=fupback[idx(i, j, Nx, Ny, 2)]/fupback[idx(i, j, Nx, Ny, 0)];
	vyxdown[idx_v(i,j,Nx)]=fdownback[idx(i, j, Nx, Ny, 2)]/fdownback[idx(i, j, Nx, Ny, 0)];
	//printf("yxbackdata[idx(i, j, (Nx+1), Ny, 2)]=%f,yxbackdata[idx(i, j, (Nx+1), Ny, 0)]=%f\n",yxbackdata[idx(i, j, (Nx+1), Ny, 2)],yxbackdata[idx(i, j, (Nx+1), Ny, 0)]);
	//printf("i = %li, j=%li, vxx[idx_v(i,j,(Nx+1))]=%f,vyx[idx_v(i,j,(Nx+1))]=%f\n",i,j,vxx[idx_v(i,j,(Nx+1))],vyx[idx_v(i,j,(Nx+1))]);
        }
    }
    
    /* Set values into vx and vy on y direction*/
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	vxyup[idx_v(i,j,Nx)]=gupback[idx(i, j, Nx, Ny, 1)]/gupback[idx(i, j, Nx, Ny, 0)];
	vxydown[idx_v(i,j,Nx)]=gdownback[idx(i, j, Nx, Ny, 1)]/gdownback[idx(i, j, Nx, Ny, 0)];
	//printf("i = %li, j=%li, yybackdata[idx(i, j, Nx, (Ny+1), 1)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",i,j,yybackdata[idx(i, j, Nx, (Ny+1), 1)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
        vyyup[idx_v(i,j,Nx)]=gupback[idx(i, j, Nx, Ny, 2)]/gupback[idx(i, j, Nx, Ny, 0)];
	vyydown[idx_v(i,j,Nx)]=gdownback[idx(i, j, Nx, Ny, 2)]/gdownback[idx(i, j, Nx, Ny, 0)];
	//printf("yybackdata[idx(i, j, Nx, (Ny+1), 2)]=%f,yybackdata[idx(i, j, Nx, (Ny+1), 0)]=%f\n",yybackdata[idx(i, j, Nx, (Ny+1), 2)],yybackdata[idx(i, j, Nx, (Ny+1), 0)]);
        }
    }
    
     /* Compute the values of tao in the whole domain on x direction*/
    for(j=0;j<Ny;j++){
      for(i=0;i<Nx;i++){
	pxup[idx_v(i,j,Nx)] =  fupback[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(fupback[idx(i, j, Nx, Ny, 3)]/fupback[idx(i, j, Nx, Ny, 0)]-0.5*(vxxup[idx_v(i,j,Nx)]*vxxup[idx_v(i,j,Nx)]+vyxup[idx_v(i,j,Nx)]*vyxup[idx_v(i,j,Nx)]));		   
	Cjup[idx(i, j, Nx, Ny, 0)] = (fupback[idx(i, j, Nx, Ny, 3)]+pxup[idx_v(i,j,Nx)])*vxxup[idx_v(i,j,Nx)];

	pxdown[idx_v(i,j,Nx)] =  fdownback[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(fdownback[idx(i, j, Nx, Ny, 3)]/fdownback[idx(i, j, Nx, Ny, 0)]-0.5*(vxxdown[idx_v(i,j,Nx)]*vxxdown[idx_v(i,j,Nx)]+vyxdown[idx_v(i,j,Nx)]*vyxdown[idx_v(i,j,Nx)]));       
	Cjdown[idx(i, j, Nx, Ny, 0)] = (fdownback[idx(i, j, Nx, Ny, 3)]+pxdown[idx_v(i,j,Nx)])*vxxdown[idx_v(i,j,Nx)];
	     } 
	}

    /* Compute the values of tao in the whole domain on y direction*/
    for (i=0;i<Nx;i++){
      for (j=0;j<Ny;j++){
	pyup[idx_v(i,j,Nx)] =  gupback[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(gupback[idx(i, j, Nx, Ny, 3)]/gupback[idx(i, j, Nx, Ny, 0)]-0.5*(vyxup[idx_v(i,j,Nx)]*vyxup[idx_v(i,j,Nx)]+vyyup[idx_v(i,j,Nx)]*vyyup[idx_v(i,j,Nx)]));	   
	Cjup[idx(i, j, Nx, Ny, 1)] = (gupback[idx(i, j, Nx, Ny, 3)]+pyup[idx_v(i,j,Nx)])*vyyup[idx_v(i,j,Nx)];

	pydown[idx_v(i,j,Nx)] =  gdownback[idx(i, j, Nx, Ny, 0)]*(gama-1.0)*(gdownback[idx(i, j, Nx, Ny, 3)]/gdownback[idx(i, j, Nx, Ny, 0)]-0.5*(vyxdown[idx_v(i,j,Nx)]*vyxdown[idx_v(i,j,Nx)]+vyydown[idx_v(i,j,Nx)]*vyydown[idx_v(i,j,Nx)]));	   
	Cjdown[idx(i, j, Nx, Ny, 1)] = (gdownback[idx(i, j, Nx, Ny, 3)]+pydown[idx_v(i,j,Nx)])*vyydown[idx_v(i,j,Nx)];
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
