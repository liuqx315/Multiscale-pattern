
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arkode/arkode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>


#define idx(i,j,Nx,Ny,k) ((k)*(Nx)*(Ny)+(j)*(Nx)+i)
#define idx_v(i,j,Nx) ((j)*(Nx)+i)

#define PI RCONST(3.1415926535897932)
#define ONE RCONST(1.0)
#define TWO RCONST(2.0)


/* user data structure */
typedef struct {
    long int Nx;    /* number of intervals     */
    long int Ny;
    realtype dx;   /* mesh spacing            */
    realtype dy;   /* mesh spacing            */
    realtype Lx;
    realtype Ly
    int k;
} *UserData;

/* User-supplied Functions Called by the Solver */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
/* Private function to check function return values */
static int check_flag(void *flagvalue, const string funcname, int opt);
/* get value of tao */
static int Gettao(N_Vector y, N_Vector tao, long int Nx, long int Ny);
/* get value of J */
static int GetCj(N_Vector y, N_Vector Cj, long int Nx, long int Ny);

// what is the input u exactly? I do it as below, is it ok?
int main(int argc, const char * argv[])
{
    /* general problem parameters */
    realtype T0 = RCONST(0.0);
    realtype Tf = RCONST(10.0);
    int Nt = 10;
    int Nvar = 4;
    UserData udata = NULL;
    realtype *data;
    long int N, NEQ, i, j;
    
    /* declare solver parameters */
    int flag;
    
    N_Vector y = NULL;
 //   N_Vector rou = NULL;
 //   N_Vector qx = NULL;
 //   N_Vector qy = NULL;
 //   N_Vector E = NULL;
    
    void *arkode_mem = NULL;
    
    /* allocate udata structure */
    udata = (UserData) malloc(sizeof(*udata));
    if (check_flag((void *) udata, "malloc", 2)) return 1;
    
    /* read problem parameter and tolerances from input file:*/
    double dx, dy;
    FILE *FID;
    FID=fopen("input_WENO2D.txt","r");
    flag = fscanf(FID," Nx = %li\n", &Nx);
    flag = fscanf(FID," Ny = %li\n", &Ny);
    flag = fscanf(FID," Lx = %lf\n", &Lx);
    flag = fscanf(FID," Ly = %lf\n", &Ly);
    fclose(FID);
    
    /* store the inputs in the UserData structure */
    udata->Nx = Nx;
    udata->Ny = Ny;
    udata->Lx = Lx;
    udata->Ly = Ly;
    
    /* open solver diagnostics output file for writing */
    FILE *DFID;
    DFID=fopen("diags_ark_WENO2D.txt","w");
    
    /* set total allocated vector length */
    NEQ = Nvar*udata->N;
    
    /* Initial problem output */
    printf("\n2D gas dynamic test problem:\n");
    printf("    Nx = %li,  Ny = %li, NEQ = %li\n", udata->Nx, udata->Ny, NEQ);
    printf("    problem parameters:  Lx = %g,  Ly = %g\n", udata->Lx, udata->Ly);
    
    /* Create serial vector of length NEQ for initial condition */
    y = N_VNew_Serial(NEQ);
    if (check_flag((void *) y, "N_VNew_Serial", 0)) return 1;
    
    /* set spatial mesh spacing */
    udata->dx = Lx/Nx;
    udata->dy = Ly/Ny;
    
    /* output mesh to disk */
    FID=fopen("WENO2D_mesh.txt","w");
    for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->dx*(i+0.5));
    fclose(FID);
    
    /* Access data array for new NVector y */
    data = N_VGetArrayPointer(y);
    if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
    
    for(j=0;j<Ny;j++){
        for (i=0; i<Nx; i++) {
            data[idx(i,j,Nx,Ny,0)] =  sin(PI*(i+0.5)*udata->dx)*sin(PI*(j+0.5)*udata->dy);  /* rou */
            data[idx(i,j,Nx,Ny,1)] =  sin(PI*(i+0.5)*udata->dx)*sin(PI*(j+0.5)*udata->dy);  /* qx */
            data[idx(i,j,Nx,Ny,2)] =  sin(PI*(i+0.5)*udata->dx)*sin(PI*(j+0.5)*udata->dy);  /* qy */
            data[idx(i,j,Nx,Ny,3)] =  sin(PI*(i+0.5)*udata->dx)*sin(PI*(j+0.5)*udata->dy);  /* E */
        }   
    }
    
    /* Call ARKodeCreate to create the solver memory */
    arkode_mem = ARKodeCreate();
    if (check_flag((void *) arkode_mem, "ARKodeCreate", 0)) return 1;
    
    /* Call init_from_file helper routine to read and set solver parameters */
    realtype rtol, atol;
    rtol = 1.e-3;
    atol = 1.e-5;
    realtype reltol  = rtol;
    realtype abstol  = atol;
    
    /* Reference solution uses default implicit method */
    flag = ARKodeInit(arkode_mem, f, NULL, T0, y);
    if (check_flag(&flag, "ARKodeInit", 1)) return 1;
    
    /* Call ARKodeSetUserData to pass rdata to user functions */
    flag = ARKodeSetUserData(arkode_mem, (void *) udata);
    if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;
    
    /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
    flag = ARKodeSetDiagnostics(arkode_mem, DFID);
    if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return 1;
    
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
    printf("\n");
    flag = ARKodeWriteParameters(arkode_mem, stdout);
    if (check_flag(&flag, "ARKodeWriteParameters", 1)) return 1;
    
    /* In loop, call ARKode, print results, and test for error.
     Break out of loop when the final output time has been reached */
    realtype t  = T0;
    realtype dTout = Tf/Nt;
    realtype tout = T0+dTout;
    
    int iout;
    for (iout=0; iout<Nt; iout++) {
        
        flag = ARKodeSetStopTime(arkode_mem, tout);
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
    long int NEQ, NEQS, i, j, k;
    int flag;
    
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
    
    realtype Epsilon;
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
    
    /* shortcuts to number of intervals, background values */
    long int Nx = udata->Nx;
    long int Ny = udata->Ny;
    realtype dx = udata->dx;
    realtype dy = udata->dy;
    realtype Lx = udata->Lx;
    realtype Ly = udata->Ly;
    
    NEQ = 4*Nx*Ny;
    NEQS = 2*Nx*Ny;
    
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
    
    flag=Gettao(y, tao, Nx, Ny);
    if (flag!=0) printf("error in Gettao function \n");
    flag=GetCj(y, Cj, Nx, Ny);
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
    
    /* iterate over domain, computing all equations */
            i=0;
    for(j=0; j<Ny; j++){
            // x direction
            IS0_px[0]=(13/12)*(ypdata[idx(Nx-2, j, Nx, Ny, 1)]-2*ypdata[idx(Nx-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(Nx-2, j, Nx, Ny, 1)]-2*ypdata[idx(Nx-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(Nx-2, j, Nx, Ny, 1)]-4*ypdata[idx(Nx-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(Nx-2, j, Nx, Ny, 1)]-4*ypdata[idx(Nx-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_px[0]=(13/12)*(ypdata[idx(Nx-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(i+1, j, Nx, Ny, 1)])*(ypdata[idx(Nx-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(i+1, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(Nx-1, j, Nx, Ny, 1)]-1*ypdata[idx(i+1, j, Nx, Ny, 1)])*(ypdata[idx(Nx-1, j, Nx, Ny, 1)]-1*ypdata[idx(i+1, j, Nx, Ny, 1)]);
            
            IS2_px[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)]);
            
            IS0_nx[0]=(13/12)*(yndata[idx(i+1, j, Nx, Ny, 1)]-2*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)])*(yndata[idx(i+1, j, Nx, Ny, 1)]-2*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)])+(1/4)*(3*yndata[idx(i+1, j, Nx, Ny, 1)]-4*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)])*(3*yndata[idx(i+1, j, Nx, Ny, 1)]-4*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)]);
            
            IS1_nx[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(i+1, j, Nx, Ny, 1)]+yndata[idx(i+2, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(i+1, j, Nx, Ny, 1)]+yndata[idx(i+2, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(i+2, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(i+2, j, Nx, Ny, 1)]);
            
            IS2_nx[0]=(13/12)*(yndata[idx(Nx-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(i+1, j, Nx, Ny, 1)])*(yndata[idx(Nx-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(i+1, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(Nx-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(i+1, j, Nx, Ny, 1)])*(yndata[idx(Nx-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(i+1, j, Nx, Ny, 1)]);
            
            
            IS0_px[1]=(13/12)*(taopdata[idx(Nx-2, j, Nx, Ny, 0)]-2*taopdata[idx(Nx-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(Nx-2, j, Nx, Ny, 0)]-2*taopdata[idx(Nx-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(Nx-2, j, Nx, Ny, 0)]-4*taopdata[idx(Nx-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(Nx-2, j, Nx, Ny, 0)]-4*taopdata[idx(Nx-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)]);
            
            IS1_px[1]=(13/12)*(taopdata[idx(Nx-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(i+1, j, Nx, Ny, 0)])*(taopdata[idx(Nx-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(Nx-1, j, Nx, Ny, 0)]-1*taopdata[idx(i+1, j, Nx, Ny, 0)])*(taopdata[idx(Nx-1, j, Nx, Ny, 0)]-1*taopdata[idx(i+1, j, Nx, Ny, 0)]);
            
            IS2_px[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS0_nx[1]=(13/12)*(taondata[idx(i+1, j, Nx, Ny, 0)]-2*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)])*(taondata[idx(i+1, j, Nx, Ny, 0)]-2*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)])+(1/4)*(3*taondata[idx(i+1, j, Nx, Ny, 0)]-4*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)])*(3*taondata[idx(i+1, j, Nx, Ny, 0)]-4*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)]);
            
            IS1_nx[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(i+1, j, Nx, Ny, 0)]+taondata[idx(i+2, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(i+1, j, Nx, Ny, 0)]+taondata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(i+2, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS2_nx[1]=(13/12)*(taondata[idx(Nx-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(i+1, j, Nx, Ny, 0)])*(taondata[idx(Nx-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(Nx-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(i+1, j, Nx, Ny, 0)])*(taondata[idx(Nx-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(i+1, j, Nx, Ny, 0)]);
            
            
            IS0_px[2]=(13/12)*(taopdata[idx(Nx-2, j, Nx, Ny, 2)]-2*taopdata[idx(Nx-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(Nx-2, j, Nx, Ny, 2)]-2*taopdata[idx(Nx-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(Nx-2, j, Nx, Ny, 2)]-4*taopdata[idx(Nx-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(Nx-2, j, Nx, Ny, 2)]-4*taopdata[idx(Nx-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_px[2]=(13/12)*(taopdata[idx(Nx-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(i+1, j, Nx, Ny, 2)])*(taopdata[idx(Nx-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(i+1, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(Nx-1, j, Nx, Ny, 2)]-1*taopdata[idx(i+1, j, Nx, Ny, 2)])*(taopdata[idx(Nx-1, j, Nx, Ny, 2)]-1*taopdata[idx(i+1, j, Nx, Ny, 2)]);
            
            IS2_px[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)]);
            
            IS0_nx[2]=(13/12)*(taondata[idx(i+1, j, Nx, Ny, 2)]-2*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)])*(taondata[idx(i+1, j, Nx, Ny, 2)]-2*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)])+(1/4)*(3*taondata[idx(i+1, j, Nx, Ny, 2)]-4*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)])*(3*taondata[idx(i+1, j, Nx, Ny, 2)]-4*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)]);
            
            IS1_nx[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(i+1, j, Nx, Ny, 2)]+taondata[idx(i+2, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(i+1, j, Nx, Ny, 2)]+taondata[idx(i+2, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(i+2, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(i+2, j, Nx, Ny, 2)]);
            
            IS2_nx[2]=(13/12)*(taondata[idx(Nx-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(i+1, j, Nx, Ny, 2)])*(taondata[idx(Nx-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(i+1, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(Nx-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(i+1, j, Nx, Ny, 2)])*(taondata[idx(Nx-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(i+1, j, Nx, Ny, 2)]);
            
            
            IS0_px[3]=(13/12)*(Cjpdata[idx(Nx-2, j, Nx, Ny, 0)]-2*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(Nx-2, j, Nx, Ny, 0)]-2*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(Nx-2, j, Nx, Ny, 0)]-4*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(Nx-2, j, Nx, Ny, 0)]-4*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)]);
            
            IS1_px[3]=(13/12)*(Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(i+1, j, Nx, Ny, 0)])*(Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(i+1, j, Nx, Ny, 0)])*(Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(i+1, j, Nx, Ny, 0)]);
            
            IS2_px[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS0_nx[3]=(13/12)*(Cjndata[idx(i+1, j, Nx, Ny, 0)]-2*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)])*(Cjndata[idx(i+1, j, Nx, Ny, 0)]-2*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)])+(1/4)*(3*Cjndata[idx(i+1, j, Nx, Ny, 0)]-4*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)])*(3*Cjndata[idx(i+1, j, Nx, Ny, 0)]-4*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)]);
            
            IS1_nx[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(i+1, j, Nx, Ny, 0)]+Cjndata[idx(i+2, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(i+1, j, Nx, Ny, 0)]+Cjndata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(i+2, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS2_nx[3]=(13/12)*(Cjndata[idx(Nx-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(i+1, j, Nx, Ny, 0)])*(Cjndata[idx(Nx-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(Nx-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(i+1, j, Nx, Ny, 0)])*(Cjndata[idx(Nx-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(i+1, j, Nx, Ny, 0)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0px[k]=(1/10)*(1/(Epsilon+IS0_px[k]))*(1/(Epsilon+IS0_px[k]));
                alpha_1px[k]=(6/10)*(1/(Epsilon+IS1_px[k]))*(1/(Epsilon+IS1_px[k]));
                alpha_2px[k]=(3/10)*(1/(Epsilon+IS2_px[k]))*(1/(Epsilon+IS2_px[k]));
                alpha_0nx[k]=(1/10)*(1/(Epsilon+IS0_nx[k]))*(1/(Epsilon+IS0_nx[k]));
                alpha_1nx[k]=(6/10)*(1/(Epsilon+IS1_nx[k]))*(1/(Epsilon+IS1_nx[k]));
                alpha_2nx[k]=(3/10)*(1/(Epsilon+IS2_nx[k]))*(1/(Epsilon+IS2_nx[k]));
                
                w0_px[k]=alpha_0px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w1_px[k]=alpha_1px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w2_px[k]=alpha_2px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w0_nx[k]=alpha_0nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
                w1_nx[k]=alpha_1nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
                w2_nx[k]=alpha_2nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
            }
            
            u_tpphx[0]=w0_px[0]*((2/6)*ypdata[idx(Nx-2, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(Nx-1, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(Nx-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i+1, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i+1, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(i+2, j, Nx, Ny, 1)]);
            
            u_tnphx[0]=w2_nx[0]*((-1/6)*yndata[idx(Nx-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+1, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i+1, j, Nx, Ny, 1)]-(1/6)*yndata[idx(i+2, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i+1, j, Nx, Ny, 1)]-(7/6)*yndata[idx(i+2, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+3, j, Nx, Ny, 1)]);
            
            u_tpnhx[0]=w0_px[0]*((2/6)*ypdata[idx(Nx-3, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(Nx-2, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(Nx-1, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(Nx-2, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(Nx-1, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(Nx-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(i+1, j, Nx, Ny, 1)]);
            
            u_tnnhx[0]=w2_nx[0]*((-1/6)*yndata[idx(Nx-2, j, Nx, Ny, 1)]+(5/6)*yndata[idx(Nx-1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(Nx-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]-(1/6)*yndata[idx(i+1, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 1)]-(7/6)*yndata[idx(i+1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+2, j, Nx, Ny, 1)]);
            
            
            u_tpphx[1]=w0_px[1]*((2/6)*taopdata[idx(Nx-2, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(Nx-1, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(Nx-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i+1, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(i+2, j, Nx, Ny, 0)]);
            
            u_tnphx[1]=w2_nx[1]*((-1/6)*taondata[idx(Nx-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+1, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*taondata[idx(i+2, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i+1, j, Nx, Ny, 0)]-(7/6)*taondata[idx(i+2, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+3, j, Nx, Ny, 0)]);
            
            u_tpnhx[1]=w0_px[1]*((2/6)*taopdata[idx(Nx-3, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(Nx-2, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(Nx-1, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(Nx-2, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(Nx-1, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(Nx-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(i+1, j, Nx, Ny, 0)]);
            
            u_tnnhx[1]=w2_nx[1]*((-1/6)*taondata[idx(Nx-2, j, Nx, Ny, 0)]+(5/6)*taondata[idx(Nx-1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(Nx-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]-(1/6)*taondata[idx(i+1, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 0)]-(7/6)*taondata[idx(i+1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+2, j, Nx, Ny, 0)]);
            
            
            u_tpphx[2]=w0_px[2]*((2/6)*taopdata[idx(Nx-2, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(Nx-1, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(Nx-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i+1, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i+1, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(i+2, j, Nx, Ny, 2)]);
            
            u_tnphx[2]=w2_nx[2]*((-1/6)*taondata[idx(Nx-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+1, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i+1, j, Nx, Ny, 2)]-(1/6)*taondata[idx(i+2, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i+1, j, Nx, Ny, 2)]-(7/6)*taondata[idx(i+2, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+3, j, Nx, Ny, 2)]);
            
            u_tpnhx[2]=w0_px[2]*((2/6)*taopdata[idx(Nx-3, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(Nx-2, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(Nx-1, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(Nx-2, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(Nx-1, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(Nx-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(i+1, j, Nx, Ny, 2)]);
            
            u_tnnhx[2]=w2_nx[2]*((-1/6)*taondata[idx(Nx-2, j, Nx, Ny, 2)]+(5/6)*taondata[idx(Nx-1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(Nx-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]-(1/6)*taondata[idx(i+1, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 2)]-(7/6)*taondata[idx(i+1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+2, j, Nx, Ny, 2)]);
            
            
            u_tpphx[3]=w0_px[3]*((2/6)*Cjpdata[idx(Nx-2, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(i+2, j, Nx, Ny, 0)]);
            
            u_tnphx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(Nx-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+3, j, Nx, Ny, 0)]);
            
            u_tpnhx[3]=w0_px[3]*((2/6)*Cjpdata[idx(Nx-3, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(Nx-2, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(Nx-2, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)]);
            
            u_tnnhx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(Nx-2, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(Nx-1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(Nx-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)]);
            
            for(k=0;k<4;k++){
                yxdata[idx(i, j, Nx, Ny, k)]=-(1/dx)*((u_tpphx[k]-u_tpnhx[k])+(u_tnphx[k]-u_tnnhx[k]));
            }
    }
    
        j=0;
        for (i=0; i<Nx; i++){
            // x direction
            IS0_py[0]=(13/12)*(ypdata[idx(i, Ny-2, Nx, Ny, 2)]-2*ypdata[idx(i, Ny-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, Ny-2, Nx, Ny, 2)]-2*ypdata[idx(i, Ny-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, Ny-2, Nx, Ny, 2)]-4*ypdata[idx(i, Ny-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, Ny-2, Nx, Ny, 2)]-4*ypdata[idx(i, Ny-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_py[0]=(13/12)*(ypdata[idx(i, Ny-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, j+1, Nx, Ny, 2)])*(ypdata[idx(i, Ny-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, j+1, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, Ny-1, Nx, Ny, 2)]-1*ypdata[idx(i, j+1, Nx, Ny, 2)])*(ypdata[idx(i, Ny-1, Nx, Ny, 2)]-1*ypdata[idx(i, j+1, Nx, Ny, 2)]);
            
            IS2_py[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)]);
            
            IS0_ny[0]=(13/12)*(yndata[idx(i, j+1, Nx, Ny, 2)]-2*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)])*(yndata[idx(i, j+1, Nx, Ny, 2)]-2*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)])+(1/4)*(3*yndata[idx(i, j+1, Nx, Ny, 2)]-4*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)])*(3*yndata[idx(i, j+1, Nx, Ny, 2)]-4*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)]);
            
            IS1_ny[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, j+1, Nx, Ny, 2)]+yndata[idx(i, j+2, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, j+1, Nx, Ny, 2)]+yndata[idx(i, j+2, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, j+2, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, j+2, Nx, Ny, 2)]);
            
            IS2_ny[0]=(13/12)*(yndata[idx(i, Ny-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, j+1, Nx, Ny, 2)])*(yndata[idx(i, Ny-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, j+1, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, Ny-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, j+1, Nx, Ny, 2)])*(yndata[idx(i, Ny-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, j+1, Nx, Ny, 2)]);
            
            
            IS0_py[1]=(13/12)*(taopdata[idx(i, Ny-2, Nx, Ny, 1)]-2*taopdata[idx(i, Ny-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, Ny-2, Nx, Ny, 1)]-2*taopdata[idx(i, Ny-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, Ny-2, Nx, Ny, 1)]-4*taopdata[idx(i, Ny-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, Ny-2, Nx, Ny, 1)]-4*taopdata[idx(i, Ny-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[1]=(13/12)*(taopdata[idx(i, Ny-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, j+1, Nx, Ny, 1)])*(taopdata[idx(i, Ny-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, Ny-1, Nx, Ny, 1)]-1*taopdata[idx(i, j+1, Nx, Ny, 1)])*(taopdata[idx(i, Ny-1, Nx, Ny, 1)]-1*taopdata[idx(i, j+1, Nx, Ny, 1)]);
            
            IS2_py[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS0_ny[1]=(13/12)*(taondata[idx(i, j+1, Nx, Ny, 1)]-2*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)])*(taondata[idx(i, j+1, Nx, Ny, 1)]-2*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)])+(1/4)*(3*taondata[idx(i, j+1, Nx, Ny, 1)]-4*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)])*(3*taondata[idx(i, j+1, Nx, Ny, 1)]-4*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)]);
            
            IS1_ny[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, j+1, Nx, Ny, 1)]+taondata[idx(i, j+2, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, j+1, Nx, Ny, 1)]+taondata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, j+2, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS2_ny[1]=(13/12)*(taondata[idx(i, Ny-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, j+1, Nx, Ny, 1)])*(taondata[idx(i, Ny-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, Ny-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, j+1, Nx, Ny, 1)])*(taondata[idx(i, Ny-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, j+1, Nx, Ny, 1)]);
            
            
            IS0_py[2]=(13/12)*(taopdata[idx(i, Ny-2, Nx, Ny, 3)]-2*taopdata[idx(i, Ny-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, Ny-2, Nx, Ny, 3)]-2*taopdata[idx(i, Ny-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, Ny-2, Nx, Ny, 3)]-4*taopdata[idx(i, Ny-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, Ny-2, Nx, Ny, 3)]-4*taopdata[idx(i, Ny-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)]);
            
            IS1_py[2]=(13/12)*(taopdata[idx(i, Ny-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, j+1, Nx, Ny, 3)])*(taopdata[idx(i, Ny-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, j+1, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, Ny-1, Nx, Ny, 3)]-1*taopdata[idx(i, j+1, Nx, Ny, 3)])*(taopdata[idx(i, Ny-1, Nx, Ny, 3)]-1*taopdata[idx(i, j+1, Nx, Ny, 3)]);
            
            IS2_py[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)]);
            
            IS0_ny[2]=(13/12)*(taondata[idx(i, j+1, Nx, Ny, 3)]-2*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)])*(taondata[idx(i, j+1, Nx, Ny, 3)]-2*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)])+(1/4)*(3*taondata[idx(i, j+1, Nx, Ny, 3)]-4*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)])*(3*taondata[idx(i, j+1, Nx, Ny, 3)]-4*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)]);
            
            IS1_ny[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, j+1, Nx, Ny, 3)]+taondata[idx(i, j+2, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, j+1, Nx, Ny, 3)]+taondata[idx(i, j+2, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, j+2, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, j+2, Nx, Ny, 3)]);
            
            IS2_ny[2]=(13/12)*(taondata[idx(i, Ny-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, j+1, Nx, Ny, 3)])*(taondata[idx(i, Ny-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, j+1, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, Ny-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, j+1, Nx, Ny, 3)])*(taondata[idx(i, Ny-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, j+1, Nx, Ny, 3)]);
            
            
            IS0_py[3]=(13/12)*(Cjpdata[idx(i, Ny-2, Nx, Ny, 1)]-2*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, Ny-2, Nx, Ny, 1)]-2*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, Ny-2, Nx, Ny, 1)]-4*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, Ny-2, Nx, Ny, 1)]-4*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[3]=(13/12)*(Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, j+1, Nx, Ny, 1)])*(Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, j+1, Nx, Ny, 1)])*(Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, j+1, Nx, Ny, 1)]);
            
            IS2_py[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS0_ny[3]=(13/12)*(Cjndata[idx(i, j+1, Nx, Ny, 1)]-2*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)])*(Cjndata[idx(i, j+1, Nx, Ny, 1)]-2*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)])+(1/4)*(3*Cjndata[idx(i, j+1, Nx, Ny, 1)]-4*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)])*(3*Cjndata[idx(i, j+1, Nx, Ny, 1)]-4*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)]);
            
            IS1_ny[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, j+1, Nx, Ny, 1)]+Cjndata[idx(i, j+2, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, j+1, Nx, Ny, 1)]+Cjndata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, j+2, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS2_ny[3]=(13/12)*(Cjndata[idx(i, Ny-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, j+1, Nx, Ny, 1)])*(Cjndata[idx(i, Ny-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, Ny-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, j+1, Nx, Ny, 1)])*(Cjndata[idx(i, Ny-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, j+1, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0py[k]=(1/10)*(1/(Epsilon+IS0_py[k]))*(1/(Epsilon+IS0_py[k]));
                alpha_1py[k]=(6/10)*(1/(Epsilon+IS1_py[k]))*(1/(Epsilon+IS1_py[k]));
                alpha_2py[k]=(3/10)*(1/(Epsilon+IS2_py[k]))*(1/(Epsilon+IS2_py[k]));
                alpha_0ny[k]=(1/10)*(1/(Epsilon+IS0_ny[k]))*(1/(Epsilon+IS0_ny[k]));
                alpha_1ny[k]=(6/10)*(1/(Epsilon+IS1_ny[k]))*(1/(Epsilon+IS1_ny[k]));
                alpha_2ny[k]=(3/10)*(1/(Epsilon+IS2_ny[k]))*(1/(Epsilon+IS2_ny[k]));
                
                w0_py[k]=alpha_0py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w1_py[k]=alpha_1py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w2_py[k]=alpha_2py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w0_ny[k]=alpha_0ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w1_ny[k]=alpha_1ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w2_ny[k]=alpha_2ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
            }
            
            u_tpphy[0]=w0_py[0]*((2/6)*ypdata[idx(i, Ny-2, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, Ny-1, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, Ny-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j+1, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j+1, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, j+2, Nx, Ny, 2)]);
            
            u_tnphy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, Ny-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+1, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j+1, Nx, Ny, 2)]-(1/6)*yndata[idx(i, j+2, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j+1, Nx, Ny, 2)]-(7/6)*yndata[idx(i, j+2, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+3, Nx, Ny, 2)]);
            
            u_tpnhy[0]=w0_py[0]*((2/6)*ypdata[idx(i, Ny-3, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, Ny-2, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, Ny-1, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, Ny-2, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, Ny-1, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, Ny-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, j+1, Nx, Ny, 2)]);
            
            u_tnnhy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, Ny-2, Nx, Ny, 2)]+(5/6)*yndata[idx(i, Ny-1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, Ny-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]-(1/6)*yndata[idx(i, j+1, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 2)]-(7/6)*yndata[idx(i, j+1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+2, Nx, Ny, 2)]);
            
            
            u_tpphy[1]=w0_py[1]*((2/6)*taopdata[idx(i, Ny-2, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, Ny-1, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, Ny-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j+1, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, j+2, Nx, Ny, 1)]);
            
            u_tnphy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, Ny-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+1, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*taondata[idx(i, j+2, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j+1, Nx, Ny, 1)]-(7/6)*taondata[idx(i, j+2, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+3, Nx, Ny, 1)]);
            
            u_tpnhy[1]=w0_py[1]*((2/6)*taopdata[idx(i, Ny-3, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, Ny-2, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, Ny-1, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, Ny-2, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, Ny-1, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, Ny-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, j+1, Nx, Ny, 1)]);
            
            u_tnnhy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, Ny-2, Nx, Ny, 1)]+(5/6)*taondata[idx(i, Ny-1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, Ny-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]-(1/6)*taondata[idx(i, j+1, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 1)]-(7/6)*taondata[idx(i, j+1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+2, Nx, Ny, 1)]);
            
            
            u_tpphy[2]=w0_py[2]*((2/6)*taopdata[idx(i, Ny-2, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, Ny-1, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, Ny-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j+1, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j+1, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, j+2, Nx, Ny, 3)]);
            
            u_tnphy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, Ny-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+1, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j+1, Nx, Ny, 3)]-(1/6)*taondata[idx(i, j+2, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j+1, Nx, Ny, 3)]-(7/6)*taondata[idx(i, j+2, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+3, Nx, Ny, 3)]);
            
            u_tpnhy[2]=w0_py[2]*((2/6)*taopdata[idx(i, Ny-3, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, Ny-2, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, Ny-1, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, Ny-2, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, Ny-1, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, Ny-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, j+1, Nx, Ny, 3)]);
            
            u_tnnhy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, Ny-2, Nx, Ny, 3)]+(5/6)*taondata[idx(i, Ny-1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, Ny-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]-(1/6)*taondata[idx(i, j+1, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 3)]-(7/6)*taondata[idx(i, j+1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+2, Nx, Ny, 3)]);
            
            
            u_tpphy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, Ny-2, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, j+2, Nx, Ny, 1)]);
            
            u_tnphy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, Ny-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+3, Nx, Ny, 1)]);
            
            u_tpnhy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, Ny-3, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, Ny-2, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, Ny-2, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)]);
            
            u_tnnhy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, Ny-2, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, Ny-1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, Ny-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++){
                yydata[idx(i, j, Nx, Ny, k)]=-(1/dy)*((u_tpphy[k]-u_tpnhy[k])+(u_tnphy[k]-u_tnnhy[k]));
            }
        }

            i=1;
    for(j=0; j<Ny; j++){
            // x direction
            IS0_px[0]=(13/12)*(ypdata[idx(Nx-1, j, Nx, Ny, 1)]-2*ypdata[idx(i-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(Nx-1, j, Nx, Ny, 1)]-2*ypdata[idx(i-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(Nx-1, j, Nx, Ny, 1)]-4*ypdata[idx(i-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(Nx-1, j, Nx, Ny, 1)]-4*ypdata[idx(i-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_px[0]=(13/12)*(ypdata[idx(i-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(i+1, j, Nx, Ny, 1)])*(ypdata[idx(i-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(i+1, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(i-1, j, Nx, Ny, 1)]-1*ypdata[idx(i+1, j, Nx, Ny, 1)])*(ypdata[idx(i-1, j, Nx, Ny, 1)]-1*ypdata[idx(i+1, j, Nx, Ny, 1)]);
            
            IS2_px[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)]);
            
            IS0_nx[0]=(13/12)*(yndata[idx(i+1, j, Nx, Ny, 1)]-2*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)])*(yndata[idx(i+1, j, Nx, Ny, 1)]-2*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)])+(1/4)*(3*yndata[idx(i+1, j, Nx, Ny, 1)]-4*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)])*(3*yndata[idx(i+1, j, Nx, Ny, 1)]-4*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)]);
            
            IS1_nx[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(i+1, j, Nx, Ny, 1)]+yndata[idx(i+2, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(i+1, j, Nx, Ny, 1)]+yndata[idx(i+2, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(i+2, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(i+2, j, Nx, Ny, 1)]);
            
            IS2_nx[0]=(13/12)*(yndata[idx(i-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(i+1, j, Nx, Ny, 1)])*(yndata[idx(i-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(i+1, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(i+1, j, Nx, Ny, 1)])*(yndata[idx(i-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(i+1, j, Nx, Ny, 1)]);
            
            
            IS0_px[1]=(13/12)*(taopdata[idx(Nx-1, j, Nx, Ny, 0)]-2*taopdata[idx(i-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(Nx-1, j, Nx, Ny, 0)]-2*taopdata[idx(i-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(Nx-1, j, Nx, Ny, 0)]-4*taopdata[idx(i-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(Nx-1, j, Nx, Ny, 0)]-4*taopdata[idx(i-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)]);
            
            IS1_px[1]=(13/12)*(taopdata[idx(i-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(i+1, j, Nx, Ny, 0)])*(taopdata[idx(i-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(i-1, j, Nx, Ny, 0)]-1*taopdata[idx(i+1, j, Nx, Ny, 0)])*(taopdata[idx(i-1, j, Nx, Ny, 0)]-1*taopdata[idx(i+1, j, Nx, Ny, 0)]);
            
            IS2_px[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS0_nx[1]=(13/12)*(taondata[idx(i+1, j, Nx, Ny, 0)]-2*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)])*(taondata[idx(i+1, j, Nx, Ny, 0)]-2*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)])+(1/4)*(3*taondata[idx(i+1, j, Nx, Ny, 0)]-4*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)])*(3*taondata[idx(i+1, j, Nx, Ny, 0)]-4*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)]);
            
            IS1_nx[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(i+1, j, Nx, Ny, 0)]+taondata[idx(i+2, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(i+1, j, Nx, Ny, 0)]+taondata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(i+2, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS2_nx[1]=(13/12)*(taondata[idx(i-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(i+1, j, Nx, Ny, 0)])*(taondata[idx(i-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(i+1, j, Nx, Ny, 0)])*(taondata[idx(i-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(i+1, j, Nx, Ny, 0)]);
            
            
            IS0_px[2]=(13/12)*(taopdata[idx(Nx-1, j, Nx, Ny, 2)]-2*taopdata[idx(i-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(Nx-1, j, Nx, Ny, 2)]-2*taopdata[idx(i-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(Nx-1, j, Nx, Ny, 2)]-4*taopdata[idx(i-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(Nx-1, j, Nx, Ny, 2)]-4*taopdata[idx(i-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_px[2]=(13/12)*(taopdata[idx(i-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(i+1, j, Nx, Ny, 2)])*(taopdata[idx(i-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(i+1, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(i-1, j, Nx, Ny, 2)]-1*taopdata[idx(i+1, j, Nx, Ny, 2)])*(taopdata[idx(i-1, j, Nx, Ny, 2)]-1*taopdata[idx(i+1, j, Nx, Ny, 2)]);
            
            IS2_px[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)]);
            
            IS0_nx[2]=(13/12)*(taondata[idx(i+1, j, Nx, Ny, 2)]-2*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)])*(taondata[idx(i+1, j, Nx, Ny, 2)]-2*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)])+(1/4)*(3*taondata[idx(i+1, j, Nx, Ny, 2)]-4*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)])*(3*taondata[idx(i+1, j, Nx, Ny, 2)]-4*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)]);
            
            IS1_nx[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(i+1, j, Nx, Ny, 2)]+taondata[idx(i+2, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(i+1, j, Nx, Ny, 2)]+taondata[idx(i+2, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(i+2, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(i+2, j, Nx, Ny, 2)]);
            
            IS2_nx[2]=(13/12)*(taondata[idx(i-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(i+1, j, Nx, Ny, 2)])*(taondata[idx(i-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(i+1, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(i+1, j, Nx, Ny, 2)])*(taondata[idx(i-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(i+1, j, Nx, Ny, 2)]);
            
            
            IS0_px[3]=(13/12)*(Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]-4*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]-4*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)]);
            
            IS1_px[3]=(13/12)*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(i+1, j, Nx, Ny, 0)])*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(i+1, j, Nx, Ny, 0)])*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(i+1, j, Nx, Ny, 0)]);
            
            IS2_px[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS0_nx[3]=(13/12)*(Cjndata[idx(i+1, j, Nx, Ny, 0)]-2*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)])*(Cjndata[idx(i+1, j, Nx, Ny, 0)]-2*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)])+(1/4)*(3*Cjndata[idx(i+1, j, Nx, Ny, 0)]-4*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)])*(3*Cjndata[idx(i+1, j, Nx, Ny, 0)]-4*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)]);
            
            IS1_nx[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(i+1, j, Nx, Ny, 0)]+Cjndata[idx(i+2, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(i+1, j, Nx, Ny, 0)]+Cjndata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(i+2, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS2_nx[3]=(13/12)*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(i+1, j, Nx, Ny, 0)])*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(i+1, j, Nx, Ny, 0)])*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(i+1, j, Nx, Ny, 0)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0px[k]=(1/10)*(1/(Epsilon+IS0_px[k]))*(1/(Epsilon+IS0_px[k]));
                alpha_1px[k]=(6/10)*(1/(Epsilon+IS1_px[k]))*(1/(Epsilon+IS1_px[k]));
                alpha_2px[k]=(3/10)*(1/(Epsilon+IS2_px[k]))*(1/(Epsilon+IS2_px[k]));
                alpha_0nx[k]=(1/10)*(1/(Epsilon+IS0_nx[k]))*(1/(Epsilon+IS0_nx[k]));
                alpha_1nx[k]=(6/10)*(1/(Epsilon+IS1_nx[k]))*(1/(Epsilon+IS1_nx[k]));
                alpha_2nx[k]=(3/10)*(1/(Epsilon+IS2_nx[k]))*(1/(Epsilon+IS2_nx[k]));
                
                w0_px[k]=alpha_0px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w1_px[k]=alpha_1px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w2_px[k]=alpha_2px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w0_nx[k]=alpha_0nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
                w1_nx[k]=alpha_1nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
                w2_nx[k]=alpha_2nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
            }
            
            u_tpphx[0]=w0_px[0]*((2/6)*ypdata[idx(Nx-1, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i+1, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i+1, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(i+2, j, Nx, Ny, 1)]);
            
            u_tnphx[0]=w2_nx[0]*((-1/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+1, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i+1, j, Nx, Ny, 1)]-(1/6)*yndata[idx(i+2, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i+1, j, Nx, Ny, 1)]-(7/6)*yndata[idx(i+2, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+3, j, Nx, Ny, 1)]);
            
            u_tpnhx[0]=w0_px[0]*((2/6)*ypdata[idx(Nx-2, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(Nx-1, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i-1, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(Nx-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(i+1, j, Nx, Ny, 1)]);
            
            u_tnnhx[0]=w2_nx[0]*((-1/6)*yndata[idx(Nx-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]-(1/6)*yndata[idx(i+1, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 1)]-(7/6)*yndata[idx(i+1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+2, j, Nx, Ny, 1)]);
            
            
            u_tpphx[1]=w0_px[1]*((2/6)*taopdata[idx(Nx-1, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i+1, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(i+2, j, Nx, Ny, 0)]);
            
            u_tnphx[1]=w2_nx[1]*((-1/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+1, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*taondata[idx(i+2, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i+1, j, Nx, Ny, 0)]-(7/6)*taondata[idx(i+2, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+3, j, Nx, Ny, 0)]);
            
            u_tpnhx[1]=w0_px[1]*((2/6)*taopdata[idx(Nx-2, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(Nx-1, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i-1, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(Nx-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(i+1, j, Nx, Ny, 0)]);
            
            u_tnnhx[1]=w2_nx[1]*((-1/6)*taondata[idx(Nx-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]-(1/6)*taondata[idx(i+1, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 0)]-(7/6)*taondata[idx(i+1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+2, j, Nx, Ny, 0)]);
            
            
            u_tpphx[2]=w0_px[2]*((2/6)*taopdata[idx(Nx-1, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i+1, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i+1, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(i+2, j, Nx, Ny, 2)]);
            
            u_tnphx[2]=w2_nx[2]*((-1/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+1, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i+1, j, Nx, Ny, 2)]-(1/6)*taondata[idx(i+2, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i+1, j, Nx, Ny, 2)]-(7/6)*taondata[idx(i+2, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+3, j, Nx, Ny, 2)]);
            
            u_tpnhx[2]=w0_px[2]*((2/6)*taopdata[idx(Nx-2, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(Nx-1, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i-1, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(Nx-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(i+1, j, Nx, Ny, 2)]);
            
            u_tnnhx[2]=w2_nx[2]*((-1/6)*taondata[idx(Nx-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]-(1/6)*taondata[idx(i+1, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 2)]-(7/6)*taondata[idx(i+1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+2, j, Nx, Ny, 2)]);
            
            
            u_tpphx[3]=w0_px[3]*((2/6)*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(i+2, j, Nx, Ny, 0)]);
            
            u_tnphx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+3, j, Nx, Ny, 0)]);
            
            u_tpnhx[3]=w0_px[3]*((2/6)*Cjpdata[idx(Nx-2, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)]);
            
            u_tnnhx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(Nx-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)]);
            
            for(k=0;k<4;k++){
                yxdata[idx(i, j, Nx, Ny, k)]=-(1/dx)*((u_tpphx[k]-u_tpnhx[k])+(u_tnphx[k]-u_tnnhx[k]));
            }
    }
    
            j=1;
        for (i=0; i<Nx; i++){
            // x direction
            IS0_py[0]=(13/12)*(ypdata[idx(i, Ny-1, Nx, Ny, 2)]-2*ypdata[idx(i, j-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, Ny-1, Nx, Ny, 2)]-2*ypdata[idx(i, j-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, Ny-1, Nx, Ny, 2)]-4*ypdata[idx(i, j-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, Ny-1, Nx, Ny, 2)]-4*ypdata[idx(i, j-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_py[0]=(13/12)*(ypdata[idx(i, j-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, j+1, Nx, Ny, 2)])*(ypdata[idx(i, j-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, j+1, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, j-1, Nx, Ny, 2)]-1*ypdata[idx(i, j+1, Nx, Ny, 2)])*(ypdata[idx(i, j-1, Nx, Ny, 2)]-1*ypdata[idx(i, j+1, Nx, Ny, 2)]);
            
            IS2_py[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)]);
            
            IS0_ny[0]=(13/12)*(yndata[idx(i, j+1, Nx, Ny, 2)]-2*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)])*(yndata[idx(i, j+1, Nx, Ny, 2)]-2*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)])+(1/4)*(3*yndata[idx(i, j+1, Nx, Ny, 2)]-4*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)])*(3*yndata[idx(i, j+1, Nx, Ny, 2)]-4*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)]);
            
            IS1_ny[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, j+1, Nx, Ny, 2)]+yndata[idx(i, j+2, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, j+1, Nx, Ny, 2)]+yndata[idx(i, j+2, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, j+2, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, j+2, Nx, Ny, 2)]);
            
            IS2_ny[0]=(13/12)*(yndata[idx(i, j-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, j+1, Nx, Ny, 2)])*(yndata[idx(i, j-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, j+1, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, j+1, Nx, Ny, 2)])*(yndata[idx(i, j-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, j+1, Nx, Ny, 2)]);
            
            
            IS0_py[1]=(13/12)*(taopdata[idx(i, Ny-1, Nx, Ny, 1)]-2*taopdata[idx(i, j-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, Ny-1, Nx, Ny, 1)]-2*taopdata[idx(i, j-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, Ny-1, Nx, Ny, 1)]-4*taopdata[idx(i, j-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, Ny-1, Nx, Ny, 1)]-4*taopdata[idx(i, j-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[1]=(13/12)*(taopdata[idx(i, j-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, j+1, Nx, Ny, 1)])*(taopdata[idx(i, j-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, j-1, Nx, Ny, 1)]-1*taopdata[idx(i, j+1, Nx, Ny, 1)])*(taopdata[idx(i, j-1, Nx, Ny, 1)]-1*taopdata[idx(i, j+1, Nx, Ny, 1)]);
            
            IS2_py[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS0_ny[1]=(13/12)*(taondata[idx(i, j+1, Nx, Ny, 1)]-2*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)])*(taondata[idx(i, j+1, Nx, Ny, 1)]-2*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)])+(1/4)*(3*taondata[idx(i, j+1, Nx, Ny, 1)]-4*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)])*(3*taondata[idx(i, j+1, Nx, Ny, 1)]-4*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)]);
            
            IS1_ny[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, j+1, Nx, Ny, 1)]+taondata[idx(i, j+2, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, j+1, Nx, Ny, 1)]+taondata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, j+2, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS2_ny[1]=(13/12)*(taondata[idx(i, j-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, j+1, Nx, Ny, 1)])*(taondata[idx(i, j-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, j+1, Nx, Ny, 1)])*(taondata[idx(i, j-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, j+1, Nx, Ny, 1)]);
            
            
            IS0_py[2]=(13/12)*(taopdata[idx(i, Ny-1, Nx, Ny, 3)]-2*taopdata[idx(i, j-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, Ny-1, Nx, Ny, 3)]-2*taopdata[idx(i, j-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, Ny-1, Nx, Ny, 3)]-4*taopdata[idx(i, j-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, Ny-1, Nx, Ny, 3)]-4*taopdata[idx(i, j-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)]);
            
            IS1_py[2]=(13/12)*(taopdata[idx(i, j-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, j+1, Nx, Ny, 3)])*(taopdata[idx(i, j-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, j+1, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, j-1, Nx, Ny, 3)]-1*taopdata[idx(i, j+1, Nx, Ny, 3)])*(taopdata[idx(i, j-1, Nx, Ny, 3)]-1*taopdata[idx(i, j+1, Nx, Ny, 3)]);
            
            IS2_py[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)]);
            
            IS0_ny[2]=(13/12)*(taondata[idx(i, j+1, Nx, Ny, 3)]-2*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)])*(taondata[idx(i, j+1, Nx, Ny, 3)]-2*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)])+(1/4)*(3*taondata[idx(i, j+1, Nx, Ny, 3)]-4*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)])*(3*taondata[idx(i, j+1, Nx, Ny, 3)]-4*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)]);
            
            IS1_ny[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, j+1, Nx, Ny, 3)]+taondata[idx(i, j+2, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, j+1, Nx, Ny, 3)]+taondata[idx(i, j+2, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, j+2, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, j+2, Nx, Ny, 3)]);
            
            IS2_ny[2]=(13/12)*(taondata[idx(i, j-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, j+1, Nx, Ny, 3)])*(taondata[idx(i, j-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, j+1, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, j+1, Nx, Ny, 3)])*(taondata[idx(i, j-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, j+1, Nx, Ny, 3)]);
            
            
            IS0_py[3]=(13/12)*(Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]-4*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]-4*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[3]=(13/12)*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, j+1, Nx, Ny, 1)])*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, j+1, Nx, Ny, 1)])*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, j+1, Nx, Ny, 1)]);
            
            IS2_py[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS0_ny[3]=(13/12)*(Cjndata[idx(i, j+1, Nx, Ny, 1)]-2*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)])*(Cjndata[idx(i, j+1, Nx, Ny, 1)]-2*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)])+(1/4)*(3*Cjndata[idx(i, j+1, Nx, Ny, 1)]-4*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)])*(3*Cjndata[idx(i, j+1, Nx, Ny, 1)]-4*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)]);
            
            IS1_ny[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, j+1, Nx, Ny, 1)]+Cjndata[idx(i, j+2, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, j+1, Nx, Ny, 1)]+Cjndata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, j+2, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS2_ny[3]=(13/12)*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, j+1, Nx, Ny, 1)])*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, j+1, Nx, Ny, 1)])*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, j+1, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0py[k]=(1/10)*(1/(Epsilon+IS0_py[k]))*(1/(Epsilon+IS0_py[k]));
                alpha_1py[k]=(6/10)*(1/(Epsilon+IS1_py[k]))*(1/(Epsilon+IS1_py[k]));
                alpha_2py[k]=(3/10)*(1/(Epsilon+IS2_py[k]))*(1/(Epsilon+IS2_py[k]));
                alpha_0ny[k]=(1/10)*(1/(Epsilon+IS0_ny[k]))*(1/(Epsilon+IS0_ny[k]));
                alpha_1ny[k]=(6/10)*(1/(Epsilon+IS1_ny[k]))*(1/(Epsilon+IS1_ny[k]));
                alpha_2ny[k]=(3/10)*(1/(Epsilon+IS2_ny[k]))*(1/(Epsilon+IS2_ny[k]));
                
                w0_py[k]=alpha_0py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w1_py[k]=alpha_1py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w2_py[k]=alpha_2py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w0_ny[k]=alpha_0ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w1_ny[k]=alpha_1ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w2_ny[k]=alpha_2ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
            }
            
            u_tpphy[0]=w0_py[0]*((2/6)*ypdata[idx(i, Ny-1, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j+1, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j+1, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, j+2, Nx, Ny, 2)]);
            
            u_tnphy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+1, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j+1, Nx, Ny, 2)]-(1/6)*yndata[idx(i, j+2, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j+1, Nx, Ny, 2)]-(7/6)*yndata[idx(i, j+2, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+3, Nx, Ny, 2)]);
            
            u_tpnhy[0]=w0_py[0]*((2/6)*ypdata[idx(i, Ny-2, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, Ny-1, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j-1, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, Ny-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, j+1, Nx, Ny, 2)]);
            
            u_tnnhy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, Ny-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]-(1/6)*yndata[idx(i, j+1, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 2)]-(7/6)*yndata[idx(i, j+1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+2, Nx, Ny, 2)]);
            
            
            u_tpphy[1]=w0_py[1]*((2/6)*taopdata[idx(i, Ny-1, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j+1, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, j+2, Nx, Ny, 1)]);
            
            u_tnphy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+1, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*taondata[idx(i, j+2, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j+1, Nx, Ny, 1)]-(7/6)*taondata[idx(i, j+2, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+3, Nx, Ny, 1)]);
            
            u_tpnhy[1]=w0_py[1]*((2/6)*taopdata[idx(i, Ny-2, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, Ny-1, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j-1, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, Ny-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, j+1, Nx, Ny, 1)]);
            
            u_tnnhy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, Ny-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]-(1/6)*taondata[idx(i, j+1, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 1)]-(7/6)*taondata[idx(i, j+1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+2, Nx, Ny, 1)]);
            
            
            u_tpphy[2]=w0_py[2]*((2/6)*taopdata[idx(i, Ny-1, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j+1, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j+1, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, j+2, Nx, Ny, 3)]);
            
            u_tnphy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+1, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j+1, Nx, Ny, 3)]-(1/6)*taondata[idx(i, j+2, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j+1, Nx, Ny, 3)]-(7/6)*taondata[idx(i, j+2, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+3, Nx, Ny, 3)]);
            
            u_tpnhy[2]=w0_py[2]*((2/6)*taopdata[idx(i, Ny-2, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, Ny-1, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j-1, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, Ny-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, j+1, Nx, Ny, 3)]);
            
            u_tnnhy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, Ny-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]-(1/6)*taondata[idx(i, j+1, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 3)]-(7/6)*taondata[idx(i, j+1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+2, Nx, Ny, 3)]);
            
            
            u_tpphy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, j+2, Nx, Ny, 1)]);
            
            u_tnphy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+3, Nx, Ny, 1)]);
            
            u_tpnhy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, Ny-2, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)]);
            
            u_tnnhy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, Ny-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++){
                yydata[idx(i, j, Nx, Ny, k)]=-(1/dy)*((u_tpphy[k]-u_tpnhy[k])+(u_tnphy[k]-u_tnnhy[k]));
            }
        }

            i=2;
    for(j=0; j<Ny; j++){
            // x direction
            IS0_px[0]=(13/12)*(ypdata[idx(i-2, j, Nx, Ny, 1)]-2*ypdata[idx(i-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(i-2, j, Nx, Ny, 1)]-2*ypdata[idx(i-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(i-2, j, Nx, Ny, 1)]-4*ypdata[idx(i-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(i-2, j, Nx, Ny, 1)]-4*ypdata[idx(i-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_px[0]=(13/12)*(ypdata[idx(i-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(i+1, j, Nx, Ny, 1)])*(ypdata[idx(i-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(i+1, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(i-1, j, Nx, Ny, 1)]-1*ypdata[idx(i+1, j, Nx, Ny, 1)])*(ypdata[idx(i-1, j, Nx, Ny, 1)]-1*ypdata[idx(i+1, j, Nx, Ny, 1)]);
            
            IS2_px[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)]);
            
            IS0_nx[0]=(13/12)*(yndata[idx(i+1, j, Nx, Ny, 1)]-2*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)])*(yndata[idx(i+1, j, Nx, Ny, 1)]-2*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)])+(1/4)*(3*yndata[idx(i+1, j, Nx, Ny, 1)]-4*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)])*(3*yndata[idx(i+1, j, Nx, Ny, 1)]-4*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)]);
            
            IS1_nx[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(i+1, j, Nx, Ny, 1)]+yndata[idx(i+2, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(i+1, j, Nx, Ny, 1)]+yndata[idx(i+2, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(i+2, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(i+2, j, Nx, Ny, 1)]);
            
            IS2_nx[0]=(13/12)*(yndata[idx(i-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(i+1, j, Nx, Ny, 1)])*(yndata[idx(i-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(i+1, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(i+1, j, Nx, Ny, 1)])*(yndata[idx(i-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(i+1, j, Nx, Ny, 1)]);
            
            
            IS0_px[1]=(13/12)*(taopdata[idx(i-2, j, Nx, Ny, 0)]-2*taopdata[idx(i-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(i-2, j, Nx, Ny, 0)]-2*taopdata[idx(i-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(i-2, j, Nx, Ny, 0)]-4*taopdata[idx(i-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(i-2, j, Nx, Ny, 0)]-4*taopdata[idx(i-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)]);
            
            IS1_px[1]=(13/12)*(taopdata[idx(i-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(i+1, j, Nx, Ny, 0)])*(taopdata[idx(i-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(i-1, j, Nx, Ny, 0)]-1*taopdata[idx(i+1, j, Nx, Ny, 0)])*(taopdata[idx(i-1, j, Nx, Ny, 0)]-1*taopdata[idx(i+1, j, Nx, Ny, 0)]);
            
            IS2_px[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS0_nx[1]=(13/12)*(taondata[idx(i+1, j, Nx, Ny, 0)]-2*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)])*(taondata[idx(i+1, j, Nx, Ny, 0)]-2*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)])+(1/4)*(3*taondata[idx(i+1, j, Nx, Ny, 0)]-4*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)])*(3*taondata[idx(i+1, j, Nx, Ny, 0)]-4*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)]);
            
            IS1_nx[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(i+1, j, Nx, Ny, 0)]+taondata[idx(i+2, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(i+1, j, Nx, Ny, 0)]+taondata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(i+2, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS2_nx[1]=(13/12)*(taondata[idx(i-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(i+1, j, Nx, Ny, 0)])*(taondata[idx(i-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(i+1, j, Nx, Ny, 0)])*(taondata[idx(i-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(i+1, j, Nx, Ny, 0)]);
            
            
            IS0_px[2]=(13/12)*(taopdata[idx(i-2, j, Nx, Ny, 2)]-2*taopdata[idx(i-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(i-2, j, Nx, Ny, 2)]-2*taopdata[idx(i-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(i-2, j, Nx, Ny, 2)]-4*taopdata[idx(i-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(i-2, j, Nx, Ny, 2)]-4*taopdata[idx(i-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_px[2]=(13/12)*(taopdata[idx(i-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(i+1, j, Nx, Ny, 2)])*(taopdata[idx(i-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(i+1, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(i-1, j, Nx, Ny, 2)]-1*taopdata[idx(i+1, j, Nx, Ny, 2)])*(taopdata[idx(i-1, j, Nx, Ny, 2)]-1*taopdata[idx(i+1, j, Nx, Ny, 2)]);
            
            IS2_px[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)]);
            
            IS0_nx[2]=(13/12)*(taondata[idx(i+1, j, Nx, Ny, 2)]-2*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)])*(taondata[idx(i+1, j, Nx, Ny, 2)]-2*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)])+(1/4)*(3*taondata[idx(i+1, j, Nx, Ny, 2)]-4*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)])*(3*taondata[idx(i+1, j, Nx, Ny, 2)]-4*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)]);
            
            IS1_nx[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(i+1, j, Nx, Ny, 2)]+taondata[idx(i+2, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(i+1, j, Nx, Ny, 2)]+taondata[idx(i+2, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(i+2, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(i+2, j, Nx, Ny, 2)]);
            
            IS2_nx[2]=(13/12)*(taondata[idx(i-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(i+1, j, Nx, Ny, 2)])*(taondata[idx(i-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(i+1, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(i+1, j, Nx, Ny, 2)])*(taondata[idx(i-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(i+1, j, Nx, Ny, 2)]);
            
            
            IS0_px[3]=(13/12)*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-2*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-2*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-4*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-4*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)]);
            
            IS1_px[3]=(13/12)*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(i+1, j, Nx, Ny, 0)])*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(i+1, j, Nx, Ny, 0)])*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(i+1, j, Nx, Ny, 0)]);
            
            IS2_px[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS0_nx[3]=(13/12)*(Cjndata[idx(i+1, j, Nx, Ny, 0)]-2*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)])*(Cjndata[idx(i+1, j, Nx, Ny, 0)]-2*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)])+(1/4)*(3*Cjndata[idx(i+1, j, Nx, Ny, 0)]-4*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)])*(3*Cjndata[idx(i+1, j, Nx, Ny, 0)]-4*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)]);
            
            IS1_nx[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(i+1, j, Nx, Ny, 0)]+Cjndata[idx(i+2, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(i+1, j, Nx, Ny, 0)]+Cjndata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(i+2, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS2_nx[3]=(13/12)*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(i+1, j, Nx, Ny, 0)])*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(i+1, j, Nx, Ny, 0)])*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(i+1, j, Nx, Ny, 0)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0px[k]=(1/10)*(1/(Epsilon+IS0_px[k]))*(1/(Epsilon+IS0_px[k]));
                alpha_1px[k]=(6/10)*(1/(Epsilon+IS1_px[k]))*(1/(Epsilon+IS1_px[k]));
                alpha_2px[k]=(3/10)*(1/(Epsilon+IS2_px[k]))*(1/(Epsilon+IS2_px[k]));
                alpha_0nx[k]=(1/10)*(1/(Epsilon+IS0_nx[k]))*(1/(Epsilon+IS0_nx[k]));
                alpha_1nx[k]=(6/10)*(1/(Epsilon+IS1_nx[k]))*(1/(Epsilon+IS1_nx[k]));
                alpha_2nx[k]=(3/10)*(1/(Epsilon+IS2_nx[k]))*(1/(Epsilon+IS2_nx[k]));
                
                w0_px[k]=alpha_0px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w1_px[k]=alpha_1px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w2_px[k]=alpha_2px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w0_nx[k]=alpha_0nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
                w1_nx[k]=alpha_1nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
                w2_nx[k]=alpha_2nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
            }
            
            u_tpphx[0]=w0_px[0]*((2/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i+1, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i+1, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(i+2, j, Nx, Ny, 1)]);
            
            u_tnphx[0]=w2_nx[0]*((-1/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+1, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i+1, j, Nx, Ny, 1)]-(1/6)*yndata[idx(i+2, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i+1, j, Nx, Ny, 1)]-(7/6)*yndata[idx(i+2, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+3, j, Nx, Ny, 1)]);
            
            u_tpnhx[0]=w0_px[0]*((2/6)*ypdata[idx(Nx-1, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i-1, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(i+1, j, Nx, Ny, 1)]);
            
            u_tnnhx[0]=w2_nx[0]*((-1/6)*yndata[idx(i-2, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]-(1/6)*yndata[idx(i+1, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 1)]-(7/6)*yndata[idx(i+1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+2, j, Nx, Ny, 1)]);
            
            
            u_tpphx[1]=w0_px[1]*((2/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i+1, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(i+2, j, Nx, Ny, 0)]);
            
            u_tnphx[1]=w2_nx[1]*((-1/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+1, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*taondata[idx(i+2, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i+1, j, Nx, Ny, 0)]-(7/6)*taondata[idx(i+2, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+3, j, Nx, Ny, 0)]);
            
            u_tpnhx[1]=w0_px[1]*((2/6)*taopdata[idx(Nx-1, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i-1, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(i+1, j, Nx, Ny, 0)]);
            
            u_tnnhx[1]=w2_nx[1]*((-1/6)*taondata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]-(1/6)*taondata[idx(i+1, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 0)]-(7/6)*taondata[idx(i+1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+2, j, Nx, Ny, 0)]);
            
            
            u_tpphx[2]=w0_px[2]*((2/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i+1, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i+1, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(i+2, j, Nx, Ny, 2)]);
            
            u_tnphx[2]=w2_nx[2]*((-1/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+1, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i+1, j, Nx, Ny, 2)]-(1/6)*taondata[idx(i+2, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i+1, j, Nx, Ny, 2)]-(7/6)*taondata[idx(i+2, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+3, j, Nx, Ny, 2)]);
            
            u_tpnhx[2]=w0_px[2]*((2/6)*taopdata[idx(Nx-1, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i-1, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(i+1, j, Nx, Ny, 2)]);
            
            u_tnnhx[2]=w2_nx[2]*((-1/6)*taondata[idx(i-2, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]-(1/6)*taondata[idx(i+1, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 2)]-(7/6)*taondata[idx(i+1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+2, j, Nx, Ny, 2)]);
            
            
            u_tpphx[3]=w0_px[3]*((2/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(i+2, j, Nx, Ny, 0)]);
            
            u_tnphx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+3, j, Nx, Ny, 0)]);
            
            u_tpnhx[3]=w0_px[3]*((2/6)*Cjpdata[idx(Nx-1, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)]);
            
            u_tnnhx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)]);
            
            for(k=0;k<4;k++){
                yxdata[idx(i, j, Nx, Ny, k)]=-(1/dx)*((u_tpphx[k]-u_tpnhx[k])+(u_tnphx[k]-u_tnnhx[k]));
            }
    }
    
            j=2;
        for (i=0; i<Nx; i++){
            // x direction
            IS0_py[0]=(13/12)*(ypdata[idx(i, j-2, Nx, Ny, 2)]-2*ypdata[idx(i, j-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, j-2, Nx, Ny, 2)]-2*ypdata[idx(i, j-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, j-2, Nx, Ny, 2)]-4*ypdata[idx(i, j-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, j-2, Nx, Ny, 2)]-4*ypdata[idx(i, j-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_py[0]=(13/12)*(ypdata[idx(i, j-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, j+1, Nx, Ny, 2)])*(ypdata[idx(i, j-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, j+1, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, j-1, Nx, Ny, 2)]-1*ypdata[idx(i, j+1, Nx, Ny, 2)])*(ypdata[idx(i, j-1, Nx, Ny, 2)]-1*ypdata[idx(i, j+1, Nx, Ny, 2)]);
            
            IS2_py[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)]);
            
            IS0_ny[0]=(13/12)*(yndata[idx(i, j+1, Nx, Ny, 2)]-2*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)])*(yndata[idx(i, j+1, Nx, Ny, 2)]-2*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)])+(1/4)*(3*yndata[idx(i, j+1, Nx, Ny, 2)]-4*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)])*(3*yndata[idx(i, j+1, Nx, Ny, 2)]-4*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)]);
            
            IS1_ny[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, j+1, Nx, Ny, 2)]+yndata[idx(i, j+2, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, j+1, Nx, Ny, 2)]+yndata[idx(i, j+2, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, j+2, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, j+2, Nx, Ny, 2)]);
            
            IS2_ny[0]=(13/12)*(yndata[idx(i, j-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, j+1, Nx, Ny, 2)])*(yndata[idx(i, j-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, j+1, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, j+1, Nx, Ny, 2)])*(yndata[idx(i, j-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, j+1, Nx, Ny, 2)]);
            
            
            IS0_py[1]=(13/12)*(taopdata[idx(i, j-2, Nx, Ny, 1)]-2*taopdata[idx(i, j-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, j-2, Nx, Ny, 1)]-2*taopdata[idx(i, j-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, j-2, Nx, Ny, 1)]-4*taopdata[idx(i, j-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, j-2, Nx, Ny, 1)]-4*taopdata[idx(i, j-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[1]=(13/12)*(taopdata[idx(i, j-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, j+1, Nx, Ny, 1)])*(taopdata[idx(i, j-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, j-1, Nx, Ny, 1)]-1*taopdata[idx(i, j+1, Nx, Ny, 1)])*(taopdata[idx(i, j-1, Nx, Ny, 1)]-1*taopdata[idx(i, j+1, Nx, Ny, 1)]);
            
            IS2_py[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS0_ny[1]=(13/12)*(taondata[idx(i, j+1, Nx, Ny, 1)]-2*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)])*(taondata[idx(i, j+1, Nx, Ny, 1)]-2*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)])+(1/4)*(3*taondata[idx(i, j+1, Nx, Ny, 1)]-4*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)])*(3*taondata[idx(i, j+1, Nx, Ny, 1)]-4*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)]);
            
            IS1_ny[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, j+1, Nx, Ny, 1)]+taondata[idx(i, j+2, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, j+1, Nx, Ny, 1)]+taondata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, j+2, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS2_ny[1]=(13/12)*(taondata[idx(i, j-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, j+1, Nx, Ny, 1)])*(taondata[idx(i, j-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, j+1, Nx, Ny, 1)])*(taondata[idx(i, j-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, j+1, Nx, Ny, 1)]);
            
            
            IS0_py[2]=(13/12)*(taopdata[idx(i, j-2, Nx, Ny, 3)]-2*taopdata[idx(i, j-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, j-2, Nx, Ny, 3)]-2*taopdata[idx(i, j-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, j-2, Nx, Ny, 3)]-4*taopdata[idx(i, j-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, j-2, Nx, Ny, 3)]-4*taopdata[idx(i, j-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)]);
            
            IS1_py[2]=(13/12)*(taopdata[idx(i, j-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, j+1, Nx, Ny, 3)])*(taopdata[idx(i, j-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, j+1, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, j-1, Nx, Ny, 3)]-1*taopdata[idx(i, j+1, Nx, Ny, 3)])*(taopdata[idx(i, j-1, Nx, Ny, 3)]-1*taopdata[idx(i, j+1, Nx, Ny, 3)]);
            
            IS2_py[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)]);
            
            IS0_ny[2]=(13/12)*(taondata[idx(i, j+1, Nx, Ny, 3)]-2*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)])*(taondata[idx(i, j+1, Nx, Ny, 3)]-2*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)])+(1/4)*(3*taondata[idx(i, j+1, Nx, Ny, 3)]-4*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)])*(3*taondata[idx(i, j+1, Nx, Ny, 3)]-4*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)]);
            
            IS1_ny[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, j+1, Nx, Ny, 3)]+taondata[idx(i, j+2, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, j+1, Nx, Ny, 3)]+taondata[idx(i, j+2, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, j+2, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, j+2, Nx, Ny, 3)]);
            
            IS2_ny[2]=(13/12)*(taondata[idx(i, j-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, j+1, Nx, Ny, 3)])*(taondata[idx(i, j-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, j+1, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, j+1, Nx, Ny, 3)])*(taondata[idx(i, j-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, j+1, Nx, Ny, 3)]);
            
            
            IS0_py[3]=(13/12)*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-2*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-2*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-4*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-4*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[3]=(13/12)*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, j+1, Nx, Ny, 1)])*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, j+1, Nx, Ny, 1)])*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, j+1, Nx, Ny, 1)]);
            
            IS2_py[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS0_ny[3]=(13/12)*(Cjndata[idx(i, j+1, Nx, Ny, 1)]-2*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)])*(Cjndata[idx(i, j+1, Nx, Ny, 1)]-2*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)])+(1/4)*(3*Cjndata[idx(i, j+1, Nx, Ny, 1)]-4*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)])*(3*Cjndata[idx(i, j+1, Nx, Ny, 1)]-4*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)]);
            
            IS1_ny[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, j+1, Nx, Ny, 1)]+Cjndata[idx(i, j+2, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, j+1, Nx, Ny, 1)]+Cjndata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, j+2, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS2_ny[3]=(13/12)*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, j+1, Nx, Ny, 1)])*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, j+1, Nx, Ny, 1)])*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, j+1, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0py[k]=(1/10)*(1/(Epsilon+IS0_py[k]))*(1/(Epsilon+IS0_py[k]));
                alpha_1py[k]=(6/10)*(1/(Epsilon+IS1_py[k]))*(1/(Epsilon+IS1_py[k]));
                alpha_2py[k]=(3/10)*(1/(Epsilon+IS2_py[k]))*(1/(Epsilon+IS2_py[k]));
                alpha_0ny[k]=(1/10)*(1/(Epsilon+IS0_ny[k]))*(1/(Epsilon+IS0_ny[k]));
                alpha_1ny[k]=(6/10)*(1/(Epsilon+IS1_ny[k]))*(1/(Epsilon+IS1_ny[k]));
                alpha_2ny[k]=(3/10)*(1/(Epsilon+IS2_ny[k]))*(1/(Epsilon+IS2_ny[k]));
                
                w0_py[k]=alpha_0py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w1_py[k]=alpha_1py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w2_py[k]=alpha_2py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w0_ny[k]=alpha_0ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w1_ny[k]=alpha_1ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w2_ny[k]=alpha_2ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
            }
            
            u_tpphy[0]=w0_py[0]*((2/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j+1, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j+1, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, j+2, Nx, Ny, 2)]);
            
            u_tnphy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+1, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j+1, Nx, Ny, 2)]-(1/6)*yndata[idx(i, j+2, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j+1, Nx, Ny, 2)]-(7/6)*yndata[idx(i, j+2, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+3, Nx, Ny, 2)]);
            
            u_tpnhy[0]=w0_py[0]*((2/6)*ypdata[idx(i, Ny-1, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j-1, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, j+1, Nx, Ny, 2)]);
            
            u_tnnhy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, j-2, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]-(1/6)*yndata[idx(i, j+1, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 2)]-(7/6)*yndata[idx(i, j+1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+2, Nx, Ny, 2)]);
            
            
            u_tpphy[1]=w0_py[1]*((2/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j+1, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, j+2, Nx, Ny, 1)]);
            
            u_tnphy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+1, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*taondata[idx(i, j+2, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j+1, Nx, Ny, 1)]-(7/6)*taondata[idx(i, j+2, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+3, Nx, Ny, 1)]);
            
            u_tpnhy[1]=w0_py[1]*((2/6)*taopdata[idx(i, Ny-1, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j-1, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, j+1, Nx, Ny, 1)]);
            
            u_tnnhy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]-(1/6)*taondata[idx(i, j+1, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 1)]-(7/6)*taondata[idx(i, j+1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+2, Nx, Ny, 1)]);
            
            
            u_tpphy[2]=w0_py[2]*((2/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j+1, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j+1, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, j+2, Nx, Ny, 3)]);
            
            u_tnphy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+1, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j+1, Nx, Ny, 3)]-(1/6)*taondata[idx(i, j+2, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j+1, Nx, Ny, 3)]-(7/6)*taondata[idx(i, j+2, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+3, Nx, Ny, 3)]);
            
            u_tpnhy[2]=w0_py[2]*((2/6)*taopdata[idx(i, Ny-1, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j-1, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, j+1, Nx, Ny, 3)]);
            
            u_tnnhy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, j-2, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]-(1/6)*taondata[idx(i, j+1, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 3)]-(7/6)*taondata[idx(i, j+1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+2, Nx, Ny, 3)]);
            
            
            u_tpphy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, j+2, Nx, Ny, 1)]);
            
            u_tnphy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+3, Nx, Ny, 1)]);
            
            u_tpnhy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, Ny-1, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)]);
            
            u_tnnhy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++){
                yydata[idx(i, j, Nx, Ny, k)]=-(1/dy)*((u_tpphy[k]-u_tpnhy[k])+(u_tnphy[k]-u_tnnhy[k]));
            }
        }

    
    
    for(j=0; j<Ny; j++){
        for (i=3; i<Nx-3; i++){
        // x direction
        IS0_px[0]=(13/12)*(ypdata[idx(i-2, j, Nx, Ny, 1)]-2*ypdata[idx(i-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(i-2, j, Nx, Ny, 1)]-2*ypdata[idx(i-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(i-2, j, Nx, Ny, 1)]-4*ypdata[idx(i-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(i-2, j, Nx, Ny, 1)]-4*ypdata[idx(i-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)]);
        
        IS1_px[0]=(13/12)*(ypdata[idx(i-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(i+1, j, Nx, Ny, 1)])*(ypdata[idx(i-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(i+1, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(i-1, j, Nx, Ny, 1)]-1*ypdata[idx(i+1, j, Nx, Ny, 1)])*(ypdata[idx(i-1, j, Nx, Ny, 1)]-1*ypdata[idx(i+1, j, Nx, Ny, 1)]);
        
        IS2_px[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)]);
        
        IS0_nx[0]=(13/12)*(yndata[idx(i+1, j, Nx, Ny, 1)]-2*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)])*(yndata[idx(i+1, j, Nx, Ny, 1)]-2*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)])+(1/4)*(3*yndata[idx(i+1, j, Nx, Ny, 1)]-4*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)])*(3*yndata[idx(i+1, j, Nx, Ny, 1)]-4*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(i+3, j, Nx, Ny, 1)]);
        
        IS1_nx[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(i+1, j, Nx, Ny, 1)]+yndata[idx(i+2, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(i+1, j, Nx, Ny, 1)]+yndata[idx(i+2, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(i+2, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(i+2, j, Nx, Ny, 1)]);
        
        IS2_nx[0]=(13/12)*(yndata[idx(i-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(i+1, j, Nx, Ny, 1)])*(yndata[idx(i-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(i+1, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(i+1, j, Nx, Ny, 1)])*(yndata[idx(i-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(i+1, j, Nx, Ny, 1)]);
        
                
        IS0_px[1]=(13/12)*(taopdata[idx(i-2, j, Nx, Ny, 0)]-2*taopdata[idx(i-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(i-2, j, Nx, Ny, 0)]-2*taopdata[idx(i-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(i-2, j, Nx, Ny, 0)]-4*taopdata[idx(i-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(i-2, j, Nx, Ny, 0)]-4*taopdata[idx(i-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)]);
    
        IS1_px[1]=(13/12)*(taopdata[idx(i-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(i+1, j, Nx, Ny, 0)])*(taopdata[idx(i-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(i-1, j, Nx, Ny, 0)]-1*taopdata[idx(i+1, j, Nx, Ny, 0)])*(taopdata[idx(i-1, j, Nx, Ny, 0)]-1*taopdata[idx(i+1, j, Nx, Ny, 0)]);
                
        IS2_px[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)]);
                
        IS0_nx[1]=(13/12)*(taondata[idx(i+1, j, Nx, Ny, 0)]-2*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)])*(taondata[idx(i+1, j, Nx, Ny, 0)]-2*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)])+(1/4)*(3*taondata[idx(i+1, j, Nx, Ny, 0)]-4*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)])*(3*taondata[idx(i+1, j, Nx, Ny, 0)]-4*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(i+3, j, Nx, Ny, 0)]);
                
        IS1_nx[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(i+1, j, Nx, Ny, 0)]+taondata[idx(i+2, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(i+1, j, Nx, Ny, 0)]+taondata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(i+2, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(i+2, j, Nx, Ny, 0)]);
                
        IS2_nx[1]=(13/12)*(taondata[idx(i-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(i+1, j, Nx, Ny, 0)])*(taondata[idx(i-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(i+1, j, Nx, Ny, 0)])*(taondata[idx(i-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(i+1, j, Nx, Ny, 0)]);
                
                
        IS0_px[2]=(13/12)*(taopdata[idx(i-2, j, Nx, Ny, 2)]-2*taopdata[idx(i-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(i-2, j, Nx, Ny, 2)]-2*taopdata[idx(i-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(i-2, j, Nx, Ny, 2)]-4*taopdata[idx(i-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(i-2, j, Nx, Ny, 2)]-4*taopdata[idx(i-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)]);
                
        IS1_px[2]=(13/12)*(taopdata[idx(i-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(i+1, j, Nx, Ny, 2)])*(taopdata[idx(i-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(i+1, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(i-1, j, Nx, Ny, 2)]-1*taopdata[idx(i+1, j, Nx, Ny, 2)])*(taopdata[idx(i-1, j, Nx, Ny, 2)]-1*taopdata[idx(i+1, j, Nx, Ny, 2)]);
                
        IS2_px[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)]);
                
        IS0_nx[2]=(13/12)*(taondata[idx(i+1, j, Nx, Ny, 2)]-2*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)])*(taondata[idx(i+1, j, Nx, Ny, 2)]-2*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)])+(1/4)*(3*taondata[idx(i+1, j, Nx, Ny, 2)]-4*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)])*(3*taondata[idx(i+1, j, Nx, Ny, 2)]-4*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(i+3, j, Nx, Ny, 2)]);
                
        IS1_nx[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(i+1, j, Nx, Ny, 2)]+taondata[idx(i+2, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(i+1, j, Nx, Ny, 2)]+taondata[idx(i+2, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(i+2, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(i+2, j, Nx, Ny, 2)]);
                
        IS2_nx[2]=(13/12)*(taondata[idx(i-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(i+1, j, Nx, Ny, 2)])*(taondata[idx(i-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(i+1, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(i+1, j, Nx, Ny, 2)])*(taondata[idx(i-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(i+1, j, Nx, Ny, 2)]);
        
                
        IS0_px[3]=(13/12)*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-2*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-2*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-4*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-4*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)]);
                
        IS1_px[3]=(13/12)*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(i+1, j, Nx, Ny, 0)])*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(i+1, j, Nx, Ny, 0)])*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(i+1, j, Nx, Ny, 0)]);
                
        IS2_px[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)]);
                
        IS0_nx[3]=(13/12)*(Cjndata[idx(i+1, j, Nx, Ny, 0)]-2*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)])*(Cjndata[idx(i+1, j, Nx, Ny, 0)]-2*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)])+(1/4)*(3*Cjndata[idx(i+1, j, Nx, Ny, 0)]-4*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)])*(3*Cjndata[idx(i+1, j, Nx, Ny, 0)]-4*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(i+3, j, Nx, Ny, 0)]);
                
        IS1_nx[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(i+1, j, Nx, Ny, 0)]+Cjndata[idx(i+2, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(i+1, j, Nx, Ny, 0)]+Cjndata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(i+2, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(i+2, j, Nx, Ny, 0)]);
                
        IS2_nx[3]=(13/12)*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(i+1, j, Nx, Ny, 0)])*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(i+1, j, Nx, Ny, 0)])*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(i+1, j, Nx, Ny, 0)]);
                
    for(k=0;k<4;k++)
    {
        Epsilon=0.000001;
        alpha_0px[k]=(1/10)*(1/(Epsilon+IS0_px[k]))*(1/(Epsilon+IS0_px[k]));
        alpha_1px[k]=(6/10)*(1/(Epsilon+IS1_px[k]))*(1/(Epsilon+IS1_px[k]));
        alpha_2px[k]=(3/10)*(1/(Epsilon+IS2_px[k]))*(1/(Epsilon+IS2_px[k]));
        alpha_0nx[k]=(1/10)*(1/(Epsilon+IS0_nx[k]))*(1/(Epsilon+IS0_nx[k]));
        alpha_1nx[k]=(6/10)*(1/(Epsilon+IS1_nx[k]))*(1/(Epsilon+IS1_nx[k]));
        alpha_2nx[k]=(3/10)*(1/(Epsilon+IS2_nx[k]))*(1/(Epsilon+IS2_nx[k]));
        
        w0_px[k]=alpha_0px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
        w1_px[k]=alpha_1px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
        w2_px[k]=alpha_2px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
        w0_nx[k]=alpha_0nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
        w1_nx[k]=alpha_1nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
        w2_nx[k]=alpha_2nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
    }
                
        u_tpphx[0]=w0_px[0]*((2/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i+1, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i+1, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(i+2, j, Nx, Ny, 1)]);
        
        u_tnphx[0]=w2_nx[0]*((-1/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+1, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i+1, j, Nx, Ny, 1)]-(1/6)*yndata[idx(i+2, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i+1, j, Nx, Ny, 1)]-(7/6)*yndata[idx(i+2, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+3, j, Nx, Ny, 1)]);
        
        u_tpnhx[0]=w0_px[0]*((2/6)*ypdata[idx(i-3, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i-1, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(i+1, j, Nx, Ny, 1)]);
        
        u_tnnhx[0]=w2_nx[0]*((-1/6)*yndata[idx(i-2, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]-(1/6)*yndata[idx(i+1, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 1)]-(7/6)*yndata[idx(i+1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+2, j, Nx, Ny, 1)]);
        
                
        u_tpphx[1]=w0_px[1]*((2/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i+1, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(i+2, j, Nx, Ny, 0)]);
                
        u_tnphx[1]=w2_nx[1]*((-1/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+1, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*taondata[idx(i+2, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i+1, j, Nx, Ny, 0)]-(7/6)*taondata[idx(i+2, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+3, j, Nx, Ny, 0)]);
                
        u_tpnhx[1]=w0_px[1]*((2/6)*taopdata[idx(i-3, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i-1, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(i+1, j, Nx, Ny, 0)]);
                
        u_tnnhx[1]=w2_nx[1]*((-1/6)*taondata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]-(1/6)*taondata[idx(i+1, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 0)]-(7/6)*taondata[idx(i+1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+2, j, Nx, Ny, 0)]);
                
                
        u_tpphx[2]=w0_px[2]*((2/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i+1, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i+1, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(i+2, j, Nx, Ny, 2)]);
                
        u_tnphx[2]=w2_nx[2]*((-1/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+1, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i+1, j, Nx, Ny, 2)]-(1/6)*taondata[idx(i+2, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i+1, j, Nx, Ny, 2)]-(7/6)*taondata[idx(i+2, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+3, j, Nx, Ny, 2)]);
                
        u_tpnhx[2]=w0_px[2]*((2/6)*taopdata[idx(i-3, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i-1, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(i+1, j, Nx, Ny, 2)]);
                
        u_tnnhx[2]=w2_nx[2]*((-1/6)*taondata[idx(i-2, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]-(1/6)*taondata[idx(i+1, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 2)]-(7/6)*taondata[idx(i+1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+2, j, Nx, Ny, 2)]);
               
                
        u_tpphx[3]=w0_px[3]*((2/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(i+2, j, Nx, Ny, 0)]);
                
        u_tnphx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+3, j, Nx, Ny, 0)]);
                
        u_tpnhx[3]=w0_px[3]*((2/6)*Cjpdata[idx(i-3, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)]);
                
        u_tnnhx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)]);
                
                for(k=0;k<4;k++){
                yxdata[idx(i, j, Nx, Ny, k)]=-(1/dx)*((u_tpphx[k]-u_tpnhx[k])+(u_tnphx[k]-u_tnnhx[k]));
                }
        }
    }
    
    for(j=3; j<Ny-3; j++){
        for (i=0; i<Nx; i++){
            // x direction
            IS0_py[0]=(13/12)*(ypdata[idx(i, j-2, Nx, Ny, 2)]-2*ypdata[idx(i, j-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, j-2, Nx, Ny, 2)]-2*ypdata[idx(i, j-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, j-2, Nx, Ny, 2)]-4*ypdata[idx(i, j-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, j-2, Nx, Ny, 2)]-4*ypdata[idx(i, j-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_py[0]=(13/12)*(ypdata[idx(i, j-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, j+1, Nx, Ny, 2)])*(ypdata[idx(i, j-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, j+1, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, j-1, Nx, Ny, 2)]-1*ypdata[idx(i, j+1, Nx, Ny, 2)])*(ypdata[idx(i, j-1, Nx, Ny, 2)]-1*ypdata[idx(i, j+1, Nx, Ny, 2)]);
            
            IS2_py[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)]);
            
            IS0_ny[0]=(13/12)*(yndata[idx(i, j+1, Nx, Ny, 2)]-2*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)])*(yndata[idx(i, j+1, Nx, Ny, 2)]-2*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)])+(1/4)*(3*yndata[idx(i, j+1, Nx, Ny, 2)]-4*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)])*(3*yndata[idx(i, j+1, Nx, Ny, 2)]-4*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, j+3, Nx, Ny, 2)]);
            
            IS1_ny[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, j+1, Nx, Ny, 2)]+yndata[idx(i, j+2, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, j+1, Nx, Ny, 2)]+yndata[idx(i, j+2, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, j+2, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, j+2, Nx, Ny, 2)]);
            
            IS2_ny[0]=(13/12)*(yndata[idx(i, j-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, j+1, Nx, Ny, 2)])*(yndata[idx(i, j-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, j+1, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, j+1, Nx, Ny, 2)])*(yndata[idx(i, j-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, j+1, Nx, Ny, 2)]);
            
            
            IS0_py[1]=(13/12)*(taopdata[idx(i, j-2, Nx, Ny, 1)]-2*taopdata[idx(i, j-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, j-2, Nx, Ny, 1)]-2*taopdata[idx(i, j-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, j-2, Nx, Ny, 1)]-4*taopdata[idx(i, j-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, j-2, Nx, Ny, 1)]-4*taopdata[idx(i, j-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[1]=(13/12)*(taopdata[idx(i, j-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, j+1, Nx, Ny, 1)])*(taopdata[idx(i, j-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, j-1, Nx, Ny, 1)]-1*taopdata[idx(i, j+1, Nx, Ny, 1)])*(taopdata[idx(i, j-1, Nx, Ny, 1)]-1*taopdata[idx(i, j+1, Nx, Ny, 1)]);
            
            IS2_py[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS0_ny[1]=(13/12)*(taondata[idx(i, j+1, Nx, Ny, 1)]-2*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)])*(taondata[idx(i, j+1, Nx, Ny, 1)]-2*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)])+(1/4)*(3*taondata[idx(i, j+1, Nx, Ny, 1)]-4*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)])*(3*taondata[idx(i, j+1, Nx, Ny, 1)]-4*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, j+3, Nx, Ny, 1)]);
            
            IS1_ny[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, j+1, Nx, Ny, 1)]+taondata[idx(i, j+2, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, j+1, Nx, Ny, 1)]+taondata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, j+2, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS2_ny[1]=(13/12)*(taondata[idx(i, j-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, j+1, Nx, Ny, 1)])*(taondata[idx(i, j-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, j+1, Nx, Ny, 1)])*(taondata[idx(i, j-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, j+1, Nx, Ny, 1)]);
            
            
            IS0_py[2]=(13/12)*(taopdata[idx(i, j-2, Nx, Ny, 3)]-2*taopdata[idx(i, j-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, j-2, Nx, Ny, 3)]-2*taopdata[idx(i, j-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, j-2, Nx, Ny, 3)]-4*taopdata[idx(i, j-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, j-2, Nx, Ny, 3)]-4*taopdata[idx(i, j-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)]);
            
            IS1_py[2]=(13/12)*(taopdata[idx(i, j-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, j+1, Nx, Ny, 3)])*(taopdata[idx(i, j-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, j+1, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, j-1, Nx, Ny, 3)]-1*taopdata[idx(i, j+1, Nx, Ny, 3)])*(taopdata[idx(i, j-1, Nx, Ny, 3)]-1*taopdata[idx(i, j+1, Nx, Ny, 3)]);
            
            IS2_py[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)]);
            
            IS0_ny[2]=(13/12)*(taondata[idx(i, j+1, Nx, Ny, 3)]-2*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)])*(taondata[idx(i, j+1, Nx, Ny, 3)]-2*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)])+(1/4)*(3*taondata[idx(i, j+1, Nx, Ny, 3)]-4*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)])*(3*taondata[idx(i, j+1, Nx, Ny, 3)]-4*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, j+3, Nx, Ny, 3)]);
            
            IS1_ny[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, j+1, Nx, Ny, 3)]+taondata[idx(i, j+2, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, j+1, Nx, Ny, 3)]+taondata[idx(i, j+2, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, j+2, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, j+2, Nx, Ny, 3)]);
            
            IS2_ny[2]=(13/12)*(taondata[idx(i, j-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, j+1, Nx, Ny, 3)])*(taondata[idx(i, j-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, j+1, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, j+1, Nx, Ny, 3)])*(taondata[idx(i, j-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, j+1, Nx, Ny, 3)]);
            
            
            IS0_py[3]=(13/12)*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-2*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-2*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-4*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-4*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[3]=(13/12)*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, j+1, Nx, Ny, 1)])*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, j+1, Nx, Ny, 1)])*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, j+1, Nx, Ny, 1)]);
            
            IS2_py[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS0_ny[3]=(13/12)*(Cjndata[idx(i, j+1, Nx, Ny, 1)]-2*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)])*(Cjndata[idx(i, j+1, Nx, Ny, 1)]-2*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)])+(1/4)*(3*Cjndata[idx(i, j+1, Nx, Ny, 1)]-4*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)])*(3*Cjndata[idx(i, j+1, Nx, Ny, 1)]-4*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, j+3, Nx, Ny, 1)]);
            
            IS1_ny[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, j+1, Nx, Ny, 1)]+Cjndata[idx(i, j+2, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, j+1, Nx, Ny, 1)]+Cjndata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, j+2, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS2_ny[3]=(13/12)*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, j+1, Nx, Ny, 1)])*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, j+1, Nx, Ny, 1)])*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, j+1, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0py[k]=(1/10)*(1/(Epsilon+IS0_py[k]))*(1/(Epsilon+IS0_py[k]));
                alpha_1py[k]=(6/10)*(1/(Epsilon+IS1_py[k]))*(1/(Epsilon+IS1_py[k]));
                alpha_2py[k]=(3/10)*(1/(Epsilon+IS2_py[k]))*(1/(Epsilon+IS2_py[k]));
                alpha_0ny[k]=(1/10)*(1/(Epsilon+IS0_ny[k]))*(1/(Epsilon+IS0_ny[k]));
                alpha_1ny[k]=(6/10)*(1/(Epsilon+IS1_ny[k]))*(1/(Epsilon+IS1_ny[k]));
                alpha_2ny[k]=(3/10)*(1/(Epsilon+IS2_ny[k]))*(1/(Epsilon+IS2_ny[k]));
                
                w0_py[k]=alpha_0py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w1_py[k]=alpha_1py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w2_py[k]=alpha_2py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w0_ny[k]=alpha_0ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w1_ny[k]=alpha_1ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w2_ny[k]=alpha_2ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
            }
            
            u_tpphy[0]=w0_py[0]*((2/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j+1, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j+1, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, j+2, Nx, Ny, 2)]);
            
            u_tnphy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+1, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j+1, Nx, Ny, 2)]-(1/6)*yndata[idx(i, j+2, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j+1, Nx, Ny, 2)]-(7/6)*yndata[idx(i, j+2, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+3, Nx, Ny, 2)]);
            
            u_tpnhy[0]=w0_py[0]*((2/6)*ypdata[idx(i, j-3, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j-1, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, j+1, Nx, Ny, 2)]);
            
            u_tnnhy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, j-2, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]-(1/6)*yndata[idx(i, j+1, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 2)]-(7/6)*yndata[idx(i, j+1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+2, Nx, Ny, 2)]);
            
            
            u_tpphy[1]=w0_py[1]*((2/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j+1, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, j+2, Nx, Ny, 1)]);
            
            u_tnphy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+1, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*taondata[idx(i, j+2, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j+1, Nx, Ny, 1)]-(7/6)*taondata[idx(i, j+2, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+3, Nx, Ny, 1)]);
            
            u_tpnhy[1]=w0_py[1]*((2/6)*taopdata[idx(i, j-3, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j-1, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, j+1, Nx, Ny, 1)]);
            
            u_tnnhy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]-(1/6)*taondata[idx(i, j+1, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 1)]-(7/6)*taondata[idx(i, j+1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+2, Nx, Ny, 1)]);
            
            
            u_tpphy[2]=w0_py[2]*((2/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j+1, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j+1, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, j+2, Nx, Ny, 3)]);
            
            u_tnphy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+1, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j+1, Nx, Ny, 3)]-(1/6)*taondata[idx(i, j+2, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j+1, Nx, Ny, 3)]-(7/6)*taondata[idx(i, j+2, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+3, Nx, Ny, 3)]);
            
            u_tpnhy[2]=w0_py[2]*((2/6)*taopdata[idx(i, j-3, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j-1, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, j+1, Nx, Ny, 3)]);
            
            u_tnnhy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, j-2, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]-(1/6)*taondata[idx(i, j+1, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 3)]-(7/6)*taondata[idx(i, j+1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+2, Nx, Ny, 3)]);
            
            
            u_tpphy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, j+2, Nx, Ny, 1)]);
            
            u_tnphy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+3, Nx, Ny, 1)]);
            
            u_tpnhy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, j-3, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)]);
            
            u_tnnhy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++){
                yydata[idx(i, j, Nx, Ny, k)]=-(1/dy)*((u_tpphy[k]-u_tpnhy[k])+(u_tnphy[k]-u_tnnhy[k]));
            }
        }
    }
    
            i=Nx-3;
    for(j=0; j<Ny; j++){
            // x direction
            IS0_px[0]=(13/12)*(ypdata[idx(i-2, j, Nx, Ny, 1)]-2*ypdata[idx(i-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(i-2, j, Nx, Ny, 1)]-2*ypdata[idx(i-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(i-2, j, Nx, Ny, 1)]-4*ypdata[idx(i-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(i-2, j, Nx, Ny, 1)]-4*ypdata[idx(i-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_px[0]=(13/12)*(ypdata[idx(i-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(i+1, j, Nx, Ny, 1)])*(ypdata[idx(i-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(i+1, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(i-1, j, Nx, Ny, 1)]-1*ypdata[idx(i+1, j, Nx, Ny, 1)])*(ypdata[idx(i-1, j, Nx, Ny, 1)]-1*ypdata[idx(i+1, j, Nx, Ny, 1)]);
            
            IS2_px[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)])*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(i+2, j, Nx, Ny, 1)]);
            
            IS0_nx[0]=(13/12)*(yndata[idx(i+1, j, Nx, Ny, 1)]-2*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(0, j, Nx, Ny, 1)])*(yndata[idx(i+1, j, Nx, Ny, 1)]-2*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(0, j, Nx, Ny, 1)])+(1/4)*(3*yndata[idx(i+1, j, Nx, Ny, 1)]-4*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(0, j, Nx, Ny, 1)])*(3*yndata[idx(i+1, j, Nx, Ny, 1)]-4*yndata[idx(i+2, j, Nx, Ny, 1)]+yndata[idx(0, j, Nx, Ny, 1)]);
            
            IS1_nx[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(i+1, j, Nx, Ny, 1)]+yndata[idx(i+2, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(i+1, j, Nx, Ny, 1)]+yndata[idx(i+2, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(i+2, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(i+2, j, Nx, Ny, 1)]);
            
            IS2_nx[0]=(13/12)*(yndata[idx(i-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(i+1, j, Nx, Ny, 1)])*(yndata[idx(i-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(i+1, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(i+1, j, Nx, Ny, 1)])*(yndata[idx(i-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(i+1, j, Nx, Ny, 1)]);
            
            
            IS0_px[1]=(13/12)*(taopdata[idx(i-2, j, Nx, Ny, 0)]-2*taopdata[idx(i-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(i-2, j, Nx, Ny, 0)]-2*taopdata[idx(i-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(i-2, j, Nx, Ny, 0)]-4*taopdata[idx(i-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(i-2, j, Nx, Ny, 0)]-4*taopdata[idx(i-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)]);
            
            IS1_px[1]=(13/12)*(taopdata[idx(i-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(i+1, j, Nx, Ny, 0)])*(taopdata[idx(i-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(i-1, j, Nx, Ny, 0)]-1*taopdata[idx(i+1, j, Nx, Ny, 0)])*(taopdata[idx(i-1, j, Nx, Ny, 0)]-1*taopdata[idx(i+1, j, Nx, Ny, 0)]);
            
            IS2_px[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)])*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS0_nx[1]=(13/12)*(taondata[idx(i+1, j, Nx, Ny, 0)]-2*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(0, j, Nx, Ny, 0)])*(taondata[idx(i+1, j, Nx, Ny, 0)]-2*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(0, j, Nx, Ny, 0)])+(1/4)*(3*taondata[idx(i+1, j, Nx, Ny, 0)]-4*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(0, j, Nx, Ny, 0)])*(3*taondata[idx(i+1, j, Nx, Ny, 0)]-4*taondata[idx(i+2, j, Nx, Ny, 0)]+taondata[idx(0, j, Nx, Ny, 0)]);
            
            IS1_nx[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(i+1, j, Nx, Ny, 0)]+taondata[idx(i+2, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(i+1, j, Nx, Ny, 0)]+taondata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(i+2, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS2_nx[1]=(13/12)*(taondata[idx(i-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(i+1, j, Nx, Ny, 0)])*(taondata[idx(i-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(i+1, j, Nx, Ny, 0)])*(taondata[idx(i-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(i+1, j, Nx, Ny, 0)]);
            
            
            IS0_px[2]=(13/12)*(taopdata[idx(i-2, j, Nx, Ny, 2)]-2*taopdata[idx(i-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(i-2, j, Nx, Ny, 2)]-2*taopdata[idx(i-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(i-2, j, Nx, Ny, 2)]-4*taopdata[idx(i-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(i-2, j, Nx, Ny, 2)]-4*taopdata[idx(i-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_px[2]=(13/12)*(taopdata[idx(i-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(i+1, j, Nx, Ny, 2)])*(taopdata[idx(i-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(i+1, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(i-1, j, Nx, Ny, 2)]-1*taopdata[idx(i+1, j, Nx, Ny, 2)])*(taopdata[idx(i-1, j, Nx, Ny, 2)]-1*taopdata[idx(i+1, j, Nx, Ny, 2)]);
            
            IS2_px[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)])*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(i+2, j, Nx, Ny, 2)]);
            
            IS0_nx[2]=(13/12)*(taondata[idx(i+1, j, Nx, Ny, 2)]-2*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(0, j, Nx, Ny, 2)])*(taondata[idx(i+1, j, Nx, Ny, 2)]-2*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(0, j, Nx, Ny, 2)])+(1/4)*(3*taondata[idx(i+1, j, Nx, Ny, 2)]-4*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(0, j, Nx, Ny, 2)])*(3*taondata[idx(i+1, j, Nx, Ny, 2)]-4*taondata[idx(i+2, j, Nx, Ny, 2)]+taondata[idx(0, j, Nx, Ny, 2)]);
            
            IS1_nx[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(i+1, j, Nx, Ny, 2)]+taondata[idx(i+2, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(i+1, j, Nx, Ny, 2)]+taondata[idx(i+2, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(i+2, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(i+2, j, Nx, Ny, 2)]);
            
            IS2_nx[2]=(13/12)*(taondata[idx(i-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(i+1, j, Nx, Ny, 2)])*(taondata[idx(i-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(i+1, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(i+1, j, Nx, Ny, 2)])*(taondata[idx(i-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(i+1, j, Nx, Ny, 2)]);
            
            
            IS0_px[3]=(13/12)*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-2*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-2*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-4*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-4*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)]);
            
            IS1_px[3]=(13/12)*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(i+1, j, Nx, Ny, 0)])*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(i+1, j, Nx, Ny, 0)])*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(i+1, j, Nx, Ny, 0)]);
            
            IS2_px[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)])*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS0_nx[3]=(13/12)*(Cjndata[idx(i+1, j, Nx, Ny, 0)]-2*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(0, j, Nx, Ny, 0)])*(Cjndata[idx(i+1, j, Nx, Ny, 0)]-2*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(0, j, Nx, Ny, 0)])+(1/4)*(3*Cjndata[idx(i+1, j, Nx, Ny, 0)]-4*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(0, j, Nx, Ny, 0)])*(3*Cjndata[idx(i+1, j, Nx, Ny, 0)]-4*Cjndata[idx(i+2, j, Nx, Ny, 0)]+Cjndata[idx(0, j, Nx, Ny, 0)]);
            
            IS1_nx[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(i+1, j, Nx, Ny, 0)]+Cjndata[idx(i+2, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(i+1, j, Nx, Ny, 0)]+Cjndata[idx(i+2, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(i+2, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(i+2, j, Nx, Ny, 0)]);
            
            IS2_nx[3]=(13/12)*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(i+1, j, Nx, Ny, 0)])*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(i+1, j, Nx, Ny, 0)])*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(i+1, j, Nx, Ny, 0)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0px[k]=(1/10)*(1/(Epsilon+IS0_px[k]))*(1/(Epsilon+IS0_px[k]));
                alpha_1px[k]=(6/10)*(1/(Epsilon+IS1_px[k]))*(1/(Epsilon+IS1_px[k]));
                alpha_2px[k]=(3/10)*(1/(Epsilon+IS2_px[k]))*(1/(Epsilon+IS2_px[k]));
                alpha_0nx[k]=(1/10)*(1/(Epsilon+IS0_nx[k]))*(1/(Epsilon+IS0_nx[k]));
                alpha_1nx[k]=(6/10)*(1/(Epsilon+IS1_nx[k]))*(1/(Epsilon+IS1_nx[k]));
                alpha_2nx[k]=(3/10)*(1/(Epsilon+IS2_nx[k]))*(1/(Epsilon+IS2_nx[k]));
                
                w0_px[k]=alpha_0px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w1_px[k]=alpha_1px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w2_px[k]=alpha_2px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w0_nx[k]=alpha_0nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
                w1_nx[k]=alpha_1nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
                w2_nx[k]=alpha_2nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
            }
            
            u_tpphx[0]=w0_px[0]*((2/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i+1, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i+1, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(i+2, j, Nx, Ny, 1)]);
            
            u_tnphx[0]=w2_nx[0]*((-1/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+1, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i+1, j, Nx, Ny, 1)]-(1/6)*yndata[idx(i+2, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i+1, j, Nx, Ny, 1)]-(7/6)*yndata[idx(i+2, j, Nx, Ny, 1)]+(2/6)*yndata[idx(0, j, Nx, Ny, 1)]);
            
            u_tpnhx[0]=w0_px[0]*((2/6)*ypdata[idx(i-3, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i-1, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(i+1, j, Nx, Ny, 1)]);
            
            u_tnnhx[0]=w2_nx[0]*((-1/6)*yndata[idx(i-2, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]-(1/6)*yndata[idx(i+1, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 1)]-(7/6)*yndata[idx(i+1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+2, j, Nx, Ny, 1)]);
            
            
            u_tpphx[1]=w0_px[1]*((2/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i+1, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(i+2, j, Nx, Ny, 0)]);
            
            u_tnphx[1]=w2_nx[1]*((-1/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+1, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*taondata[idx(i+2, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i+1, j, Nx, Ny, 0)]-(7/6)*taondata[idx(i+2, j, Nx, Ny, 0)]+(2/6)*taondata[idx(0, j, Nx, Ny, 0)]);
            
            u_tpnhx[1]=w0_px[1]*((2/6)*taopdata[idx(i-3, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i-1, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(i+1, j, Nx, Ny, 0)]);
            
            u_tnnhx[1]=w2_nx[1]*((-1/6)*taondata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]-(1/6)*taondata[idx(i+1, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 0)]-(7/6)*taondata[idx(i+1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+2, j, Nx, Ny, 0)]);
            
            
            u_tpphx[2]=w0_px[2]*((2/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i+1, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i+1, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(i+2, j, Nx, Ny, 2)]);
            
            u_tnphx[2]=w2_nx[2]*((-1/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+1, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i+1, j, Nx, Ny, 2)]-(1/6)*taondata[idx(i+2, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i+1, j, Nx, Ny, 2)]-(7/6)*taondata[idx(i+2, j, Nx, Ny, 2)]+(2/6)*taondata[idx(0, j, Nx, Ny, 2)]);
            
            u_tpnhx[2]=w0_px[2]*((2/6)*taopdata[idx(i-3, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i-1, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(i+1, j, Nx, Ny, 2)]);
            
            u_tnnhx[2]=w2_nx[2]*((-1/6)*taondata[idx(i-2, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]-(1/6)*taondata[idx(i+1, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 2)]-(7/6)*taondata[idx(i+1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+2, j, Nx, Ny, 2)]);
            
            
            u_tpphx[3]=w0_px[3]*((2/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(i+2, j, Nx, Ny, 0)]);
            
            u_tnphx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(0, j, Nx, Ny, 0)]);
            
            u_tpnhx[3]=w0_px[3]*((2/6)*Cjpdata[idx(i-3, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)]);
            
            u_tnnhx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+2, j, Nx, Ny, 0)]);
            
            for(k=0;k<4;k++){
                yxdata[idx(i, j, Nx, Ny, k)]=-(1/dx)*((u_tpphx[k]-u_tpnhx[k])+(u_tnphx[k]-u_tnnhx[k]));
            }
    }
    
            j=Ny-3;
        for (i=0; i<Nx; i++){
            // x direction
            IS0_py[0]=(13/12)*(ypdata[idx(i, j-2, Nx, Ny, 2)]-2*ypdata[idx(i, j-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, j-2, Nx, Ny, 2)]-2*ypdata[idx(i, j-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, j-2, Nx, Ny, 2)]-4*ypdata[idx(i, j-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, j-2, Nx, Ny, 2)]-4*ypdata[idx(i, j-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_py[0]=(13/12)*(ypdata[idx(i, j-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, j+1, Nx, Ny, 2)])*(ypdata[idx(i, j-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, j+1, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, j-1, Nx, Ny, 2)]-1*ypdata[idx(i, j+1, Nx, Ny, 2)])*(ypdata[idx(i, j-1, Nx, Ny, 2)]-1*ypdata[idx(i, j+1, Nx, Ny, 2)]);
            
            IS2_py[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)])*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, j+2, Nx, Ny, 2)]);
            
            IS0_ny[0]=(13/12)*(yndata[idx(i, j+1, Nx, Ny, 2)]-2*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, 0, Nx, Ny, 2)])*(yndata[idx(i, j+1, Nx, Ny, 2)]-2*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, 0, Nx, Ny, 2)])+(1/4)*(3*yndata[idx(i, j+1, Nx, Ny, 2)]-4*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, 0, Nx, Ny, 2)])*(3*yndata[idx(i, j+1, Nx, Ny, 2)]-4*yndata[idx(i, j+2, Nx, Ny, 2)]+yndata[idx(i, 0, Nx, Ny, 2)]);
            
            IS1_ny[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, j+1, Nx, Ny, 2)]+yndata[idx(i, j+2, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, j+1, Nx, Ny, 2)]+yndata[idx(i, j+2, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, j+2, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, j+2, Nx, Ny, 2)]);
            
            IS2_ny[0]=(13/12)*(yndata[idx(i, j-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, j+1, Nx, Ny, 2)])*(yndata[idx(i, j-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, j+1, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, j+1, Nx, Ny, 2)])*(yndata[idx(i, j-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, j+1, Nx, Ny, 2)]);
            
            
            IS0_py[1]=(13/12)*(taopdata[idx(i, j-2, Nx, Ny, 1)]-2*taopdata[idx(i, j-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, j-2, Nx, Ny, 1)]-2*taopdata[idx(i, j-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, j-2, Nx, Ny, 1)]-4*taopdata[idx(i, j-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, j-2, Nx, Ny, 1)]-4*taopdata[idx(i, j-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[1]=(13/12)*(taopdata[idx(i, j-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, j+1, Nx, Ny, 1)])*(taopdata[idx(i, j-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, j-1, Nx, Ny, 1)]-1*taopdata[idx(i, j+1, Nx, Ny, 1)])*(taopdata[idx(i, j-1, Nx, Ny, 1)]-1*taopdata[idx(i, j+1, Nx, Ny, 1)]);
            
            IS2_py[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)])*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS0_ny[1]=(13/12)*(taondata[idx(i, j+1, Nx, Ny, 1)]-2*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, 0, Nx, Ny, 1)])*(taondata[idx(i, j+1, Nx, Ny, 1)]-2*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, 0, Nx, Ny, 1)])+(1/4)*(3*taondata[idx(i, j+1, Nx, Ny, 1)]-4*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, 0, Nx, Ny, 1)])*(3*taondata[idx(i, j+1, Nx, Ny, 1)]-4*taondata[idx(i, j+2, Nx, Ny, 1)]+taondata[idx(i, 0, Nx, Ny, 1)]);
            
            IS1_ny[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, j+1, Nx, Ny, 1)]+taondata[idx(i, j+2, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, j+1, Nx, Ny, 1)]+taondata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, j+2, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS2_ny[1]=(13/12)*(taondata[idx(i, j-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, j+1, Nx, Ny, 1)])*(taondata[idx(i, j-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, j+1, Nx, Ny, 1)])*(taondata[idx(i, j-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, j+1, Nx, Ny, 1)]);
            
            
            IS0_py[2]=(13/12)*(taopdata[idx(i, j-2, Nx, Ny, 3)]-2*taopdata[idx(i, j-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, j-2, Nx, Ny, 3)]-2*taopdata[idx(i, j-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, j-2, Nx, Ny, 3)]-4*taopdata[idx(i, j-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, j-2, Nx, Ny, 3)]-4*taopdata[idx(i, j-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)]);
            
            IS1_py[2]=(13/12)*(taopdata[idx(i, j-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, j+1, Nx, Ny, 3)])*(taopdata[idx(i, j-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, j+1, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, j-1, Nx, Ny, 3)]-1*taopdata[idx(i, j+1, Nx, Ny, 3)])*(taopdata[idx(i, j-1, Nx, Ny, 3)]-1*taopdata[idx(i, j+1, Nx, Ny, 3)]);
            
            IS2_py[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)])*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, j+2, Nx, Ny, 3)]);
            
            IS0_ny[2]=(13/12)*(taondata[idx(i, j+1, Nx, Ny, 3)]-2*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, 0, Nx, Ny, 3)])*(taondata[idx(i, j+1, Nx, Ny, 3)]-2*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, 0, Nx, Ny, 3)])+(1/4)*(3*taondata[idx(i, j+1, Nx, Ny, 3)]-4*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, 0, Nx, Ny, 3)])*(3*taondata[idx(i, j+1, Nx, Ny, 3)]-4*taondata[idx(i, j+2, Nx, Ny, 3)]+taondata[idx(i, 0, Nx, Ny, 3)]);
            
            IS1_ny[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, j+1, Nx, Ny, 3)]+taondata[idx(i, j+2, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, j+1, Nx, Ny, 3)]+taondata[idx(i, j+2, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, j+2, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, j+2, Nx, Ny, 3)]);
            
            IS2_ny[2]=(13/12)*(taondata[idx(i, j-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, j+1, Nx, Ny, 3)])*(taondata[idx(i, j-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, j+1, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, j+1, Nx, Ny, 3)])*(taondata[idx(i, j-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, j+1, Nx, Ny, 3)]);
            
            
            IS0_py[3]=(13/12)*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-2*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-2*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-4*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-4*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[3]=(13/12)*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, j+1, Nx, Ny, 1)])*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, j+1, Nx, Ny, 1)])*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, j+1, Nx, Ny, 1)]);
            
            IS2_py[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)])*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS0_ny[3]=(13/12)*(Cjndata[idx(i, j+1, Nx, Ny, 1)]-2*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, 0, Nx, Ny, 1)])*(Cjndata[idx(i, j+1, Nx, Ny, 1)]-2*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, 0, Nx, Ny, 1)])+(1/4)*(3*Cjndata[idx(i, j+1, Nx, Ny, 1)]-4*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, 0, Nx, Ny, 1)])*(3*Cjndata[idx(i, j+1, Nx, Ny, 1)]-4*Cjndata[idx(i, j+2, Nx, Ny, 1)]+Cjndata[idx(i, 0, Nx, Ny, 1)]);
            
            IS1_ny[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, j+1, Nx, Ny, 1)]+Cjndata[idx(i, j+2, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, j+1, Nx, Ny, 1)]+Cjndata[idx(i, j+2, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, j+2, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, j+2, Nx, Ny, 1)]);
            
            IS2_ny[3]=(13/12)*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, j+1, Nx, Ny, 1)])*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, j+1, Nx, Ny, 1)])*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, j+1, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0py[k]=(1/10)*(1/(Epsilon+IS0_py[k]))*(1/(Epsilon+IS0_py[k]));
                alpha_1py[k]=(6/10)*(1/(Epsilon+IS1_py[k]))*(1/(Epsilon+IS1_py[k]));
                alpha_2py[k]=(3/10)*(1/(Epsilon+IS2_py[k]))*(1/(Epsilon+IS2_py[k]));
                alpha_0ny[k]=(1/10)*(1/(Epsilon+IS0_ny[k]))*(1/(Epsilon+IS0_ny[k]));
                alpha_1ny[k]=(6/10)*(1/(Epsilon+IS1_ny[k]))*(1/(Epsilon+IS1_ny[k]));
                alpha_2ny[k]=(3/10)*(1/(Epsilon+IS2_ny[k]))*(1/(Epsilon+IS2_ny[k]));
                
                w0_py[k]=alpha_0py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w1_py[k]=alpha_1py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w2_py[k]=alpha_2py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w0_ny[k]=alpha_0ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w1_ny[k]=alpha_1ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w2_ny[k]=alpha_2ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
            }
            
            u_tpphy[0]=w0_py[0]*((2/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j+1, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j+1, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, j+2, Nx, Ny, 2)]);
            
            u_tnphy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+1, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j+1, Nx, Ny, 2)]-(1/6)*yndata[idx(i, j+2, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j+1, Nx, Ny, 2)]-(7/6)*yndata[idx(i, j+2, Nx, Ny, 2)]+(2/6)*yndata[idx(i, 0, Nx, Ny, 2)]);
            
            u_tpnhy[0]=w0_py[0]*((2/6)*ypdata[idx(i, j-3, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j-1, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, j+1, Nx, Ny, 2)]);
            
            u_tnnhy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, j-2, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]-(1/6)*yndata[idx(i, j+1, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 2)]-(7/6)*yndata[idx(i, j+1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+2, Nx, Ny, 2)]);
            
            
            u_tpphy[1]=w0_py[1]*((2/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j+1, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, j+2, Nx, Ny, 1)]);
            
            u_tnphy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+1, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*taondata[idx(i, j+2, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j+1, Nx, Ny, 1)]-(7/6)*taondata[idx(i, j+2, Nx, Ny, 1)]+(2/6)*taondata[idx(i, 0, Nx, Ny, 1)]);
            
            u_tpnhy[1]=w0_py[1]*((2/6)*taopdata[idx(i, j-3, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j-1, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, j+1, Nx, Ny, 1)]);
            
            u_tnnhy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]-(1/6)*taondata[idx(i, j+1, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 1)]-(7/6)*taondata[idx(i, j+1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+2, Nx, Ny, 1)]);
            
            
            u_tpphy[2]=w0_py[2]*((2/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j+1, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j+1, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, j+2, Nx, Ny, 3)]);
            
            u_tnphy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+1, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j+1, Nx, Ny, 3)]-(1/6)*taondata[idx(i, j+2, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j+1, Nx, Ny, 3)]-(7/6)*taondata[idx(i, j+2, Nx, Ny, 3)]+(2/6)*taondata[idx(i, 0, Nx, Ny, 3)]);
            
            u_tpnhy[2]=w0_py[2]*((2/6)*taopdata[idx(i, j-3, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j-1, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, j+1, Nx, Ny, 3)]);
            
            u_tnnhy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, j-2, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]-(1/6)*taondata[idx(i, j+1, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 3)]-(7/6)*taondata[idx(i, j+1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+2, Nx, Ny, 3)]);
            
            
            u_tpphy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, j+2, Nx, Ny, 1)]);
            
            u_tnphy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, 0, Nx, Ny, 1)]);
            
            u_tpnhy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, j-3, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)]);
            
            u_tnnhy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+2, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++){
                yydata[idx(i, j, Nx, Ny, k)]=-(1/dy)*((u_tpphy[k]-u_tpnhy[k])+(u_tnphy[k]-u_tnnhy[k]));
            }
    }

            i=Nx-2;
    for(j=0; j<Ny; j++){
            // x direction
            IS0_px[0]=(13/12)*(ypdata[idx(i-2, j, Nx, Ny, 1)]-2*ypdata[idx(i-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(i-2, j, Nx, Ny, 1)]-2*ypdata[idx(i-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(i-2, j, Nx, Ny, 1)]-4*ypdata[idx(i-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(i-2, j, Nx, Ny, 1)]-4*ypdata[idx(i-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_px[0]=(13/12)*(ypdata[idx(i-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(i+1, j, Nx, Ny, 1)])*(ypdata[idx(i-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(i+1, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(i-1, j, Nx, Ny, 1)]-1*ypdata[idx(i+1, j, Nx, Ny, 1)])*(ypdata[idx(i-1, j, Nx, Ny, 1)]-1*ypdata[idx(i+1, j, Nx, Ny, 1)]);
            
            IS2_px[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(0, j, Nx, Ny, 1)])*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(0, j, Nx, Ny, 1)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(0, j, Nx, Ny, 1)])*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(i+1, j, Nx, Ny, 1)]+ypdata[idx(0, j, Nx, Ny, 1)]);
            
            IS0_nx[0]=(13/12)*(yndata[idx(i+1, j, Nx, Ny, 1)]-2*yndata[idx(0, j, Nx, Ny, 1)]+yndata[idx(1, j, Nx, Ny, 1)])*(yndata[idx(i+1, j, Nx, Ny, 1)]-2*yndata[idx(0, j, Nx, Ny, 1)]+yndata[idx(1, j, Nx, Ny, 1)])+(1/4)*(3*yndata[idx(i+1, j, Nx, Ny, 1)]-4*yndata[idx(0, j, Nx, Ny, 1)]+yndata[idx(1, j, Nx, Ny, 1)])*(3*yndata[idx(i+1, j, Nx, Ny, 1)]-4*yndata[idx(0, j, Nx, Ny, 1)]+yndata[idx(1, j, Nx, Ny, 1)]);
            
            IS1_nx[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(i+1, j, Nx, Ny, 1)]+yndata[idx(0, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(i+1, j, Nx, Ny, 1)]+yndata[idx(0, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(0, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(0, j, Nx, Ny, 1)]);
            
            IS2_nx[0]=(13/12)*(yndata[idx(i-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(i+1, j, Nx, Ny, 1)])*(yndata[idx(i-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(i+1, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(i+1, j, Nx, Ny, 1)])*(yndata[idx(i-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(i+1, j, Nx, Ny, 1)]);
            
            
            IS0_px[1]=(13/12)*(taopdata[idx(i-2, j, Nx, Ny, 0)]-2*taopdata[idx(i-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(i-2, j, Nx, Ny, 0)]-2*taopdata[idx(i-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(i-2, j, Nx, Ny, 0)]-4*taopdata[idx(i-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(i-2, j, Nx, Ny, 0)]-4*taopdata[idx(i-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)]);
            
            IS1_px[1]=(13/12)*(taopdata[idx(i-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(i+1, j, Nx, Ny, 0)])*(taopdata[idx(i-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(i-1, j, Nx, Ny, 0)]-1*taopdata[idx(i+1, j, Nx, Ny, 0)])*(taopdata[idx(i-1, j, Nx, Ny, 0)]-1*taopdata[idx(i+1, j, Nx, Ny, 0)]);
            
            IS2_px[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(0, j, Nx, Ny, 0)])*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(0, j, Nx, Ny, 0)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(0, j, Nx, Ny, 0)])*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(i+1, j, Nx, Ny, 0)]+taopdata[idx(0, j, Nx, Ny, 0)]);
            
            IS0_nx[1]=(13/12)*(taondata[idx(i+1, j, Nx, Ny, 0)]-2*taondata[idx(0, j, Nx, Ny, 0)]+taondata[idx(1, j, Nx, Ny, 0)])*(taondata[idx(i+1, j, Nx, Ny, 0)]-2*taondata[idx(0, j, Nx, Ny, 0)]+taondata[idx(1, j, Nx, Ny, 0)])+(1/4)*(3*taondata[idx(i+1, j, Nx, Ny, 0)]-4*taondata[idx(0, j, Nx, Ny, 0)]+taondata[idx(1, j, Nx, Ny, 0)])*(3*taondata[idx(i+1, j, Nx, Ny, 0)]-4*taondata[idx(0, j, Nx, Ny, 0)]+taondata[idx(1, j, Nx, Ny, 0)]);
            
            IS1_nx[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(i+1, j, Nx, Ny, 0)]+taondata[idx(0, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(i+1, j, Nx, Ny, 0)]+taondata[idx(0, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(0, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(0, j, Nx, Ny, 0)]);
            
            IS2_nx[1]=(13/12)*(taondata[idx(i-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(i+1, j, Nx, Ny, 0)])*(taondata[idx(i-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(i+1, j, Nx, Ny, 0)])*(taondata[idx(i-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(i+1, j, Nx, Ny, 0)]);
            
            
            IS0_px[2]=(13/12)*(taopdata[idx(i-2, j, Nx, Ny, 2)]-2*taopdata[idx(i-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(i-2, j, Nx, Ny, 2)]-2*taopdata[idx(i-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(i-2, j, Nx, Ny, 2)]-4*taopdata[idx(i-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(i-2, j, Nx, Ny, 2)]-4*taopdata[idx(i-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_px[2]=(13/12)*(taopdata[idx(i-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(i+1, j, Nx, Ny, 2)])*(taopdata[idx(i-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(i+1, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(i-1, j, Nx, Ny, 2)]-1*taopdata[idx(i+1, j, Nx, Ny, 2)])*(taopdata[idx(i-1, j, Nx, Ny, 2)]-1*taopdata[idx(i+1, j, Nx, Ny, 2)]);
            
            IS2_px[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(0, j, Nx, Ny, 2)])*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(0, j, Nx, Ny, 2)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(0, j, Nx, Ny, 2)])*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(i+1, j, Nx, Ny, 2)]+taopdata[idx(0, j, Nx, Ny, 2)]);
            
            IS0_nx[2]=(13/12)*(taondata[idx(i+1, j, Nx, Ny, 2)]-2*taondata[idx(0, j, Nx, Ny, 2)]+taondata[idx(1, j, Nx, Ny, 2)])*(taondata[idx(i+1, j, Nx, Ny, 2)]-2*taondata[idx(0, j, Nx, Ny, 2)]+taondata[idx(1, j, Nx, Ny, 2)])+(1/4)*(3*taondata[idx(i+1, j, Nx, Ny, 2)]-4*taondata[idx(0, j, Nx, Ny, 2)]+taondata[idx(1, j, Nx, Ny, 2)])*(3*taondata[idx(i+1, j, Nx, Ny, 2)]-4*taondata[idx(0, j, Nx, Ny, 2)]+taondata[idx(1, j, Nx, Ny, 2)]);
            
            IS1_nx[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(i+1, j, Nx, Ny, 2)]+taondata[idx(0, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(i+1, j, Nx, Ny, 2)]+taondata[idx(0, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(0, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(0, j, Nx, Ny, 2)]);
            
            IS2_nx[2]=(13/12)*(taondata[idx(i-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(i+1, j, Nx, Ny, 2)])*(taondata[idx(i-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(i+1, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(i+1, j, Nx, Ny, 2)])*(taondata[idx(i-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(i+1, j, Nx, Ny, 2)]);
            
            
            IS0_px[3]=(13/12)*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-2*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-2*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-4*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-4*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)]);
            
            IS1_px[3]=(13/12)*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(i+1, j, Nx, Ny, 0)])*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(i+1, j, Nx, Ny, 0)])*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(i+1, j, Nx, Ny, 0)]);
            
            IS2_px[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(0, j, Nx, Ny, 0)])*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(0, j, Nx, Ny, 0)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(0, j, Nx, Ny, 0)])*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(i+1, j, Nx, Ny, 0)]+Cjpdata[idx(0, j, Nx, Ny, 0)]);
            
            IS0_nx[3]=(13/12)*(Cjndata[idx(i+1, j, Nx, Ny, 0)]-2*Cjndata[idx(0, j, Nx, Ny, 0)]+Cjndata[idx(1, j, Nx, Ny, 0)])*(Cjndata[idx(i+1, j, Nx, Ny, 0)]-2*Cjndata[idx(0, j, Nx, Ny, 0)]+Cjndata[idx(1, j, Nx, Ny, 0)])+(1/4)*(3*Cjndata[idx(i+1, j, Nx, Ny, 0)]-4*Cjndata[idx(0, j, Nx, Ny, 0)]+Cjndata[idx(1, j, Nx, Ny, 0)])*(3*Cjndata[idx(i+1, j, Nx, Ny, 0)]-4*Cjndata[idx(0, j, Nx, Ny, 0)]+Cjndata[idx(1, j, Nx, Ny, 0)]);
            
            IS1_nx[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(i+1, j, Nx, Ny, 0)]+Cjndata[idx(0, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(i+1, j, Nx, Ny, 0)]+Cjndata[idx(0, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(0, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(0, j, Nx, Ny, 0)]);
            
            IS2_nx[3]=(13/12)*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(i+1, j, Nx, Ny, 0)])*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(i+1, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(i+1, j, Nx, Ny, 0)])*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(i+1, j, Nx, Ny, 0)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0px[k]=(1/10)*(1/(Epsilon+IS0_px[k]))*(1/(Epsilon+IS0_px[k]));
                alpha_1px[k]=(6/10)*(1/(Epsilon+IS1_px[k]))*(1/(Epsilon+IS1_px[k]));
                alpha_2px[k]=(3/10)*(1/(Epsilon+IS2_px[k]))*(1/(Epsilon+IS2_px[k]));
                alpha_0nx[k]=(1/10)*(1/(Epsilon+IS0_nx[k]))*(1/(Epsilon+IS0_nx[k]));
                alpha_1nx[k]=(6/10)*(1/(Epsilon+IS1_nx[k]))*(1/(Epsilon+IS1_nx[k]));
                alpha_2nx[k]=(3/10)*(1/(Epsilon+IS2_nx[k]))*(1/(Epsilon+IS2_nx[k]));
                
                w0_px[k]=alpha_0px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w1_px[k]=alpha_1px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w2_px[k]=alpha_2px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w0_nx[k]=alpha_0nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
                w1_nx[k]=alpha_1nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
                w2_nx[k]=alpha_2nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
            }
            
            u_tpphx[0]=w0_px[0]*((2/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i+1, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i+1, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(0, j, Nx, Ny, 1)]);
            
            u_tnphx[0]=w2_nx[0]*((-1/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i+1, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i+1, j, Nx, Ny, 1)]-(1/6)*yndata[idx(0, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i+1, j, Nx, Ny, 1)]-(7/6)*yndata[idx(i+2, j, Nx, Ny, 1)]+(2/6)*yndata[idx(0, j, Nx, Ny, 1)]);
            
            u_tpnhx[0]=w0_px[0]*((2/6)*ypdata[idx(i-3, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i-1, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(i+1, j, Nx, Ny, 1)]);
            
            u_tnnhx[0]=w2_nx[0]*((-1/6)*yndata[idx(i-2, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]-(1/6)*yndata[idx(i+1, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 1)]-(7/6)*yndata[idx(i+1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(0, j, Nx, Ny, 1)]);
            
            
            u_tpphx[1]=w0_px[1]*((2/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i+1, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(0, j, Nx, Ny, 0)]);
            
            u_tnphx[1]=w2_nx[1]*((-1/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i+1, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*taondata[idx(0, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i+1, j, Nx, Ny, 0)]-(7/6)*taondata[idx(0, j, Nx, Ny, 0)]+(2/6)*taondata[idx(0, j, Nx, Ny, 0)]);
            
            u_tpnhx[1]=w0_px[1]*((2/6)*taopdata[idx(i-3, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i-1, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(i+1, j, Nx, Ny, 0)]);
            
            u_tnnhx[1]=w2_nx[1]*((-1/6)*taondata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]-(1/6)*taondata[idx(i+1, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 0)]-(7/6)*taondata[idx(i+1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(0, j, Nx, Ny, 0)]);
            
            
            u_tpphx[2]=w0_px[2]*((2/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i+1, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i+1, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(0, j, Nx, Ny, 2)]);
            
            u_tnphx[2]=w2_nx[2]*((-1/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i+1, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i+1, j, Nx, Ny, 2)]-(1/6)*taondata[idx(0, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i+1, j, Nx, Ny, 2)]-(7/6)*taondata[idx(0, j, Nx, Ny, 2)]+(2/6)*taondata[idx(1, j, Nx, Ny, 2)]);
            
            u_tpnhx[2]=w0_px[2]*((2/6)*taopdata[idx(i-3, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i-1, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(i+1, j, Nx, Ny, 2)]);
            
            u_tnnhx[2]=w2_nx[2]*((-1/6)*taondata[idx(i-2, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]-(1/6)*taondata[idx(i+1, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 2)]-(7/6)*taondata[idx(i+1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(0, j, Nx, Ny, 2)]);
            
            
            u_tpphx[3]=w0_px[3]*((2/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(0, j, Nx, Ny, 0)]);
            
            u_tnphx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(0, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(0, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(1, j, Nx, Ny, 0)]);
            
            u_tpnhx[3]=w0_px[3]*((2/6)*Cjpdata[idx(i-3, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(i+1, j, Nx, Ny, 0)]);
            
            u_tnnhx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(i+1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(0, j, Nx, Ny, 0)]);
            
            for(k=0;k<4;k++){
                yxdata[idx(i, j, Nx, Ny, k)]=-(1/dx)*((u_tpphx[k]-u_tpnhx[k])+(u_tnphx[k]-u_tnnhx[k]));
            }
    }
    
            j=Ny-2;
        for (i=0; i<Nx; i++){
            // x direction
            IS0_py[0]=(13/12)*(ypdata[idx(i, j-2, Nx, Ny, 2)]-2*ypdata[idx(i, j-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, j-2, Nx, Ny, 2)]-2*ypdata[idx(i, j-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, j-2, Nx, Ny, 2)]-4*ypdata[idx(i, j-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, j-2, Nx, Ny, 2)]-4*ypdata[idx(i, j-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_py[0]=(13/12)*(ypdata[idx(i, j-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, j+1, Nx, Ny, 2)])*(ypdata[idx(i, j-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, j+1, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, j-1, Nx, Ny, 2)]-1*ypdata[idx(i, j+1, Nx, Ny, 2)])*(ypdata[idx(i, j-1, Nx, Ny, 2)]-1*ypdata[idx(i, j+1, Nx, Ny, 2)]);
            
            IS2_py[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, 0, Nx, Ny, 2)])*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, 0, Nx, Ny, 2)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, 0, Nx, Ny, 2)])*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, j+1, Nx, Ny, 2)]+ypdata[idx(i, 0, Nx, Ny, 2)]);
            
            IS0_ny[0]=(13/12)*(yndata[idx(i, j+1, Nx, Ny, 2)]-2*yndata[idx(i, 0, Nx, Ny, 2)]+yndata[idx(i, 1, Nx, Ny, 2)])*(yndata[idx(i, j+1, Nx, Ny, 2)]-2*yndata[idx(i, 0, Nx, Ny, 2)]+yndata[idx(i, 1, Nx, Ny, 2)])+(1/4)*(3*yndata[idx(i, j+1, Nx, Ny, 2)]-4*yndata[idx(i, 0, Nx, Ny, 2)]+yndata[idx(i, 1, Nx, Ny, 2)])*(3*yndata[idx(i, j+1, Nx, Ny, 2)]-4*yndata[idx(i, 0, Nx, Ny, 2)]+yndata[idx(i, 1, Nx, Ny, 2)]);
            
            IS1_ny[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, j+1, Nx, Ny, 2)]+yndata[idx(i, 0, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, j+1, Nx, Ny, 2)]+yndata[idx(i, 0, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, 0, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, 0, Nx, Ny, 2)]);
            
            IS2_ny[0]=(13/12)*(yndata[idx(i, j-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, j+1, Nx, Ny, 2)])*(yndata[idx(i, j-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, j+1, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, j+1, Nx, Ny, 2)])*(yndata[idx(i, j-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, j+1, Nx, Ny, 2)]);
            
            
            IS0_py[1]=(13/12)*(taopdata[idx(i, j-2, Nx, Ny, 1)]-2*taopdata[idx(i, j-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, j-2, Nx, Ny, 1)]-2*taopdata[idx(i, j-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, j-2, Nx, Ny, 1)]-4*taopdata[idx(i, j-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, j-2, Nx, Ny, 1)]-4*taopdata[idx(i, j-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[1]=(13/12)*(taopdata[idx(i, j-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, j+1, Nx, Ny, 1)])*(taopdata[idx(i, j-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, j-1, Nx, Ny, 1)]-1*taopdata[idx(i, j+1, Nx, Ny, 1)])*(taopdata[idx(i, j-1, Nx, Ny, 1)]-1*taopdata[idx(i, j+1, Nx, Ny, 1)]);
            
            IS2_py[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, 0, Nx, Ny, 1)])*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, 0, Nx, Ny, 1)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, 0, Nx, Ny, 1)])*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, j+1, Nx, Ny, 1)]+taopdata[idx(i, 0, Nx, Ny, 1)]);
            
            IS0_ny[1]=(13/12)*(taondata[idx(i, j+1, Nx, Ny, 1)]-2*taondata[idx(i, 0, Nx, Ny, 1)]+taondata[idx(i, 1, Nx, Ny, 1)])*(taondata[idx(i, j+1, Nx, Ny, 1)]-2*taondata[idx(i, 0, Nx, Ny, 1)]+taondata[idx(i, 1, Nx, Ny, 1)])+(1/4)*(3*taondata[idx(i, j+1, Nx, Ny, 1)]-4*taondata[idx(i, 0, Nx, Ny, 1)]+taondata[idx(i, 1, Nx, Ny, 1)])*(3*taondata[idx(i, j+1, Nx, Ny, 1)]-4*taondata[idx(i, 0, Nx, Ny, 1)]+taondata[idx(i, 1, Nx, Ny, 1)]);
            
            IS1_ny[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, j+1, Nx, Ny, 1)]+taondata[idx(i, 0, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, j+1, Nx, Ny, 1)]+taondata[idx(i, 0, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, 0, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, 0, Nx, Ny, 1)]);
            
            IS2_ny[1]=(13/12)*(taondata[idx(i, j-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, j+1, Nx, Ny, 1)])*(taondata[idx(i, j-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, j+1, Nx, Ny, 1)])*(taondata[idx(i, j-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, j+1, Nx, Ny, 1)]);
            
            
            IS0_py[2]=(13/12)*(taopdata[idx(i, j-2, Nx, Ny, 3)]-2*taopdata[idx(i, j-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, j-2, Nx, Ny, 3)]-2*taopdata[idx(i, j-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, j-2, Nx, Ny, 3)]-4*taopdata[idx(i, j-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, j-2, Nx, Ny, 3)]-4*taopdata[idx(i, j-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)]);
            
            IS1_py[2]=(13/12)*(taopdata[idx(i, j-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, j+1, Nx, Ny, 3)])*(taopdata[idx(i, j-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, j+1, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, j-1, Nx, Ny, 3)]-1*taopdata[idx(i, j+1, Nx, Ny, 3)])*(taopdata[idx(i, j-1, Nx, Ny, 3)]-1*taopdata[idx(i, j+1, Nx, Ny, 3)]);
            
            IS2_py[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, 0, Nx, Ny, 3)])*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, 0, Nx, Ny, 3)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, 0, Nx, Ny, 3)])*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, j+1, Nx, Ny, 3)]+taopdata[idx(i, 0, Nx, Ny, 3)]);
            
            IS0_ny[2]=(13/12)*(taondata[idx(i, j+1, Nx, Ny, 3)]-2*taondata[idx(i, 0, Nx, Ny, 3)]+taondata[idx(i, 1, Nx, Ny, 3)])*(taondata[idx(i, j+1, Nx, Ny, 3)]-2*taondata[idx(i, 0, Nx, Ny, 3)]+taondata[idx(i, 1, Nx, Ny, 3)])+(1/4)*(3*taondata[idx(i, j+1, Nx, Ny, 3)]-4*taondata[idx(i, 0, Nx, Ny, 3)]+taondata[idx(i, 1, Nx, Ny, 3)])*(3*taondata[idx(i, j+1, Nx, Ny, 3)]-4*taondata[idx(i, 0, Nx, Ny, 3)]+taondata[idx(i, 1, Nx, Ny, 3)]);
            
            IS1_ny[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, j+1, Nx, Ny, 3)]+taondata[idx(i, 0, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, j+1, Nx, Ny, 3)]+taondata[idx(i, 0, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, 0, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, 0, Nx, Ny, 3)]);
            
            IS2_ny[2]=(13/12)*(taondata[idx(i, j-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, j+1, Nx, Ny, 3)])*(taondata[idx(i, j-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, j+1, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, j+1, Nx, Ny, 3)])*(taondata[idx(i, j-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, j+1, Nx, Ny, 3)]);
            
            
            IS0_py[3]=(13/12)*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-2*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-2*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-4*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-4*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[3]=(13/12)*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, j+1, Nx, Ny, 1)])*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, j+1, Nx, Ny, 1)])*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, j+1, Nx, Ny, 1)]);
            
            IS2_py[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, 0, Nx, Ny, 1)])*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, 0, Nx, Ny, 1)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, 0, Nx, Ny, 1)])*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, j+1, Nx, Ny, 1)]+Cjpdata[idx(i, 0, Nx, Ny, 1)]);
            
            IS0_ny[3]=(13/12)*(Cjndata[idx(i, j+1, Nx, Ny, 1)]-2*Cjndata[idx(i, 0, Nx, Ny, 1)]+Cjndata[idx(i, 1, Nx, Ny, 1)])*(Cjndata[idx(i, j+1, Nx, Ny, 1)]-2*Cjndata[idx(i, 0, Nx, Ny, 1)]+Cjndata[idx(i, 1, Nx, Ny, 1)])+(1/4)*(3*Cjndata[idx(i, j+1, Nx, Ny, 1)]-4*Cjndata[idx(i, 0, Nx, Ny, 1)]+Cjndata[idx(i, 1, Nx, Ny, 1)])*(3*Cjndata[idx(i, j+1, Nx, Ny, 1)]-4*Cjndata[idx(i, 0, Nx, Ny, 1)]+Cjndata[idx(i, 1, Nx, Ny, 1)]);
            
            IS1_ny[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, j+1, Nx, Ny, 1)]+Cjndata[idx(i, 0, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, j+1, Nx, Ny, 1)]+Cjndata[idx(i, 0, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, 0, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, 0, Nx, Ny, 1)]);
            
            IS2_ny[3]=(13/12)*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, j+1, Nx, Ny, 1)])*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, j+1, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, j+1, Nx, Ny, 1)])*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, j+1, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0py[k]=(1/10)*(1/(Epsilon+IS0_py[k]))*(1/(Epsilon+IS0_py[k]));
                alpha_1py[k]=(6/10)*(1/(Epsilon+IS1_py[k]))*(1/(Epsilon+IS1_py[k]));
                alpha_2py[k]=(3/10)*(1/(Epsilon+IS2_py[k]))*(1/(Epsilon+IS2_py[k]));
                alpha_0ny[k]=(1/10)*(1/(Epsilon+IS0_ny[k]))*(1/(Epsilon+IS0_ny[k]));
                alpha_1ny[k]=(6/10)*(1/(Epsilon+IS1_ny[k]))*(1/(Epsilon+IS1_ny[k]));
                alpha_2ny[k]=(3/10)*(1/(Epsilon+IS2_ny[k]))*(1/(Epsilon+IS2_ny[k]));
                
                w0_py[k]=alpha_0py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w1_py[k]=alpha_1py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w2_py[k]=alpha_2py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w0_ny[k]=alpha_0ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w1_ny[k]=alpha_1ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w2_ny[k]=alpha_2ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
            }
            
            u_tpphy[0]=w0_py[0]*((2/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j+1, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j+1, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, 0, Nx, Ny, 2)]);
            
            u_tnphy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j+1, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j+1, Nx, Ny, 2)]-(1/6)*yndata[idx(i, 0, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j+1, Nx, Ny, 2)]-(7/6)*yndata[idx(i, 0, Nx, Ny, 2)]+(2/6)*yndata[idx(i, 1, Nx, Ny, 2)]);
            
            u_tpnhy[0]=w0_py[0]*((2/6)*ypdata[idx(i, j-3, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j-1, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, j+1, Nx, Ny, 2)]);
            
            u_tnnhy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, j-2, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]-(1/6)*yndata[idx(i, j+1, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 2)]-(7/6)*yndata[idx(i, j+1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, 0, Nx, Ny, 2)]);
            
            
            u_tpphy[1]=w0_py[1]*((2/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j+1, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, 0, Nx, Ny, 1)]);
            
            u_tnphy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j+1, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*taondata[idx(i, 0, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j+1, Nx, Ny, 1)]-(7/6)*taondata[idx(i, 0, Nx, Ny, 1)]+(2/6)*taondata[idx(i, 1, Nx, Ny, 1)]);
            
            u_tpnhy[1]=w0_py[1]*((2/6)*taopdata[idx(i, j-3, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j-1, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, j+1, Nx, Ny, 1)]);
            
            u_tnnhy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]-(1/6)*taondata[idx(i, j+1, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 1)]-(7/6)*taondata[idx(i, j+1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, 0, Nx, Ny, 1)]);
            
            
            u_tpphy[2]=w0_py[2]*((2/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j+1, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j+1, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, 0, Nx, Ny, 3)]);
            
            u_tnphy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j+1, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j+1, Nx, Ny, 3)]-(1/6)*taondata[idx(i, 0, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j+1, Nx, Ny, 3)]-(7/6)*taondata[idx(i, 0, Nx, Ny, 3)]+(2/6)*taondata[idx(i, 1, Nx, Ny, 3)]);
            
            u_tpnhy[2]=w0_py[2]*((2/6)*taopdata[idx(i, j-3, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j-1, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, j+1, Nx, Ny, 3)]);
            
            u_tnnhy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, j-2, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]-(1/6)*taondata[idx(i, j+1, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 3)]-(7/6)*taondata[idx(i, j+1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, 0, Nx, Ny, 3)]);
            
            
            u_tpphy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, 0, Nx, Ny, 1)]);
            
            u_tnphy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, 0, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, 0, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, 1, Nx, Ny, 1)]);
            
            u_tpnhy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, j-3, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, j+1, Nx, Ny, 1)]);
            
            u_tnnhy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, j+1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, 0, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++){
                yydata[idx(i, j, Nx, Ny, k)]=-(1/dy)*((u_tpphy[k]-u_tpnhy[k])+(u_tnphy[k]-u_tnnhy[k]));
            }
    }

            i=Nx-1;
    for(j=0; j<Ny; j++){
            // x direction
            IS0_px[0]=(13/12)*(ypdata[idx(i-2, j, Nx, Ny, 1)]-2*ypdata[idx(i-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(i-2, j, Nx, Ny, 1)]-2*ypdata[idx(i-1, j, Nx, Ny, 1)]+ypdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(i-2, j, Nx, Ny, 1)]-4*ypdata[idx(i-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)])*(ypdata[idx(i-2, j, Nx, Ny, 1)]-4*ypdata[idx(i-1, j, Nx, Ny, 1)]+3*ypdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_px[0]=(13/12)*(ypdata[idx(i-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(0, j, Nx, Ny, 1)])*(ypdata[idx(i-1, j, Nx, Ny, 1)]-2*ypdata[idx(i, j, Nx, Ny, 1)]+ypdata[idx(0, j, Nx, Ny, 1)])+(1/4)*(ypdata[idx(i-1, j, Nx, Ny, 1)]-1*ypdata[idx(0, j, Nx, Ny, 1)])*(ypdata[idx(i-1, j, Nx, Ny, 1)]-1*ypdata[idx(0, j, Nx, Ny, 1)]);
            
            IS2_px[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(0, j, Nx, Ny, 1)]+ypdata[idx(1, j, Nx, Ny, 1)])*(ypdata[idx(i, j, Nx, Ny, 1)]-2*ypdata[idx(0, j, Nx, Ny, 1)]+ypdata[idx(1, j, Nx, Ny, 1)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(0, j, Nx, Ny, 1)]+ypdata[idx(1, j, Nx, Ny, 1)])*(3*ypdata[idx(i, j, Nx, Ny, 1)]-4*ypdata[idx(0, j, Nx, Ny, 1)]+ypdata[idx(1, j, Nx, Ny, 1)]);
            
            IS0_nx[0]=(13/12)*(yndata[idx(0, j, Nx, Ny, 1)]-2*yndata[idx(1, j, Nx, Ny, 1)]+yndata[idx(2, j, Nx, Ny, 1)])*(yndata[idx(0, j, Nx, Ny, 1)]-2*yndata[idx(1, j, Nx, Ny, 1)]+yndata[idx(2, j, Nx, Ny, 1)])+(1/4)*(3*yndata[idx(0, j, Nx, Ny, 1)]-4*yndata[idx(1, j, Nx, Ny, 1)]+yndata[idx(2, j, Nx, Ny, 1)])*(3*yndata[idx(0, j, Nx, Ny, 1)]-4*yndata[idx(1, j, Nx, Ny, 1)]+yndata[idx(2, j, Nx, Ny, 1)]);
            
            IS1_nx[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(0, j, Nx, Ny, 1)]+yndata[idx(1, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-2*yndata[idx(0, j, Nx, Ny, 1)]+yndata[idx(1, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(1, j, Nx, Ny, 1)])*(yndata[idx(i, j, Nx, Ny, 1)]-1*yndata[idx(1, j, Nx, Ny, 1)]);
            
            IS2_nx[0]=(13/12)*(yndata[idx(i-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(0, j, Nx, Ny, 1)])*(yndata[idx(i-1, j, Nx, Ny, 1)]-2*yndata[idx(i, j, Nx, Ny, 1)]+yndata[idx(0, j, Nx, Ny, 1)])+(1/4)*(yndata[idx(i-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(0, j, Nx, Ny, 1)])*(yndata[idx(i-1, j, Nx, Ny, 1)]-4*yndata[idx(i, j, Nx, Ny, 1)]+3*yndata[idx(0, j, Nx, Ny, 1)]);
            
            
            IS0_px[1]=(13/12)*(taopdata[idx(i-2, j, Nx, Ny, 0)]-2*taopdata[idx(i-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(i-2, j, Nx, Ny, 0)]-2*taopdata[idx(i-1, j, Nx, Ny, 0)]+taopdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(i-2, j, Nx, Ny, 0)]-4*taopdata[idx(i-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)])*(taopdata[idx(i-2, j, Nx, Ny, 0)]-4*taopdata[idx(i-1, j, Nx, Ny, 0)]+3*taopdata[idx(i, j, Nx, Ny, 0)]);
            
            IS1_px[1]=(13/12)*(taopdata[idx(i-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(0, j, Nx, Ny, 0)])*(taopdata[idx(i-1, j, Nx, Ny, 0)]-2*taopdata[idx(i, j, Nx, Ny, 0)]+taopdata[idx(0, j, Nx, Ny, 0)])+(1/4)*(taopdata[idx(i-1, j, Nx, Ny, 0)]-1*taopdata[idx(0, j, Nx, Ny, 0)])*(taopdata[idx(i-1, j, Nx, Ny, 0)]-1*taopdata[idx(0, j, Nx, Ny, 0)]);
            
            IS2_px[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(0, j, Nx, Ny, 0)]+taopdata[idx(1, j, Nx, Ny, 0)])*(taopdata[idx(i, j, Nx, Ny, 0)]-2*taopdata[idx(0, j, Nx, Ny, 0)]+taopdata[idx(1, j, Nx, Ny, 0)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(0, j, Nx, Ny, 0)]+taopdata[idx(1, j, Nx, Ny, 0)])*(3*taopdata[idx(i, j, Nx, Ny, 0)]-4*taopdata[idx(0, j, Nx, Ny, 0)]+taopdata[idx(1, j, Nx, Ny, 0)]);
            
            IS0_nx[1]=(13/12)*(taondata[idx(0, j, Nx, Ny, 0)]-2*taondata[idx(1, j, Nx, Ny, 0)]+taondata[idx(2, j, Nx, Ny, 0)])*(taondata[idx(0, j, Nx, Ny, 0)]-2*taondata[idx(1, j, Nx, Ny, 0)]+taondata[idx(2, j, Nx, Ny, 0)])+(1/4)*(3*taondata[idx(0, j, Nx, Ny, 0)]-4*taondata[idx(1, j, Nx, Ny, 0)]+taondata[idx(2, j, Nx, Ny, 0)])*(3*taondata[idx(0, j, Nx, Ny, 0)]-4*taondata[idx(1, j, Nx, Ny, 0)]+taondata[idx(2, j, Nx, Ny, 0)]);
            
            IS1_nx[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(0, j, Nx, Ny, 0)]+taondata[idx(1, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-2*taondata[idx(0, j, Nx, Ny, 0)]+taondata[idx(1, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(1, j, Nx, Ny, 0)])*(taondata[idx(i, j, Nx, Ny, 0)]-1*taondata[idx(1, j, Nx, Ny, 0)]);
            
            IS2_nx[1]=(13/12)*(taondata[idx(i-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(0, j, Nx, Ny, 0)])*(taondata[idx(i-1, j, Nx, Ny, 0)]-2*taondata[idx(i, j, Nx, Ny, 0)]+taondata[idx(0, j, Nx, Ny, 0)])+(1/4)*(taondata[idx(i-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(0, j, Nx, Ny, 0)])*(taondata[idx(i-1, j, Nx, Ny, 0)]-4*taondata[idx(i, j, Nx, Ny, 0)]+3*taondata[idx(0, j, Nx, Ny, 0)]);
            
            
            IS0_px[2]=(13/12)*(taopdata[idx(i-2, j, Nx, Ny, 2)]-2*taopdata[idx(i-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(i-2, j, Nx, Ny, 2)]-2*taopdata[idx(i-1, j, Nx, Ny, 2)]+taopdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(i-2, j, Nx, Ny, 2)]-4*taopdata[idx(i-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)])*(taopdata[idx(i-2, j, Nx, Ny, 2)]-4*taopdata[idx(i-1, j, Nx, Ny, 2)]+3*taopdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_px[2]=(13/12)*(taopdata[idx(i-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(0, j, Nx, Ny, 2)])*(taopdata[idx(i-1, j, Nx, Ny, 2)]-2*taopdata[idx(i, j, Nx, Ny, 2)]+taopdata[idx(0, j, Nx, Ny, 2)])+(1/4)*(taopdata[idx(i-1, j, Nx, Ny, 2)]-1*taopdata[idx(0, j, Nx, Ny, 2)])*(taopdata[idx(i-1, j, Nx, Ny, 2)]-1*taopdata[idx(0, j, Nx, Ny, 2)]);
            
            IS2_px[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(0, j, Nx, Ny, 2)]+taopdata[idx(1, j, Nx, Ny, 2)])*(taopdata[idx(i, j, Nx, Ny, 2)]-2*taopdata[idx(0, j, Nx, Ny, 2)]+taopdata[idx(1, j, Nx, Ny, 2)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(0, j, Nx, Ny, 2)]+taopdata[idx(1, j, Nx, Ny, 2)])*(3*taopdata[idx(i, j, Nx, Ny, 2)]-4*taopdata[idx(0, j, Nx, Ny, 2)]+taopdata[idx(1, j, Nx, Ny, 2)]);
            
            IS0_nx[2]=(13/12)*(taondata[idx(0, j, Nx, Ny, 2)]-2*taondata[idx(1, j, Nx, Ny, 2)]+taondata[idx(2, j, Nx, Ny, 2)])*(taondata[idx(0, j, Nx, Ny, 2)]-2*taondata[idx(1, j, Nx, Ny, 2)]+taondata[idx(2, j, Nx, Ny, 2)])+(1/4)*(3*taondata[idx(0, j, Nx, Ny, 2)]-4*taondata[idx(1, j, Nx, Ny, 2)]+taondata[idx(2, j, Nx, Ny, 2)])*(3*taondata[idx(0, j, Nx, Ny, 2)]-4*taondata[idx(1, j, Nx, Ny, 2)]+taondata[idx(2, j, Nx, Ny, 2)]);
            
            IS1_nx[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(0, j, Nx, Ny, 2)]+taondata[idx(1, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-2*taondata[idx(0, j, Nx, Ny, 2)]+taondata[idx(1, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(1, j, Nx, Ny, 2)])*(taondata[idx(i, j, Nx, Ny, 2)]-1*taondata[idx(1, j, Nx, Ny, 2)]);
            
            IS2_nx[2]=(13/12)*(taondata[idx(i-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(0, j, Nx, Ny, 2)])*(taondata[idx(i-1, j, Nx, Ny, 2)]-2*taondata[idx(i, j, Nx, Ny, 2)]+taondata[idx(0, j, Nx, Ny, 2)])+(1/4)*(taondata[idx(i-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(0, j, Nx, Ny, 2)])*(taondata[idx(i-1, j, Nx, Ny, 2)]-4*taondata[idx(i, j, Nx, Ny, 2)]+3*taondata[idx(0, j, Nx, Ny, 2)]);
            
            
            IS0_px[3]=(13/12)*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-2*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-2*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+Cjpdata[idx(i, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-4*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)])*(Cjpdata[idx(i-2, j, Nx, Ny, 0)]-4*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+3*Cjpdata[idx(i, j, Nx, Ny, 0)]);
            
            IS1_px[3]=(13/12)*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(0, j, Nx, Ny, 0)])*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-2*Cjpdata[idx(i, j, Nx, Ny, 0)]+Cjpdata[idx(0, j, Nx, Ny, 0)])+(1/4)*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(0, j, Nx, Ny, 0)])*(Cjpdata[idx(i-1, j, Nx, Ny, 0)]-1*Cjpdata[idx(0, j, Nx, Ny, 0)]);
            
            IS2_px[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(0, j, Nx, Ny, 0)]+Cjpdata[idx(1, j, Nx, Ny, 0)])*(Cjpdata[idx(i, j, Nx, Ny, 0)]-2*Cjpdata[idx(0, j, Nx, Ny, 0)]+Cjpdata[idx(1, j, Nx, Ny, 0)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(0, j, Nx, Ny, 0)]+Cjpdata[idx(1, j, Nx, Ny, 0)])*(3*Cjpdata[idx(i, j, Nx, Ny, 0)]-4*Cjpdata[idx(0, j, Nx, Ny, 0)]+Cjpdata[idx(1, j, Nx, Ny, 0)]);
            
            IS0_nx[3]=(13/12)*(Cjndata[idx(0, j, Nx, Ny, 0)]-2*Cjndata[idx(1, j, Nx, Ny, 0)]+Cjndata[idx(2, j, Nx, Ny, 0)])*(Cjndata[idx(0, j, Nx, Ny, 0)]-2*Cjndata[idx(1, j, Nx, Ny, 0)]+Cjndata[idx(2, j, Nx, Ny, 0)])+(1/4)*(3*Cjndata[idx(0, j, Nx, Ny, 0)]-4*Cjndata[idx(1, j, Nx, Ny, 0)]+Cjndata[idx(2, j, Nx, Ny, 0)])*(3*Cjndata[idx(0, j, Nx, Ny, 0)]-4*Cjndata[idx(1, j, Nx, Ny, 0)]+Cjndata[idx(2, j, Nx, Ny, 0)]);
            
            IS1_nx[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(0, j, Nx, Ny, 0)]+Cjndata[idx(1, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-2*Cjndata[idx(0, j, Nx, Ny, 0)]+Cjndata[idx(1, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(1, j, Nx, Ny, 0)])*(Cjndata[idx(i, j, Nx, Ny, 0)]-1*Cjndata[idx(1, j, Nx, Ny, 0)]);
            
            IS2_nx[3]=(13/12)*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(0, j, Nx, Ny, 0)])*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-2*Cjndata[idx(i, j, Nx, Ny, 0)]+Cjndata[idx(0, j, Nx, Ny, 0)])+(1/4)*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(0, j, Nx, Ny, 0)])*(Cjndata[idx(i-1, j, Nx, Ny, 0)]-4*Cjndata[idx(i, j, Nx, Ny, 0)]+3*Cjndata[idx(0, j, Nx, Ny, 0)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0px[k]=(1/10)*(1/(Epsilon+IS0_px[k]))*(1/(Epsilon+IS0_px[k]));
                alpha_1px[k]=(6/10)*(1/(Epsilon+IS1_px[k]))*(1/(Epsilon+IS1_px[k]));
                alpha_2px[k]=(3/10)*(1/(Epsilon+IS2_px[k]))*(1/(Epsilon+IS2_px[k]));
                alpha_0nx[k]=(1/10)*(1/(Epsilon+IS0_nx[k]))*(1/(Epsilon+IS0_nx[k]));
                alpha_1nx[k]=(6/10)*(1/(Epsilon+IS1_nx[k]))*(1/(Epsilon+IS1_nx[k]));
                alpha_2nx[k]=(3/10)*(1/(Epsilon+IS2_nx[k]))*(1/(Epsilon+IS2_nx[k]));
                
                w0_px[k]=alpha_0px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w1_px[k]=alpha_1px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w2_px[k]=alpha_2px[k]/(alpha_0px[k]+alpha_1px[k]+alpha_2px[k]);
                w0_nx[k]=alpha_0nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
                w1_nx[k]=alpha_1nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
                w2_nx[k]=alpha_2nx[k]/(alpha_0nx[k]+alpha_1nx[k]+alpha_2nx[k]);
            }
            
            u_tpphx[0]=w0_px[0]*((2/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i+1, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(0, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(1, j, Nx, Ny, 1)]);
            
            u_tnphx[0]=w2_nx[0]*((-1/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]+(2/6)*yndata[idx(0, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 1)]+(5/6)*yndata[idx(0, j, Nx, Ny, 1)]-(1/6)*yndata[idx(1, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(0, j, Nx, Ny, 1)]-(7/6)*yndata[idx(1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(2, j, Nx, Ny, 1)]);
            
            u_tpnhx[0]=w0_px[0]*((2/6)*ypdata[idx(i-3, j, Nx, Ny, 1)]-(7/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]+(11/6)*ypdata[idx(i-1, j, Nx, Ny, 1)])+w1_px[0]*((-1/6)*ypdata[idx(i-2, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 1)])+w2_px[0]*((2/6)*ypdata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 1)]-(1/6)*ypdata[idx(0, j, Nx, Ny, 1)]);
            
            u_tnnhx[0]=w2_nx[0]*((-1/6)*yndata[idx(i-2, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(2/6)*yndata[idx(i, j, Nx, Ny, 1)])+w1_nx[0]*((2/6)*yndata[idx(i-1, j, Nx, Ny, 1)]+(5/6)*yndata[idx(i, j, Nx, Ny, 1)]-(1/6)*yndata[idx(0, j, Nx, Ny, 1)])+w0_nx[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 1)]-(7/6)*yndata[idx(0, j, Nx, Ny, 1)]+(2/6)*yndata[idx(1, j, Nx, Ny, 1)]);
            
            
            u_tpphx[1]=w0_px[1]*((2/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(0, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(0, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(1, j, Nx, Ny, 0)]);
            
            u_tnphx[1]=w2_nx[1]*((-1/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]+(2/6)*taondata[idx(0, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 0)]+(5/6)*taondata[idx(0, j, Nx, Ny, 0)]-(1/6)*taondata[idx(1, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(0, j, Nx, Ny, 0)]-(7/6)*taondata[idx(1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(2, j, Nx, Ny, 0)]);
            
            u_tpnhx[1]=w0_px[1]*((2/6)*taopdata[idx(i-3, j, Nx, Ny, 0)]-(7/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]+(11/6)*taopdata[idx(i-1, j, Nx, Ny, 0)])+w1_px[1]*((-1/6)*taopdata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 0)])+w2_px[1]*((2/6)*taopdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 0)]-(1/6)*taopdata[idx(0, j, Nx, Ny, 0)]);
            
            u_tnnhx[1]=w2_nx[1]*((-1/6)*taondata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*taondata[idx(i, j, Nx, Ny, 0)])+w1_nx[1]*((2/6)*taondata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*taondata[idx(i, j, Nx, Ny, 0)]-(1/6)*taondata[idx(0, j, Nx, Ny, 0)])+w0_nx[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 0)]-(7/6)*taondata[idx(0, j, Nx, Ny, 0)]+(2/6)*taondata[idx(1, j, Nx, Ny, 0)]);
            
            
            u_tpphx[2]=w0_px[2]*((2/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(0, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(0, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(1, j, Nx, Ny, 2)]);
            
            u_tnphx[2]=w2_nx[2]*((-1/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]+(2/6)*taondata[idx(0, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 2)]+(5/6)*taondata[idx(0, j, Nx, Ny, 2)]-(1/6)*taondata[idx(1, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(0, j, Nx, Ny, 2)]-(7/6)*taondata[idx(1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(2, j, Nx, Ny, 2)]);
            
            u_tpnhx[2]=w0_px[2]*((2/6)*taopdata[idx(i-3, j, Nx, Ny, 2)]-(7/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]+(11/6)*taopdata[idx(i-1, j, Nx, Ny, 2)])+w1_px[2]*((-1/6)*taopdata[idx(i-2, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 2)])+w2_px[2]*((2/6)*taopdata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 2)]-(1/6)*taopdata[idx(0, j, Nx, Ny, 2)]);
            
            u_tnnhx[2]=w2_nx[2]*((-1/6)*taondata[idx(i-2, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(2/6)*taondata[idx(i, j, Nx, Ny, 2)])+w1_nx[2]*((2/6)*taondata[idx(i-1, j, Nx, Ny, 2)]+(5/6)*taondata[idx(i, j, Nx, Ny, 2)]-(1/6)*taondata[idx(0, j, Nx, Ny, 2)])+w0_nx[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 2)]-(7/6)*taondata[idx(0, j, Nx, Ny, 2)]+(2/6)*taondata[idx(1, j, Nx, Ny, 2)]);
            
            
            u_tpphx[3]=w0_px[3]*((2/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(0, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(0, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(1, j, Nx, Ny, 0)]);
            
            u_tnphx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(0, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(0, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(1, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(0, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(2, j, Nx, Ny, 0)]);
            
            u_tpnhx[3]=w0_px[3]*((2/6)*Cjpdata[idx(i-3, j, Nx, Ny, 0)]-(7/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]+(11/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)])+w1_px[3]*((-1/6)*Cjpdata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 0)])+w2_px[3]*((2/6)*Cjpdata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjpdata[idx(0, j, Nx, Ny, 0)]);
            
            u_tnnhx[3]=w2_nx[3]*((-1/6)*Cjndata[idx(i-2, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 0)])+w1_nx[3]*((2/6)*Cjndata[idx(i-1, j, Nx, Ny, 0)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(1/6)*Cjndata[idx(0, j, Nx, Ny, 0)])+w0_nx[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 0)]-(7/6)*Cjndata[idx(0, j, Nx, Ny, 0)]+(2/6)*Cjndata[idx(1, j, Nx, Ny, 0)]);
            
            for(k=0;k<4;k++){
                yxdata[idx(i, j, Nx, Ny, k)]=-(1/dx)*((u_tpphx[k]-u_tpnhx[k])+(u_tnphx[k]-u_tnnhx[k]));
            }
    }
    
            j=Ny-1;
        for (i=0; i<Nx; i++){
            // x direction
            IS0_py[0]=(13/12)*(ypdata[idx(i, j-2, Nx, Ny, 2)]-2*ypdata[idx(i, j-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, j-2, Nx, Ny, 2)]-2*ypdata[idx(i, j-1, Nx, Ny, 2)]+ypdata[idx(i, j, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, j-2, Nx, Ny, 2)]-4*ypdata[idx(i, j-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)])*(ypdata[idx(i, j-2, Nx, Ny, 2)]-4*ypdata[idx(i, j-1, Nx, Ny, 2)]+3*ypdata[idx(i, j, Nx, Ny, 2)]);
            
            IS1_py[0]=(13/12)*(ypdata[idx(i, j-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, 0, Nx, Ny, 2)])*(ypdata[idx(i, j-1, Nx, Ny, 2)]-2*ypdata[idx(i, j, Nx, Ny, 2)]+ypdata[idx(i, 0, Nx, Ny, 2)])+(1/4)*(ypdata[idx(i, j-1, Nx, Ny, 2)]-1*ypdata[idx(i, 0, Nx, Ny, 2)])*(ypdata[idx(i, j-1, Nx, Ny, 2)]-1*ypdata[idx(i, 0, Nx, Ny, 2)]);
            
            IS2_py[0]=(13/12)*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, 0, Nx, Ny, 2)]+ypdata[idx(i, 1, Nx, Ny, 2)])*(ypdata[idx(i, j, Nx, Ny, 2)]-2*ypdata[idx(i, 0, Nx, Ny, 2)]+ypdata[idx(i, 1, Nx, Ny, 2)])+(1/4)*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, 0, Nx, Ny, 2)]+ypdata[idx(i, 1, Nx, Ny, 2)])*(3*ypdata[idx(i, j, Nx, Ny, 2)]-4*ypdata[idx(i, 0, Nx, Ny, 2)]+ypdata[idx(i, 1, Nx, Ny, 2)]);
            
            IS0_ny[0]=(13/12)*(yndata[idx(i, 0, Nx, Ny, 2)]-2*yndata[idx(i, 1, Nx, Ny, 2)]+yndata[idx(i, 2, Nx, Ny, 2)])*(yndata[idx(i, 0, Nx, Ny, 2)]-2*yndata[idx(i, 1, Nx, Ny, 2)]+yndata[idx(i, 2, Nx, Ny, 2)])+(1/4)*(3*yndata[idx(i, 0, Nx, Ny, 2)]-4*yndata[idx(i, 1, Nx, Ny, 2)]+yndata[idx(i, 2, Nx, Ny, 2)])*(3*yndata[idx(i, 0, Nx, Ny, 2)]-4*yndata[idx(i, 1, Nx, Ny, 2)]+yndata[idx(i, 2, Nx, Ny, 2)]);
            
            IS1_ny[0]=(13/12)*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, 0, Nx, Ny, 2)]+yndata[idx(i, 1, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-2*yndata[idx(i, 0, Nx, Ny, 2)]+yndata[idx(i, 1, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, 1, Nx, Ny, 2)])*(yndata[idx(i, j, Nx, Ny, 2)]-1*yndata[idx(i, 1, Nx, Ny, 2)]);
            
            IS2_ny[0]=(13/12)*(yndata[idx(i, j-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, 0, Nx, Ny, 2)])*(yndata[idx(i, j-1, Nx, Ny, 2)]-2*yndata[idx(i, j, Nx, Ny, 2)]+yndata[idx(i, 0, Nx, Ny, 2)])+(1/4)*(yndata[idx(i, j-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, 0, Nx, Ny, 2)])*(yndata[idx(i, j-1, Nx, Ny, 2)]-4*yndata[idx(i, j, Nx, Ny, 2)]+3*yndata[idx(i, 0, Nx, Ny, 2)]);
            
            
            IS0_py[1]=(13/12)*(taopdata[idx(i, j-2, Nx, Ny, 1)]-2*taopdata[idx(i, j-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, j-2, Nx, Ny, 1)]-2*taopdata[idx(i, j-1, Nx, Ny, 1)]+taopdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, j-2, Nx, Ny, 1)]-4*taopdata[idx(i, j-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)])*(taopdata[idx(i, j-2, Nx, Ny, 1)]-4*taopdata[idx(i, j-1, Nx, Ny, 1)]+3*taopdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[1]=(13/12)*(taopdata[idx(i, j-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, 0, Nx, Ny, 1)])*(taopdata[idx(i, j-1, Nx, Ny, 1)]-2*taopdata[idx(i, j, Nx, Ny, 1)]+taopdata[idx(i, 0, Nx, Ny, 1)])+(1/4)*(taopdata[idx(i, j-1, Nx, Ny, 1)]-1*taopdata[idx(i, 0, Nx, Ny, 1)])*(taopdata[idx(i, j-1, Nx, Ny, 1)]-1*taopdata[idx(i, 0, Nx, Ny, 1)]);
            
            IS2_py[1]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, 0, Nx, Ny, 1)]+taopdata[idx(i, 1, Nx, Ny, 1)])*(taopdata[idx(i, j, Nx, Ny, 1)]-2*taopdata[idx(i, 0, Nx, Ny, 1)]+taopdata[idx(i, 1, Nx, Ny, 1)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, 0, Nx, Ny, 1)]+taopdata[idx(i, 1, Nx, Ny, 1)])*(3*taopdata[idx(i, j, Nx, Ny, 1)]-4*taopdata[idx(i, 0, Nx, Ny, 1)]+taopdata[idx(i, 1, Nx, Ny, 1)]);
            
            IS0_ny[1]=(13/12)*(taondata[idx(i, 0, Nx, Ny, 1)]-2*taondata[idx(i, 1, Nx, Ny, 1)]+taondata[idx(i, 2, Nx, Ny, 1)])*(taondata[idx(i, 0, Nx, Ny, 1)]-2*taondata[idx(i, 1, Nx, Ny, 1)]+taondata[idx(i, 2, Nx, Ny, 1)])+(1/4)*(3*taondata[idx(i, 0, Nx, Ny, 1)]-4*taondata[idx(i, 1, Nx, Ny, 1)]+taondata[idx(i, 2, Nx, Ny, 1)])*(3*taondata[idx(i, 0, Nx, Ny, 1)]-4*taondata[idx(i, 1, Nx, Ny, 1)]+taondata[idx(i, 2, Nx, Ny, 1)]);
            
            IS1_ny[1]=(13/12)*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, 0, Nx, Ny, 1)]+taondata[idx(i, 1, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-2*taondata[idx(i, 0, Nx, Ny, 1)]+taondata[idx(i, 1, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, 1, Nx, Ny, 1)])*(taondata[idx(i, j, Nx, Ny, 1)]-1*taondata[idx(i, 1, Nx, Ny, 1)]);
            
            IS2_ny[1]=(13/12)*(taondata[idx(i, j-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, 0, Nx, Ny, 1)])*(taondata[idx(i, j-1, Nx, Ny, 1)]-2*taondata[idx(i, j, Nx, Ny, 1)]+taondata[idx(i, 0, Nx, Ny, 1)])+(1/4)*(taondata[idx(i, j-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, 0, Nx, Ny, 1)])*(taondata[idx(i, j-1, Nx, Ny, 1)]-4*taondata[idx(i, j, Nx, Ny, 1)]+3*taondata[idx(i, 0, Nx, Ny, 1)]);
            
            
            IS0_py[2]=(13/12)*(taopdata[idx(i, j-2, Nx, Ny, 3)]-2*taopdata[idx(i, j-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, j-2, Nx, Ny, 3)]-2*taopdata[idx(i, j-1, Nx, Ny, 3)]+taopdata[idx(i, j, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, j-2, Nx, Ny, 3)]-4*taopdata[idx(i, j-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)])*(taopdata[idx(i, j-2, Nx, Ny, 3)]-4*taopdata[idx(i, j-1, Nx, Ny, 3)]+3*taopdata[idx(i, j, Nx, Ny, 3)]);
            
            IS1_py[2]=(13/12)*(taopdata[idx(i, j-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, 0, Nx, Ny, 3)])*(taopdata[idx(i, j-1, Nx, Ny, 3)]-2*taopdata[idx(i, j, Nx, Ny, 3)]+taopdata[idx(i, 0, Nx, Ny, 3)])+(1/4)*(taopdata[idx(i, j-1, Nx, Ny, 3)]-1*taopdata[idx(i, 0, Nx, Ny, 3)])*(taopdata[idx(i, j-1, Nx, Ny, 3)]-1*taopdata[idx(i, 0, Nx, Ny, 3)]);
            
            IS2_py[2]=(13/12)*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, 0, Nx, Ny, 3)]+taopdata[idx(i, 1, Nx, Ny, 3)])*(taopdata[idx(i, j, Nx, Ny, 3)]-2*taopdata[idx(i, 0, Nx, Ny, 3)]+taopdata[idx(i, 1, Nx, Ny, 3)])+(1/4)*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, 0, Nx, Ny, 3)]+taopdata[idx(i, 1, Nx, Ny, 3)])*(3*taopdata[idx(i, j, Nx, Ny, 3)]-4*taopdata[idx(i, 0, Nx, Ny, 3)]+taopdata[idx(i, 1, Nx, Ny, 3)]);
            
            IS0_ny[2]=(13/12)*(taondata[idx(i, 0, Nx, Ny, 3)]-2*taondata[idx(i, 1, Nx, Ny, 3)]+taondata[idx(i, 2, Nx, Ny, 3)])*(taondata[idx(i, 0, Nx, Ny, 3)]-2*taondata[idx(i, 1, Nx, Ny, 3)]+taondata[idx(i, 2, Nx, Ny, 3)])+(1/4)*(3*taondata[idx(i, 0, Nx, Ny, 3)]-4*taondata[idx(i, 1, Nx, Ny, 3)]+taondata[idx(i, 2, Nx, Ny, 3)])*(3*taondata[idx(i, 0, Nx, Ny, 3)]-4*taondata[idx(i, 1, Nx, Ny, 3)]+taondata[idx(i, 2, Nx, Ny, 3)]);
            
            IS1_ny[2]=(13/12)*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, 0, Nx, Ny, 3)]+taondata[idx(i, 1, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-2*taondata[idx(i, 0, Nx, Ny, 3)]+taondata[idx(i, 1, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, 1, Nx, Ny, 3)])*(taondata[idx(i, j, Nx, Ny, 3)]-1*taondata[idx(i, 1, Nx, Ny, 3)]);
            
            IS2_ny[2]=(13/12)*(taondata[idx(i, j-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, 0, Nx, Ny, 3)])*(taondata[idx(i, j-1, Nx, Ny, 3)]-2*taondata[idx(i, j, Nx, Ny, 3)]+taondata[idx(i, 0, Nx, Ny, 3)])+(1/4)*(taondata[idx(i, j-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, 0, Nx, Ny, 3)])*(taondata[idx(i, j-1, Nx, Ny, 3)]-4*taondata[idx(i, j, Nx, Ny, 3)]+3*taondata[idx(i, 0, Nx, Ny, 3)]);
            
            
            IS0_py[3]=(13/12)*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-2*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-2*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+Cjpdata[idx(i, j, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-4*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)])*(Cjpdata[idx(i, j-2, Nx, Ny, 1)]-4*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+3*Cjpdata[idx(i, j, Nx, Ny, 1)]);
            
            IS1_py[3]=(13/12)*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, 0, Nx, Ny, 1)])*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-2*Cjpdata[idx(i, j, Nx, Ny, 1)]+Cjpdata[idx(i, 0, Nx, Ny, 1)])+(1/4)*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, 0, Nx, Ny, 1)])*(Cjpdata[idx(i, j-1, Nx, Ny, 1)]-1*Cjpdata[idx(i, 0, Nx, Ny, 1)]);
            
            IS2_py[3]=(13/12)*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, 0, Nx, Ny, 1)]+Cjpdata[idx(i, 1, Nx, Ny, 1)])*(Cjpdata[idx(i, j, Nx, Ny, 1)]-2*Cjpdata[idx(i, 0, Nx, Ny, 1)]+Cjpdata[idx(i, 1, Nx, Ny, 1)])+(1/4)*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, 0, Nx, Ny, 1)]+Cjpdata[idx(i, 1, Nx, Ny, 1)])*(3*Cjpdata[idx(i, j, Nx, Ny, 1)]-4*Cjpdata[idx(i, 0, Nx, Ny, 1)]+Cjpdata[idx(i, 1, Nx, Ny, 1)]);
            
            IS0_ny[3]=(13/12)*(Cjndata[idx(i, 0, Nx, Ny, 1)]-2*Cjndata[idx(i, 1, Nx, Ny, 1)]+Cjndata[idx(i, 2, Nx, Ny, 1)])*(Cjndata[idx(i, 0, Nx, Ny, 1)]-2*Cjndata[idx(i, 1, Nx, Ny, 1)]+Cjndata[idx(i, 2, Nx, Ny, 1)])+(1/4)*(3*Cjndata[idx(i, 0, Nx, Ny, 1)]-4*Cjndata[idx(i, 1, Nx, Ny, 1)]+Cjndata[idx(i, 2, Nx, Ny, 1)])*(3*Cjndata[idx(i, 0, Nx, Ny, 1)]-4*Cjndata[idx(i, 1, Nx, Ny, 1)]+Cjndata[idx(i, 2, Nx, Ny, 1)]);
            
            IS1_ny[3]=(13/12)*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, 0, Nx, Ny, 1)]+Cjndata[idx(i, 1, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-2*Cjndata[idx(i, 0, Nx, Ny, 1)]+Cjndata[idx(i, 1, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, 1, Nx, Ny, 1)])*(Cjndata[idx(i, j, Nx, Ny, 1)]-1*Cjndata[idx(i, 1, Nx, Ny, 1)]);
            
            IS2_ny[3]=(13/12)*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, 0, Nx, Ny, 1)])*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-2*Cjndata[idx(i, j, Nx, Ny, 1)]+Cjndata[idx(i, 0, Nx, Ny, 1)])+(1/4)*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, 0, Nx, Ny, 1)])*(Cjndata[idx(i, j-1, Nx, Ny, 1)]-4*Cjndata[idx(i, j, Nx, Ny, 1)]+3*Cjndata[idx(i, 0, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++)
            {
                Epsilon=0.000001;
                alpha_0py[k]=(1/10)*(1/(Epsilon+IS0_py[k]))*(1/(Epsilon+IS0_py[k]));
                alpha_1py[k]=(6/10)*(1/(Epsilon+IS1_py[k]))*(1/(Epsilon+IS1_py[k]));
                alpha_2py[k]=(3/10)*(1/(Epsilon+IS2_py[k]))*(1/(Epsilon+IS2_py[k]));
                alpha_0ny[k]=(1/10)*(1/(Epsilon+IS0_ny[k]))*(1/(Epsilon+IS0_ny[k]));
                alpha_1ny[k]=(6/10)*(1/(Epsilon+IS1_ny[k]))*(1/(Epsilon+IS1_ny[k]));
                alpha_2ny[k]=(3/10)*(1/(Epsilon+IS2_ny[k]))*(1/(Epsilon+IS2_ny[k]));
                
                w0_py[k]=alpha_0py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w1_py[k]=alpha_1py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w2_py[k]=alpha_2py[k]/(alpha_0py[k]+alpha_1py[k]+alpha_2py[k]);
                w0_ny[k]=alpha_0ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w1_ny[k]=alpha_1ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
                w2_ny[k]=alpha_2ny[k]/(alpha_0ny[k]+alpha_1ny[k]+alpha_2ny[k]);
            }
            
            u_tpphy[0]=w0_py[0]*((2/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, 0, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, 0, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, 1, Nx, Ny, 2)]);
            
            u_tnphy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]+(2/6)*yndata[idx(i, 0, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j, Nx, Ny, 2)]+(5/6)*yndata[idx(i, 0, Nx, Ny, 2)]-(1/6)*yndata[idx(i, 1, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, 0, Nx, Ny, 2)]-(7/6)*yndata[idx(i, 1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, 2, Nx, Ny, 2)]);
            
            u_tpnhy[0]=w0_py[0]*((2/6)*ypdata[idx(i, j-3, Nx, Ny, 2)]-(7/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]+(11/6)*ypdata[idx(i, j-1, Nx, Ny, 2)])+w1_py[0]*((-1/6)*ypdata[idx(i, j-2, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(2/6)*ypdata[idx(i, j, Nx, Ny, 2)])+w2_py[0]*((2/6)*ypdata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*ypdata[idx(i, j, Nx, Ny, 2)]-(1/6)*ypdata[idx(i, 0, Nx, Ny, 2)]);
            
            u_tnnhy[0]=w2_ny[0]*((-1/6)*yndata[idx(i, j-2, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(2/6)*yndata[idx(i, j, Nx, Ny, 2)])+w1_ny[0]*((2/6)*yndata[idx(i, j-1, Nx, Ny, 2)]+(5/6)*yndata[idx(i, j, Nx, Ny, 2)]-(1/6)*yndata[idx(i, 0, Nx, Ny, 2)])+w0_ny[0]*((11/6)*yndata[idx(i, j, Nx, Ny, 2)]-(7/6)*yndata[idx(i, 0, Nx, Ny, 2)]+(2/6)*yndata[idx(i, 1, Nx, Ny, 2)]);
            
            
            u_tpphy[1]=w0_py[1]*((2/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, 0, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, 0, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, 1, Nx, Ny, 1)]);
            
            u_tnphy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]+(2/6)*taondata[idx(i, 0, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j, Nx, Ny, 1)]+(5/6)*taondata[idx(i, 0, Nx, Ny, 1)]-(1/6)*taondata[idx(i, 1, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, 0, Nx, Ny, 1)]-(7/6)*taondata[idx(i, 1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, 2, Nx, Ny, 1)]);
            
            u_tpnhy[1]=w0_py[1]*((2/6)*taopdata[idx(i, j-3, Nx, Ny, 1)]-(7/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]+(11/6)*taopdata[idx(i, j-1, Nx, Ny, 1)])+w1_py[1]*((-1/6)*taopdata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 1)])+w2_py[1]*((2/6)*taopdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 1)]-(1/6)*taopdata[idx(i, 0, Nx, Ny, 1)]);
            
            u_tnnhy[1]=w2_ny[1]*((-1/6)*taondata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*taondata[idx(i, j, Nx, Ny, 1)])+w1_ny[1]*((2/6)*taondata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*taondata[idx(i, j, Nx, Ny, 1)]-(1/6)*taondata[idx(i, 0, Nx, Ny, 1)])+w0_ny[1]*((11/6)*taondata[idx(i, j, Nx, Ny, 1)]-(7/6)*taondata[idx(i, 0, Nx, Ny, 1)]+(2/6)*taondata[idx(i, 1, Nx, Ny, 1)]);
            
            
            u_tpphy[2]=w0_py[2]*((2/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, 0, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, 0, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, 1, Nx, Ny, 3)]);
            
            u_tnphy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]+(2/6)*taondata[idx(i, 0, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j, Nx, Ny, 3)]+(5/6)*taondata[idx(i, 0, Nx, Ny, 3)]-(1/6)*taondata[idx(i, 1, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, 0, Nx, Ny, 3)]-(7/6)*taondata[idx(i, 1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, 2, Nx, Ny, 3)]);
            
            u_tpnhy[2]=w0_py[2]*((2/6)*taopdata[idx(i, j-3, Nx, Ny, 3)]-(7/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]+(11/6)*taopdata[idx(i, j-1, Nx, Ny, 3)])+w1_py[2]*((-1/6)*taopdata[idx(i, j-2, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(2/6)*taopdata[idx(i, j, Nx, Ny, 3)])+w2_py[2]*((2/6)*taopdata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taopdata[idx(i, j, Nx, Ny, 3)]-(1/6)*taopdata[idx(i, 0, Nx, Ny, 3)]);
            
            u_tnnhy[2]=w2_ny[2]*((-1/6)*taondata[idx(i, j-2, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(2/6)*taondata[idx(i, j, Nx, Ny, 3)])+w1_ny[2]*((2/6)*taondata[idx(i, j-1, Nx, Ny, 3)]+(5/6)*taondata[idx(i, j, Nx, Ny, 3)]-(1/6)*taondata[idx(i, 0, Nx, Ny, 3)])+w0_ny[2]*((11/6)*taondata[idx(i, j, Nx, Ny, 3)]-(7/6)*taondata[idx(i, 0, Nx, Ny, 3)]+(2/6)*taondata[idx(i, 1, Nx, Ny, 3)]);
            
            
            u_tpphy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, 0, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, 0, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, 1, Nx, Ny, 1)]);
            
            u_tnphy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, 0, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, 0, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, 1, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, 0, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, 1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, 2, Nx, Ny, 1)]);
            
            u_tpnhy[3]=w0_py[3]*((2/6)*Cjpdata[idx(i, j-3, Nx, Ny, 1)]-(7/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]+(11/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)])+w1_py[3]*((-1/6)*Cjpdata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*Cjpdata[idx(i, j, Nx, Ny, 1)])+w2_py[3]*((2/6)*Cjpdata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjpdata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjpdata[idx(i, 0, Nx, Ny, 1)]);
            
            u_tnnhy[3]=w2_ny[3]*((-1/6)*Cjndata[idx(i, j-2, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, j, Nx, Ny, 1)])+w1_ny[3]*((2/6)*Cjndata[idx(i, j-1, Nx, Ny, 1)]+(5/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(1/6)*Cjndata[idx(i, 0, Nx, Ny, 1)])+w0_ny[3]*((11/6)*Cjndata[idx(i, j, Nx, Ny, 1)]-(7/6)*Cjndata[idx(i, 0, Nx, Ny, 1)]+(2/6)*Cjndata[idx(i, 1, Nx, Ny, 1)]);
            
            for(k=0;k<4;k++){
                yydata[idx(i, j, Nx, Ny, k)]=-(1/dy)*((u_tpphy[k]-u_tpnhy[k])+(u_tnphy[k]-u_tnnhy[k]));
            }
        }

    
    
    
    for(j=0; j<Ny; j++){
        for (i=0; i<Nx; i++){
            for (k=0; k<4; k++){
                dYdata[idx(i, j, Nx, Ny, k)]=yxdata[idx(i, j, Nx, Ny, k)]+yydata[idx(i, j, Nx, Ny, k)];
            }
        }
    }
    
   
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


static int Gettao(N_Vector y, N_Vector tao, long int Nx, long int Ny)
{
    long int k,i,j, NEQ;
    NEQ = Nx*Ny;
    N_Vector vx = NULL;
    N_Vector vy = NULL;
    
    vx = N_VNew_Serial(NEQ);
    if (check_flag((void *) vx, "N_VNew_Serial", 0)) return 1;
    vy = N_VNew_Serial(NEQ);
    if (check_flag((void *) vy, "N_VNew_Serial", 0)) return 1;
    
    vx_data = N_VGetArrayPointer(vx);
    if (check_flag((void *) vx_data, "N_VGetArrayPointer", 0)) return 1;
    vy_data = N_VGetArrayPointer(vy);
    if (check_flag((void *) vy_data, "N_VGetArrayPointer", 0)) return 1;
    data = N_VGetArrayPointer(y);
    if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
    tao_data = N_VGetArrayPointer(tao);
    if (check_flag((void *) tao_data, "N_VGetArrayPointer", 0)) return 1;
    
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            vx_data[idx_v(i,j,Nx)]=data[idx(i, j, Nx, Ny, 1)]/data[idx(i, j, Nx, Ny, 0)];
            vy_data[idx_v(i,j,Nx)]=data[idx(i, j, Nx, Ny, 2)]/data[idx(i, j, Nx, Ny, 0)];
        }
    }
    //    for (k=0;k<4;k++){
    
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            tao_data[idx(i, j, Nx, Ny, 0)] = data[idx(i, j, Nx, Ny, 1)]*vx_data[idx_v(i,j,Nx)]+data[idx(i, j, Nx, Ny, 0)];
            tao_data[idx(i, j, Nx, Ny, 1)] = data[idx(i, j, Nx, Ny, 1)]*vy_data[idx_v(i,j,Nx)];
            tao_data[idx(i, j, Nx, Ny, 2)] = data[idx(i, j, Nx, Ny, 2)]*vx_data[idx_v(i,j,Nx)];
            tao_data[idx(i, j, Nx, Ny, 3)] = data[idx(i, j, Nx, Ny, 2)]*vy_data[idx_v(i,j,Nx)]+data[idx(i, j, Nx, Ny, 0)];
            //tao_xx[idx_v(i,j,Nx)]=qx[idx_v(i,j,Nx)]*vx[idx_v(i,j,Nx)]+rou[idx_v(i, j, Nx)];
            //tao_xy[idx_v(i,j,Nx)]=qx[idx_v(i,j,Nx)]*vy[idx_v(i,j,Nx)];
            //tao_yx[idx_v(i,j,Nx)]=qy[idx_v(i,j,Nx)]*vx[idx_v(i,j,Nx)];
            //tao_yy[idx_v(i,j,Nx)]=qy[idx_v(i,j,Nx)]*vy[idx_v(i,j,Nx)]+rou[idx_v(i, j, Nx)];
        }
    }
    
    /* Free vectors */
    N_VDestroy_Serial(vx);
    N_VDestroy_Serial(vy);
    
    return 0;
    //    }

}

static int GetCj(N_Vector y, N_Vector Cj, long int Nx, long int Ny)
{
    long int k,i,j, NEQ;
    NEQ = Nx*Ny;
    N_Vector vx = NULL;
    N_Vector vy = NULL;
    
    vx = N_VNew_Serial(NEQ);
    if (check_flag((void *) vx, "N_VNew_Serial", 0)) return 1;
    vy = N_VNew_Serial(NEQ);
    if (check_flag((void *) vy, "N_VNew_Serial", 0)) return 1;
    
    vx_data = N_VGetArrayPointer(vx);
    if (check_flag((void *) vx_data, "N_VGetArrayPointer", 0)) return 1;
    vy_data = N_VGetArrayPointer(vy);
    if (check_flag((void *) vy_data, "N_VGetArrayPointer", 0)) return 1;
    data = N_VGetArrayPointer(y);
    if (check_flag((void *) data, "N_VGetArrayPointer", 0)) return 1;
    Cj_data = N_VGetArrayPointer(Cj);
    if (check_flag((void *) Cj_data, "N_VGetArrayPointer", 0)) return 1;

    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            vx_data[idx_v(i,j,Nx)]=data[idx(i, j, Nx, Ny, 1)]/data[idx(i, j, Nx, Ny, 0)];
            vy_data[idx_v(i,j,Nx)]=data[idx(i, j, Nx, Ny, 2)]/data[idx(i, j, Nx, Ny, 0)];
        }
    }    //    for (k=0;k<4;k++){
    
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            Cj_data[idx(i, j, Nx, Ny, 0)] = data[idx(i, j, Nx, Ny, 1)]*data[idx(i, j, Nx, Ny, 1)]*vx_data[idx_v(i,j,Nx)]+data[idx(i, j, Nx, Ny, 0)]*data[idx(i, j, Nx, Ny, 1)]+data[idx(i, j, Nx, Ny, 1)]*data[idx(i, j, Nx, Ny, 2)]*vy_data[idx_v(i,j,Nx)];
            Cj_data[idx(i, j, Nx, Ny, 1)] = data[idx(i, j, Nx, Ny, 2)]*data[idx(i, j, Nx, Ny, 2)]*vy_data[idx_v(i,j,Nx)]+data[idx(i, j, Nx, Ny, 0)]*data[idx(i, j, Nx, Ny, 2)]+data[idx(i, j, Nx, Ny, 1)]*data[idx(i, j, Nx, Ny, 2)]*vx_data[idx_v(i,j,Nx)];
            
        }
    }

    //for(j=0;j<Ny;j++){
      //  for(i=0;i<Nx;i++){
        //    Cj_x[idx_v(i,j,Nx)]=qx[idx_v(i,j,Nx)]*qx[idx_v(i,j,Nx)]*vx[idx_v(i,j,Nx)]+rou[idx_v(i, j, Nx)]*qx[idx_v(i,j,Nx)]+qx[idx_v(i,j,Nx)]*qy[idx_v(i,j,Nx)]*vy[idx_v(i,j,Nx)];
         //   Cj_y[idx_v(i,j,Nx)]=qy[idx_v(i,j,Nx)]*qy[idx_v(i,j,Nx)]*vy[idx_v(i,j,Nx)]+rou[idx_v(i, j, Nx)]*qy[idx_v(i,j,Nx)]+qx[idx_v(i,j,Nx)]*qy[idx_v(i,j,Nx)]*vx[idx_v(i,j,Nx)];
       // }
  //  }
    /* Free vectors */
    N_VDestroy_Serial(vx);
    N_VDestroy_Serial(vy);
    
    return 0;
    //    }
}


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








