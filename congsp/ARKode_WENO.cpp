// should I do it like heat2D.cpp? I mean use N_Vector and ark_mem, something like this? function init_from_file? how to test my code? Do you have an example? how to get the max absolute eigenvalue?
// can't compile and hg pull has to merge branches，type emacs in terminal X11 doesn't pop up

#include <iostream>
#include <sundials/sundials_types.h>


#define idx(i,j,Nx,Ny,k) ((k)*(Nx)*(Ny)+(j)*(Nx)+i)
#define idx_v(i,j,Nx) ((j)*(Nx)+i)

// what is the input u exactly? I do it as below, is it ok?
int main(int argc, const char * argv[])
{
    int k,i,j;
    // Can I use realtype like this?
    realtype *rou = new realtype[Nx*Ny];
    realtype *qx = new realtype[Nx*Ny];
    realtype *qy = new realtype[Nx*Ny];
    realtype *E = new realtype[Nx*Ny];
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            rou[idx_v(i, j, Nx)]=u[idx(i, j, Nx, Ny, 0)];
            qx[idx_v(i, j, Nx)]=u[idx(i, j, Nx, Ny, 1)];
            qy[idx_v(i, j, Nx)]=u[idx(i, j, Nx, Ny, 2)];
            E[idx_v(i, j, Nx)]=u[idx(i, j, Nx, Ny, 3)];
        }
    }
// how to get Nx and Ny? divide tao to 4 different short vectors, is it ok?
    realtype *tao_xx = new realtype[Nx*Ny];
    realtype *tao_xy = new realtype[Nx*Ny];
    realtype *tao_yx = new realtype[Nx*Ny];
    realtype *tao_yy = new realtype[Nx*Ny];
    realtype *Cj_x = new realtype[Nx*Ny];
    realtype *Cj_y = new realtype[Nx*Ny];
    
    Gettao(rou, qx, qy, E, Nx, Ny, tao_xx, tao_xy, tao_yx, tao_yy);
    GetCj(rou, qx, qy, E, Nx, Ny, Cj_x, Cj_y);
    
    realtype *F = new realtype[4*Nx*Ny];
    
//    for (k=0;k<4;k++){
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            F[idx(i,j,Nx,Ny,0)]=WENO(i,j,0,Nx,Ny,dx, qx,u,0)+WENO(i,j,0,Nx,Ny,dy, qy,u,1);
        }
    }
//    }
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            F[idx(i,j,Nx,Ny,1)]=WENO(i,j,1,Nx,Ny,dx, tao_xx,u,0)+WENO(i,j,1,Nx,Ny,dy, tao_xy,u,1);
        }
    }

    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            F[idx(i,j,Nx,Ny,2)]=WENO(i,j,2,Nx,Ny,dx, tao_yx,u,0)+WENO(i,j,2,Nx,Ny,dy, tao_yy,u,1);
        }
    }

    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            F[idx(i,j,Nx,Ny,3)]=WENO(i,j,3,Nx,Ny,dx, Cj_x,u,0)+WENO(i,j,3,Nx,Ny,dy, Cj_y,u,1);
        }
    }
    
    delete []rou;
    delete []qx;
    delete []qy;
    delete []E;
    delete []tao_xx;
    delete []tao_xy;
    delete []tao_yx;
    delete []tao_yy;
    delete []Cj_x;
    delete []Cj_y;
    delete []F;
    
}

static realtype WENO(int i, int j, int k, int Nx, int Ny, realtype dx, realtype *f, realtype *u, int flag){
    realtype result;
    realtype IS0_p,IS1_p,IS2_p,IS0_n,IS1_n,IS2_n,Epsilon,alpha_0p,alpha_1p,alpha_2p,alpha_0n,alpha_1n,alpha_2n,w0_p,w1_p,w2_p,w0_n,w1_n,w2_n,u_tpph,u_tpnh,u_tnph,u_tnnh;
    int m,n;
    realtype f_p = new realtype [Nx*Ny];
    realtype f_n = new realtype [Nx*Ny];
    for (n=0;n<Ny;n++){
        for (m=0;m<Nx;m++){
            f_p[idx_v(m, n, Nx)]=0.5*(f[idx_v(m, n, Nx)]+abs(max(d_f[idx_v(m, n, Nx)]))*u[idx(m, n, Nx, Ny, k)]);
            f_n[idx_v(m, n, Nx)]=0.5*(f[idx_v(m, n, Nx)]-abs(max(d_f[idx_v(m, n, Nx)]))*u[idx(m, n, Nx, Ny, k)]);
        }
    }
    if (flag==0){
        if(i==0){
            IS0_p=(13/12)*(f_p[idx_v(Nx-2,j,Nx)]-2*f_p[idx_v(Nx-1,j,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(Nx-2,j,Nx)]-4*f_p[idx_v(Nx-1,j,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(Nx-1,j,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i+1,j,Nx)])^2+(1/4)*(f_p[idx_v(Nx-1,j,Nx)]-1*f_p[idx_v(i+1,j,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i+1,j,Nx)]+f_p[idx_v(i+2,j,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i+1,j,Nx)]+f_p[idx_v(i+2,j,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(i+1,j,Nx)]-2*f_n[idx_v(i+2,j,Nx)]+f_n[idx_v(i+3,j,Nx)])^2+(1/4)*(3*f_n[idx_v(i+1,j,Nx)]-4*f_n[idx_v(i+2,j,Nx)]+f_n[idx_v(i+3,j,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i+1,j,Nx)]+f_n[idx_v(i+2,j,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(i+2,j,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(Nx-1,j,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i+1,j,Nx)])^2+(1/4)*(f_n[idx_v(Nx-1,j,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i+1,j,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(Nx-2,j,Nx)]-(7/6)*f_p[idx_v(Nx-1,j,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(Nx-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i+1,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i+1,j,Nx)]-(1/6)*f_p[idx_v(i+2,j,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(Nx-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i+1,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i+1,j,Nx)]-(1/6)*f_n[idx_v(i+2,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i+1,j,Nx)]-(7/6)*f_n[idx_v(i+2,j,Nx)]+(2/6)*f_n[idx_v(i+3,j,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(Nx-3,j,Nx)]-(7/6)*f_p[idx_v(Nx-2,j,Nx)]+(11/6)*f_p[idx_v(Nx-1,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(Nx-2,j,Nx)]+(5/6)*f_p[idx_v(Nx-1,j,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(Nx-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i+1,j,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(Nx-2,j,Nx)]+(5/6)*f_n[idx_v(Nx-1,j,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(Nx-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i+1,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i+1,j,Nx)]+(2/6)*f_n[idx_v(i+2,j,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }
        
        if(i==1){
            IS0_p=(13/12)*(f_p[idx_v(Nx-1,j,Nx)]-2*f_p[idx_v(i-1,j,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(Nx-1,j,Nx)]-4*f_p[idx_v(i-1,j,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(i-1,j,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i+1,j,Nx)])^2+(1/4)*(f_p[idx_v(i-1,j,Nx)]-1*f_p[idx_v(i+1,j,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i+1,j,Nx)]+f_p[idx_v(i+2,j,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i+1,j,Nx)]+f_p[idx_v(i+2,j,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(i+1,j,Nx)]-2*f_n[idx_v(i+2,j,Nx)]+f_n[idx_v(i+3,j,Nx)])^2+(1/4)*(3*f_n[idx_v(i+1,j,Nx)]-4*f_n[idx_v(i+2,j,Nx)]+f_n[idx_v(i+3,j,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i+1,j,Nx)]+f_n[idx_v(i+2,j,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(i+2,j,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(i-1,j,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i+1,j,Nx)])^2+(1/4)*(f_n[idx_v(i-1,j,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i+1,j,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(Nx-1,j,Nx)]-(7/6)*f_p[idx_v(i-1,j,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i+1,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i+1,j,Nx)]-(1/6)*f_p[idx_v(i+2,j,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i+1,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i+1,j,Nx)]-(1/6)*f_n[idx_v(i+2,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i+1,j,Nx)]-(7/6)*f_n[idx_v(i+2,j,Nx)]+(2/6)*f_n[idx_v(i+3,j,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(Nx-2,j,Nx)]-(7/6)*f_p[idx_v(Nx-1,j,Nx)]+(11/6)*f_p[idx_v(i-1,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(Nx-1,j,Nx)]+(5/6)*f_p[idx_v(i-1,j,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i+1,j,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(Nx-1,j,Nx)]+(5/6)*f_n[idx_v(i-1,j,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i+1,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i+1,j,Nx)]+(2/6)*f_n[idx_v(i+2,j,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }

        if(i==2){
            IS0_p=(13/12)*(f_p[idx_v(i-2,j,Nx)]-2*f_p[idx_v(i-1,j,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(i-2,j,Nx)]-4*f_p[idx_v(i-1,j,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(i-1,j,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i+1,j,Nx)])^2+(1/4)*(f_p[idx_v(i-1,j,Nx)]-1*f_p[idx_v(i+1,j,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i+1,j,Nx)]+f_p[idx_v(i+2,j,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i+1,j,Nx)]+f_p[idx_v(i+2,j,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(i+1,j,Nx)]-2*f_n[idx_v(i+2,j,Nx)]+f_n[idx_v(i+3,j,Nx)])^2+(1/4)*(3*f_n[idx_v(i+1,j,Nx)]-4*f_n[idx_v(i+2,j,Nx)]+f_n[idx_v(i+3,j,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i+1,j,Nx)]+f_n[idx_v(i+2,j,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(i+2,j,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(i-1,j,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i+1,j,Nx)])^2+(1/4)*(f_n[idx_v(i-1,j,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i+1,j,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(i-2,j,Nx)]-(7/6)*f_p[idx_v(i-1,j,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i+1,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i+1,j,Nx)]-(1/6)*f_p[idx_v(i+2,j,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i+1,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i+1,j,Nx)]-(1/6)*f_n[idx_v(i+2,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i+1,j,Nx)]-(7/6)*f_n[idx_v(i+2,j,Nx)]+(2/6)*f_n[idx_v(i+3,j,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(Nx-1,j,Nx)]-(7/6)*f_p[idx_v(i-2,j,Nx)]+(11/6)*f_p[idx_v(i-1,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i-2,j,Nx)]+(5/6)*f_p[idx_v(i-1,j,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i+1,j,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(i-2,j,Nx)]+(5/6)*f_n[idx_v(i-1,j,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i+1,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i+1,j,Nx)]+(2/6)*f_n[idx_v(i+2,j,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }

        if(i==Nx-1){
            IS0_p=(13/12)*(f_p[idx_v(i-2,j,Nx)]-2*f_p[idx_v(i-1,j,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(i-2,j,Nx)]-4*f_p[idx_v(i-1,j,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(i-1,j,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(0,j,Nx)])^2+(1/4)*(f_p[idx_v(i-1,j,Nx)]-1*f_p[idx_v(0,j,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(0,j,Nx)]+f_p[idx_v(1,j,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(0,j,Nx)]+f_p[idx_v(1,j,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(0,j,Nx)]-2*f_n[idx_v(1,j,Nx)]+f_n[idx_v(2,j,Nx)])^2+(1/4)*(3*f_n[idx_v(0,j,Nx)]-4*f_n[idx_v(1,j,Nx)]+f_n[idx_v(2,j,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(0,j,Nx)]+f_n[idx_v(1,j,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(1,j,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(i-1,j,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(0,j,Nx)])^2+(1/4)*(f_n[idx_v(i-1,j,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(0,j,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(i-2,j,Nx)]-(7/6)*f_p[idx_v(i-1,j,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(0,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(0,j,Nx)]-(1/6)*f_p[idx_v(1,j,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(0,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(0,j,Nx)]-(1/6)*f_n[idx_v(1,j,Nx)])+w0_n*((11/6)*f_n[idx_v(0,j,Nx)]-(7/6)*f_n[idx_v(1,j,Nx)]+(2/6)*f_n[idx_v(2,j,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(i-3,j,Nx)]-(7/6)*f_p[idx_v(i-2,j,Nx)]+(11/6)*f_p[idx_v(i-1,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i-2,j,Nx)]+(5/6)*f_p[idx_v(i-1,j,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(0,j,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(i-2,j,Nx)]+(5/6)*f_n[idx_v(i-1,j,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(0,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(0,j,Nx)]+(2/6)*f_n[idx_v(1,j,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }

        if(i==Nx-2){
            IS0_p=(13/12)*(f_p[idx_v(i-2,j,Nx)]-2*f_p[idx_v(i-1,j,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(i-2,j,Nx)]-4*f_p[idx_v(i-1,j,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(i-1,j,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i+1,j,Nx)])^2+(1/4)*(f_p[idx_v(i-1,j,Nx)]-1*f_p[idx_v(i+1,j,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i+1,j,Nx)]+f_p[idx_v(0,j,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i+1,j,Nx)]+f_p[idx_v(0,j,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(i+1,j,Nx)]-2*f_n[idx_v(0,j,Nx)]+f_n[idx_v(1,j,Nx)])^2+(1/4)*(3*f_n[idx_v(i+1,j,Nx)]-4*f_n[idx_v(0,j,Nx)]+f_n[idx_v(1,j,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i+1,j,Nx)]+f_n[idx_v(0,j,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(0,j,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(i-1,j,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i+1,j,Nx)])^2+(1/4)*(f_n[idx_v(i-1,j,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i+1,j,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(i-2,j,Nx)]-(7/6)*f_p[idx_v(i-1,j,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i+1,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i+1,j,Nx)]-(1/6)*f_p[idx_v(0,j,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i+1,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i+1,j,Nx)]-(1/6)*f_n[idx_v(0,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i+1,j,Nx)]-(7/6)*f_n[idx_v(0,j,Nx)]+(2/6)*f_n[idx_v(1,j,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(i-3,j,Nx)]-(7/6)*f_p[idx_v(i-2,j,Nx)]+(11/6)*f_p[idx_v(i-1,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i-2,j,Nx)]+(5/6)*f_p[idx_v(i-1,j,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i+1,j,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(i-2,j,Nx)]+(5/6)*f_n[idx_v(i-1,j,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i+1,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i+1,j,Nx)]+(2/6)*f_n[idx_v(0,j,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }

        if(i==Nx-3){
            IS0_p=(13/12)*(f_p[idx_v(i-2,j,Nx)]-2*f_p[idx_v(i-1,j,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(i-2,j,Nx)]-4*f_p[idx_v(i-1,j,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(i-1,j,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i+1,j,Nx)])^2+(1/4)*(f_p[idx_v(i-1,j,Nx)]-1*f_p[idx_v(i+1,j,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i+1,j,Nx)]+f_p[idx_v(i+2,j,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i+1,j,Nx)]+f_p[idx_v(i+2,j,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(i+1,j,Nx)]-2*f_n[idx_v(i+2,j,Nx)]+f_n[idx_v(0,j,Nx)])^2+(1/4)*(3*f_n[idx_v(i+1,j,Nx)]-4*f_n[idx_v(i+2,j,Nx)]+f_n[idx_v(0,j,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i+1,j,Nx)]+f_n[idx_v(i+2,j,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(i+2,j,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(i-1,j,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i+1,j,Nx)])^2+(1/4)*(f_n[idx_v(i-1,j,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i+1,j,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(i-2,j,Nx)]-(7/6)*f_p[idx_v(i-1,j,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i+1,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i+1,j,Nx)]-(1/6)*f_p[idx_v(i+2,j,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i+1,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i+1,j,Nx)]-(1/6)*f_n[idx_v(i+2,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i+1,j,Nx)]-(7/6)*f_n[idx_v(i+2,j,Nx)]+(2/6)*f_n[idx_v(0,j,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(i-3,j,Nx)]-(7/6)*f_p[idx_v(i-2,j,Nx)]+(11/6)*f_p[idx_v(i-1,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i-2,j,Nx)]+(5/6)*f_p[idx_v(i-1,j,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i+1,j,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(i-2,j,Nx)]+(5/6)*f_n[idx_v(i-1,j,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i+1,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i+1,j,Nx)]+(2/6)*f_n[idx_v(i+2,j,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }

        if(i>2&&i<Nx-3){
            IS0_p=(13/12)*(f_p[idx_v(i-2,j,Nx)]-2*f_p[idx_v(i-1,j,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(i-2,j,Nx)]-4*f_p[idx_v(i-1,j,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
    
            IS1_p=(13/12)*(f_p[idx_v(i-1,j,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i+1,j,Nx)])^2+(1/4)*(f_p[idx_v(i-1,j,Nx)]-1*f_p[idx_v(i+1,j,Nx)])^2;
    
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i+1,j,Nx)]+f_p[idx_v(i+2,j,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i+1,j,Nx)]+f_p[idx_v(i+2,j,Nx)])^2;
    
            IS0_n=(13/12)*(f_n[idx_v(i+1,j,Nx)]-2*f_n[idx_v(i+2,j,Nx)]+f_n[idx_v(i+3,j,Nx)])^2+(1/4)*(3*f_n[idx_v(i+1,j,Nx)]-4*f_n[idx_v(i+2,j,Nx)]+f_n[idx_v(i+3,j,Nx)])^2;
    
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i+1,j,Nx)]+f_n[idx_v(i+2,j,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(i+2,j,Nx)])^2;
    
            IS2_n=(13/12)*(f_n[idx_v(i-1,j,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i+1,j,Nx)])^2+(1/4)*(f_n[idx_v(i-1,j,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i+1,j,Nx)])^2;
    
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
    
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
    
            u_tpph=w0_p*((2/6)*f_p[idx_v(i-2,j,Nx)]-(7/6)*f_p[idx_v(i-1,j,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i+1,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i+1,j,Nx)]-(1/6)*f_p[idx_v(i+2,j,Nx)]);
    
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i+1,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i+1,j,Nx)]-(1/6)*f_n[idx_v(i+2,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i+1,j,Nx)]-(7/6)*f_n[idx_v(i+2,j,Nx)]+(2/6)*f_n[idx_v(i+3,j,Nx)]);
    
            u_tpnh=w0_p*((2/6)*f_p[idx_v(i-3,j,Nx)]-(7/6)*f_p[idx_v(i-2,j,Nx)]+(11/6)*f_p[idx_v(i-1,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i-2,j,Nx)]+(5/6)*f_p[idx_v(i-1,j,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i-1,j,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i+1,j,Nx)]);
    
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(i-2,j,Nx)]+(5/6)*f_n[idx_v(i-1,j,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i-1,j,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i+1,j,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i+1,j,Nx)]+(2/6)*f_n[idx_v(i+2,j,Nx)]);
    
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }
    }

    if (flag==1){
        
        if(j==0){
            IS0_p=(13/12)*(f_p[idx_v(i,Ny-2,Nx)]-2*f_p[idx_v(i,Ny-1,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(i,Ny-2,Nx)]-4*f_p[idx_v(i,Ny-1,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(i,Ny-1,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i,j+1,Nx)])^2+(1/4)*(f_p[idx_v(i,Ny-1,Nx)]-1*f_p[idx_v(i,j+1,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i,j+1,Nx)]+f_p[idx_v(i,j+2,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i,j+1,Nx)]+f_p[idx_v(i,j+2,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(i,j+1,Nx)]-2*f_n[idx_v(i,j+2,Nx)]+f_n[idx_v(i,j+3,Nx)])^2+(1/4)*(3*f_n[idx_v(i,j+1,Nx)]-4*f_n[idx_v(i,j+2,Nx)]+f_n[idx_v(i,j+3,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i,j+1,Nx)]+f_n[idx_v(i,j+2,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(i,j+2,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(i,Ny-1,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i,j+1,Nx)])^2+(1/4)*(f_n[idx_v(i,Ny-1,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i,j+1,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(i,Ny-2,Nx)]-(7/6)*f_p[idx_v(i,Ny-1,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,Ny-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i,j+1,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i,j+1,Nx)]-(1/6)*f_p[idx_v(i,j+2,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i,Ny-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i,j+1,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i,j+1,Nx)]-(1/6)*f_n[idx_v(i,j+2,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j+1,Nx)]-(7/6)*f_n[idx_v(i,j+2,Nx)]+(2/6)*f_n[idx_v(i,j+3,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(i,Ny-3,Nx)]-(7/6)*f_p[idx_v(i,Ny-2,Nx)]+(11/6)*f_p[idx_v(i,Ny-1,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,Ny-2,Nx)]+(5/6)*f_p[idx_v(i,Ny-1,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,Ny-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i,j+1,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(i,Ny-2,Nx)]+(5/6)*f_n[idx_v(i,Ny-1,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,Ny-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i,j+1,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i,j+1,Nx)]+(2/6)*f_n[idx_v(i,j+2,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }
        
        if(j==1){
            IS0_p=(13/12)*(f_p[idx_v(i,Ny-1,Nx)]-2*f_p[idx_v(i,j-1,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(i,Ny-1,Nx)]-4*f_p[idx_v(i,j-1,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(i,j-1,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i,j+1,Nx)])^2+(1/4)*(f_p[idx_v(i,j-1,Nx)]-1*f_p[idx_v(i,j+1,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i,j+1,Nx)]+f_p[idx_v(i,j+2,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i,j+1,Nx)]+f_p[idx_v(i,j+2,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(i,j+1,Nx)]-2*f_n[idx_v(i,j+2,Nx)]+f_n[idx_v(i,j+3,Nx)])^2+(1/4)*(3*f_n[idx_v(i,j+1,Nx)]-4*f_n[idx_v(i,j+2,Nx)]+f_n[idx_v(i,j+3,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i,j+1,Nx)]+f_n[idx_v(i,j+2,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(i,j+2,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(i,j-1,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i,j+1,Nx)])^2+(1/4)*(f_n[idx_v(i,j-1,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i,j+1,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(i,Ny-1,Nx)]-(7/6)*f_p[idx_v(i,j-1,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,j-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i,j+1,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i,j+1,Nx)]-(1/6)*f_p[idx_v(i,j+2,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i,j-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i,j+1,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i,j+1,Nx)]-(1/6)*f_n[idx_v(i,j+2,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j+1,Nx)]-(7/6)*f_n[idx_v(i,j+2,Nx)]+(2/6)*f_n[idx_v(i,j+3,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(i,Ny-2,Nx)]-(7/6)*f_p[idx_v(i,Ny-1,Nx)]+(11/6)*f_p[idx_v(i,j-1,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,Ny-1,Nx)]+(5/6)*f_p[idx_v(i,j-1,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i,j+1,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(i,Ny-1,Nx)]+(5/6)*f_n[idx_v(i,j-1,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i,j+1,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i,j+1,Nx)]+(2/6)*f_n[idx_v(i,j+2,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }

        if(i==2){
            IS0_p=(13/12)*(f_p[idx_v(i,j-2,Nx)]-2*f_p[idx_v(i,j-1,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(i,j-2,Nx)]-4*f_p[idx_v(i,j-1,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(i,j-1,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i,j+1,Nx)])^2+(1/4)*(f_p[idx_v(i,j-1,Nx)]-1*f_p[idx_v(i,j+1,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i,j+1,Nx)]+f_p[idx_v(i,j+2,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i,j+1,Nx)]+f_p[idx_v(i,j+2,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(i,j+1,Nx)]-2*f_n[idx_v(i,j+2,Nx)]+f_n[idx_v(i,j+3,Nx)])^2+(1/4)*(3*f_n[idx_v(i,j+1,Nx)]-4*f_n[idx_v(i,j+2,Nx)]+f_n[idx_v(i,j+3,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i,j+1,Nx)]+f_n[idx_v(i,j+2,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(i,j+2,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(i,j-1,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i,j+1,Nx)])^2+(1/4)*(f_n[idx_v(i,j-1,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i,j+1,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(i,j-2,Nx)]-(7/6)*f_p[idx_v(i,j-1,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,j-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i,j+1,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i,j+1,Nx)]-(1/6)*f_p[idx_v(i,j+2,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i,j-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i,j+1,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i,j+1,Nx)]-(1/6)*f_n[idx_v(i,j+2,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j+1,Nx)]-(7/6)*f_n[idx_v(i,j+2,Nx)]+(2/6)*f_n[idx_v(i,j+3,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(i,Ny-1,Nx)]-(7/6)*f_p[idx_v(i,j-2,Nx)]+(11/6)*f_p[idx_v(i,j-1,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,j-2,Nx)]+(5/6)*f_p[idx_v(i,j-1,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i,j+1,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(i,j-2,Nx)]+(5/6)*f_n[idx_v(i,j-1,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i,j+1,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i,j+1,Nx)]+(2/6)*f_n[idx_v(i,j+2,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }
        
        if(j==Ny-1){
            IS0_p=(13/12)*(f_p[idx_v(i,j-2,Nx)]-2*f_p[idx_v(i,j-1,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(i,j-2,Nx)]-4*f_p[idx_v(i,j-1,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(i,j-1,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i,0,Nx)])^2+(1/4)*(f_p[idx_v(i,j-1,Nx)]-1*f_p[idx_v(i,0,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i,0,Nx)]+f_p[idx_v(i,1,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i,0,Nx)]+f_p[idx_v(i,1,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(i,0,Nx)]-2*f_n[idx_v(i,1,Nx)]+f_n[idx_v(i,2,Nx)])^2+(1/4)*(3*f_n[idx_v(i,0,Nx)]-4*f_n[idx_v(i,1,Nx)]+f_n[idx_v(i,2,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i,0,Nx)]+f_n[idx_v(i,1,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(i,1,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(i,j-1,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i,0,Nx)])^2+(1/4)*(f_n[idx_v(i,j-1,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i,0,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(i,j-2,Nx)]-(7/6)*f_p[idx_v(i,j-1,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,j-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i,0,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i,0,Nx)]-(1/6)*f_p[idx_v(i,1,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i,j-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i,0,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i,0,Nx)]-(1/6)*f_n[idx_v(i,1,Nx)])+w0_n*((11/6)*f_n[idx_v(i,0,Nx)]-(7/6)*f_n[idx_v(i,1,Nx)]+(2/6)*f_n[idx_v(i,2,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(i,j-3,Nx)]-(7/6)*f_p[idx_v(i,j-2,Nx)]+(11/6)*f_p[idx_v(i,j-1,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,j-2,Nx)]+(5/6)*f_p[idx_v(i,j-1,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i,0,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(i,j-2,Nx)]+(5/6)*f_n[idx_v(i,j-1,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i,0,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i,0,Nx)]+(2/6)*f_n[idx_v(i,1,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }
        
        if(i==Nx-2){
            IS0_p=(13/12)*(f_p[idx_v(i,j-2,Nx)]-2*f_p[idx_v(i,j-1,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(i,j-2,Nx)]-4*f_p[idx_v(i,j-1,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(i,j-1,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i,j+1,Nx)])^2+(1/4)*(f_p[idx_v(i,j-1,Nx)]-1*f_p[idx_v(i,j+1,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i,j+1,Nx)]+f_p[idx_v(i,0,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i,j+1,Nx)]+f_p[idx_v(i,0,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(i,j+1,Nx)]-2*f_n[idx_v(i,0,Nx)]+f_n[idx_v(i,1,Nx)])^2+(1/4)*(3*f_n[idx_v(i,j+1,Nx)]-4*f_n[idx_v(i,0,Nx)]+f_n[idx_v(i,1,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i,j+1,Nx)]+f_n[idx_v(i,0,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(i,0,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(i,j-1,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i,j+1,Nx)])^2+(1/4)*(f_n[idx_v(i,j-1,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i,j+1,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(i,j-2,Nx)]-(7/6)*f_p[idx_v(i,j-1,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,j-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i,j+1,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i,j+1,Nx)]-(1/6)*f_p[idx_v(i,0,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i,j-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i,j+1,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i,j+1,Nx)]-(1/6)*f_n[idx_v(i,0,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j+1,Nx)]-(7/6)*f_n[idx_v(i,0,Nx)]+(2/6)*f_n[idx_v(i,1,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(i,j-3,Nx)]-(7/6)*f_p[idx_v(i,j-2,Nx)]+(11/6)*f_p[idx_v(i,j-1,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,j-2,Nx)]+(5/6)*f_p[idx_v(i,j-1,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i,j+1,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(i,j-2,Nx)]+(5/6)*f_n[idx_v(i,j-1,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i,j+1,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i,j+1,Nx)]+(2/6)*f_n[idx_v(i,0,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }

        if(i==Nx-3){
            IS0_p=(13/12)*(f_p[idx_v(i,j-2,Nx)]-2*f_p[idx_v(i,j-1,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(i,j-2,Nx)]-4*f_p[idx_v(i,j-1,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(i,j-1,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i,j+1,Nx)])^2+(1/4)*(f_p[idx_v(i,j+1,Nx)]-1*f_p[idx_v(i,j+1,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i,j+1,Nx)]+f_p[idx_v(i,j+2,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i,j+1,Nx)]+f_p[idx_v(i,j+2,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(i,j+1,Nx)]-2*f_n[idx_v(i,j+2,Nx)]+f_n[idx_v(i,0,Nx)])^2+(1/4)*(3*f_n[idx_v(i,j+1,Nx)]-4*f_n[idx_v(i,j+2,Nx)]+f_n[idx_v(i,0,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i,j+1,Nx)]+f_n[idx_v(i,j+2,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(i,j+2,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(i,j-1,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i,j+1,Nx)])^2+(1/4)*(f_n[idx_v(i,j-1,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i,j+1,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(i,j-2,Nx)]-(7/6)*f_p[idx_v(i,j-1,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,j-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i,j+1,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i,j+1,Nx)]-(1/6)*f_p[idx_v(i,j+2,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i,j-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i,j+1,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i,j+1,Nx)]-(1/6)*f_n[idx_v(i,j+2,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j+1,Nx)]-(7/6)*f_n[idx_v(i,j+2,Nx)]+(2/6)*f_n[idx_v(i,0,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(i,j-3,Nx)]-(7/6)*f_p[idx_v(i,j-2,Nx)]+(11/6)*f_p[idx_v(i,j-1,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,j-2,Nx)]+(5/6)*f_p[idx_v(i,j-1,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i,j+1,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(i,j-2,Nx)]+(5/6)*f_n[idx_v(i,j-1,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i,j+1,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i,j+1,Nx)]+(2/6)*f_n[idx_v(i,j+2,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }

        if(j>2&&i<Ny-3){
            IS0_p=(13/12)*(f_p[idx_v(i,j-2,Nx)]-2*f_p[idx_v(i,j-1,Nx)]+f_p[idx_v(i,j,Nx)])^2+(1/4)*(f_p[idx_v(i,j-2,Nx)]-4*f_p[idx_v(i,j-1,Nx)]+3*f_p[idx_v(i,j,Nx)])^2;
            
            IS1_p=(13/12)*(f_p[idx_v(i,j-1,Nx)]-2*f_p[idx_v(i,j,Nx)]+f_p[idx_v(i,j+1,Nx)])^2+(1/4)*(f_p[idx_v(i,j-1,Nx)]-1*f_p[idx_v(i,j+1,Nx)])^2;
            
            IS2_p=(13/12)*(f_p[idx_v(i,j,Nx)]-2*f_p[idx_v(i,j+1,Nx)]+f_p[idx_v(i,j+2,Nx)])^2+(1/4)*(3*f_p[idx_v(i,j,Nx)]-4*f_p[idx_v(i,j+1,Nx)]+f_p[idx_v(i,j+2,Nx)])^2;
            
            IS0_n=(13/12)*(f_n[idx_v(i,j+1,Nx)]-2*f_n[idx_v(i,j+2,Nx)]+f_n[idx_v(i,j+3,Nx)])^2+(1/4)*(3*f_n[idx_v(i,j+1,Nx)]-4*f_n[idx_v(i,j+2,Nx)]+f_n[idx_v(i,j+3,Nx)])^2;
            
            IS1_n=(13/12)*(f_n[idx_v(i,j,Nx)]-2*f_n[idx_v(i,j+1,Nx)]+f_n[idx_v(i,j+2,Nx)])^2+(1/4)*(f_n[idx_v(i,j,Nx)]-1*f_n[idx_v(i,j+2,Nx)])^2;
            
            IS2_n=(13/12)*(f_n[idx_v(i,j-1,Nx)]-2*f_n[idx_v(i,j,Nx)]+f_n[idx_v(i,j+1,Nx)])^2+(1/4)*(f_n[idx_v(i,j-1,Nx)]-4*f_n[idx_v(i,j,Nx)]+3*f_n[idx_v(i,j+1,Nx)])^2;
            
            Epsilon=10e-6;
            alpha_0p=(1/10)*(1/(Epsilon+IS0_p))^2;
            alpha_1p=(6/10)*(1/(Epsilon+IS1_p))^2;
            alpha_2p=(3/10)*(1/(Epsilon+IS2_p))^2;
            alpha_0n=(1/10)*(1/(Epsilon+IS0_n))^2;
            alpha_1n=(6/10)*(1/(Epsilon+IS1_n))^2;
            alpha_2n=(3/10)*(1/(Epsilon+IS2_n))^2;
            
            w0_p=alpha_0p/(alpha_0p+alpha_1p+alpha_2p);
            w1_p=alpha_1p/(alpha_0p+alpha_1p+alpha_2p);
            w2_p=alpha_2p/(alpha_0p+alpha_1p+alpha_2p);
            w0_n=alpha_0n/(alpha_0n+alpha_1n+alpha_2n);
            w1_n=alpha_1n/(alpha_0n+alpha_1n+alpha_2n);
            w2_n=alpha_2n/(alpha_0n+alpha_1n+alpha_2n);
            
            u_tpph=w0_p*((2/6)*f_p[idx_v(i,j-2,Nx)]-(7/6)*f_p[idx_v(i,j-1,Nx)]+(11/6)*f_p[idx_v(i,j,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,j-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]+(2/6)*f_p[idx_v(i,j+1,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j,Nx)]+(5/6)*f_p[idx_v(i,j+1,Nx)]-(1/6)*f_p[idx_v(i,j+2,Nx)]);
            
            u_tnph=w2_n*((-1/6)*f_n[idx_v(i,j-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]+(2/6)*f_n[idx_v(i,j+1,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j,Nx)]+(5/6)*f_n[idx_v(i,j+1,Nx)]-(1/6)*f_n[idx_v(i,j+2,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j+1,Nx)]-(7/6)*f_n[idx_v(i,j+2,Nx)]+(2/6)*f_n[idx_v(i,j+3,Nx)]);
            
            u_tpnh=w0_p*((2/6)*f_p[idx_v(i,j-3,Nx)]-(7/6)*f_p[idx_v(i,j-2,Nx)]+(11/6)*f_p[idx_v(i,j-1,Nx)])+w1_p*((-1/6)*f_p[idx_v(i,j-2,Nx)]+(5/6)*f_p[idx_v(i,j-1,Nx)]+(2/6)*f_p[idx_v(i,j,Nx)])+w2_p*((2/6)*f_p[idx_v(i,j-1,Nx)]+(5/6)*f_p[idx_v(i,j,Nx)]-(1/6)*f_p[idx_v(i,j+1,Nx)]);
            
            u_tnnh=w2_n*((-1/6)*f_n[idx_v(i,j-2,Nx)]+(5/6)*f_n[idx_v(i,j-1,Nx)]+(2/6)*f_n[idx_v(i,j,Nx)])+w1_n*((2/6)*f_n[idx_v(i,j-1,Nx)]+(5/6)*f_n[idx_v(i,j,Nx)]-(1/6)*f_n[idx_v(i,j+1,Nx)])+w0_n*((11/6)*f_n[idx_v(i,j,Nx)]-(7/6)*f_n[idx_v(i,j+1,Nx)]+(2/6)*f_n[idx_v(i,j+2,Nx)]);
            
            result=-(1/dx)*((u_tpph-u_tpnh)+(u_tnph-u_tnnh));
        }

    }
    delete []f_p;
    delete []f_n;
    return result;
}


static int Gettao(realtype *rou, realtype *qx, realtype *qy, realtype *E, int Nx, int Ny, realtype *tao_xx, realtype *tao_xy, realtype *tao_yx, realtype *tao_yy){
    int k,i,j;
    realtype *vx = new realtype [Nx*Ny];
    realtype *vy = new realtype [Nx*Ny];
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            vx[idx_v(i,j,Nx)]=qx[idx_v(i, j, Nx)]/rou[idx_v(i, j, Nx)];
            vy[idx_v(i,j,Nx)]=qy[idx_v(i, j, Nx)]/rou[idx_v(i, j, Nx)];
         }
    }
//    for (k=0;k<4;k++){
    
        for(j=0;j<Ny;j++){
            for(i=0;i<Nx;i++){
                tao_xx[idx_v(i,j,Nx)]=qx[idx_v(i,j,Nx)]*vx[idx_v(i,j,Nx)]+rou[idx_v(i, j, Nx)];
                tao_xy[idx_v(i,j,Nx)]=qx[idx_v(i,j,Nx)]*vy[idx_v(i,j,Nx)];
                tao_yx[idx_v(i,j,Nx)]=qy[idx_v(i,j,Nx)]*vx[idx_v(i,j,Nx)];
                tao_yy[idx_v(i,j,Nx)]=qy[idx_v(i,j,Nx)]*vy[idx_v(i,j,Nx)]+rou[idx_v(i, j, Nx)];
            }
        }
    delete *vx;
    delete *vy;
    return 0;
//    }
}


static int GetCj(realtype *rou, realtype *qx, realtype *qy, realtype *E, int Nx, int Ny, realtype *Cj_x, realtype *Cj_y){
    int k,i,j;
    realtype *vx = new realtype [Nx*Ny];
    realtype *vy = new realtype [Nx*Ny];
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            vx[idx_v(i,j,Nx)]=qx[idx_v(i, j, Nx)]/rou[idx_v(i, j, Nx)];
            vy[idx_v(i,j,Nx)]=qy[idx_v(i, j, Nx)]/rou[idx_v(i, j, Nx)];
        }
    }
    //    for (k=0;k<4;k++){
    
    for(j=0;j<Ny;j++){
        for(i=0;i<Nx;i++){
            Cj_x[idx_v(i,j,Nx)]=qx[idx_v(i,j,Nx)]*qx[idx_v(i,j,Nx)]*vx[idx_v(i,j,Nx)]+rou[idx_v(i, j, Nx)]*qx[idx_v(i,j,Nx)]+qx[idx_v(i,j,Nx)]*qy[idx_v(i,j,Nx)]*vy[idx_v(i,j,Nx)];
            Cj_y[idx_v(i,j,Nx)]=qy[idx_v(i,j,Nx)]*qy[idx_v(i,j,Nx)]*vy[idx_v(i,j,Nx)]+rou[idx_v(i, j, Nx)]*qy[idx_v(i,j,Nx)]+qx[idx_v(i,j,Nx)]*qy[idx_v(i,j,Nx)]*vx[idx_v(i,j,Nx)];
        }
    }
    delete *vx;
    delete *vy;
    return 0;
    //    }
}
