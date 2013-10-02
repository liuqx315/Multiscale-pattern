/*---------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2013, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 ----------------------------------------------------------------
 This testing routine is designed to exercise the parallel 
 NVector routines.
---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_parallel.h>
#include <sundials/sundials_math.h>

/* floating point "equality" comparison, failure update macro */
#define FNEQ(a,b) ( fabs(a-b)/fabs(b) > 1.0e-15 )

/* vector length on each process */
#define VEC_LEN 100


/* Outer testing routine -- note that all tests should have 
   data arrays with 100 total entries */
int main(int argc, char *argv[]) {

  /* initialize MPI */
  MPI_Init(&argc, &argv);
  int nprocs, my_id;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &my_id);

  /* initialize success/failure flag */
  int failure = 0;

  /* determine global vector size */
  long int loc_len = VEC_LEN;
  long int glob_len = 0;
  MPI_Allreduce(&loc_len, &glob_len, 1, MPI_LONG, MPI_SUM, comm);
  
  /* Make NVector Uvec out of user data */
  realtype *Udata = malloc(loc_len * sizeof(realtype));
  N_Vector Uvec = N_VMake_Parallel(comm, loc_len, glob_len, Udata);

  /* Make NVector Xvec out of user data (approach 2) */
  realtype *Xdata = malloc(loc_len * sizeof(realtype));
  N_Vector Xvec = N_VNewEmpty_Parallel(comm, loc_len, glob_len);
  N_VSetArrayPointer(Xdata,Xvec);

  /* Make NVector Yvec from scratch */
  N_Vector Yvec = N_VNew_Parallel(comm, loc_len, glob_len);

  /* Make NVector Zvec by cloning Uvec */
  N_Vector Zvec = N_VClone(Uvec);

  /* Test 1 -- N_VConst and N_VMaxNorm */
  realtype a;
  N_VConst( 1.0, Uvec);  a = N_VMaxNorm(Uvec);
  if (FNEQ(a, 1.0)) {
    printf("  Failed test -- N_VConst and N_VMaxNorm,  U=1, X=2, Y=3, Z=-4\n");
    printf("     ||U||_inf = %g  (should be 1.0)\n",a);
    failure++;
  }

  N_VConst( 2.0, Xvec);  a = N_VMaxNorm(Xvec);
  if (FNEQ(a, 2.0)) {
    printf("  Failed test -- N_VConst and N_VMaxNorm,  U=1, X=2, Y=3, Z=-4\n");
    printf("     ||X||_inf = %g  (should be 2.0)\n",a);
    failure++;
  }

  N_VConst( 3.0, Yvec);  a = N_VMaxNorm(Yvec);
  if (FNEQ(a, 3.0)) {
    printf("  Failed test -- N_VConst and N_VMaxNorm,  U=1, X=2, Y=3, Z=-4\n");
    printf("     ||Y||_inf = %g  (should be 3.0)\n",a);
    failure++;
  }

  N_VConst(-4.0, Zvec);  a = N_VMaxNorm(Zvec);
  if (FNEQ(a, 4.0)) {
    printf("  Failed test -- N_VConst and N_VMaxNorm,  U=1, X=2, Y=3, Z=-4\n");
    printf("     ||Z||_inf = %g  (should be 4.0)\n",a);
    failure++;
  }

  /* Test 2 -- N_VSpace */
  long int nfloat, nint;
  N_VSpace(Uvec, &nfloat, &nint);
  if ((nfloat != glob_len) || (nint != 2*nprocs)) {
    printf("  Failed test -- N_VSpace\n");
    printf("     the N_Vector requires %li reals and %li ints (should be %li, %i)\n",
	   nfloat, nint, glob_len, 2*nprocs);
    failure++;
  }

  /* Test 3 -- N_VScale: first copy U into X, then put 3*U into Y and output */
  N_VScale( 1.0, Uvec, Xvec );  a = N_VMaxNorm(Xvec);
  if (FNEQ(a, 1.0)) {
    printf("  Failed test -- N_VScale,  X=U,  Y=3*U\n");
    printf("     ||X||_inf = %g  (should be 1.0)\n",a);
    failure++;
  }

  N_VScale( 3.0, Xvec, Yvec );  a = N_VMaxNorm(Yvec);
  if (FNEQ(a, 3.0)) {
    printf("  Failed test -- N_VScale,  X=U,  Y=3*U\n");
    printf("     ||Y||_inf = %g  (should be 3.0)\n",a);
    failure++;
  }

  /* Test 4 -- N_VLinearSum: X=Xvec==Uvec, Y=Yvec==3*Uvec, a=0.5, b=-0.5 */
  N_VLinearSum( 0.5, Xvec, -0.5, Yvec, Zvec );
  a = N_VMaxNorm(Zvec);
  if (FNEQ(a, 1.0)) {
    printf("  Failed test -- N_VLinearSum,  Z=0.5*X-0.5*Y\n");
    printf("     ||Z||_inf = %g  (should be 1.0)\n",a);
    failure++;
  }

  /* Test 5 -- N_VLinearSum: X=Xvec==Uvec, Y=Yvec==3*Uvec, a=-1.0, b=1.0 */
  N_VLinearSum( -1.0, Xvec, 1.0, Yvec, Zvec );
  a = N_VMaxNorm(Zvec);
  if (FNEQ(a, 2.0)) {
    printf("  Failed test -- N_VLinearSum,  Z=-1.0*X+1.0*Y\n");
    printf("     ||Z||_inf = %g  (should be 2.0)\n",a);
    failure++;
  }

  /* Test 6 -- N_VLinearSum: X=Xvec==Uvec, Y=Yvec==3*Uvec, a=1.0, b=0.5 */
  N_VLinearSum( 1.0, Xvec, 0.5, Yvec, Zvec );
  a = N_VMaxNorm(Zvec);
  if (FNEQ(a, 2.5)) {
    printf("  Failed test -- N_VLinearSum,  Z=X+0.5*Y\n");
    printf("     ||Z||_inf = %g  (should be 2.5)\n",a);
    failure++;
  }

  /* Test 7 -- N_VLinearSum: X=Xvec==Uvec, Y=Yvec==3*Uvec, a=5.0, b=1.0, Z=Yvec */
  N_VLinearSum( 0.5, Xvec, 1.0, Yvec, Yvec );
  a = N_VMaxNorm(Yvec);
  if (FNEQ(a, 3.5)) {
    printf("  Failed test -- N_VLinearSum,  Y=0.5*X+1.0*Y\n");
    printf("     ||Y||_inf = %g  (should be 3.5)\n",a);
    failure++;
  }

  /* Test 8 -- N_VProd: X=Xvec==Uvec, Y=Yvec==3*Uvec */
  N_VScale( 3.0, Uvec, Yvec );
  N_VProd( Xvec, Yvec, Zvec );
  a = N_VMaxNorm(Zvec);
  if (FNEQ(a, 3.0)) {
    printf("  Failed test -- N_VScale and N_VProd,  Z=X.*(3*U)\n");
    printf("     ||Z||_inf = %g  (should be 3.0)\n",a);
    failure++;
  }

  /* Test 9 -- N_VDiv: X=Xvec=Uvec, Y=Yvec==3*Uvec */
  N_VDiv( Xvec, Yvec, Zvec );
  a = N_VMaxNorm(Zvec);
  if (FNEQ(a, 1.0/3.0)) {
    printf("  Failed test -- N_VDiv,  Z=X./Y\n");
    printf("     ||Z||_inf = %g  (should be 1/3)\n",a);
    failure++;
  }

  /* Test 10 -- N_VAddConst: X=Xvec==Uvec, b=-20000, Z=Yvec */
  N_VAddConst( Xvec, -20000.0, Zvec );
  a = N_VMaxNorm(Zvec);
  if (FNEQ(a, 19999.0)) {
    printf("  Failed test -- N_VAddConst,  Z=X-20000\n");
    printf("     ||Z||_inf = %g  (should be 19999.0)\n",a);
    failure++;
  }

  /* Test 11 -- N_VAbs: X=Yvec */
  N_VAbs( Zvec, Yvec );
  a = N_VMaxNorm(Yvec);
  if (FNEQ(a, 19999.0)) {
    printf("  Failed test -- N_VAbs,  Y=|Z|\n");
    printf("     ||Y||_inf = %g  (should be 19999.0)\n",a);
    failure++;
  }

  /* Test 12 -- N_VInv: X=Yvec */
  N_VInv( Yvec, Zvec );
  a = N_VMaxNorm(Zvec);
  if (FNEQ(a, 1.0/19999.0)) {
    printf("  Failed test -- N_VInv,  Z=1./Y\n");
    printf("     ||Z||_inf = %g  (should be 5.00025e-5)\n",a);
    failure++;
  }

  /* Test 13 -- N_VDotProd: X=Xvec==Uvec, Y=Yvec==3*Uvec */
  N_VScale( 3.0, Uvec, Yvec );
  a = N_VDotProd( Xvec, Yvec );
  if (FNEQ(a, 3.0*glob_len)) {
    printf("  Failed test -- N_VScale and N_VDotProd\n");
    printf("     (X,3*U) = %g  (should be %g)\n",a,3.0*glob_len);
    failure++;
  }

  /* Test 14 -- N_VWrmsNorm: X=Xvec==Uvec, W=Yvec==1.5 */
  N_VConst( 1.5, Yvec );
  a = N_VWrmsNorm( Xvec, Yvec );
  if (FNEQ(a, 1.5)) {
    printf("  Failed test -- N_VConst and N_VWrmsNorm, Y=1.5\n");
    printf("     ||X||_Y = %g  (should be 1.5)\n",a);
    failure++;
  }

  /* Test 15 -- N_VWrmsNormMask: X=Xvec==Uvec, W=Yvec==1.5, ID=Zvec==1 */
  N_VConst( 1.0, Zvec );
  a = N_VWrmsNormMask( Xvec, Yvec, Zvec );
  if (FNEQ(a, 1.5)) {
    printf("  Failed test -- N_VConst and N_VWrmsNormMask\n");
    printf("     ||X||_(Y,1) = %g  (should be 1.5)\n",a);
    failure++;
  }

  /* Test 16 -- N_VMin: X=Xvec==Uvec */
  a = N_VMin( Xvec );
  if (FNEQ(a, 1.0)) {
    printf("  Failed test -- N_VMin\n");
    printf("     min(X) = %g  (should be 1.0)\n",a);
    failure++;
  }

  /* Test 17 -- N_VWL2Norm: X=Xvec==Uvec, W=Yvec==2*Uvec */
  N_VScale( 2.0, Uvec, Yvec );
  a = N_VWL2Norm( Xvec, Yvec );
  if (FNEQ(a, sqrt(4.0*glob_len))) {
    printf("  Failed test -- N_VWL2Norm,  Y=2*U\n");
    printf("     ||X||_Y = %g  (should be %g)\n",a,sqrt(4.0*glob_len));
    failure++;
  }

  /* Test 18 -- N_VL1Norm: X=Xvec==Uvec */
  a = N_VL1Norm( Xvec );
  if (FNEQ(a, 1.0*glob_len)) {
    printf("  Failed test -- N_VL1Norm\n");
    printf("     ||X||_1 = %g  (should be %g)\n",a,1.0*glob_len);
    failure++;
  }

  /* Test 19 -- N_VCompare: X=Xvec==Uvec, c=a==20000.0 */
  N_VCompare( 2000.0, Xvec, Zvec );
  a = N_VMaxNorm(Zvec);
  if (FNEQ(a, 0.0)) {
    printf("  Failed test -- N_VCompare,  Z = 2000 vs X\n");
    printf("     ||Z||_inf = %g  (should be 0.0)\n",a);
    failure++;
  }

  /* Test 20 -- N_VInvTest: X=Xvec==Uvec */
  a = N_VMaxNorm(Zvec);
  booleantype booltest = N_VInvTest( Xvec, Zvec );
  a = N_VMaxNorm(Zvec);
  if ((booltest != 1) || FNEQ(a, 1.0)) {
    printf("  Failed test -- N_VInvTest,  Z=1./X\n");
    printf("     test = %i, ||Z||_inf = %g  (should be 1, 1.0)\n", booltest, a);
    failure++;
  }

  /* Test 21 -- N_VMinQuotient: num=Xvec==Uvec, denom=Yvec==4d-5*Uvec */
  N_VScale( 0.00004, Uvec, Yvec );
  a = N_VMinQuotient( Xvec, Yvec );
  if (FNEQ(a, 25000.0)) {
    printf("  Failed test -- N_VMinQuotient,  Y=0.00004*X\n");
    printf("     min(X./Y) = %g  (should be 25000.0)\n",a);
    failure++;
  }

  
  /* Accumulate final SUCCESS/FAIL result */
  int max_fails = 0;
  MPI_Allreduce(&failure, &max_fails, 1, MPI_INT, MPI_MAX, comm);

  /* Root node prints final SUCCESS/FAIL result */
  if (my_id == 0) {
    if (max_fails) 
      printf("  FAIL: Parallel NVector module failed %i tests\n",max_fails);
    else 
      printf("  SUCCESS: Parallel NVector module passed all tests\n");
  }


  /* Free all test vectors */
  N_VDestroy(Uvec);
  N_VDestroy(Xvec);
  N_VDestroy(Yvec);
  N_VDestroy(Zvec);  

  /* Free locally allocated vector data */
  free(Udata);
  free(Xdata);
  
  /* finalize MPI and return */
  MPI_Finalize();

  return 0;
} 

/* end of testing routine ************************************************/
