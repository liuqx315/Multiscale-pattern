/*----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
 Copyright (c) 2013, Southern Methodist University.
 All rights reserved.
 For details, see the LICENSE file.
 ----------------------------------------------------------------
 sparse matrix test routine:
 
 The following program tests a variety of aspects of the sparse 
 matrix module in SUNDIALS.
----------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_sparse.h>
#include <sundials/sundials_direct.h>


/* Main Program */
int main() {

  /* local variables */
  int i, j, nnz, iret;
  SlsMat A, B, C, D, E, F;
  DlsMat G, H;
  realtype *x, *y, *z;

  /* set matrix dimensions, density */
  const int m  = 15;     /* rows         */
  const int n  = 20;     /* columns      */
  const int nz = 5*n;    /* max nonzeros */

  /* initialize all pointers to NULL */
  A = B = C = D = E = F = NULL;
  G = H = NULL;
  x = y = z = NULL;

  /* print current sparse matrices */
  printf("Just-created sparse matrices:\n");
  printf("\nA (%ix%i all zeros, %i nonzeros):\n",m,m,nz);
  A = NewSparseMat(m,m,nz);   /* square matrix */
  PrintSparseMat(A);

  printf("\nB (%ix%i all zeros, %i nonzeros):\n",m,m,nz);
  B = NewSparseMat(m,m,nz);   /* square matrix */
  PrintSparseMat(B);

  printf("\nC (%ix%i all zeros, %i nonzeros):\n",m,n,nz);
  C = NewSparseMat(m,n,nz);   /* rectangular matrix */
  PrintSparseMat(C);

  printf("\nD (%ix%i all zeros, %i nonzeros):\n",m,n,nz);
  D = NewSparseMat(m,n,nz);   /* rectangular matrix */
  PrintSparseMat(D);

  /* set values to zero and print */
  printf("Manually zeroed-out sparse matrices:\n");
  printf("\nA = 0:\n");
  SlsSetToZero(A);
  PrintSparseMat(A);

  printf("\nB = 0:\n");
  SlsSetToZero(B);
  PrintSparseMat(B);

  printf("\nC = 0:\n");
  SlsSetToZero(C);
  PrintSparseMat(C);

  SlsSetToZero(D);
  printf("\nD = 0:\n");
  PrintSparseMat(D);

  /* fill some matrix values */
  printf("\nA = 2I:\n");
  AddIdentitySparseMat(A);     /* A = I */
  AddIdentitySparseMat(A);     /* A = 2I */
  PrintSparseMat(A);

  printf("\nB = 5I:\n");
  CopySparseMat(A,B);          /* B = I */
  ScaleSparseMat(2.5,B);       /* B = 5I */
  PrintSparseMat(B);

  printf("\nC has nontrivial shape/values:\n");
  nnz = 0;
  for (j=0; j<n; j++) {        /* C has some 'random' stuff */
    C->colptrs[j] = nnz;       /* point to start of jth column */
    i = j-5;                   /* set a value in the 5th superdiag */
    if (i>=0 && i<m) {
      C->rowvals[nnz] = i;
      C->data[nnz] = -1.0*(j-5);
      nnz++;
    }
    i = j-2;                   /* set a value in the 2nd superdiag */
    if (i>=0 && i<m) {
      C->rowvals[nnz] = i;
      C->data[nnz] = 2.0*(j+1);
      nnz++;
    }
    i = j;                     /* set a diagonal value (not in row 3) */
    if (i>=0 && i<m && j!=3) {
      C->rowvals[nnz] = i;
      C->data[nnz] = 1.0*(j+4);
      nnz++;
    }
    i = j+1;                   /* set a value in 1st subdiag */
    if (i>=0 && i<m) {
      C->rowvals[nnz] = i;
      C->data[nnz] = 1.0*(j+4);
      nnz++;
    }      
    i = j+3;                   /* set a value in 3rd subdiag */
    if (i>=0 && i<m) {
      C->rowvals[nnz] = i;
      C->data[nnz] = -3.0*(j+2);
      nnz++;
    }      
  }
  C->colptrs[j] = nnz;         /* indicate end of data */
  PrintSparseMat(C);

  printf("\nD = I:\n");
  AddIdentitySparseMat(D);     /* D = 2*I + 2*C */
  PrintSparseMat(D);
  printf("\nD = 2I:\n");
  ScaleSparseMat(2.0,D);
  PrintSparseMat(D);
  printf("\nD = 2I+C:\n");
  SlsAddMat(D,C);
  PrintSparseMat(D);
  printf("\nD = 2I+2C:\n");
  SlsAddMat(D,C);
  PrintSparseMat(D);

  /* create some vectors */
  x = (double *) malloc(m * sizeof(realtype));
  y = (double *) malloc(n * sizeof(realtype));
  z = (double *) malloc(m * sizeof(realtype));

  /* try some matrix-vector products */
  printf("\nx is a vector of ones:\n");
  for (i=0; i<m; i++)  x[i] = 1.0;
  for (i=0; i<m; i++)
    printf("  %16g\n", x[i]);

  printf("\nB is:\n");
  PrintSparseMat(B);
  printf("\nz = B*x: (matvec returned %i)\n",iret);
  iret = SlsMatvec(B, x, z);
  for (i=0; i<m; i++)
    printf("  %16g\n", z[i]);

  printf("\ny is a vector of ones:\n");
  for (i=0; i<n; i++)  y[i] = 1.0;
  for (i=0; i<n; i++)
    printf("  %16g\n", y[i]);

  printf("\nC is:\n");
  PrintSparseMat(C);
  printf("\nz = C*y: (matvec returned %i)\n",iret);
  iret = SlsMatvec(C, y, z);
  for (i=0; i<m; i++)
    printf("  %16g\n", z[i]);

  printf("\nnontrivial x vector:\n");
  for (i=0; i<m; i++)  x[i] = 1.0*(m-i+1);
  for (i=0; i<m; i++)
    printf("  %16g\n", x[i]);

  printf("\nA is:\n");
  PrintSparseMat(A);
  printf("\nz = A*x: (matvec returned %i)\n",iret);
  iret = SlsMatvec(A, x, z);
  for (i=0; i<m; i++)
    printf("  %16g\n", z[i]);

  printf("\nnontrivial y vector:\n");
  for (i=0; i<n; i++)  y[i] = 1.0*(n-i+1);
  for (i=0; i<n; i++)
    printf("  %16g\n", y[i]);

  printf("\nD is:\n");
  PrintSparseMat(D);
  printf("\nz = D*y: (matvec returned %i)\n",iret);
  iret = SlsMatvec(D, y, z);
  for (i=0; i<m; i++)
    printf("  %16g\n", z[i]);


  /* create/fill dense matrix, G */
  G = NewDenseMat(m, n);
  AddIdentity(G);   /* G = I */
  DENSE_ELEM(G, 0, 1)   =  2.0;
  DENSE_ELEM(G, 1, 5)   = -6.0;
  DENSE_ELEM(G, 2, 0)   =  3.0;
  DENSE_ELEM(G, 3, 7)   =  4.0;
  DENSE_ELEM(G, 4, 2)   = -7.0;
  DENSE_ELEM(G, 5, 0)   = -2.0;
  DENSE_ELEM(G, 6, 8)   =  9.0;
  DENSE_ELEM(G, 7, 19)  =  1.0;
  DENSE_ELEM(G, 8, 12)  = 12.0;
  DENSE_ELEM(G, 9, 3)   = -2.0;
  DENSE_ELEM(G, 10, 7)  =  5.0;
  DENSE_ELEM(G, 11, 4)  =  7.0;
  DENSE_ELEM(G, 12, 15) =  6.0;
  DENSE_ELEM(G, 13, 17) =  9.0;
  DENSE_ELEM(G, 14, 19) = -8.0;
  printf("\nDense matrix G is:\n");
  PrintMat(G);
  E = SlsConvertDls(G);
  printf("\nE = sparse(G):\n");
  PrintSparseMat(E);
  
  
  /* create/fill band matrix, H */
  H = NewBandMat(n, 2, 3, 5);
  AddIdentity(H);   /* H = I */
  BAND_ELEM(H, 0, 1)   =  2.0;
  BAND_ELEM(H, 1, 4)   = -6.0;
  BAND_ELEM(H, 2, 0)   =  3.0;
  BAND_ELEM(H, 3, 4)   =  4.0;
  BAND_ELEM(H, 4, 2)   = -7.0;
  BAND_ELEM(H, 5, 3)   = -2.0;
  BAND_ELEM(H, 6, 8)   =  9.0;
  BAND_ELEM(H, 7, 10)  =  1.0;
  BAND_ELEM(H, 8, 10)  = 12.0;
  BAND_ELEM(H, 9, 8)   = -2.0;
  BAND_ELEM(H, 10, 8)  =  5.0;
  BAND_ELEM(H, 11, 12) =  7.0;
  BAND_ELEM(H, 12, 12) =  6.0;
  BAND_ELEM(H, 13, 15) =  9.0;
  BAND_ELEM(H, 14, 17) = -8.0;
  BAND_ELEM(H, 15, 13) =  4.0;
  BAND_ELEM(H, 16, 17) = -7.0;
  BAND_ELEM(H, 17, 15) =  2.0;
  BAND_ELEM(H, 18, 19) =  1.0;
  BAND_ELEM(H, 19, 18) = -3.0;
  printf("\nBand matrix H is:\n");
  PrintMat(H);
  F = SlsConvertDls(H);
  printf("\nF = sparse(H):\n");
  PrintSparseMat(F);
  

  
  /* clean up */
  free(x);
  free(y);
  free(z);
  DestroySparseMat(A);
  DestroySparseMat(B);
  DestroySparseMat(C);
  DestroySparseMat(D);
  DestroySparseMat(E);
  DestroySparseMat(F);
  DestroyMat(G);
  DestroyMat(H);

  return 0;
}



/*---- end of file ----*/
