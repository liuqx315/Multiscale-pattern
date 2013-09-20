/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006-07-05 15:32:37 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Steven Smith @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2013, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This file (companion of nvector_openmp.h) contains the
 * implementation needed for the Fortran initialization of openmp
 * vector operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fnvector_openmp.h"

/* Define global vector variables */

N_Vector F2C_CVODE_vec;
N_Vector F2C_CVODE_vecQ;
N_Vector *F2C_CVODE_vecS;
N_Vector F2C_CVODE_vecB;
N_Vector F2C_CVODE_vecQB;

N_Vector F2C_IDA_vec;
N_Vector F2C_IDA_vecQ;
N_Vector *F2C_IDA_vecS;
N_Vector F2C_IDA_vecB;
N_Vector F2C_IDA_vecQB;

N_Vector F2C_KINSOL_vec;

N_Vector F2C_ARKODE_vec;

/* Fortran callable interfaces */

void FNV_INITS(int *code, long int *N, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = NULL;
    F2C_CVODE_vec = N_VNewEmpty_OpenMP(*N);
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vec = NULL;
    F2C_IDA_vec = N_VNewEmpty_OpenMP(*N);
    if (F2C_IDA_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = NULL;
    F2C_KINSOL_vec = N_VNewEmpty_OpenMP(*N);
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    F2C_ARKODE_vec = NULL;
    F2C_ARKODE_vec = N_VNewEmpty_OpenMP(*N);
    if (F2C_ARKODE_vec == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITS_Q(int *code, long int *Nq, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecQ = NULL;
    F2C_CVODE_vecQ = N_VNewEmpty_OpenMP(*Nq);
    if (F2C_CVODE_vecQ == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecQ = NULL;
    F2C_IDA_vecQ = N_VNewEmpty_OpenMP(*Nq);
    if (F2C_IDA_vecQ == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITS_B(int *code, long int *NB, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecB = NULL;
    F2C_CVODE_vecB = N_VNewEmpty_OpenMP(*NB);
    if (F2C_CVODE_vecB == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecB = NULL;
    F2C_IDA_vecB = N_VNewEmpty_OpenMP(*NB);
    if (F2C_IDA_vecB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITS_QB(int *code, long int *NqB, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecQB = NULL;
    F2C_CVODE_vecQB = N_VNewEmpty_OpenMP(*NqB);
    if (F2C_CVODE_vecQB == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecQB = NULL;
    F2C_IDA_vecQB = N_VNewEmpty_OpenMP(*NqB);
    if (F2C_IDA_vecQB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITS_S(int *code, int *Ns, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecS = NULL;
    F2C_CVODE_vecS = (N_Vector *) N_VCloneVectorArrayEmpty_OpenMP(*Ns, F2C_CVODE_vec);
    if (F2C_CVODE_vecS == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecS = NULL;
    F2C_IDA_vecS = (N_Vector *) N_VCloneVectorArrayEmpty_OpenMP(*Ns, F2C_IDA_vec);
    if (F2C_IDA_vecS == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


