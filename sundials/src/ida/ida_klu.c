/*
 * -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date: $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2013, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the IDAKLU linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "ida_impl.h"
#include "ida_sparse_impl.h"
#include "sundials/sundials_klu_impl.h"
#include "sundials/sundials_math.h"

/* Constants */

#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* IDAKLU linit, lsetup, lsolve, and lfree routines */
 
static int IDAKLUInit(IDAMem IDA_mem);

static int IDAKLUSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
		       N_Vector rrp, N_Vector tmp1,
		       N_Vector tmp2, N_Vector tmp3);

static int IDAKLUSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
			     N_Vector ycur, N_Vector ypcur, N_Vector rrcur);

static int IDAKLUFree(IDAMem IDA_mem);

/* 
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
 */

/*
 * -----------------------------------------------------------------
 * IDAKLU
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the IDA / KLU linear solver module.  
 * IDAKLU first calls the existing lfree routine if this is not NULL.
 * Then it sets the ida_linit, ida_lsetup, ida_lsolve, ida_lperf, and
 * ida_lfree fields in (*IDA_mem) to be IDAKLUInit, IDAKLUSetup,
 * IDAKLUSolve, NULL, and IDAKLUFree, respectively.
 * It allocates memory for a structure of type IDAkluMemRec and sets
 * the ida_lmem field in (*IDA_mem) to the address of this structure.
 * It sets setupNonNull in (*IDA_mem) to TRUE.
 * Finally, it allocates memory for KLU.
 * The return value is IDASLS_SUCCESS = 0, IDASLS_LMEM_FAIL = -1,
 * or IDASLS_ILL_INPUT = -2.
 *
 * NOTE: The KLU linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, IDAKLU will first 
 *       test for a compatible N_Vector internal
 *       representation by checking that the functions N_VGetArrayPointer
 *       and N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int IDAKLU(void *ida_mem, int n, int nnz)
{
  IDAMem IDA_mem;
  IDASlsMem idasls_mem;
  KLUData klu_data;
  int flag;

  /* Return immediately if ida_mem is NULL. */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASSLS", "IDAKLU", 
		    MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test if the NVECTOR package is compatible with the Direct solver */
  if(IDA_mem->ida_tempv1->ops->nvgetarraypointer == NULL ||
     IDA_mem->ida_tempv1->ops->nvsetarraypointer == NULL) {
    IDAProcessError(IDA_mem, IDASLS_ILL_INPUT, "IDASSLS", "IDAKLU", 
		    MSGSP_BAD_NVECTOR);
    return(IDASLS_ILL_INPUT);
  }

  if (IDA_mem->ida_lfree != NULL) flag = IDA_mem->ida_lfree(IDA_mem);

  /* Set five main function fields in IDA_mem. */
  IDA_mem->ida_linit  = IDAKLUInit;
  IDA_mem->ida_lsetup = IDAKLUSetup;
  IDA_mem->ida_lsolve = IDAKLUSolve;
  IDA_mem->ida_lperf  = NULL;
  IDA_mem->ida_lfree  = IDAKLUFree;

  /* Get memory for IDASlsMemRec. */
  idasls_mem = (IDASlsMem) malloc(sizeof(struct IDASlsMemRec));
  if (idasls_mem == NULL) {
    IDAProcessError(IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", "IDAKLU", 
		    MSGSP_MEM_FAIL);
    return(IDASLS_MEM_FAIL);
  }

  /* Get memory for KLUData. */
  klu_data = (KLUData)malloc(sizeof(struct KLUDataRec));
  if (klu_data == NULL) {
    IDAProcessError(IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", "IDAKLU", 
		    MSGSP_MEM_FAIL);
    return(IDASLS_MEM_FAIL);
  }

  IDA_mem->ida_setupNonNull = TRUE;

  /* Set default Jacobian routine and Jacobian data */
  idasls_mem->s_jaceval = NULL;
  idasls_mem->s_jacdata = IDA_mem->ida_user_data;

  /* Allocate memory for the sparse Jacobian */
  idasls_mem->s_JacMat = NewSparseMat(n, n, nnz);
  if (idasls_mem->s_JacMat == NULL) {
    IDAProcessError(IDA_mem, IDASLS_MEM_FAIL, "IDASSLS", "IDAKLU", 
		    MSGSP_MEM_FAIL);
    return(IDASLS_MEM_FAIL);
  }

  /* Allocate structures for KLU */

  /* DO I ALLOCATE COMMON????*/

  klu_data->s_Symbolic = (klu_symbolic *)malloc(sizeof(klu_symbolic));
  klu_data->s_Numeric = (klu_numeric *)malloc(sizeof(klu_numeric));

  /* Set default parameters for KLU */
  flag = klu_defaults(&klu_data->s_Common);
  if (flag == 0) {
    IDAProcessError(IDA_mem, IDASLS_PACKAGE_FAIL, "IDASSLS", "IDAKLU", 
		    MSGSP_PACKAGE_FAIL);
    return(IDASLS_PACKAGE_FAIL);
  }

  /* Set ordering to COLAMD as the idas default use.
     Users can set a different value with IDAKLUSetOrdering,
     and the user-set value is loaded before any call to klu_analyze in
     IDAKLUSetup.  */
  klu_data->s_ordering = 1;
  klu_data->s_Common.ordering = klu_data->s_ordering;

  /* Attach linear solver memory to the integrator memory */
  idasls_mem->s_solver_data = (void *) klu_data;
  IDA_mem->ida_lmem = idasls_mem;

  idasls_mem->s_last_flag = IDASLS_SUCCESS;

  return(IDASLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDAKLU interface functions
 * -----------------------------------------------------------------
 */

/*
  This routine does remaining initializations specific to the IDAKLU
  linear solver module.  
  It returns 0 if successful.
*/

static int IDAKLUInit(IDAMem IDA_mem)
{
  IDASlsMem idasls_mem;

  idasls_mem = (IDASlsMem)IDA_mem->ida_lmem;

  idasls_mem->s_nje = 0;
  idasls_mem->s_first_factorize = 1;

  idasls_mem->s_last_flag = 0;
  return(0);
}

/*
  This routine does the setup operations for the IDAKLU linear 
  solver module.  It calls the Jacobian evaluation routine,
  updates counters, and calls the LU factorization routine.
  The return value is either
     IDASLS_SUCCESS = 0  if successful,
     +1  if the jac routine failed recoverably or the
         LU factorization failed, or
     -1  if the jac routine failed unrecoverably.
*/

static int IDAKLUSetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp,
		       N_Vector rrp, N_Vector tmp1, N_Vector tmp2,
		       N_Vector tmp3)
{
  int retval, last_flag;
  realtype tn, cj;
  IDASlsMem idasls_mem;
  IDASlsSparseJacFn jaceval;
  KLUData klu_data;
  SlsMat JacMat;
  void *jacdata;
  
  idasls_mem = (IDASlsMem) (IDA_mem->ida_lmem);
  tn = IDA_mem->ida_tn; 
  cj = IDA_mem->ida_cj;

  klu_data = (KLUData) idasls_mem->s_solver_data;

  last_flag = idasls_mem->s_last_flag;
  jaceval = idasls_mem->s_jaceval;
  jacdata = idasls_mem->s_jacdata;
  JacMat = idasls_mem->s_JacMat;

  /* Check that Jacobian eval routine is set */
  if (jaceval == NULL) {
    IDAProcessError(IDA_mem, IDASLS_JAC_NOSET, "IDASSLS", "IDAKLUSetup", 
		    MSGSP_JAC_NOSET);
    free(idasls_mem); idasls_mem = NULL;
    return(IDASLS_JAC_NOSET);
  }

  /* Increment nje counter and call Jacobian eval routine. */
  idasls_mem->s_nje++;
  retval = jaceval(tn, cj, yyp, ypp, rrp, JacMat, jacdata, 
		   tmp1, tmp2, tmp3);

  if (retval < 0) {
    IDAProcessError(IDA_mem, IDASLS_JACFUNC_UNRECVR, "IDASSLS", 
		    "IDAKLUSetup", MSGSP_JACFUNC_FAILED);
    last_flag = IDASLS_JACFUNC_UNRECVR;
    return(IDASLS_JACFUNC_UNRECVR);
  }
  if (retval > 0) {
    last_flag = IDASLS_JACFUNC_RECVR;
    return(+1);
  }

  if (idasls_mem->s_first_factorize) {
    /* ------------------------------------------------------------
       Get the symbolic factorization
       ------------------------------------------------------------*/ 
    /* Update the ordering option with any user-updated values from 
       calls to IDAKLUSetOrdering */
    klu_data->s_Common.ordering = klu_data->s_ordering;

    klu_data->s_Symbolic = klu_analyze(JacMat->N, JacMat->colptrs, 
				       JacMat->rowvals, &(klu_data->s_Common));
    if (klu_data->s_Symbolic == NULL) {
      IDAProcessError(IDA_mem, IDASLS_PACKAGE_FAIL, "IDASSLS", "IDAKLUSetup", 
		      MSGSP_PACKAGE_FAIL);
      return(IDASLS_PACKAGE_FAIL);
    }

    idasls_mem->s_first_factorize = 0;
  }

  /* ------------------------------------------------------------
     Compute the LU factorization of  the Jacobian.
     ------------------------------------------------------------*/
  klu_data->s_Numeric = klu_factor(JacMat->colptrs, JacMat->rowvals, JacMat->data, 
				   klu_data->s_Symbolic, &(klu_data->s_Common));

  if (klu_data->s_Numeric == NULL) {
    IDAProcessError(IDA_mem, IDASLS_PACKAGE_FAIL, "IDASSLS", "IDAKLUSetup", 
		    MSGSP_PACKAGE_FAIL);
    return(IDASLS_PACKAGE_FAIL);
  }


  last_flag = IDASLS_SUCCESS;

  return(0);
}

/*
  This routine handles the solve operation for the IDAKLU linear
  solver module.  It calls the KLU solve routine, scales the
  solution vector according to cjratio, then returns IDASLS_SUCCESS = 0.
*/

static int IDAKLUSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
		       N_Vector ycur, N_Vector ypcur, N_Vector rrcur)
{
  int last_flag, flag;
  realtype cjratio;
  IDASlsMem idasls_mem;
  KLUData klu_data;
  SlsMat JacMat;
  realtype *bd;
  
  idasls_mem = (IDASlsMem) IDA_mem->ida_lmem;
  JacMat = idasls_mem->s_JacMat;
  cjratio = IDA_mem->ida_cjratio;

  klu_data = (KLUData) idasls_mem->s_solver_data;
  last_flag = idasls_mem->s_last_flag;

  bd = N_VGetArrayPointer(b);

  /* Call KLU to solve the linear system */
  flag = klu_solve(klu_data->s_Symbolic, klu_data->s_Numeric, JacMat->N, 1, bd, 
	    &(klu_data->s_Common));
  if (flag == 0) {
    IDAProcessError(IDA_mem, IDASLS_PACKAGE_FAIL, "IDASSLS", "IDAKLUSolve", 
		    MSGSP_PACKAGE_FAIL);
    return(IDASLS_PACKAGE_FAIL);
  }

  /* Scale the correction to account for change in cj. */
  if (cjratio != ONE) N_VScale(TWO/(ONE + cjratio), b, b);

  last_flag = IDASLS_SUCCESS;
  return(IDASLS_SUCCESS);
}

/*
  This routine frees memory specific to the IDAKLU linear solver.
*/

static int IDAKLUFree(IDAMem IDA_mem)
{
  IDASlsMem idasls_mem;
  KLUData klu_data;
  
  idasls_mem = (IDASlsMem) IDA_mem->ida_lmem;
  klu_data = (KLUData) idasls_mem->s_solver_data;

  klu_free_numeric(&(klu_data->s_Numeric), &(klu_data->s_Common));
  klu_free_symbolic(&(klu_data->s_Symbolic), &(klu_data->s_Common));

  free(klu_data->s_Symbolic);

  if (idasls_mem->s_JacMat) {
    DestroySparseMat(idasls_mem->s_JacMat);
    idasls_mem->s_JacMat = NULL;
  }

  free(klu_data); 
  free(IDA_mem->ida_lmem); 

  return(IDASLS_SUCCESS);
}

/* 
 * -----------------------------------------------------------------
 * Optional Input Specification Functions
 * -----------------------------------------------------------------
 *
 * IDAKLUSetOrdering sets the ordering used by KLU for reducing fill.
 * Options are: 0 for AMD, 1 for COLAMD, and 2 for the natural ordering.
 * The default used in IDA is 1 for COLAMD.
 * -----------------------------------------------------------------
 */

int IDAKLUSetOrdering(void *ida_mem_v, int ordering_choice)
{
  IDAMem ida_mem;
  IDASlsMem idasls_mem;
  KLUData klu_data;

 /* Return immediately if ida_mem is NULL */
  if (ida_mem_v == NULL) {
    IDAProcessError(NULL, IDASLS_MEM_NULL, "IDASLS", "IDAKLUSetOrdering",
		    MSGSP_IDAMEM_NULL);
    return(IDASLS_MEM_NULL);
  }
  ida_mem = (IDAMem) ida_mem_v;

 /* Return if ordering choice argument is not valid */
  if ( (ordering_choice != 0) && (ordering_choice != 1) && 
       (ordering_choice != 2) ) {
    IDAProcessError(NULL, IDASLS_ILL_INPUT, "IDASLS", "IDAKLUSetOrdering",
		    MSGSP_ILL_INPUT);
    return(IDASLS_ILL_INPUT);
  }

  idasls_mem = (IDASlsMem) ida_mem->ida_lmem;
  klu_data = (KLUData) idasls_mem->s_solver_data;

  klu_data->s_ordering = ordering_choice;

  return(IDASLS_SUCCESS);
}
