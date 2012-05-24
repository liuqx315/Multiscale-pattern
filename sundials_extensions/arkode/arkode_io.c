/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
----------------------------------------------------------------- 
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 This is the implementation file for the optional input and 
 output functions for the ARKODE solver.
---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include <sundials/sundials_types.h>


#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)



/*===============================================================
 ARKODE optional input functions
===============================================================*/

/*---------------------------------------------------------------
 ARKodeSetDefaults:

 Resets all optional inputs to ARKode default values.  Does not 
 change problem-defining function pointers fe and fi or 
 user_data pointer.  Also leaves alone any data 
 structures/options related to root-finding (those can be reset 
 using ARKodeRootInit).
---------------------------------------------------------------*/
int ARKodeSetDefaults(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetDefaults", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Set default values for integrator optional inputs */
  ark_mem->ark_expstab     = ARKExpStab;
  ark_mem->ark_estab_data  = ark_mem;
  ark_mem->ark_hadapt      = ARKAdapt;
  ark_mem->ark_hadapt_data = ark_mem;
  ark_mem->ark_itol        = ARK_NN;
  ark_mem->ark_user_efun   = FALSE;
  ark_mem->ark_linear      = FALSE;
  ark_mem->ark_explicit    = FALSE;
  ark_mem->ark_implicit    = FALSE;
  ark_mem->ark_user_Ae     = FALSE;
  ark_mem->ark_user_Ai     = FALSE;
  ark_mem->ark_efun        = NULL;
  ark_mem->ark_e_data      = NULL;
  ark_mem->ark_ehfun       = ARKErrHandler;
  ark_mem->ark_eh_data     = ark_mem;
  ark_mem->ark_errfp       = stderr;
  ark_mem->ark_qmax        = Q_MAX;
  ark_mem->ark_mxstep      = MXSTEP_DEFAULT;
  ark_mem->ark_mxhnil      = 10;
  ark_mem->ark_hin         = ZERO;
  ark_mem->ark_hmin        = RCONST(0.0);
  ark_mem->ark_hmax_inv    = RCONST(0.0);
  ark_mem->ark_tstopset    = FALSE;
  ark_mem->ark_maxcor      = 3;
  ark_mem->ark_maxnef      = 7;
  ark_mem->ark_maxncf      = 10;
  ark_mem->ark_nlscoef     = RCONST(0.1);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetErrHandlerFn:

 Specifies the error handler function
---------------------------------------------------------------*/
int ARKodeSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun, void *eh_data)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetErrHandlerFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_ehfun = ehfun;
  ark_mem->ark_eh_data = eh_data;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetErrFile:

 Specifies the FILE pointer for output (NULL means no messages)
---------------------------------------------------------------*/
int ARKodeSetErrFile(void *arkode_mem, FILE *errfp)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetErrFile", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_errfp = errfp;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetUserData:

 Specifies the user data pointer for f
---------------------------------------------------------------*/
int ARKodeSetUserData(void *arkode_mem, void *user_data)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetUserData", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_user_data = user_data;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetOrd:

 Specifies the method order
---------------------------------------------------------------*/
int ARKodeSetOrd(void *arkode_mem, int ord)
{
  ARKodeMem ark_mem;
  int qmax_alloc;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxOrd", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  if (ord <= 0) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetMaxOrd", MSGARK_NEG_MAXORD);
    return(ARK_ILL_INPUT);
  }
  
  /* Cannot increase maximum order beyond the value that
     was used when allocating memory */
  qmax_alloc = ark_mem->ark_qmax_alloc;

  if (ord > qmax_alloc) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetMaxOrd", MSGARK_BAD_MAXORD);
    return(ARK_ILL_INPUT);
  }

  /***** CHANGE THIS TO ARK_MEM->ARK_Q = ORD WHEN WE TRANSITION TO ARK METHOD!!! *****/
  ark_mem->ark_qmax = ord;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetLinear:

 Specifies that the implicit portion of the problem is linear, 
 and to tighten the linear solver tolerances while taking only 
 one Newton iteration.
---------------------------------------------------------------*/
int ARKodeSetLinear(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxOrd", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->ark_linear = TRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetERK:

 Specifies that the implicit portion of the problem is disabled, 
 and to use an explicit RK method.
---------------------------------------------------------------*/
int ARKodeSetERK(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetERK", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* ensure that fe is defined */
  if (ark_mem->ark_fe == NULL) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetERK", MSGARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  ark_mem->ark_explicit = TRUE;
  ark_mem->ark_implicit = FALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetIRK:

 Specifies that the explicit portion of the problem is disabled, 
 and to use an implicit RK method.
---------------------------------------------------------------*/
int ARKodeSetIRK(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetIRK", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* ensure that fi is defined */
  if (ark_mem->ark_fi == NULL) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetIRK", MSGARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  ark_mem->ark_implicit = TRUE;
  ark_mem->ark_explicit = FALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetERKTable:

 Specifies to use a customized Butcher table for the explicit 
 portion of the system.  Of these, the only optional argument is 
 bdense (which can be NULL).
---------------------------------------------------------------*/
int ARKodeSetERKTable(void *arkode_mem, int s, int q, int p,
		      realtype *c, realtype **A, realtype *b, 
		      realtype *bembed, realtype **bdense)
{
  int i, j;
  ARKodeMem ark_mem;
  booleantype match;
  realtype tol = RCONST(1.0e-12);
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetERKTable", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for legal inputs */
  if (s > S_MAX) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetERKTable", "s exceeds S_MAX");
    return(ARK_ILL_INPUT);
  }
  if ((c == NULL) || (A == NULL) || (b == NULL) || (bembed == NULL)) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetERKTable", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* if IRK table already set, ensure that shared coeffs match */
  if (ark_mem->ark_user_Ai) {
    match = TRUE;
    if (ark_mem->ark_stages != s)  match = FALSE;
    if (!match) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKodeSetERKTable", "s does not match Ai table");
      return(ARK_ILL_INPUT);
    }
    if (ark_mem->ark_q != q)  match = FALSE;
    if (ark_mem->ark_p != p)  match = FALSE;
    for (i=0; i<s; i++) {
      if (ABS(ark_mem->ark_c[i]  - c[i])      > tol)  match = FALSE;
      if (ABS(ark_mem->ark_b[i]  - b[i])      > tol)  match = FALSE;
      if (ABS(ark_mem->ark_b2[i] - bembed[i]) > tol)  match = FALSE;
    }
    if ((ark_mem->ark_bd[0][0] != ZERO) && (bdense != NULL)) {
      for (i=0; i<s; i++) 
	for (j=0; j<s; j++) 
	  if (ABS(ark_mem->ark_bd[i][j] - bdense[i][j]) > tol)  match = FALSE;
    }
    if (!match) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKodeSetERKTable", "shared Butcher coeffs don't match");
      return(ARK_ILL_INPUT);
    }
  }

  /* set the relevant parameters */
  ark_mem->ark_stages = s;
  ark_mem->ark_q = q;
  ark_mem->ark_p = p;
  for (i=0; i<s; i++) {
    ark_mem->ark_c[i]  = c[i];
    ark_mem->ark_b[i]  = b[i];
    ark_mem->ark_b2[i] = bembed[i];
    for (j=0; j<s; j++) {
      ark_mem->ark_Ae[i][j] = A[i][j];
    }
  }

  /* set the dense coefficients (if supplied) */
  if (bdense != NULL) {
    for (i=0; i<s; i++) {
      for (j=0; j<s; j++) {
	ark_mem->ark_bd[i][j] = bdense[i][j];
      }
    }
  }

  /* remark that this table was supplied by the user */
  ark_mem->ark_user_Ae = TRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetIRKTable:

 Specifies to use a customized Butcher table for the implicit 
 portion of the system.
---------------------------------------------------------------*/
int ARKodeSetIRKTable(void *arkode_mem, int s, int q, int p,
		      realtype *c, realtype **A, realtype *b, 
		      realtype *bembed, realtype **bdense)
{
  int i, j;
  ARKodeMem ark_mem;
  booleantype match;
  realtype tol = RCONST(1.0e-12);
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetIRKTable", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for legal inputs */
  if (s > S_MAX) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetIRKTable", "s exceeds S_MAX");
    return(ARK_ILL_INPUT);
  }
  if ((c == NULL) || (A == NULL) || (b == NULL) || (bembed == NULL)) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetIRKTable", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* if ERK table already set, ensure that shared coeffs match */
  if (ark_mem->ark_user_Ae) {
    match = TRUE;
    if (ark_mem->ark_stages != s)  match = FALSE;
    if (!match) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKodeSetIRKTable", "s does not match Ae table");
      return(ARK_ILL_INPUT);
    }
    if (ark_mem->ark_q != q)  match = FALSE;
    if (ark_mem->ark_p != p)  match = FALSE;
    for (i=0; i<s; i++) {
      if (ABS(ark_mem->ark_c[i]  - c[i])      > tol)  match = FALSE;
      if (ABS(ark_mem->ark_b[i]  - b[i])      > tol)  match = FALSE;
      if (ABS(ark_mem->ark_b2[i] - bembed[i]) > tol)  match = FALSE;
    }
    if ((ark_mem->ark_bd[0][0] != ZERO) && (bdense != NULL)) {
      for (i=0; i<s; i++) 
	for (j=0; j<s; j++) 
	  if (ABS(ark_mem->ark_bd[i][j] - bdense[i][j]) > tol)  match = FALSE;
    }
    if (!match) {
      ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		      "ARKodeSetIRKTable", "shared Butcher coeffs don't match");
      return(ARK_ILL_INPUT);
    }
  }

  /* set the relevant parameters */
  ark_mem->ark_stages = s;
  ark_mem->ark_q = q;
  ark_mem->ark_p = p;
  for (i=0; i<s; i++) {
    ark_mem->ark_c[i]  = c[i];
    ark_mem->ark_b[i]  = b[i];
    ark_mem->ark_b2[i] = bembed[i];
    for (j=0; j<s; j++) {
      ark_mem->ark_Ai[i][j] = A[i][j];
    }
  }

  /* set the dense coefficients (if supplied) */
  if (bdense != NULL) {
    for (i=0; i<s; i++) {
      for (j=0; j<s; j++) {
	ark_mem->ark_bd[i][j] = bdense[i][j];
      }
    }
  }

  /* remark that this table was supplied by the user */
  ark_mem->ark_user_Ai = TRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxNumSteps:

 Specifies the maximum number of integration steps
---------------------------------------------------------------*/
int ARKodeSetMaxNumSteps(void *arkode_mem, long int mxsteps)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxNumSteps", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  /* Passing mxsteps=0 sets the default. Passing mxsteps<0 disables the test. */
  if (mxsteps == 0)
    ark_mem->ark_mxstep = MXSTEP_DEFAULT;
  else
    ark_mem->ark_mxstep = mxsteps;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxHnilWarns:

 Specifies the maximum number of warnings for small h
---------------------------------------------------------------*/
int ARKodeSetMaxHnilWarns(void *arkode_mem, int mxhnil)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxHnilWarns", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_mxhnil = mxhnil;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetInitStep:

 Specifies the initial step size
---------------------------------------------------------------*/
int ARKodeSetInitStep(void *arkode_mem, realtype hin)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetInitStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_hin = hin;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMinStep:

 Specifies the minimum step size
---------------------------------------------------------------*/
int ARKodeSetMinStep(void *arkode_mem, realtype hmin)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMinStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  if (hmin<0) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetMinStep", MSGARK_NEG_HMIN);
    return(ARK_ILL_INPUT);
  }

  /* Passing 0 sets hmin = zero */
  if (hmin == ZERO) {
    ark_mem->ark_hmin = RCONST(0.0);
    return(ARK_SUCCESS);
  }

  if (hmin * ark_mem->ark_hmax_inv > ONE) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetMinStep", MSGARK_BAD_HMIN_HMAX);
    return(ARK_ILL_INPUT);
  }

  ark_mem->ark_hmin = hmin;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxStep:

 Specifies the maximum step size
---------------------------------------------------------------*/
int ARKodeSetMaxStep(void *arkode_mem, realtype hmax)
{
  realtype hmax_inv;
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxStep", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  if (hmax < 0) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetMaxStep", MSGARK_NEG_HMAX);
    return(ARK_ILL_INPUT);
  }

  /* Passing 0 sets hmax = infinity */
  if (hmax == ZERO) {
    ark_mem->ark_hmax_inv = RCONST(0.0);
    return(ARK_SUCCESS);
  }

  hmax_inv = ONE/hmax;
  if (hmax_inv * ark_mem->ark_hmin > ONE) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetMaxStep", MSGARK_BAD_HMIN_HMAX);
    return(ARK_ILL_INPUT);
  }

  ark_mem->ark_hmax_inv = hmax_inv;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetStopTime:

 Specifies the time beyond which the integration is not to proceed.
---------------------------------------------------------------*/
int ARKodeSetStopTime(void *arkode_mem, realtype tstop)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetStopTime", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* If ARKode was called at least once, test if tstop is legal
   * (i.e. if it was not already passed).
   * If ARKodeSetStopTime is called before the first call to ARKode,
   * tstop will be checked in ARKode. */
  if (ark_mem->ark_nst > 0) {

    if ( (tstop - ark_mem->ark_tn) * ark_mem->ark_h < ZERO ) {
      ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "ARKodeSetStopTime", MSGARK_BAD_TSTOP, ark_mem->ark_tn);
      return(ARK_ILL_INPUT);
    }

  }

  ark_mem->ark_tstop = tstop;
  ark_mem->ark_tstopset = TRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetAdaptMethod:

 Specifies the time step adaptivity algorithm (and associated 
 parameters) to use.
---------------------------------------------------------------*/
int ARKodeSetAdaptMethod(void *arkode_mem, int imethod, 
			 realtype *adapt_params)
{
  /* FILL THIS IN!!!! */
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetAdaptivityFn:

 Specifies the user-provided time step adaptivity function to use.
---------------------------------------------------------------*/
int ARKodeSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun,
			  void *h_data)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetAdaptivityFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_hadapt = hfun;
  ark_mem->ark_hadapt_data = h_data;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetStabilityFn:

 Specifies the user-provided explicit time step stability 
 function to use.
---------------------------------------------------------------*/
int ARKodeSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab,
			 void *estab_data)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetStabilityFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_expstab = EStab;
  ark_mem->ark_estab_data = estab_data;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxErrTestFails:

 Specifies the maximum number of error test failures during one
 step try.
---------------------------------------------------------------*/
int ARKodeSetMaxErrTestFails(void *arkode_mem, int maxnef)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxErrTestFails", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_maxnef = maxnef;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxConvFails:

 Specifies the maximum number of nonlinear convergence failures 
 during one step try.
---------------------------------------------------------------*/
int ARKodeSetMaxConvFails(void *arkode_mem, int maxncf)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxConvFails", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_maxncf = maxncf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxNonlinIters:

 Specifies the maximum number of nonlinear iterations during
 one solve.
---------------------------------------------------------------*/
int ARKodeSetMaxNonlinIters(void *arkode_mem, int maxcor)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxNonlinIters", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_maxcor = maxcor;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetNonlinConvCoef:

 Specifies the coeficient in the nonlinear solver convergence
 test
---------------------------------------------------------------*/
int ARKodeSetNonlinConvCoef(void *arkode_mem, realtype nlscoef)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetNonlinConvCoef", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_nlscoef = nlscoef;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetRootDirection:

 Specifies the direction of zero-crossings to be monitored.
 The default is to monitor both crossings.
---------------------------------------------------------------*/
int ARKodeSetRootDirection(void *arkode_mem, int *rootdir)
{
  ARKodeMem ark_mem;
  int i, nrt;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetRootDirection", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  nrt = ark_mem->ark_nrtfn;
  if (nrt==0) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetRootDirection", MSGARK_NO_ROOT);
    return(ARK_ILL_INPUT);    
  }

  for(i=0; i<nrt; i++) ark_mem->ark_rootdir[i] = rootdir[i];

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetNoInactiveRootWarn:

 Disables issuing a warning if some root function appears
 to be identically zero at the beginning of the integration
---------------------------------------------------------------*/
int ARKodeSetNoInactiveRootWarn(void *arkode_mem)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetNoInactiveRootWarn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_mxgnull = 0;
  
  return(ARK_SUCCESS);
}


/*===============================================================
 ARKODE optional output functions
===============================================================*/

/*---------------------------------------------------------------
 ARKodeGetNumSteps:

 Returns the current number of integration steps
---------------------------------------------------------------*/
int ARKodeGetNumSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumSteps", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->ark_nst;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumRhsEvals:

 Returns the current number of calls to f
---------------------------------------------------------------*/
int ARKodeGetNumRhsEvals(void *arkode_mem, long int *nfevals)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumRhsEvals", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *nfevals = ark_mem->ark_nfe;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumLinSolvSetups:

 Returns the current number of calls to the linear solver setup routine
---------------------------------------------------------------*/
int ARKodeGetNumLinSolvSetups(void *arkode_mem, long int *nlinsetups)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumLinSolvSetups", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *nlinsetups = ark_mem->ark_nsetups;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumErrTestFails:

 Returns the current number of error test failures
---------------------------------------------------------------*/
int ARKodeGetNumErrTestFails(void *arkode_mem, long int *netfails)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumErrTestFails", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *netfails = ark_mem->ark_netf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetLastOrder:

 Returns the order on the last succesful step
---------------------------------------------------------------*/
int ARKodeGetLastOrder(void *arkode_mem, int *qlast)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetLastOrder", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *qlast = ark_mem->ark_qu;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetCurrentOrder:

 Returns the order to be attempted on the next step
---------------------------------------------------------------*/
int ARKodeGetCurrentOrder(void *arkode_mem, int *qcur)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetCurrentOrder", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *qcur = ark_mem->ark_next_q;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetActualInitStep:

 Returns the step size used on the first step
---------------------------------------------------------------*/
int ARKodeGetActualInitStep(void *arkode_mem, realtype *hinused)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetActualInitStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *hinused = ark_mem->ark_h0u;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetLastStep:

 Returns the step size used on the last successful step
---------------------------------------------------------------*/
int ARKodeGetLastStep(void *arkode_mem, realtype *hlast)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetLastStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *hlast = ark_mem->ark_hu;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetCurrentStep:

 Returns the step size to be attempted on the next step
---------------------------------------------------------------*/
int ARKodeGetCurrentStep(void *arkode_mem, realtype *hcur)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetCurrentStep", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;
  
  *hcur = ark_mem->ark_next_h;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetCurrentTime:

 Returns the current value of the independent variable
---------------------------------------------------------------*/
int ARKodeGetCurrentTime(void *arkode_mem, realtype *tcur)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetCurrentTime", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *tcur = ark_mem->ark_tn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetTolScaleFactor:

 Returns a suggested factor for scaling tolerances
---------------------------------------------------------------*/
int ARKodeGetTolScaleFactor(void *arkode_mem, realtype *tolsfact)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetTolScaleFactor", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *tolsfact = ark_mem->ark_tolsf;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetErrWeights:

 This routine returns the current weight vector.
---------------------------------------------------------------*/
int ARKodeGetErrWeights(void *arkode_mem, N_Vector eweight)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetErrWeights", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  N_VScale(ONE, ark_mem->ark_ewt, eweight);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetEstLocalErrors:

 Returns an estimate of the local error
---------------------------------------------------------------*/
int ARKodeGetEstLocalErrors(void *arkode_mem, N_Vector ele)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetEstLocalErrors", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  N_VScale(ONE, ark_mem->ark_acor, ele);

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetWorkSpace:

 Returns integrator work space requirements
---------------------------------------------------------------*/
int ARKodeGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetWorkSpace", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *leniw = ark_mem->ark_liw;
  *lenrw = ark_mem->ark_lrw;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetIntegratorStats:

 Returns integrator statistics
---------------------------------------------------------------*/
int ARKodeGetIntegratorStats(void *arkode_mem, long int *nsteps, 
			     long int *expsteps, long int *accsteps, 
			     long int *convsteps, long int *nfevals, 
			     long int *nlinsetups, long int *netfails, 
			     realtype *hinused, realtype *hlast, 
			     realtype *hcur, realtype *tcur)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetIntegratorStats", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->ark_nst;
  *nfevals = ark_mem->ark_nfe;
  *nlinsetups = ark_mem->ark_nsetups;
  *netfails = ark_mem->ark_netf;
  *hinused = ark_mem->ark_h0u;
  *hlast = ark_mem->ark_hu;
  *hcur = ark_mem->ark_next_h;
  *tcur = ark_mem->ark_tn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumGEvals:

 Returns the current number of calls to g (for rootfinding)
---------------------------------------------------------------*/
int ARKodeGetNumGEvals(void *arkode_mem, long int *ngevals)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumGEvals", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *ngevals = ark_mem->ark_nge;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetRootInfo:

 Returns pointer to array rootsfound showing roots found
---------------------------------------------------------------*/
int ARKodeGetRootInfo(void *arkode_mem, int *rootsfound)
{
  ARKodeMem ark_mem;
  int i, nrt;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetRootInfo", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  nrt = ark_mem->ark_nrtfn;

  for (i=0; i<nrt; i++) rootsfound[i] = ark_mem->ark_iroots[i];

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumNonlinSolvIters:

 Returns the current number of iterations in the nonlinear solver
---------------------------------------------------------------*/
int ARKodeGetNumNonlinSolvIters(void *arkode_mem, long int *nniters)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumNonlinSolvIters", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *nniters = ark_mem->ark_nni;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumNonlinSolvConvFails:

 Returns the current number of convergence failures in the
 nonlinear solver
---------------------------------------------------------------*/
int ARKodeGetNumNonlinSolvConvFails(void *arkode_mem, long int *nncfails)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumNonlinSolvConvFails", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *nncfails = ark_mem->ark_ncfn;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNonlinSolvStats:

 Returns nonlinear solver statistics
---------------------------------------------------------------*/
int ARKodeGetNonlinSolvStats(void *arkode_mem, long int *nniters, 
                            long int *nncfails)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNonlinSolvStats", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;

  *nniters = ark_mem->ark_nni;
  *nncfails = ark_mem->ark_ncfn;

  return(ARK_SUCCESS);
}

/*-----------------------------------------------------------------*/

char *ARKodeGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(24*sizeof(char));

  switch(flag) {
  case ARK_SUCCESS:
    sprintf(name,"ARK_SUCCESS");
    break;
  case ARK_TSTOP_RETURN:
    sprintf(name,"ARK_TSTOP_RETURN");
    break;
  case ARK_ROOT_RETURN:
    sprintf(name,"ARK_ROOT_RETURN");
    break;
  case ARK_TOO_MUCH_WORK:
    sprintf(name,"ARK_TOO_MUCH_WORK");
    break;
  case ARK_TOO_MUCH_ACC:
    sprintf(name,"ARK_TOO_MUCH_ACC");
    break;
  case ARK_ERR_FAILURE:
    sprintf(name,"ARK_ERR_FAILURE");
    break;
  case ARK_CONV_FAILURE:
    sprintf(name,"ARK_CONV_FAILURE");
    break;
  case ARK_LINIT_FAIL:
    sprintf(name,"ARK_LINIT_FAIL");
    break;
  case ARK_LSETUP_FAIL:
    sprintf(name,"ARK_LSETUP_FAIL");
    break;
  case ARK_LSOLVE_FAIL:
    sprintf(name,"ARK_LSOLVE_FAIL");
    break;
  case ARK_RHSFUNC_FAIL:
    sprintf(name,"ARK_RHSFUNC_FAIL");
    break;
  case ARK_FIRST_RHSFUNC_ERR:
    sprintf(name,"ARK_FIRST_RHSFUNC_ERR");
    break;
  case ARK_REPTD_RHSFUNC_ERR:
    sprintf(name,"ARK_REPTD_RHSFUNC_ERR");
    break;
  case ARK_UNREC_RHSFUNC_ERR:
    sprintf(name,"ARK_UNREC_RHSFUNC_ERR");
    break;
  case ARK_RTFUNC_FAIL:
    sprintf(name,"ARK_RTFUNC_FAIL");
    break;
  case ARK_MEM_FAIL:
    sprintf(name,"ARK_MEM_FAIL");
    break;
  case ARK_MEM_NULL:
    sprintf(name,"ARK_MEM_NULL");
    break;
  case ARK_ILL_INPUT:
    sprintf(name,"ARK_ILL_INPUT");
    break;
  case ARK_NO_MALLOC:
    sprintf(name,"ARK_NO_MALLOC");
    break;
  case ARK_BAD_K:
    sprintf(name,"ARK_BAD_K");
    break;
  case ARK_BAD_T:
    sprintf(name,"ARK_BAD_T");
    break;
  case ARK_BAD_DKY:
    sprintf(name,"ARK_BAD_DKY");
    break;
  case ARK_TOO_CLOSE:
    sprintf(name,"ARK_TOO_CLOSE");
    break;    
  default:
    sprintf(name,"NONE");
  }

  return(name);
}


/*---------------------------------------------------------------
      EOF
---------------------------------------------------------------*/
