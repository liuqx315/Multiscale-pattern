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
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>


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
  int i, j;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetDefaults", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Set default values for integrator optional inputs */
  ark_mem->ark_dense_q          = QDENSE_DEF;
  ark_mem->ark_expstab          = ARKExpStab;
  ark_mem->ark_estab_data       = ark_mem;
  ark_mem->ark_hadapt           = NULL;
  ark_mem->ark_hadapt_data      = NULL;
  ark_mem->ark_hadapt_cfl       = CFLFAC;
  ark_mem->ark_hadapt_safety    = SAFETY;
  ark_mem->ark_hadapt_bias      = BIAS;
  ark_mem->ark_hadapt_growth    = GROWTH;
  ark_mem->ark_hadapt_lbound    = HFIXED_LB;
  ark_mem->ark_hadapt_ubound    = HFIXED_UB;
  ark_mem->ark_hadapt_imethod   = 0;
  ark_mem->ark_hadapt_k1        = AD0_K1;
  ark_mem->ark_hadapt_k2        = AD0_K2;
  ark_mem->ark_hadapt_k3        = AD0_K3;
  ark_mem->ark_predictor        = 3;
  ark_mem->ark_reltol           = 1.e-4;
  ark_mem->ark_Sabstol          = 1.e-9;
  ark_mem->ark_itol             = ARK_SS;
  ark_mem->ark_user_efun        = FALSE;
  ark_mem->ark_efun             = ARKEwtSet;
  ark_mem->ark_e_data           = NULL;
  ark_mem->ark_linear           = FALSE;
  ark_mem->ark_explicit         = FALSE;
  ark_mem->ark_implicit         = FALSE;
  ark_mem->ark_ehfun            = ARKErrHandler;
  ark_mem->ark_eh_data          = ark_mem;
  ark_mem->ark_errfp            = stderr;
  ark_mem->ark_q                = Q_DEFAULT;
  ark_mem->ark_mxstep           = MXSTEP_DEFAULT;
  ark_mem->ark_mxhnil           = MXHNIL;
  ark_mem->ark_hin              = ZERO;
  ark_mem->ark_hmin             = ZERO;
  ark_mem->ark_hmax_inv         = ZERO;
  ark_mem->ark_tstopset         = ZERO;
  ark_mem->ark_tstopset         = FALSE;
  ark_mem->ark_maxcor           = MAXCOR;
  ark_mem->ark_maxnef           = MAXNEF;
  ark_mem->ark_maxncf           = MAXNCF;
  ark_mem->ark_nlscoef          = NLSCOEF;
  ark_mem->ark_etamx1           = ETAMX1;
  ark_mem->ark_etamxf           = ETAMXF;
  ark_mem->ark_small_nef        = SMALL_NEF;
  ark_mem->ark_etacf            = ETACF;
  ark_mem->ark_crdown           = CRDOWN;
  ark_mem->ark_rdiv             = RDIV;
  ark_mem->ark_dgmax            = DGMAX;
  ark_mem->ark_msbp             = MSBP;
  ark_mem->ark_diagfp           = NULL;
  ark_mem->ark_report           = FALSE;
  ark_mem->ark_stages           = 0;
  ark_mem->ark_istage           = 0;
  ark_mem->ark_p                = 0;
  for (i=0; i<ARK_S_MAX; i++) {
    for (j=0; j<ARK_S_MAX; j++) {
      ARK_A(ark_mem->ark_Ae,i,j) = ZERO;
      ARK_A(ark_mem->ark_Ai,i,j) = ZERO;
    }
    ark_mem->ark_c[i]   = ZERO;
    ark_mem->ark_b[i]   = ZERO;
    ark_mem->ark_b2[i]  = ZERO;
  }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetOptimalParams:

 Sets all adaptivity and solver parameters to our 'best guess' 
 values, for a given integration method (ERK, DIRK, ARK) and a 
 given method order.  Should only be called after the method
 order and integration method have been set.
---------------------------------------------------------------*/
int ARKodeSetOptimalParams(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetDefaults", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* Choose values based on method, order */

  /*    explicit */
  if (ark_mem->ark_explicit) {
    ark_mem->ark_hadapt_imethod = 1;
    ark_mem->ark_hadapt_safety  = RCONST(0.99);
    ark_mem->ark_hadapt_bias    = RCONST(1.2);
    ark_mem->ark_hadapt_growth  = RCONST(25.0);
    ark_mem->ark_hadapt_k1      = RCONST(0.8);
    ark_mem->ark_hadapt_k2      = RCONST(0.31);
    ark_mem->ark_etamxf         = RCONST(0.3);

  /*    implicit */
  } else if (ark_mem->ark_implicit) {
    switch (ark_mem->ark_q) {
    case 3:
      ark_mem->ark_hadapt_imethod   = 2;
      ark_mem->ark_hadapt_safety    = RCONST(0.957);
      ark_mem->ark_hadapt_bias      = RCONST(1.9);
      ark_mem->ark_hadapt_growth    = RCONST(17.6);
      ark_mem->ark_etamxf           = RCONST(0.45);
      ark_mem->ark_nlscoef          = RCONST(0.22);
      ark_mem->ark_crdown           = RCONST(0.17);
      ark_mem->ark_rdiv             = RCONST(2.3);
      ark_mem->ark_dgmax            = RCONST(0.19);
      ark_mem->ark_msbp             = 60;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    case 4:
      ark_mem->ark_hadapt_imethod   = 0;
      ark_mem->ark_hadapt_safety    = RCONST(0.988);
      ark_mem->ark_hadapt_bias      = RCONST(1.2);
      ark_mem->ark_hadapt_growth    = RCONST(31.5);
      ark_mem->ark_hadapt_k1        = RCONST(0.535);
      ark_mem->ark_hadapt_k2        = RCONST(0.209);
      ark_mem->ark_hadapt_k3        = RCONST(0.148);
      ark_mem->ark_etamxf           = RCONST(0.33);
      ark_mem->ark_nlscoef          = RCONST(0.24);
      ark_mem->ark_crdown           = RCONST(0.26);
      ark_mem->ark_rdiv             = RCONST(2.3);
      ark_mem->ark_dgmax            = RCONST(0.16);
      ark_mem->ark_msbp             = 31;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    case 5:
      ark_mem->ark_hadapt_imethod   = 0;
      ark_mem->ark_hadapt_safety    = RCONST(0.937);
      ark_mem->ark_hadapt_bias      = RCONST(3.3);
      ark_mem->ark_hadapt_growth    = RCONST(22.0);
      ark_mem->ark_hadapt_k1        = RCONST(0.56);
      ark_mem->ark_hadapt_k2        = RCONST(0.338);
      ark_mem->ark_hadapt_k3        = RCONST(0.14);
      ark_mem->ark_etamxf           = RCONST(0.44);
      ark_mem->ark_nlscoef          = RCONST(0.25);
      ark_mem->ark_crdown           = RCONST(0.4);
      ark_mem->ark_rdiv             = RCONST(2.3);
      ark_mem->ark_dgmax            = RCONST(0.32);
      ark_mem->ark_msbp             = 31;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    }

  /*    imex */
  } else {
    switch (ark_mem->ark_q) {
    case 3:
      ark_mem->ark_hadapt_imethod   = 0;
      ark_mem->ark_hadapt_safety    = RCONST(0.965);
      ark_mem->ark_hadapt_bias      = RCONST(1.42);
      ark_mem->ark_hadapt_growth    = RCONST(28.7);
      ark_mem->ark_hadapt_k1        = RCONST(0.54);
      ark_mem->ark_hadapt_k2        = RCONST(0.36);
      ark_mem->ark_hadapt_k3        = RCONST(0.14);
      ark_mem->ark_etamxf           = RCONST(0.46);
      ark_mem->ark_nlscoef          = RCONST(0.22);
      ark_mem->ark_crdown           = RCONST(0.17);
      ark_mem->ark_rdiv             = RCONST(2.3);
      ark_mem->ark_dgmax            = RCONST(0.19);
      ark_mem->ark_msbp             = 60;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    case 4:
      ark_mem->ark_hadapt_imethod   = 0;
      ark_mem->ark_hadapt_safety    = RCONST(0.97);
      ark_mem->ark_hadapt_bias      = RCONST(1.35);
      ark_mem->ark_hadapt_growth    = RCONST(25.0);
      ark_mem->ark_hadapt_k1        = RCONST(0.543);
      ark_mem->ark_hadapt_k2        = RCONST(0.297);
      ark_mem->ark_hadapt_k3        = RCONST(0.14);
      ark_mem->ark_etamxf           = RCONST(0.47);
      ark_mem->ark_nlscoef          = RCONST(0.24);
      ark_mem->ark_crdown           = RCONST(0.26);
      ark_mem->ark_rdiv             = RCONST(2.3);
      ark_mem->ark_dgmax            = RCONST(0.16);
      ark_mem->ark_msbp             = 31;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    case 5:
      ark_mem->ark_hadapt_imethod   = 1;
      ark_mem->ark_hadapt_safety    = RCONST(0.993);
      ark_mem->ark_hadapt_bias      = RCONST(1.15);
      ark_mem->ark_hadapt_growth    = RCONST(28.5);
      ark_mem->ark_hadapt_k1        = RCONST(0.8);
      ark_mem->ark_hadapt_k2        = RCONST(0.35);
      ark_mem->ark_etamxf           = RCONST(0.3);
      ark_mem->ark_nlscoef          = RCONST(0.25);
      ark_mem->ark_crdown           = RCONST(0.4);
      ark_mem->ark_rdiv             = RCONST(2.3);
      ark_mem->ark_dgmax            = RCONST(0.32);
      ark_mem->ark_msbp             = 31;
      ark_mem->ark_small_nef        = SMALL_NEF;
      ark_mem->ark_etacf            = ETACF;
      break;
    }
  }
  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetErrHandlerFn:

 Specifies the error handler function
---------------------------------------------------------------*/
int ARKodeSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun, 
			  void *eh_data)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetErrHandlerFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* set user-provided values, or defaults, depending on argument */
  if (ehfun == NULL) {
    ark_mem->ark_ehfun = ARKErrHandler;
    ark_mem->ark_eh_data = ark_mem;
  } else {
    ark_mem->ark_ehfun   = ehfun;
    ark_mem->ark_eh_data = eh_data;
  }

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
 ARKodeSetDiagnostics:

 Specifies to enable solver diagnostics, and specifies the FILE
 pointer for output (diagfp==NULL disables output)
---------------------------------------------------------------*/
int ARKodeSetDiagnostics(void *arkode_mem, FILE *diagfp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetDiagnostics", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  ark_mem->ark_diagfp = diagfp;
  if (diagfp != NULL) {
    ark_mem->ark_report = TRUE;
  } else {
    ark_mem->ark_report = FALSE;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetOrder:

 Specifies the method order

 ** Note in documentation that this should not be called along 
    with ARKodeSetERKTable, ARKodeSetIRKTable, ARKodeSetARKTable, 
    ARKodeSetERKTableNum, ARKodeSetIRKTableNum or 
    ARKodeSetARKTableNum.  This routine is used to specify a 
    desired method order using default Butcher tables, whereas 
    any user-supplied table will have their own order associated 
    with them.
---------------------------------------------------------------*/
int ARKodeSetOrder(void *arkode_mem, int ord)
{
  int i, j;

  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxOrd", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* set user-provided value, or default, depending on argument */
  if (ord <= 0) {
    ark_mem->ark_q = Q_DEFAULT;
  } else {
    ark_mem->ark_q = ord;
  }

  /* clear Butcher tables, since user is requesting a change in method
     or a reset to defaults.  Tables will be set in ARKInitialSetup. */
  ark_mem->ark_stages = 0;
  ark_mem->ark_istage = 0;
  ark_mem->ark_p = 0;
  for (i=0; i<ARK_S_MAX; i++) {
    for (j=0; j<ARK_S_MAX; j++) {
      ARK_A(ark_mem->ark_Ae,i,j) = ZERO;
      ARK_A(ark_mem->ark_Ai,i,j) = ZERO;
    }
    ark_mem->ark_c[i]   = ZERO;
    ark_mem->ark_b[i]   = ZERO;
    ark_mem->ark_b2[i]  = ZERO;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetDenseOrder:

 Specifies the polynomial order for dense output.  Allowed values
 range from 0 to min(q,5), where q is the order of the time 
 integration method.  Illegal values imply to use the default.
---------------------------------------------------------------*/
int ARKodeSetDenseOrder(void *arkode_mem, int dord)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetMaxOrd", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check input */
  if (dord > 5) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetMaxOrd", "Dense output order must be <= 5");
    return(ARK_ILL_INPUT);
  }
  /* NOTE: we check that dord < q internally, to allow for subsequent 
     changes via ARKodeSetOrder */

  /* set user-provided value, or default, depending on argument */
  /* if ((dord < 0) || (dord > 5)) { */
  if ((dord < 0) || (dord > 3)) {
    ark_mem->ark_dense_q = 3;
  } else {
    ark_mem->ark_dense_q = dord;
  }

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
		    "ARKodeSetLinear", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->ark_linear = TRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetNonlinear:

 Specifies that the implicit portion of the problem is nonlinear.
 Used to undo a previous call to ARKodeSetLinear.
---------------------------------------------------------------*/
int ARKodeSetNonlinear(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetNonlinear", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->ark_linear = FALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetExplicit:

 Specifies that the implicit portion of the problem is disabled, 
 and to use an explicit RK method.
---------------------------------------------------------------*/
int ARKodeSetExplicit(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetExplicit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* ensure that fe is defined */
  if (ark_mem->ark_fe == NULL) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetExplicit", MSGARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  ark_mem->ark_explicit = TRUE;
  ark_mem->ark_implicit = FALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetImplicit:

 Specifies that the explicit portion of the problem is disabled, 
 and to use an implicit RK method.
---------------------------------------------------------------*/
int ARKodeSetImplicit(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetImplicit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* ensure that fi is defined */
  if (ark_mem->ark_fi == NULL) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetImplicit", MSGARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  ark_mem->ark_implicit = TRUE;
  ark_mem->ark_explicit = FALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetImEx:

 Specifies that the specifies that problem has both implicit and
 explicit parts, and to use an ARK method (this is the default).
---------------------------------------------------------------*/
int ARKodeSetImEx(void *arkode_mem)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetImEx", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* ensure that fe and fi are defined */
  if (ark_mem->ark_fe == NULL) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetImEx", MSGARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }
  if (ark_mem->ark_fi == NULL) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetImEx", MSGARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  /* set the relevant parameters */
  ark_mem->ark_explicit = FALSE;
  ark_mem->ark_implicit = FALSE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetERKTable:

 Specifies to use a customized Butcher table for the explicit 
 portion of the system (automatically calls ARKodeSetExplicit).
---------------------------------------------------------------*/
int ARKodeSetERKTable(void *arkode_mem, int s, int q, int p,
		      realtype *c, realtype *A, realtype *b, 
		      realtype *bembed)
{
  int i, j;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetERKTable", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for legal inputs */
  if (s > ARK_S_MAX) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetERKTable", "s exceeds ARK_S_MAX");
    return(ARK_ILL_INPUT);
  }
  if ((c == NULL) || (A == NULL) || (b == NULL) || (bembed == NULL)) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetERKTable", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters */
  ark_mem->ark_stages = 0;
  ark_mem->ark_q = 0;
  ark_mem->ark_p = 0;
  for (i=0; i<ARK_S_MAX; i++) {
    ark_mem->ark_c[i]  = 0.0;
    ark_mem->ark_b[i]  = 0.0;
    ark_mem->ark_b2[i] = 0.0;
    for (j=0; j<ARK_S_MAX; j++) 
      ARK_A(ark_mem->ark_Ae,i,j) = 0.0;
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
      ARK_A(ark_mem->ark_Ae,i,j) = A[i*s + j];
    }
  }
  
  /* set method as purely explicit */
  if (ARKodeSetExplicit(arkode_mem) != ARK_SUCCESS) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetERKTable", MSGARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetIRKTable:

 Specifies to use a customized Butcher table for the implicit 
 portion of the system (automatically calls ARKodeSetImplicit).
---------------------------------------------------------------*/
int ARKodeSetIRKTable(void *arkode_mem, int s, int q, int p,
		      realtype *c, realtype *A, realtype *b, 
		      realtype *bembed)
{
  int i, j;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetIRKTable", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for legal inputs */
  if (s > ARK_S_MAX) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetIRKTable", "s exceeds ARK_S_MAX");
    return(ARK_ILL_INPUT);
  }
  if ((c == NULL) || (A == NULL) || (b == NULL) || (bembed == NULL)) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetIRKTable", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters */
  ark_mem->ark_stages = 0;
  ark_mem->ark_q = 0;
  ark_mem->ark_p = 0;
  for (i=0; i<ARK_S_MAX; i++) {
    ark_mem->ark_c[i]  = 0.0;
    ark_mem->ark_b[i]  = 0.0;
    ark_mem->ark_b2[i] = 0.0;
    for (j=0; j<ARK_S_MAX; j++) 
      ARK_A(ark_mem->ark_Ai,i,j) = 0.0;
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
      ARK_A(ark_mem->ark_Ai,i,j) = A[i*s + j];
    }
  }

  /* set method as purely implicit */
  if (ARKodeSetImplicit(arkode_mem) != ARK_SUCCESS) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetIRKTable", MSGARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetARKTables:

 Specifies to use customized Butcher tables for the ImEx system
 (automatically calls ARKodeSetImEx).
---------------------------------------------------------------*/
int ARKodeSetARKTables(void *arkode_mem, int s, int q, int p,
		       realtype *c, realtype *Ai, realtype *Ae, 
		       realtype *b, realtype *bembed)
{
  int i, j;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetARKTables", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for legal inputs */
  if (s > ARK_S_MAX) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetARKTables", "s exceeds ARK_S_MAX");
    return(ARK_ILL_INPUT);
  }
  if ((c == NULL) || (Ai == NULL) || (Ae == NULL) || 
      (b == NULL) || (bembed == NULL)) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetARKTables", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }

  /* clear any existing parameters */
  ark_mem->ark_stages = 0;
  ark_mem->ark_q = 0;
  ark_mem->ark_p = 0;
  for (i=0; i<ARK_S_MAX; i++) {
    ark_mem->ark_c[i]  = 0.0;
    ark_mem->ark_b[i]  = 0.0;
    ark_mem->ark_b2[i] = 0.0;
    for (j=0; j<ARK_S_MAX; j++) {
      ARK_A(ark_mem->ark_Ai,i,j) = 0.0;
      ARK_A(ark_mem->ark_Ae,i,j) = 0.0;
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
      ARK_A(ark_mem->ark_Ai,i,j) = Ai[i*s + j];
      ARK_A(ark_mem->ark_Ae,i,j) = Ae[i*s + j];
    }
  }

  /* set method as ImEx */
  if (ARKodeSetImEx(arkode_mem) != ARK_SUCCESS) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetARKTables", MSGARK_MISSING_F);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetERKTableNum:

 Specifies to use a pre-existing Butcher table for the explicit 
 portion of the problem, based on the integer flag held in 
 ARKodeLoadButcherTable() within the file arkode_butcher.c 
 (automatically calls ARKodeSetExplicit).
---------------------------------------------------------------*/
int ARKodeSetERKTableNum(void *arkode_mem, int itable)
{
  int iflag;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetERKTableNum", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* fill in table based on argument */
  iflag = ARKodeLoadButcherTable(itable, &ark_mem->ark_stages, 
				 &ark_mem->ark_q, 
				 &ark_mem->ark_p, 
				 ark_mem->ark_Ae, 
				 ark_mem->ark_b, 
				 ark_mem->ark_c, 
				 ark_mem->ark_b2);
  /* check that requested table is legal */
  if (iflag != ARK_SUCCESS) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetERKTableNum", 
		    "Illegal ERK table number");
    return(ARK_ILL_INPUT);
  }

  /* set method as purely explicit */
  if (ARKodeSetExplicit(arkode_mem) != ARK_SUCCESS) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetERKTableNum", MSGARK_MISSING_FE);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetIRKTableNum:

 Specifies to use a pre-existing Butcher table for the implicit 
 portion of the problem, based on the integer flag held in 
 ARKodeLoadButcherTable() within the file arkode_butcher.c
 (automatically calls ARKodeSetImplicit).
---------------------------------------------------------------*/
int ARKodeSetIRKTableNum(void *arkode_mem, int itable)
{
  int iflag;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetImplicit", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* fill in table based on argument */
  iflag = ARKodeLoadButcherTable(itable, &ark_mem->ark_stages, 
				 &ark_mem->ark_q, 
				 &ark_mem->ark_p, 
				 ark_mem->ark_Ai, 
				 ark_mem->ark_b, 
				 ark_mem->ark_c, 
				 ark_mem->ark_b2);
  /* check that requested table is legal */
  if (iflag != ARK_SUCCESS) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetIRKTableNum", 
		    "Illegal IRK table number");
    return(ARK_ILL_INPUT);
  }

  /* set method as purely implicit */
  if (ARKodeSetImplicit(arkode_mem) != ARK_SUCCESS) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetIRKTableNum", MSGARK_MISSING_FI);
    return(ARK_ILL_INPUT);
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetARKTableNum:

 Specifies to use pre-existing Butcher tables for the ImEx system,
 based on the integer flags held in ARKodeLoadButcherTable() 
 within the file arkode_butcher.c (automatically calls ARKodeSetImEx).
---------------------------------------------------------------*/
int ARKodeSetARKTableNum(void *arkode_mem, int itable, int etable)
{
  int iflag, eflag;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetARKTableNum", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* ensure that tables match */
  iflag = 1;
  if ((etable == 3)  && (itable == 16))  iflag = 0;
  if ((etable == 6)  && (itable == 22))  iflag = 0;
  if ((etable == 11) && (itable == 26))  iflag = 0;
  if (iflag) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetARKTableNum", 
		    "Incompatible Butcher tables for ARK method");
    return(ARK_ILL_INPUT);
  }

  /* fill in tables based on arguments */
  iflag = ARKodeLoadButcherTable(itable, &ark_mem->ark_stages, 
				 &ark_mem->ark_q, 
				 &ark_mem->ark_p, 
				 ark_mem->ark_Ai, 
				 ark_mem->ark_b, 
				 ark_mem->ark_c, 
				 ark_mem->ark_b2);
  eflag = ARKodeLoadButcherTable(etable, &ark_mem->ark_stages, 
				 &ark_mem->ark_q, 
				 &ark_mem->ark_p, 
				 ark_mem->ark_Ae, 
				 ark_mem->ark_b, 
				 ark_mem->ark_c, 
				 ark_mem->ark_b2);

  /* check that requested tables are legal */
  if (iflag != ARK_SUCCESS) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetARKTableNum", 
		    "Illegal IRK table number");
    return(ARK_ILL_INPUT);
  }
  if (eflag != ARK_SUCCESS) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetARKTableNum", 
		    "Illegal ERK table number");
    return(ARK_ILL_INPUT);
  }

  /* set method as ImEx */
  if (ARKodeSetImEx(arkode_mem) != ARK_SUCCESS) {
    ARKProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetARKTableNum", MSGARK_MISSING_F);
    return(ARK_ILL_INPUT);
  }

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

  /* Passing mxhnil=0 sets the default, otherwise use input. */
  if (mxhnil == 0) {
    ark_mem->ark_mxhnil = 10;
  } else {
    ark_mem->ark_mxhnil = mxhnil;
  }

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

  /* Passing hin<=0 sets the default, otherwise use input. */
  if (hin <= ZERO) {
    ark_mem->ark_hin = ZERO;
  } else {
    ark_mem->ark_hin = hin;
  }

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
    ark_mem->ark_hmax_inv = ZERO;
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
     (i.e. if it was not already passed).
     If ARKodeSetStopTime is called before the first call to ARKode,
     tstop will be checked in ARKode. */
  if (ark_mem->ark_nst > 0) {
    if ( (tstop - ark_mem->ark_tn) * ark_mem->ark_h < ZERO ) {
      ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		      "ARKodeSetStopTime", MSGARK_BAD_TSTOP, ark_mem->ark_tn);
      return(ARK_ILL_INPUT);
    }
  }

  ark_mem->ark_tstop    = tstop;
  ark_mem->ark_tstopset = TRUE;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetAdaptivityMethod:

 Specifies the built-in time step adaptivity algorithm (and 
 associated parameters) to use.  Any zero-valued parameter will 
 imply a reset to the default value.  Any negative parameter will 
 be left unchanged from the previous value.
---------------------------------------------------------------*/
int ARKodeSetAdaptivityMethod(void *arkode_mem, int imethod, 
			      realtype *adapt_params)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetStopTime", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* check for allowable parameters */
  if (imethod > 5) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetAdaptivityMethod", "Illegal imethod");
    return(ARK_ILL_INPUT);
  }

  if (adapt_params == NULL) {
    ARKProcessError(ark_mem, ARK_ILL_INPUT, "ARKODE", 
		    "ARKodeSetAdaptivityMethod", "Illegal adapt_params");
    return(ARK_ILL_INPUT);
  }

  /* set method if change requested */
  if (imethod >= 0)
    ark_mem->ark_hadapt_imethod = imethod;

  /* set positive-valued paramters into ark_mem, 
     otherwise set defaults for the chosen method */
  ark_mem->ark_hadapt_cfl = adapt_params[0];
  if (adapt_params[0] == ZERO) 
    ark_mem->ark_hadapt_cfl = CFLFAC;
    
  ark_mem->ark_hadapt_safety = adapt_params[1];
  if (adapt_params[1] == ZERO) 
    ark_mem->ark_hadapt_safety = SAFETY;

  ark_mem->ark_hadapt_bias = adapt_params[2];
  if (adapt_params[2] == ZERO) 
    ark_mem->ark_hadapt_bias = BIAS;

  ark_mem->ark_hadapt_growth = adapt_params[3];
  if (adapt_params[3] == ZERO) 
    ark_mem->ark_hadapt_growth = GROWTH;

  ark_mem->ark_hadapt_lbound = adapt_params[4];
  if (adapt_params[4] == ZERO) 
    ark_mem->ark_hadapt_lbound = HFIXED_LB;

  ark_mem->ark_hadapt_ubound = adapt_params[5];
  if (adapt_params[5] == ZERO) 
    ark_mem->ark_hadapt_ubound = HFIXED_UB;

  ark_mem->ark_hadapt_k1 = adapt_params[6];
  if (adapt_params[6] == ZERO) 
    switch (ark_mem->ark_hadapt_imethod) {
    case (0):
      ark_mem->ark_hadapt_k1 = AD0_K1; break;
    case (1):
      ark_mem->ark_hadapt_k1 = AD1_K1; break;
    case (2):
      ark_mem->ark_hadapt_k1 = AD2_K1; break;
    case (3):
      ark_mem->ark_hadapt_k1 = AD3_K1; break;
    case (4):
      ark_mem->ark_hadapt_k1 = AD4_K1; break;
    case (5):
      ark_mem->ark_hadapt_k1 = AD5_K1; break;
    }

  ark_mem->ark_hadapt_k2 = adapt_params[7];
  if (adapt_params[7] == ZERO) 
    switch (ark_mem->ark_hadapt_imethod) {
    case (0):
      ark_mem->ark_hadapt_k2 = AD0_K2; break;
    case (1):
      ark_mem->ark_hadapt_k2 = AD1_K2; break;
    case (3):
      ark_mem->ark_hadapt_k2 = AD3_K2; break;
    case (4):
      ark_mem->ark_hadapt_k2 = AD4_K2; break;
    case (5):
      ark_mem->ark_hadapt_k2 = AD5_K2; break;
    }

  ark_mem->ark_hadapt_k3 = adapt_params[8];
  if (adapt_params[8] == ZERO) 
    switch (ark_mem->ark_hadapt_imethod) {
    case (0):
      ark_mem->ark_hadapt_k3 = AD0_K3; break;
    case (5):
      ark_mem->ark_hadapt_k3 = AD5_K3; break;
    }

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

  /* NULL argument sets default, otherwise set inputs */
  if (hfun == NULL) {
    ark_mem->ark_hadapt         = NULL;
    ark_mem->ark_hadapt_data    = NULL;
    ark_mem->ark_hadapt_imethod = 0;
  } else {
    ark_mem->ark_hadapt         = hfun;
    ark_mem->ark_hadapt_data    = h_data;
    ark_mem->ark_hadapt_imethod = -1;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetAdaptivityConstants:

 Specifies the user-provided time step adaptivity constants
 etamx1, etamxf, etacf and small_nef.  Illegal inputs imply 
 a reset to the default value.
---------------------------------------------------------------*/
int ARKodeSetAdaptivityConstants(void *arkode_mem, realtype etamx1,
				 realtype etamxf, realtype etacf,
				 int small_nef)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetAdaptivityFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* 0 or 0.0 argument sets default, otherwise set inputs */
  if (etamx1 <= 0.0) {
    ark_mem->ark_etamx1 = ETAMX1;
  } else {
    ark_mem->ark_etamx1 = etamx1;
  }
  if (etamxf <= 0.0) {
    ark_mem->ark_etamxf = ETAMXF;
  } else {
    ark_mem->ark_etamxf = etamxf;
  }
  if (etacf <= 0.0) {
    ark_mem->ark_etacf = ETACF;
  } else {
    ark_mem->ark_etacf = etacf;
  }
  if (small_nef <= 0) {
    ark_mem->ark_small_nef = SMALL_NEF;
  } else {
    ark_mem->ark_small_nef = small_nef;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetNewtonConstants:

 Specifies the user-provided nonlinear convergence constants
 crdown and rdiv.  Zero-valued inputs imply a reset to the 
 default value.
---------------------------------------------------------------*/
int ARKodeSetNewtonConstants(void *arkode_mem, realtype crdown,
			     realtype rdiv)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetAdaptivityFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* 0.0 argument sets default, otherwise set inputs */
  if (crdown == 0.0) {
    ark_mem->ark_crdown = CRDOWN;
  } else {
    ark_mem->ark_crdown = crdown;
  }
  if (rdiv == 0.0) {
    ark_mem->ark_rdiv = RDIV;
  } else {
    ark_mem->ark_rdiv = rdiv;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetLSetupConstants:

 Specifies the user-provided linear setup decision constants
 dgmax and msbp.  Zero-valued inputs imply a reset to the 
 default value.
---------------------------------------------------------------*/
int ARKodeSetLSetupConstants(void *arkode_mem, realtype dgmax,
			     int msbp)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetAdaptivityFn", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* 0 or 0.0 argument sets default, otherwise set inputs */
  if (dgmax == 0.0) {
    ark_mem->ark_dgmax = DGMAX;
  } else {
    ark_mem->ark_dgmax = dgmax;
  }
  if (msbp == 0) {
    ark_mem->ark_msbp = MSBP;
  } else {
    ark_mem->ark_msbp = msbp;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetPredictorMethod:

 Specifies the method to use for predicting implicit solutions.  
 Non-default choices are {1,2,3}, all others will use default 
 (trivial) predictor.
---------------------------------------------------------------*/
int ARKodeSetPredictorMethod(void *arkode_mem, int pred_method)
{
  ARKodeMem ark_mem;

  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetStopTime", MSGARK_NO_MEM);
    return (ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;
  ark_mem->ark_predictor = pred_method;

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

  /* NULL argument sets default, otherwise set inputs */
  if (EStab == NULL) {
    ark_mem->ark_expstab    = ARKExpStab;
    ark_mem->ark_estab_data = ark_mem;
  } else {
    ark_mem->ark_expstab    = EStab;
    ark_mem->ark_estab_data = estab_data;
  }

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

  /* argument <= 0 sets default, otherwise set input */
  if (maxnef <= 0) {
    ark_mem->ark_maxnef = MAXNEF;
  } else {
    ark_mem->ark_maxnef = maxnef;
  }

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

  /* argument <= 0 sets default, otherwise set input */
  if (maxncf <= 0) {
    ark_mem->ark_maxncf = MAXNCF;
  } else {
    ark_mem->ark_maxncf = maxncf;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetMaxNonlinIters:

 Specifies the maximum number of nonlinear iterations during
 one solve.  Zero-valued input implies a reset to the 
 default value.
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

  /* argument <= 0 sets default, otherwise set input */
  if (maxcor <= 0) {
    ark_mem->ark_maxcor = MAXCOR;
  } else {
    ark_mem->ark_maxcor = maxcor;
  }

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeSetNonlinConvCoef:

 Specifies the coefficient in the nonlinear solver convergence
 test.  Zero-valued input implies a reset to the default value.
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

  /* argument <= 0 sets default, otherwise set input */
  if (nlscoef <= ZERO) {
    ark_mem->ark_nlscoef = NLSCOEF;
  } else {
    ark_mem->ark_nlscoef = nlscoef;
  }

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
 ARKodeGetNumExpSteps:

 Returns the current number of stability-limited steps
---------------------------------------------------------------*/
int ARKodeGetNumExpSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumSteps", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->ark_nst_exp;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumAccSteps:

 Returns the current number of accuracy-limited steps
---------------------------------------------------------------*/
int ARKodeGetNumAccSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumSteps", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->ark_nst_acc;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumConvSteps:

 Returns the current number of convergence-limited steps
---------------------------------------------------------------*/
int ARKodeGetNumConvSteps(void *arkode_mem, long int *nsteps)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumSteps", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps = ark_mem->ark_nst_con;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumRhsEvals:

 Returns the current number of calls to fe and fi
---------------------------------------------------------------*/
int ARKodeGetNumRhsEvals(void *arkode_mem, long int *fe_evals,
			 long int *fi_evals)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetNumRhsEvals", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *fe_evals = ark_mem->ark_nfe;
  *fi_evals = ark_mem->ark_nfi;

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeGetNumLinSolvSetups:

 Returns the current number of calls to the lsetup routine
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

  *hlast = ark_mem->ark_hold;

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
 ARKodeGetCurrentButcherTables:

 Returns the explicit and implicit Butcher tables currently in 
 use.  The tables should be declared statically as 
 Ae[ARK_S_MAX*ARK_S_MAX] and Ai[ARK_S_MAX*ARK_S_MAX], and the 
 arrays c, b and b2 should all have length ARK_S_MAX.
---------------------------------------------------------------*/
int ARKodeGetCurrentButcherTables(void *arkode_mem, 
				  int *s, int *q, int *p,
				  realtype *Ai, realtype *Ae, 
				  realtype *c, realtype *b,
				  realtype *b2)
{
  int i,j;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetCurrentTime", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *s = ark_mem->ark_stages;
  *q = ark_mem->ark_q;
  *p = ark_mem->ark_p;
  for (i=0; i<ARK_S_MAX; i++) {
    for (j=0; j<ARK_S_MAX; j++) {
      ARK_A(Ae,i,j) = ARK_A(ark_mem->ark_Ae,i,j);
      ARK_A(Ai,i,j) = ARK_A(ark_mem->ark_Ai,i,j);
    }
    c[i]  = ark_mem->ark_c[i];
    b[i]  = ark_mem->ark_b[i];
    b2[i] = ark_mem->ark_b2[i];
  }

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
 ARKodeGetEstLocalErrors:  (to be updated)

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
			     long int *convsteps, long int *fe_evals, 
			     long int *fi_evals, long int *nlinsetups, 
			     long int *netfails, realtype *hinused, 
			     realtype *hlast, realtype *hcur, 
			     realtype *tcur)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeGetIntegratorStats", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  *nsteps     = ark_mem->ark_nst;
  *expsteps   = ark_mem->ark_nst_exp;
  *accsteps   = ark_mem->ark_nst_acc;
  *convsteps  = ark_mem->ark_nst_con;
  *fe_evals   = ark_mem->ark_nfe;
  *fi_evals   = ark_mem->ark_nfi;
  *nlinsetups = ark_mem->ark_nsetups;
  *netfails   = ark_mem->ark_netf;
  *hinused    = ark_mem->ark_h0u;
  *hlast      = ark_mem->ark_hold;
  *hcur       = ark_mem->ark_next_h;
  *tcur       = ark_mem->ark_tn;

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
  int i, nrt;
  ARKodeMem ark_mem;
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

 Returns the current number of nonlinear solver iterations 
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

 Returns the current number of nonlinear solver convergence fails
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



/*===============================================================
  ARKODE parameter output
===============================================================*/

/*---------------------------------------------------------------
 ARKodeWriteParameters:

 Outputs all solver parameters to the provided file pointer.
---------------------------------------------------------------*/
int ARKodeWriteParameters(void *arkode_mem, FILE *fp)
{
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetDefaults", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* print integrator parameters to file */
  fprintf(fp, "ARKode solver parameters:\n");
  fprintf(fp, "  Method order %i\n",ark_mem->ark_q);
  fprintf(fp, "  Dense output order %i\n",ark_mem->ark_dense_q);
  if (ark_mem->ark_linear)  fprintf(fp, "  Linear problem\n");
  if (ark_mem->ark_explicit) {
    fprintf(fp, "  Explicit integrator\n");
  } else if (ark_mem->ark_implicit) {
    fprintf(fp, "  Implicit integrator\n");
  } else {
    fprintf(fp, "  ImEx integrator\n");
  }
  if (ark_mem->ark_hadapt == NULL) {
    fprintf(fp, "  Time step adaptivity method %i\n", ark_mem->ark_hadapt_imethod);
    fprintf(fp, "     Safety factor = %g\n", ark_mem->ark_hadapt_safety);
    fprintf(fp, "     Bias factor = %g\n", ark_mem->ark_hadapt_bias);
    fprintf(fp, "     Growth factor = %g\n", ark_mem->ark_hadapt_growth);
    fprintf(fp, "     Step growth lower bound = %g\n", ark_mem->ark_hadapt_lbound);
    fprintf(fp, "     Step growth upper bound = %g\n", ark_mem->ark_hadapt_ubound);
    fprintf(fp, "     k1 = %g\n", ark_mem->ark_hadapt_k1);
    fprintf(fp, "     k2 = %g\n", ark_mem->ark_hadapt_k2);
    fprintf(fp, "     k3 = %g\n", ark_mem->ark_hadapt_k3);
  } else {
    fprintf(fp, "  User provided time step adaptivity function\n");
  }
  if (ark_mem->ark_itol == ARK_WF) {
    fprintf(fp, "  User provided error weight function\n");
  } else {
    fprintf(fp, "  Solver relative tolerance = %g\n", ark_mem->ark_reltol);
    if (ark_mem->ark_itol == ARK_SS) {
      fprintf(fp, "  Solver absolute tolerance = %g\n", ark_mem->ark_Sabstol);
    } else {
      fprintf(fp, "  Vector-valued solver absolute tolerance\n");
    }
  }
  if (ark_mem->ark_hin != ZERO)  
    fprintf(fp, "  Initial step size = %g\n",ark_mem->ark_hin);
  if (ark_mem->ark_hmin != ZERO)  
    fprintf(fp, "  Minimum step size = %g\n",ark_mem->ark_hmin);
  if (ark_mem->ark_hmax_inv != ZERO)  
    fprintf(fp, "  Maximum step size = %g\n",ONE/ark_mem->ark_hmax_inv);
  fprintf(fp, "  Maximum number of error test failures = %i\n",ark_mem->ark_maxnef);
  fprintf(fp, "  Maximum number of convergence test failures = %i\n",ark_mem->ark_maxncf);
  fprintf(fp, "  Maximum step increase (first step) = %g\n",ark_mem->ark_etamx1);
  fprintf(fp, "  Step reduction factor on multiple error fails = %g\n",ark_mem->ark_etamxf);
  fprintf(fp, "  Minimum error fails before above factor is used = %i\n",ark_mem->ark_small_nef);
  fprintf(fp, "  Step reduction factor on Newton convergence failure = %g\n",ark_mem->ark_etacf);

  if (!ark_mem->ark_implicit) {
    if (ark_mem->ark_expstab == ARKExpStab) {
      fprintf(fp, "  Default explicit stability function\n");
    } else {
      fprintf(fp, "  User provided explicit stability function\n");
    }
    fprintf(fp, "  Explicit safety factor = %g\n",ark_mem->ark_hadapt_cfl);
  }
  if (!ark_mem->ark_explicit) {
    fprintf(fp, "  Newton predictor method = %i\n",ark_mem->ark_predictor);
    fprintf(fp, "  Newton solver tolerance coefficient = %g\n",ark_mem->ark_nlscoef);
    fprintf(fp, "  Maximum number of Newton corrections = %i\n",ark_mem->ark_maxcor);
    fprintf(fp, "  Nonlinear convergence rate constant = %g\n",ark_mem->ark_crdown);
    fprintf(fp, "  Nonlinear divergence tolerance = %g\n",ark_mem->ark_rdiv);
    fprintf(fp, "  Gamma factor LSetup tolerance = %g\n",ark_mem->ark_dgmax);
    fprintf(fp, "  Number of steps between LSetup calls = %i\n",ark_mem->ark_msbp);
  }
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}


/*---------------------------------------------------------------
 ARKodeWriteButcher:

 Outputs Butcher tables to the provided file pointer.
---------------------------------------------------------------*/
int ARKodeWriteButcher(void *arkode_mem, FILE *fp)
{
  int i, j;
  ARKodeMem ark_mem;
  if (arkode_mem==NULL) {
    ARKProcessError(NULL, ARK_MEM_NULL, "ARKODE", 
		    "ARKodeSetDefaults", MSGARK_NO_MEM);
    return(ARK_MEM_NULL);
  }
  ark_mem = (ARKodeMem) arkode_mem;

  /* print Butcher tables to file */
  fprintf(fp, "\nARKode Butcher tables (stages = %i):\n", ark_mem->ark_stages);
  if (!ark_mem->ark_implicit) {
    fprintf(fp, "  Explicit Butcher table:\n");
    for (i=0; i<ark_mem->ark_stages; i++) {
      fprintf(fp, "     %.5f",ark_mem->ark_c[i]);
      for (j=0; j<ark_mem->ark_stages; j++) 
	fprintf(fp, " %.5f",ARK_A(ark_mem->ark_Ae,i,j));
      fprintf(fp,"\n");
    }
    fprintf(fp, "            ");
    for (j=0; j<ark_mem->ark_stages; j++) 
      fprintf(fp, " %.5f",ark_mem->ark_b[j]);
    fprintf(fp,"\n");
    fprintf(fp, "            ");
    for (j=0; j<ark_mem->ark_stages; j++) 
      fprintf(fp, " %.5f",ark_mem->ark_b2[j]);
    fprintf(fp,"\n");
  }
  if (!ark_mem->ark_explicit) {
    fprintf(fp, "  Implicit Butcher table:\n");
    for (i=0; i<ark_mem->ark_stages; i++) {
      fprintf(fp, "     %.5f",ark_mem->ark_c[i]);
      for (j=0; j<ark_mem->ark_stages; j++) 
	fprintf(fp, " %.5f",ARK_A(ark_mem->ark_Ai,i,j));
      fprintf(fp,"\n");
    }
    fprintf(fp, "            ");
    for (j=0; j<ark_mem->ark_stages; j++) 
      fprintf(fp, " %.5f",ark_mem->ark_b[j]);
    fprintf(fp,"\n");
    fprintf(fp, "            ");
    for (j=0; j<ark_mem->ark_stages; j++) 
      fprintf(fp, " %.5f",ark_mem->ark_b2[j]);
    fprintf(fp,"\n");
  }
  fprintf(fp, "\n");

  return(ARK_SUCCESS);
}




/*---------------------------------------------------------------
      EOF
---------------------------------------------------------------*/
