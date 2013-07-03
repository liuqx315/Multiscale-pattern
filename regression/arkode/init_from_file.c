/*---------------------------------------------------------------
 $Revision: 1.0 $
 $Date:  $
-----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Utility routine to read input parameters from a specified
 file, and call associated "set" routines to specify options to
 to ARKode solver.
---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <arkode/arkode.h>
#include <sundials/sundials_types.h>

#define MAX_LINE_LENGTH 512


/* Main Program */
int init_from_file(void *ark_mem, char *fname, ARKRhsFn f, 
		   ARKRhsFn fe, ARKRhsFn fi, realtype T0, 
		   N_Vector y0, int *ImEx, int *dorder, 
		   realtype *RTol, realtype *ATol) {

  /* declare available solver parameters (with default values) */
  int order = 0; 
  int imex = 0;  
  int adapt_method = 0; 
  int small_nef = 0; 
  int msbp = 0; 
  int maxcor = 0; 
  int predictor = 0;
  int maxnef = 0;
  int maxncf = 0;
  int mxhnil = 0;
  int mxsteps = 0;
  int dense_order = -1;
  int btable = -1;
  double cflfac = 0.0;
  double safety = 0.0;
  double bias = 0.0;
  double growth = 0.0;
  double hfixed_lb = 0.0;
  double hfixed_ub = 0.0;
  double k1 = 0.0;
  double k2 = 0.0;
  double k3 = 0.0;
  double etamx1 = 0.0;
  double etamxf = 0.0;
  double etacf = 0.0;
  double crdown = 0.0;
  double rdiv = 0.0;
  double dgmax = 0.0;
  double nlscoef = 0.0;
  double h0 = 0.0;
  double hmin = 0.0;
  double hmax = 0.0;
  double rtol = 0.0;
  double atol = 0.0;

  /* open parameter file */
  FILE *fptr = NULL;
  fptr = fopen(fname,"r");
  if (fptr == NULL) {
    fprintf(stderr, "set_from_file error: cannot open parameter file %s\n", fname);
    return 1;
  }

  /* read solver parameters from file */
  int ret;
  char line[MAX_LINE_LENGTH];
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    /* initialize return flag for line */
    ret = 0;

    /* read parameter */
    ret += sscanf(line,"order = %i", &order);
    ret += sscanf(line,"dense_order = %i", &dense_order);
    ret += sscanf(line,"imex = %i", &imex);
    ret += sscanf(line,"btable = %i",  &btable);
    ret += sscanf(line,"adapt_method = %i", &adapt_method);
    ret += sscanf(line,"maxnef = %i", &maxnef);
    ret += sscanf(line,"maxncf = %i", &maxncf);
    ret += sscanf(line,"mxhnil = %i", &mxhnil);
    ret += sscanf(line,"mxsteps = %i", &mxsteps);
    ret += sscanf(line,"cflfac = %lf", &cflfac);
    ret += sscanf(line,"safety = %lf", &safety);
    ret += sscanf(line,"bias = %lf", &bias);
    ret += sscanf(line,"growth = %lf", &growth);
    ret += sscanf(line,"hfixed_lb = %lf", &hfixed_lb);
    ret += sscanf(line,"hfixed_ub = %lf", &hfixed_ub);
    ret += sscanf(line,"k1 = %lf", &k1);
    ret += sscanf(line,"k2 = %lf", &k2);
    ret += sscanf(line,"k3 = %lf", &k3);
    ret += sscanf(line,"etamx1 = %lf", &etamx1);
    ret += sscanf(line,"etamxf = %lf", &etamxf);
    ret += sscanf(line,"etacf = %lf", &etacf);
    ret += sscanf(line,"small_nef = %i", &small_nef);
    ret += sscanf(line,"crdown = %lf", &crdown);
    ret += sscanf(line,"rdiv = %lf", &rdiv);
    ret += sscanf(line,"dgmax = %lf", &dgmax);
    ret += sscanf(line,"predictor = %i", &predictor);
    ret += sscanf(line,"msbp = %i", &msbp);
    ret += sscanf(line,"maxcor = %i", &maxcor);
    ret += sscanf(line,"nlscoef = %lf", &nlscoef);
    ret += sscanf(line,"h0 = %lf", &h0);
    ret += sscanf(line,"hmin = %lf", &hmin);
    ret += sscanf(line,"hmax = %lf", &hmax);
    ret += sscanf(line,"rtol = %lf", &rtol);
    ret += sscanf(line,"atol = %lf", &atol);

    /* if unable to read the line (and it looks suspicious) issue a warning */
    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#')
      fprintf(stderr, "set_from_file Warning: parameter line was not interpreted:\n%s", line);

  }
  fclose(fptr);


  /*** check for allowable inputs ***/

  /* check that y0 is not NULL */
  if (y0 == NULL) {
    fprintf(stderr, "set_from_file error: cannot initialize problem with y0 == NULL!\n");
    return 1;
  }

  /* ensure that "imex" agrees with user-supplied rhs functions */
  if ((imex == 2) && (fe == NULL || fi == NULL)) {
    fprintf(stderr, "set_from_file error: imex problem but fe or fi is NULL!\n");
    return 1;
  }
  if ((imex == 0 || imex == 1) && (f == NULL)) {
    fprintf(stderr, "set_from_file error: implicit or explicit problem but f is NULL!\n");
    return 1;
  }



  /*** set outputs to be used by problem ***/
  *ImEx = imex;
  *dorder = dense_order;
  *RTol = rtol;
  *ATol = atol;



  /*** Call ARKode routines to initialize integrator and set options ***/

  /* initialize the integrator memory  */
  switch (imex) {
  case 0:         /* purely implicit */
    ret = ARKodeInit(ark_mem, NULL, f, T0, y0);  break;
  case 1:         /* purely explicit */
    ret = ARKodeInit(ark_mem, f, NULL, T0, y0);  break;
  default:        /* imex */
    ret = ARKodeInit(ark_mem, fe, fi, T0, y0);   break;
  }
  if (ret != 0) {
    fprintf(stderr, "set_from_file error in ARKodeInit = %i\n",ret);
    return 1;
  }

  /* set RK order, or specify individual Butcher table -- "order" overrides "btable" */
  if (order != 0) {     /*  */
    ret = ARKodeSetOrder(ark_mem, order);
    if (ret != 0) {
      fprintf(stderr,"set_from_file error in ARKodeSetOrder = %i\n",ret);
      return 1;
    }
  } else if (btable != -1) {
    if (imex == 1) {  
      ret = ARKodeSetERKTableNum(ark_mem, btable);
      if (ret != 0) {
	fprintf(stderr,"set_from_file error in ARKodeSetERKTableNum = %i\n",ret);
	return 1;
      }
    } else if (imex == 0) {  
      ret = ARKodeSetIRKTableNum(ark_mem, btable);
      if (ret != 0) {
	fprintf(stderr,"set_from_file error in ARKodeSetIRKTableNum = %i\n",ret);
	return 1;
      }
    } else if (imex == 2) {  
      int btable2;
      if (btable == 3)   btable2 = 16;
      if (btable == 6)   btable2 = 22;
      if (btable == 11)  btable2 = 26;
      ret = ARKodeSetARKTableNum(ark_mem, btable2, btable);
      if (ret != 0) {
	fprintf(stderr,"set_from_file error in ARKodeSetARKTableNum = %i\n",ret);
	return 1;
      }
    }
  }

  /* set dense output order */
  ret = ARKodeSetDenseOrder(ark_mem, dense_order);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetDenseOrder = %i\n",ret);
    return 1;
  }

  /* set cfl stability fraction */
  ret = ARKodeSetCFLFraction(ark_mem, cflfac);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetCFLFraction = %i\n",ret);
    return 1;
  }

  /* set safety factor */
  ret = ARKodeSetSafetyFactor(ark_mem, safety);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetSafetyFactor = %i\n",ret);
    return 1;
  }

  /* set error bias */
  ret = ARKodeSetErrorBias(ark_mem, bias);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetErrorBias = %i\n",ret);
    return 1;
  }

  /* set step growth factor */
  ret = ARKodeSetMaxGrowth(ark_mem, growth);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMaxGrowth = %i\n",ret);
    return 1;
  }

  /* set fixed step size bounds */
  ret = ARKodeSetFixedStepBounds(ark_mem, hfixed_lb, hfixed_ub);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetFixedStepBounds = %i\n",ret);
    return 1;
  }

  /* set time step adaptivity method */
  realtype adapt_params[] = {k1, k2, k3};
  int idefault = 1;
  if (fabs(k1)+fabs(k2)+fabs(k3) > 0.0)  idefault=0;
  ret = ARKodeSetAdaptivityMethod(ark_mem, adapt_method, idefault, 0, adapt_params);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetAdaptivityMethod = %i\n",ret);
    return 1;
  }

  /* set first step growth factor */
  ret = ARKodeSetMaxFirstGrowth(ark_mem, etamx1);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMaxFirstGrowth = %i\n",ret);
    return 1;
  }

  /* set error failure growth factor */
  ret = ARKodeSetMaxEFailGrowth(ark_mem, etamxf);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMaxEFailGrowth = %i\n",ret);
    return 1;
  }

  /* set number of fails before using above threshold */
  ret = ARKodeSetSmallNumEFails(ark_mem, small_nef);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetSmallNumEFails = %i\n",ret);
    return 1;
  }

  /* set convergence failure growth factor */
  ret = ARKodeSetMaxCFailGrowth(ark_mem, etacf);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMaxCFailGrowth = %i\n",ret);
    return 1;
  }

  /* set Newton method convergence rate constant */
  ret = ARKodeSetNewtonCRDown(ark_mem, crdown);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetNewtonCRDown = %i\n",ret);
    return 1;
  }

  /* set Newton method divergence constant */
  ret = ARKodeSetNewtonRDiv(ark_mem, rdiv);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetNewtonRDiv = %i\n",ret);
    return 1;
  }

  /* set linear solver setup constants */
  ret = ARKodeSetDeltaGammaMax(ark_mem, dgmax);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetDeltaGammaMax = %i\n",ret);
    return 1;
  }

  /* set linear solver setup constants */
  ret = ARKodeSetMaxStepsBetweenLSet(ark_mem, msbp);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMaxStepsBetweenLSet = %i\n",ret);
    return 1;
  }

  /* set predictor method */
  ret = ARKodeSetPredictorMethod(ark_mem, predictor);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetPredictorMethod = %i\n",ret);
    return 1;
  }

  /* set maximum Newton iterations */
  ret = ARKodeSetMaxNonlinIters(ark_mem, maxcor);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMaxNonlinIters = %i\n",ret);
    return 1;
  }

  /* set Newton solver tolerance coefficient */
  ret = ARKodeSetNonlinConvCoef(ark_mem, nlscoef);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMaxNonlinIters = %i\n",ret);
    return 1;
  }

  /* set initial time step size */
  ret = ARKodeSetInitStep(ark_mem, h0);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetInitStep = %i\n",ret);
    return 1;
  }

  /* set minimum time step size */
  ret = ARKodeSetMinStep(ark_mem, hmin);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMinStep = %i\n",ret);
    return 1;
  }

  /* set maximum time step size */
  ret = ARKodeSetMaxStep(ark_mem, hmax);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMaxStep = %i\n",ret);
    return 1;
  }

  /* set maximum allowed error test failures */
  ret = ARKodeSetMaxErrTestFails(ark_mem, maxnef);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMaxErrTestFails = %i\n",ret);
    return 1;
  }

  /* set maximum allowed convergence failures */
  ret = ARKodeSetMaxConvFails(ark_mem, maxncf);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMaxConvFails = %i\n",ret);
    return 1;
  }

  /* set maximum allowed hnil warnings */
  ret = ARKodeSetMaxHnilWarns(ark_mem, mxhnil);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMaxHnilWarns = %i\n",ret);
    return 1;
  }

  /* set maximum allowed steps */
  ret = ARKodeSetMaxNumSteps(ark_mem, mxsteps);
  if (ret != 0) {
    fprintf(stderr,"set_from_file error in ARKodeSetMaxNumSteps = %i\n",ret);
    return 1;
  }

  return 0;
}


/*---- end of file ----*/
