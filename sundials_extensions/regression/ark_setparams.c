/* -----------------------------------------------------------------
 * $Revision: $
 * $Date: $
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * Utility function to read solver parameters from the input file
 * 'solve_params.txt' and pass them to ARKode.
 * 
 * Returns 0 on success, 1 on failure.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <arkode/arkode.h>


int ark_SetParams(void *arkode_mem, int *idense)
{
  /* declare solver parameters */
  int flag, order, dense_order, btable, adapt_method, small_nef, 
    msbp, maxcor, predictor;
  realtype adapt_params[9], etamx1, etamxf, etacf, crdown, rdiv, 
    dgmax, nlscoef;

  /* open solver parameter file */
  FILE *FID;
  FID=fopen("solve_params.txt","r");
  fscanf(FID,"order = %i\n",  &order);
  fscanf(FID,"dense_order = %i\n", &dense_order);
  fscanf(FID,"btable = %i\n",  &btable);
  fscanf(FID,"adapt_method = %i\n", &adapt_method);
  fscanf(FID,"cflfac = %lf\n", &(adapt_params[0]));
  fscanf(FID,"safety = %lf\n", &(adapt_params[1]));
  fscanf(FID,"bias = %lf\n", &(adapt_params[2]));
  fscanf(FID,"growth = %lf\n", &(adapt_params[3]));
  fscanf(FID,"hfixed_lb = %lf\n", &(adapt_params[4]));
  fscanf(FID,"hfixed_ub = %lf\n", &(adapt_params[5]));
  fscanf(FID,"k1 = %lf\n", &(adapt_params[6]));
  fscanf(FID,"k2 = %lf\n", &(adapt_params[7]));
  fscanf(FID,"k3 = %lf\n", &(adapt_params[8]));
  fscanf(FID,"etamx1 = %lf\n", &etamx1);
  fscanf(FID,"etamxf = %lf\n", &etamxf);
  fscanf(FID,"etacf = %lf\n", &etacf);
  fscanf(FID,"small_nef = %i\n", &small_nef);
  fscanf(FID,"crdown = %lf\n", &crdown);
  fscanf(FID,"rdiv = %lf\n", &rdiv);
  fscanf(FID,"dgmax = %lf\n", &dgmax);
  fscanf(FID,"predictor = %i\n", &predictor);
  fscanf(FID,"msbp = %i\n", &msbp);
  fscanf(FID,"maxcor = %i\n", &maxcor);
  fscanf(FID,"nlscoef = %lf\n", &nlscoef);
  fclose(FID);

  /* Call ARKodeSet routines to insert solver parameters */
  if (order != 0) {     /* order overrides btable */
    flag = ARKodeSetOrder(arkode_mem, order);
    if (flag != 0) {
      fprintf(stderr,"Error in ARKodeSetOrder = %i\n",flag);
      return(1);
    }
  } else if (btable != 0) {
    flag = ARKodeSetIRKTableNum(arkode_mem, btable);
    if (flag != 0) {
      fprintf(stderr,"Error in ARKodeSetIRKTableNum = %i\n",flag);
      return(1);
    }
  }
  flag = ARKodeSetDenseOrder(arkode_mem, dense_order);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetDenseOrder = %i\n",flag);
    return(1);
  }
  flag = ARKodeSetAdaptivityMethod(arkode_mem, adapt_method, adapt_params);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetAdaptMethod = %i\n",flag);
    return(1);
  }
  flag = ARKodeSetAdaptivityConstants(arkode_mem, etamx1, etamxf, etacf, small_nef);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetAdaptConstants = %i\n",flag);
    return(1);
  }
  flag = ARKodeSetNewtonConstants(arkode_mem, crdown, rdiv);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetNewtonConstants = %i\n",flag);
    return(1);
  }
  flag = ARKodeSetLSetupConstants(arkode_mem, dgmax, msbp);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetLSetupConstants = %i\n",flag);
    return(1);
  }
  flag = ARKodeSetPredictorMethod(arkode_mem, predictor);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetPredictorMethod = %i\n",flag);
    return(1);
  }
  flag = ARKodeSetMaxNonlinIters(arkode_mem, maxcor);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetMaxNonlinIters = %i\n",flag);
    return(1);
  }
  flag = ARKodeSetNonlinConvCoef(arkode_mem, nlscoef);
  if (flag != 0) {
    fprintf(stderr,"Error in ARKodeSetMaxNonlinIters = %i\n",flag);
    return(1);
  }

  /* If (dense_order != -1), tell integrator to use dense output */
  if (dense_order != -1) {
    *idense = 1;
  } else {
    *idense = 0;
  }

  return(0);
}

/*---- end of file ----*/
