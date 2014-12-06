/*---------------------------------------------------------------
 $Revision: $
 $Date: $
-----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Example problem:
 
 The following is a simple example problem with analytical 
 solution,
    du/dt = 100*u - 400*v
    dv/dt = 100*u + 100*v,
 with initial conditions u(0) = v(0) = 2.  The analytical 
 solution to this problem is 
    u(t) = exp(100*t)*(2*cos(200*t) - 4*sin(200*t)),
    v(t) = exp(100*t)*(2*cos(200*t) + sin(200*t)),
 for t in the interval [0.0, 0.25].
 
 The problem itself is unstable, having eigenvalues 100 +/- 200i.
 The goal of this test problem is to investigate whether our
 built-in integration methods incorrectly obtain stable 
 solutions for this unstable problem.
---------------------------------------------------------------*/

// Header files 
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <arkode/arkode.h>
#include <nvector/nvector_serial.h>
#include <arkode/arkode_dense.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>

using namespace std;

// User-supplied Functions Called by the Solver 
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int JacI(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Private function to check function return values 
static int check_flag(void *flagvalue, const string funcname, int opt);
static int sol(realtype t, N_Vector y);

// Parameter input helper function 
int init_from_file(void *ark_mem, char *fname, ARKRhsFn f, 
		   ARKRhsFn fe, ARKRhsFn fi, realtype T0, 
		   N_Vector y0, int *ImEx, int *dorder, 
		   int *fxpt, realtype *RTol, realtype *ATol);

// Main Program 
int main()
{
  // general problem parameters 
  realtype T0 = RCONST(0.0);       // initial time
  realtype Tf = RCONST(0.25);      // final time
  realtype dTout = RCONST(0.025);  // time between outputs
  long int NEQ = 2;                // number of dependent vars.

  // declare solver parameters 
  int flag;                      // reusable error-checking flag
  int dense_order;               // order of accuracy for dense output
  int idense;                    // flag denoting tstop vs interpolated output
  int imex;                      // flag denoting integrator type
  int fixedpt;                   // flag denoting use of fixed-point nonlinear solver
  N_Vector y = NULL;             // empty vector for storing solution
  N_Vector ytrue = NULL;         // empty vector for analytical solution
  void *arkode_mem = NULL;       // empty ARKode memory structure

  // open solver diagnostics output file for writing 
  FILE *DFID;
  DFID=fopen("diags_ark_hyperstable.txt","w");
  
  // Initial problem output 
  cout << "\nHyperstability test problem:\n\n";

  // Initialize data structures
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  ytrue = N_VNew_Serial(NEQ);
  if (check_flag((void *)ytrue, "N_VNew_Serial", 0)) return 1;
  NV_Ith_S(y,0) = 2.0;
  NV_Ith_S(y,1) = 2.0;
  arkode_mem = ARKodeCreate();
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return 1;
  
  // Call init_from_file helper routine to read and set solver parameters 
  realtype rtol, atol;
  flag = init_from_file(arkode_mem, "solve_params.txt", f, fe, fi, T0,
			y, &imex, &dense_order, &fixedpt, &rtol, &atol);
  if (check_flag(&flag, "init_from_file", 1)) return 1;
  if (rtol <= 0.0)  rtol = 1.e-6;
  if (atol <= 0.0)  atol = 1.e-8;
  realtype reltol = rtol;
  realtype abstol = atol;
  
  // If (dense_order == -1), tell integrator to use tstop 
  if (dense_order == -1) {
    idense = 0;
  } else {    // otherwise tell integrator to use dense output 
    idense = 1;
  }

  // Set routines
  flag = ARKodeSetDiagnostics(arkode_mem, DFID);
  if (check_flag(&flag, "ARKodeSetDiagnostics", 1)) return 1;
  flag = ARKodeSetMaxNumSteps(arkode_mem, 10000);
  if (check_flag(&flag, "ARKodeSetMaxNumSteps", 1)) return 1;
  flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
  if (check_flag(&flag, "ARKodeSStolerances", 1)) return 1;

  // Linear solver specification
  flag = ARKDense(arkode_mem, NEQ);        // Specify dense linear solver 
  if (check_flag(&flag, "ARKDense", 1)) return 1;
  switch (imex) {                          // Set Jacobian routine
  case 0:         // purely implicit 
    flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);   break;
    if (check_flag(&flag, "ARKDlsSetDenseJacFn", 1)) return 1;
  case 1:         // purely explicit 
    break;
  default:        // imex 
    flag = ARKDlsSetDenseJacFn(arkode_mem, JacI);  break;
    if (check_flag(&flag, "ARKDlsSetDenseJacFn", 1)) return 1;
  }

  // Write all solver parameters to stdout 
  cout << "\n";
  flag = ARKodeWriteParameters(arkode_mem, stdout);
  if (check_flag(&flag, "ARKodeWriteParameters", 1)) return 1;

  /* Main time-stepping loop: calls ARKode to perform the integration, then
     prints results.  Stops when the final time has been reached */
  realtype t = T0;
  realtype tout = T0+dTout;
  realtype y0, y1, yt0, yt1;
  realtype y0err, y1err, err2=0.0, errI=0.0;
  int Nt=0;
  printf("      t        y0          y1          err0        err1\n");
  printf("   --------------------------------------------------------\n");
  while (Tf - t > 1.0e-8) {

    if (!idense)
      flag = ARKodeSetStopTime(arkode_mem, tout);
    flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);   // call integrator
    if (flag >= 0) {                                      // successful solve: update time
      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    } else {                                              // unsuccessful solve: break
      cerr << "Solver failure, stopping integration\n";
      return 1;
    }
    if (check_flag(&flag, "ARKode", 1)) break;
    y0 = NV_Ith_S(y,0);                                   // access/print solution & error
    y1 = NV_Ith_S(y,1);
    flag = sol(t, ytrue);
    yt0 = NV_Ith_S(ytrue,0);
    yt1 = NV_Ith_S(ytrue,1);
    y0err = fabs(y0-yt0)/fabs(yt0);
    y1err = fabs(y1-yt1)/fabs(yt1);
    printf("  %8.4f  %10.3e  %10.3e  %10.3e  %10.3e\n", t, y0, y1, y0err, y1err);
  }
  printf("   --------------------------------------------------------\n");

  // Print some final statistics 
  long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
  flag = ARKodeGetNumSteps(arkode_mem, &nst);
  check_flag(&flag, "ARKodeGetNumSteps", 1);
  flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  check_flag(&flag, "ARKodeGetNumStepAttempts", 1);
  flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  check_flag(&flag, "ARKodeGetNumRhsEvals", 1);
  flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
  check_flag(&flag, "ARKodeGetNumLinSolvSetups", 1);
  flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
  check_flag(&flag, "ARKodeGetNumErrTestFails", 1);
  flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  check_flag(&flag, "ARKodeGetNumNonlinSolvIters", 1);
  flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  check_flag(&flag, "ARKodeGetNumNonlinSolvConvFails", 1);
  flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
  check_flag(&flag, "ARKDlsGetNumJacEvals", 1);
  flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
  check_flag(&flag, "ARKDlsGetNumRhsEvals", 1);

  cout << "\nFinal Solver Statistics:\n";
  cout << "   Internal solver steps = " << nst << " (attempted = " << nst_a << ")\n";
  cout << "   Total RHS evals:  Fe = " << nfe << ",  Fi = " << nfi << "\n";
  cout << "   Total linear solver setups = " << nsetups << "\n";
  cout << "   Total RHS evals for setting up the linear system = " << nfeLS << "\n";
  cout << "   Total number of Jacobian evaluations = " << nje << "\n";
  cout << "   Total number of nonlinear iterations = " << nni << "\n";
  cout << "   Total number of nonlinear solver convergence failures = " << ncfn << "\n";
  cout << "   Total number of error test failures = " << netf << "\n\n";

  // Clean up and return with successful completion
  N_VDestroy_Serial(y);        // Free vectors
  N_VDestroy_Serial(ytrue);
  ARKodeFree(&arkode_mem);     // Free integrator memory
  fclose(DFID);                // close diagnostics output file
  return 0;
}

/*-------------------------------
 * Functions called by the solver
 *-------------------------------*/

// f routine to compute the ODE RHS function f(t,y). 
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y0 = NV_Ith_S(y,0);                // access current solution values
  realtype y1 = NV_Ith_S(y,1);

  // fill in the RHS function
  NV_Ith_S(ydot,0) = 100.0*y0 - 400.0*y1;
  NV_Ith_S(ydot,1) = 100.0*y0 + 100.0*y1;
 
  return 0;                                  // Return with success
}

// fe routine to compute the explicit portion of f(t,y). 
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y0 = NV_Ith_S(y,0);                // access current solution values
  realtype y1 = NV_Ith_S(y,1);

  // fill in the RHS function
  NV_Ith_S(ydot,0) = 100.0*y0;
  NV_Ith_S(ydot,1) = 100.0*y0 + 100.0*y1;
 
  return 0;                                  // Return with success
}

// fi routine to compute the implicit portion of f(t,y). 
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y0 = NV_Ith_S(y,0);                // access current solution values
  realtype y1 = NV_Ith_S(y,1);

  // fill in the RHS function
  NV_Ith_S(ydot,0) = -400.0*y1;
  NV_Ith_S(ydot,1) = 0.0;
 
  return 0;                                  // Return with success
}

// Jacobian routine to compute J(t,y) = df/dy. 
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Fill in Jacobian matrix
  DENSE_ELEM(J,0,0) =  100.0;
  DENSE_ELEM(J,0,1) = -400.0;
  DENSE_ELEM(J,1,0) =  100.0;
  DENSE_ELEM(J,1,1) =  100.0;
  return 0;                        // Return with success
}


// Jacobian routine to compute J(t,y) = dfi/dy. 
static int JacI(long int N, realtype t,
		N_Vector y, N_Vector fy, DlsMat J, void *user_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  // Fill in Jacobian matrix
  DENSE_ELEM(J,0,0) =  0.0;
  DENSE_ELEM(J,0,1) = -400.0;
  DENSE_ELEM(J,1,0) =  0.0;
  DENSE_ELEM(J,1,1) =  0.0;
  return 0;                        // Return with success
}



/*-------------------------------
 * Private helper functions
 *-------------------------------*/

/* Check function return value...
    opt == 0 means SUNDIALS function allocates memory so check if
             returned NULL pointer
    opt == 1 means SUNDIALS function returns a flag so check if
             flag >= 0
    opt == 2 means function allocates memory so check if returned
             NULL pointer  
*/
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

// sol routine to compute the ODE solution y(t). 
static int sol(realtype t, N_Vector y)
{
  NV_Ith_S(y,0) = exp(100.0*t)*(2.0*cos(200.0*t) - 4.0*sin(200.0*t));
  NV_Ith_S(y,1) = exp(100.0*t)*(2.0*cos(200.0*t) + sin(200.0*t));
  return 0;                              // Return with success
}


//---- end of file ----
