/*---------------------------------------------------------------
 $Revision: $
 $Date: $
-----------------------------------------------------------------
 Programmer(s): Daniel R. Reynolds @ SMU
-----------------------------------------------------------------
 Example problem:
 
 The following is a simple example problem with analytical 
 solution,
    dy/dt = A*y
 where A = V*D*Vi, 
      V = [1 -1 1; -1 2 1; 0 -1 2];
      Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1];
      D = [-0.5 0 0; 0 -0.1 0; 0 0 lam];
 where lam is a large negative number. The analytical solution 
 to this problem is 
   Y(t) = V*exp(D*t)*Vi*Y0
 for t in the interval [0.0, 0.05], with initial condition: 
 y(0) = [1,1,1]'.
 
 The stiffness of the problem is directly proportional to the 
 value of "lamda", which is specified through an input file, 
 along with the desired relative and absolute tolerances.  The
 value of lamda should be negative to result in a well-posed 
 ODE; for values with magnitude larger than 100 the problem 
 becomes quite stiff.

 In the example input file, we choose lamda = -100.
 
 This program solves the problem with the DIRK method,
 Newton iteration with the ARKDENSE dense linear solver, and a
 user-supplied Jacobian routine.
 Output is printed every 1.0 units of time (10 total).
 Run statistics (optional outputs) are printed at the end.
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
static int sol(realtype t, realtype lam, N_Vector y);
static int dense_MM(DlsMat A, DlsMat B, DlsMat C);

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
  realtype Tf = RCONST(0.05);      // final time
  realtype dTout = RCONST(0.005);  // time between outputs
  long int NEQ = 3;                // number of dependent vars.

  // declare solver parameters 
  int flag;                      // reusable error-checking flag
  int dense_order;               // order of accuracy for dense output
  int idense;                    // flag denoting tstop vs interpolated output
  int imex;                      // flag denoting integrator type
  int fixedpt;                   // flag denoting use of fixed-point nonlinear solver
  N_Vector y = NULL;             // empty vector for storing solution
  N_Vector ytrue = NULL;         // empty vector for analytical solution
  void *arkode_mem = NULL;       // empty ARKode memory structure

  /* read problem parameter and tolerances from input file:
     lamda  - problem stiffness parameter */
  double lamda_;
  FILE *FID;
  FID=fopen("input_analytic_sys.txt","r");
  flag = fscanf(FID,"  lamda = %lf\n",  &lamda_);
  fclose(FID);

  // convert the input to 'realtype' format 
  realtype lamda  = lamda_;

  // open solver diagnostics output file for writing 
  FILE *DFID;
  DFID=fopen("diags_ark_analytic_sys.txt","w");
  
  // Initial problem output 
  cout << "\nAnalytical ODE test problem:\n";
  cout << "    lamda = " << lamda << "\n\n";

  // Initialize data structures
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
  ytrue = N_VNew_Serial(NEQ);
  if (check_flag((void *)ytrue, "N_VNew_Serial", 0)) return 1;
  NV_Ith_S(y,0) = 1.0;
  NV_Ith_S(y,1) = 1.0;
  NV_Ith_S(y,2) = 1.0;
  arkode_mem = ARKodeCreate();
  if (check_flag((void *)arkode_mem, "ARKodeCreate", 0)) return 1;
  
  // Call init_from_file helper routine to read and set solver parameters 
  realtype rtol, atol;
  flag = init_from_file(arkode_mem, "solve_params.txt", f, fe, fi, T0,
			y, &imex, &dense_order, &fixedpt, &rtol, &atol);
  if (check_flag(&flag, "init_from_file", 1)) return 1;
  if (rtol <= 0.0)  rtol = 1.e-6;
  if (atol <= 0.0)  atol = 1.e-10;
  realtype reltol = rtol;
  realtype abstol = atol;
  
  // If (dense_order == -1), tell integrator to use tstop 
  if (dense_order == -1) {
    idense = 0;
  } else {    // otherwise tell integrator to use dense output 
    idense = 1;
  }

  // Set routines
  flag = ARKodeSetUserData(arkode_mem, (void *) &lamda);
  if (check_flag(&flag, "ARKodeSetUserData", 1)) return 1;
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
  realtype y0, y1, y2, yt0, yt1, yt2;
  realtype y0err, y1err, y2err, err2=0.0, errI=0.0;
  int Nt=0;
  printf("      t        y0        y1        y2        err0        err1        err2\n");
  printf("   --------------------------------------------------------------------------\n");
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
    y2 = NV_Ith_S(y,2);
    flag = sol(t, lamda, ytrue);
    yt0 = NV_Ith_S(ytrue,0);
    yt1 = NV_Ith_S(ytrue,1);
    yt2 = NV_Ith_S(ytrue,2);
    y0err = fabs(y0-yt0);
    y1err = fabs(y1-yt1);
    y2err = fabs(y2-yt2);
    errI = (errI > y0err) ? errI : y0err;
    errI = (errI > y1err) ? errI : y1err;
    errI = (errI > y2err) ? errI : y2err;
    err2 += y0err*y0err + y1err*y1err + y2err*y2err;
    Nt++;
    printf("  %8.4f  %8.5f  %8.5f  %8.5f  %10.3e  %10.3e  %10.3e\n", 
	   t, y0, y1, y2, y0err, y1err, y2err);
  }
  err2 = sqrt(err2 / 3.0 / Nt);
  printf("   --------------------------------------------------------------------------\n");

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
  cout << "   Error: max = " << errI << ", rms = " << err2 << "\n";
  cout << "   Oversolve = " << reltol/(err2+1.e-10*reltol) << "\n\n";

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
  realtype *rdata = (realtype *) user_data;   // cast user_data to realtype
  realtype lam = rdata[0];                    // set shortcut for stiffness parameter
  realtype y0 = NV_Ith_S(y,0);                // access current solution values
  realtype y1 = NV_Ith_S(y,1);
  realtype y2 = NV_Ith_S(y,2);
  realtype yd0, yd1, yd2;

  // fill in the RHS function: f(t,y) = V*D*Vi*y
  yd0 = 0.25*(5.0*y0 + 1.0*y1 - 3.0*y2);     // yd = Vi*y
  yd1 = 0.25*(2.0*y0 + 2.0*y1 - 2.0*y2);
  yd2 = 0.25*(1.0*y0 + 1.0*y1 + 1.0*y2);
  y0  = -0.5*yd0;                            //  y = D*yd
  y1  = -0.1*yd1;
  y2  =  lam*yd2;
  yd0 =  1.0*y0 - 1.0*y1 + 1.0*y2;           // yd = V*y
  yd1 = -1.0*y0 + 2.0*y1 + 1.0*y2;
  yd2 =  0.0*y0 - 1.0*y1 + 2.0*y2;
  NV_Ith_S(ydot,0) = yd0;
  NV_Ith_S(ydot,1) = yd1;
  NV_Ith_S(ydot,2) = yd2;
 
  return 0;                                  // Return with success
}

// fe routine to compute the explicit portion of f(t,y). 
static int fe(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;   // cast user_data to realtype
  realtype lam = rdata[0];                    // set shortcut for stiffness parameter
  realtype y0 = NV_Ith_S(y,0);                // access current solution values
  realtype y1 = NV_Ith_S(y,1);
  realtype y2 = NV_Ith_S(y,2);
  realtype yd0, yd1, yd2;

  // fill in the RHS function: f(t,y) = V*D*Vi*y
  yd0 = 0.25*(5.0*y0 + 1.0*y1 - 3.0*y2);     // yd = Vi*y
  yd1 = 0.25*(2.0*y0 + 2.0*y1 - 2.0*y2);
  yd2 = 0.25*(1.0*y0 + 1.0*y1 + 1.0*y2);
  y0  = -0.5*yd0;                            //  y = D*yd
  y1  = -0.1*yd1;
  y2  =  0.0;                                // note: lam=0 here
  yd0 =  1.0*y0 - 1.0*y1 + 1.0*y2;           // yd = V*y
  yd1 = -1.0*y0 + 2.0*y1 + 1.0*y2;
  yd2 =  0.0*y0 - 1.0*y1 + 2.0*y2;
  NV_Ith_S(ydot,0) = yd0;
  NV_Ith_S(ydot,1) = yd1;
  NV_Ith_S(ydot,2) = yd2;
 
  return 0;                                  // Return with success
}

// fi routine to compute the implicit portion of f(t,y). 
static int fi(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *rdata = (realtype *) user_data;   // cast user_data to realtype
  realtype lam = rdata[0];                    // set shortcut for stiffness parameter
  realtype y0 = NV_Ith_S(y,0);                // access current solution values
  realtype y1 = NV_Ith_S(y,1);
  realtype y2 = NV_Ith_S(y,2);
  realtype yd0, yd1, yd2;

  // fill in the RHS function: f(t,y) = V*D*Vi*y
  yd0 = 0.25*(5.0*y0 + 1.0*y1 - 3.0*y2);     // yd = Vi*y
  yd1 = 0.25*(2.0*y0 + 2.0*y1 - 2.0*y2);
  yd2 = 0.25*(1.0*y0 + 1.0*y1 + 1.0*y2);
  y0  = 0.0;                                 //  y = D*yd
  y1  = 0.0;                                 // note: all but lam are 0.0
  y2  = lam*yd2;
  yd0 =  1.0*y0 - 1.0*y1 + 1.0*y2;           // yd = V*y
  yd1 = -1.0*y0 + 2.0*y1 + 1.0*y2;
  yd2 =  0.0*y0 - 1.0*y1 + 2.0*y2;
  NV_Ith_S(ydot,0) = yd0;
  NV_Ith_S(ydot,1) = yd1;
  NV_Ith_S(ydot,2) = yd2;
 
  return 0;                                  // Return with success
}

// Jacobian routine to compute J(t,y) = df/dy. 
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *rdata = (realtype *) user_data;   // cast user_data to realtype
  realtype lam = rdata[0];                    // set shortcut for stiffness parameter
  DlsMat V  = NewDenseMat(3,3);               // create temporary DlsMat objects
  DlsMat D  = NewDenseMat(3,3);
  DlsMat Vi = NewDenseMat(3,3);

  DenseScale(0.0, V);     // initialize temporary matrices to zero
  DenseScale(0.0, D);
  DenseScale(0.0, Vi);

  // Fill in temporary matrices:
  //    V = [1 -1 1; -1 2 1; 0 -1 2]
  DENSE_ELEM(V,0,0) =  1.0;
  DENSE_ELEM(V,0,1) = -1.0;
  DENSE_ELEM(V,0,2) =  1.0;
  DENSE_ELEM(V,1,0) = -1.0;
  DENSE_ELEM(V,1,1) =  2.0;
  DENSE_ELEM(V,1,2) =  1.0;
  DENSE_ELEM(V,2,0) =  0.0;
  DENSE_ELEM(V,2,1) = -1.0;
  DENSE_ELEM(V,2,2) =  2.0;

  //    Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1]
  DENSE_ELEM(Vi,0,0) =  0.25*5.0;
  DENSE_ELEM(Vi,0,1) =  0.25*1.0;
  DENSE_ELEM(Vi,0,2) = -0.25*3.0;
  DENSE_ELEM(Vi,1,0) =  0.25*2.0;
  DENSE_ELEM(Vi,1,1) =  0.25*2.0;
  DENSE_ELEM(Vi,1,2) = -0.25*2.0;
  DENSE_ELEM(Vi,2,0) =  0.25*1.0;
  DENSE_ELEM(Vi,2,1) =  0.25*1.0;
  DENSE_ELEM(Vi,2,2) =  0.25*1.0;

  //    D = [-0.5 0 0; 0 -0.1 0; 0 0 lam]
  DENSE_ELEM(D,0,0) = -0.5;
  DENSE_ELEM(D,1,1) = -0.1;
  DENSE_ELEM(D,2,2) = lam;

  // Compute J = V*D*Vi
  if (dense_MM(D,Vi,J) != 0) {     // J = D*Vi
    cerr << "matmul error\n";
    return 1;
  }
  if (dense_MM(V,J,D) != 0) {      // D = V*J [= V*D*Vi]
    cerr << "matmul error\n";
    return 1;
  }
  DenseCopy(D, J);                 // J = D [= V*D*Vi]

  return 0;                        // Return with success
}


// Jacobian routine to compute J(t,y) = dfi/dy. 
static int JacI(long int N, realtype t,
		N_Vector y, N_Vector fy, DlsMat J, void *user_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype *rdata = (realtype *) user_data;   // cast user_data to realtype
  realtype lam = rdata[0];                    // set shortcut for stiffness parameter
  DlsMat V  = NewDenseMat(3,3);               // create temporary DlsMat objects
  DlsMat D  = NewDenseMat(3,3);
  DlsMat Vi = NewDenseMat(3,3);

  DenseScale(0.0, V);     // initialize temporary matrices to zero
  DenseScale(0.0, D);
  DenseScale(0.0, Vi);

  // Fill in temporary matrices:
  //    V = [1 -1 1; -1 2 1; 0 -1 2]
  DENSE_ELEM(V,0,0) =  1.0;
  DENSE_ELEM(V,0,1) = -1.0;
  DENSE_ELEM(V,0,2) =  1.0;
  DENSE_ELEM(V,1,0) = -1.0;
  DENSE_ELEM(V,1,1) =  2.0;
  DENSE_ELEM(V,1,2) =  1.0;
  DENSE_ELEM(V,2,0) =  0.0;
  DENSE_ELEM(V,2,1) = -1.0;
  DENSE_ELEM(V,2,2) =  2.0;

  //    Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1]
  DENSE_ELEM(Vi,0,0) =  0.25*5.0;
  DENSE_ELEM(Vi,0,1) =  0.25*1.0;
  DENSE_ELEM(Vi,0,2) = -0.25*3.0;
  DENSE_ELEM(Vi,1,0) =  0.25*2.0;
  DENSE_ELEM(Vi,1,1) =  0.25*2.0;
  DENSE_ELEM(Vi,1,2) = -0.25*2.0;
  DENSE_ELEM(Vi,2,0) =  0.25*1.0;
  DENSE_ELEM(Vi,2,1) =  0.25*1.0;
  DENSE_ELEM(Vi,2,2) =  0.25*1.0;

  //    D = [-0.5 0 0; 0 -0.1 0; 0 0 lam]
  DENSE_ELEM(D,2,2) = lam;

  // Compute J = V*D*Vi
  if (dense_MM(D,Vi,J) != 0) {     // J = D*Vi
    cerr << "matmul error\n";
    return 1;
  }
  if (dense_MM(V,J,D) != 0) {      // D = V*J [= V*D*Vi]
    cerr << "matmul error\n";
    return 1;
  }
  DenseCopy(D, J);                 // J = D [= V*D*Vi]

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
static int sol(realtype t, realtype lam, N_Vector y)
{
  realtype y0, y1, y2;
  realtype x0 = 1.0, x1 = 1.0, x2 = 1.0;
  y0 = 0.25*(5.0*x0 + 1.0*x1 - 3.0*x2);  //   y = Vi*x [= Vi*y0] 
  y1 = 0.25*(2.0*x0 + 2.0*x1 - 2.0*x2);
  y2 = 0.25*(1.0*x0 + 1.0*x1 + 1.0*x2);
  x0  = exp(-0.5*t)*y0;                  //   x = exp(D*t)*y 
  x1  = exp(-0.1*t)*y1;
  x2  = exp( lam*t)*y2;
  y0 =  1.0*x0 - 1.0*x1 + 1.0*x2;        //   y = V*x 
  y1 = -1.0*x0 + 2.0*x1 + 1.0*x2;
  y2 =  0.0*x0 - 1.0*x1 + 2.0*x2;
  NV_Ith_S(y,0) = y0;
  NV_Ith_S(y,1) = y1;
  NV_Ith_S(y,2) = y2;

  return 0;                              // Return with success
}

// DlsMat matrix-multiply utility routine: C = A*B.
static int dense_MM(DlsMat A, DlsMat B, DlsMat C)
{
  // check for legal dimensions
  if ((A->N != B->M) || (C->M != A->M) || (C->N != B->N)) {
    cerr << "\n matmul error: dimension mismatch\n\n";
    return 1;
  }

  realtype **adata = A->cols;     // access data and extents
  realtype **bdata = B->cols;
  realtype **cdata = C->cols;
  long int m = C->M;
  long int n = C->N;
  long int l = A->N;
  int i, j, k;
  DenseScale(0.0, C);             // initialize output

  // perform multiply (not optimal, but fine for 3x3 matrices)
  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
      for (k=0; k<l; k++)
     cdata[i][j] += adata[i][k] * bdata[k][j];

  return 0;                       // Return with success
}


//---- end of file ----
