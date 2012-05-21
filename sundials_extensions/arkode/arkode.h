/* -----------------------------------------------------------------
 * $Revision: 1.0 $
 * $Date:  $
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * This is the interface file for the main ARKODE integrator.
 * -----------------------------------------------------------------
 *
 * ARKODE is used to solve numerically the ordinary initial value
 * problem:
 *
 *               M(t)*y' = f(t,y),
 *                 y(t0) = y0,
 *
 * where t0, y0 in R^N, M(t)*y: R x R^N -> R^N, and 
 * f: R x R^N -> R^N are given.
 *
 * ----------------------------------------------------------------- */

#ifndef _ARKODE_H
#define _ARKODE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdio.h>

#include <sundials/sundials_nvector.h>

/* =================================================================
 *                        ARKODE CONSTANTS
 * ================================================================= */

/* -----------------------------------------------------------------
 * Enumerations for inputs to ARKodeCreate and ARKode.
 * -----------------------------------------------------------------
 * Symbolic constants for the lmm and iter parameters to ARKodeCreate
 * and the input parameter itask to ARKode, are given below.
 *
 * itask: The itask input parameter to ARKode indicates the job
 *        of the solver for the next user step. The ARK_NORMAL
 *        itask is to have the solver take internal steps until
 *        it has reached or just passed the user specified tout
 *        parameter. The solver then interpolates in order to
 *        return an approximate value of y(tout). The ARK_ONE_STEP
 *        option tells the solver to just take one internal step
 *        and return the solution at the point reached by that step.
 * ----------------------------------------------------------------- */

/* itask */
#define ARK_NORMAL         1
#define ARK_ONE_STEP       2

/* ----------------------------------------
 * ARKODE return flags
 * ---------------------------------------- */

#define ARK_SUCCESS               0
#define ARK_TSTOP_RETURN          1
#define ARK_ROOT_RETURN           2

#define ARK_WARNING              99

#define ARK_TOO_MUCH_WORK        -1
#define ARK_TOO_MUCH_ACC         -2
#define ARK_ERR_FAILURE          -3
#define ARK_CONV_FAILURE         -4

#define ARK_LINIT_FAIL           -5
#define ARK_LSETUP_FAIL          -6
#define ARK_LSOLVE_FAIL          -7
#define ARK_RHSFUNC_FAIL         -8
#define ARK_FIRST_RHSFUNC_ERR    -9
#define ARK_REPTD_RHSFUNC_ERR    -10
#define ARK_UNREC_RHSFUNC_ERR    -11
#define ARK_RTFUNC_FAIL          -12

#define ARK_MEM_FAIL             -20
#define ARK_MEM_NULL             -21
#define ARK_ILL_INPUT            -22
#define ARK_NO_MALLOC            -23
#define ARK_BAD_K                -24
#define ARK_BAD_T                -25
#define ARK_BAD_DKY              -26
#define ARK_TOO_CLOSE            -27

/* =================================================================
 *                         FUNCTION TYPES
 * ================================================================= */

/* -----------------------------------------------------------------
 * Type : ARKRhsFn
 * -----------------------------------------------------------------
 * The f functions which define the right hand side of the ODE
 * system M*y' = fe(t,y) + fi(t,y) must have type ARKRhsFn.
 * f takes as input the independent variable value t, and the
 * dependent variable vector y.  It stores the result of fe(t,y) 
 * or fi(t,y) in the vector ydot.  The y and ydot arguments are of 
 * type N_Vector.
 * (Allocation of memory for ydot is handled within ARKODE)
 * The user_data parameter is the same as the user_data
 * parameter set by the user through the ARKodeSetUserData routine.
 * This user-supplied pointer is passed to the user's fe or fi
 * function every time it is called.
 *
 * A ARKRhsFn should return 0 if successful, a negative value if
 * an unrecoverable error occured, and a positive value if a 
 * recoverable error (e.g. invalid y values) occured. 
 * If an unrecoverable occured, the integration is halted. 
 * If a recoverable error occured, then (in most cases) ARKODE
 * will try to correct and retry.
 * ----------------------------------------------------------------- */

typedef int (*ARKRhsFn)(realtype t, N_Vector y,
			N_Vector ydot, void *user_data);

/* -----------------------------------------------------------------
 * Type : ARKRootFn
 * -----------------------------------------------------------------
 * A function g, which defines a set of functions g_i(t,y) whose
 * roots are sought during the integration, must have type ARKRootFn.
 * The function g takes as input the independent variable value
 * t, and the dependent variable vector y.  It stores the nrtfn
 * values g_i(t,y) in the realtype array gout.
 * (Allocation of memory for gout is handled within ARKODE.)
 * The user_data parameter is the same as that passed by the user
 * to the ARKodeSetUserData routine.  This user-supplied pointer is
 * passed to the user's g function every time it is called.
 *
 * A ARKRootFn should return 0 if successful or a non-zero value
 * if an error occured (in which case the integration will be halted).
 * ----------------------------------------------------------------- */

typedef int (*ARKRootFn)(realtype t, N_Vector y, 
			 realtype *gout, void *user_data);

/* -----------------------------------------------------------------
 * Type : ARKEwtFn
 * -----------------------------------------------------------------
 * A function e, which sets the error weight vector ewt, must have
 * type ARKEwtFn.  The function e takes as input the current 
 * dependent variable y. It must set the vector of error weights 
 * used in the WRMS norm:
 * 
 *   ||y||_WRMS = sqrt [ 1/N * sum ( ewt_i * y_i)^2 ]
 *
 * Typically, the vector ewt has components:
 * 
 *   ewt_i = 1 / (reltol * |y_i| + abstol_i)
 *
 * The user_data parameter is the same as that passed by the user
 * to the ARKodeSetUserData routine.  This user-supplied pointer is
 * passed to the user's e function every time it is called.
 * A ARKEwtFn e must return 0 if the error weight vector has been
 * successfuly set and a non-zero value otherwise.
 * ----------------------------------------------------------------- */

typedef int (*ARKEwtFn)(N_Vector y, N_Vector ewt, void *user_data);

/* -----------------------------------------------------------------
 * Type : ARKErrHandlerFn
 * -----------------------------------------------------------------
 * A function eh, which handles error messages, must have type
 * ARKErrHandlerFn.  The function eh takes as input the error code, 
 * the name of the module reporting the error, the error message, 
 * and a pointer to user data, the same as that passed to 
 * ARKodeSetUserData.
 * 
 * All error codes are negative, except ARK_WARNING which indicates 
 * a warning (the solver continues).
 *
 * A ARKErrHandlerFn has no return value.
 * ----------------------------------------------------------------- */
  
typedef void (*ARKErrHandlerFn)(int error_code, const char *module, 
				const char *function, char *msg, 
				void *user_data); 

/* -----------------------------------------------------------------
 * Type : ARKAdaptFn
 * -----------------------------------------------------------------
 * A function which sets the new time step h, must have type 
 * ARKAdaptFn.  The function takes as input the current dependent 
 * variable y, the current time t, the current time step size h, 
 * the last 3 error estimates, the method order q, the embedding 
 * order p, and a pointer to user data. The function must set the 
 * scalar step size for the upcoming time step.  This value will 
 * subsequently be bounded by the user-supplied values for the 
 * minimum and maximum allowed time step, and the time step 
 * satisfying the explicit stability restriction.  The user_data 
 * parameter is the same as that passed by the user to the 
 * ARKodeSetUserData routine.  This user-supplied pointer is passed 
 * to the function every time it is called.
 *
 * A ARKAdaptFn must return 0 if the new time step has been
 * successfuly set and a non-zero value otherwise.
 * ----------------------------------------------------------------- */

typedef int (*ARKAdaptFn)(N_Vector y, realtype t, realtype h, 
			  realtype e1, realtype e2, realtype e3, 
			  int q, int p, realtype *hnew, 
			  void *user_data);

/* -----------------------------------------------------------------
 * Type : ARKExpStabFn
 * -----------------------------------------------------------------
 * A function which returns the time step satisfying the stability 
 * restriction for the explicit portion of the ODE.  The function 
 * takes as input the current dependent variable y, the current 
 * time t, and a pointer to user data. The function must set the 
 * scalar step size satisfying the stability restriction for the 
 * upcoming time step.  This value will subsequently be bounded by 
 * the user-supplied values for the minimum and maximum allowed 
 * time step, and the accuracy-based time step.  The user_data 
 * parameter is the same as that passed by the user to the 
 * ARKodeSetUserData routine.  This user-supplied pointer is passed
 * to the function every time it is called.
 *
 * If this function is not supplied (NULL), or if it returns a
 * negative time step size, then ARKode will assume that there is
 * no explicit stability restriction on the time step size.
 *
 * A ARKExpStabFn must return 0 if the step size limit has been
 * successfuly set and a non-zero value otherwise.
 * ----------------------------------------------------------------- */

typedef int (*ARKExpStabFn)(N_Vector y, realtype t, 
			    realtype *hstab, void *user_data);

/* =================================================================
 *                       USER-CALLABLE ROUTINES
 * ================================================================= */

/* -----------------------------------------------------------------
 * Function : ARKodeCreate
 * -----------------------------------------------------------------
 * ARKodeCreate creates an internal memory block for a problem to
 * be solved by ARKODE.
 *
 * If successful, ARKodeCreate returns a pointer to initialized
 * problem memory. This pointer should be passed to ARKodeInit.
 * If an initialization error occurs, ARKodeCreate prints an error
 * message to standard err and returns NULL.
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT void *ARKodeCreate();

/* -----------------------------------------------------------------
 * Integrator optional input specification functions
 * -----------------------------------------------------------------
 * The following functions can be called to set optional inputs
 * to values other than the defaults given below:
 *
 * Function                 |  Optional input / [ default value ]
 * -----------------------------------------------------------------
 *                          |
 * ARKodeSetErrHandlerFn    | user-provided ErrHandler function.
 *                          | [internal]
 *                          |
 * ARKodeSetErrFile         | the file pointer for an error file
 *                          | where all ARKODE warning and error
 *                          | messages will be written if the default
 *                          | internal error handling function is used. 
 *                          | This parameter can be stdout (standard 
 *                          | output), stderr (standard error), or a 
 *                          | file pointer (corresponding to a user 
 *                          | error file opened for writing) returned 
 *                          | by fopen.
 *                          | If not called, then all messages will
 *                          | be written to the standard error stream.
 *                          | [stderr]
 *                          |
 * ARKodeSetUserData        | a pointer to user data that will be
 *                          | passed to the user's f function every
 *                          | time f is called.
 *                          | [NULL]
 *                          |
 * ARKodeSetOrd             | method order to be used by the solver.
 *                          | [4]
 *                          |
 * ARKodeSetERK             | specifies that implicit portion of 
 *                          | problem is disabled, and to use an 
 *                          | explicit RK method.
 *                          | [0]
 *                          |
 * ARKodeSetIRK             | specifies that explicit portion of 
 *                          | problem is disabled, and to use an 
 *                          | implicit RK method.
 *                          | [0]
 *                          |
 * ARKodeSetERKTable        | specifies to use a customized Butcher 
 *                          | table for the explicit portion of the 
 *                          | system.
 *                          | [determined by ARKODE based on order]
 *                          |
 * ARKodeSetIRKTable        | specifies to use a customized Butcher 
 *                          | table for the implicit portion of the 
 *                          | system.
 *                          | [determined by ARKODE based on order]
 *                          |
 * ARKodeSetMaxNumSteps     | maximum number of internal steps to be
 *                          | taken by the solver in its attempt to
 *                          | reach tout.
 *                          | [500]
 *                          |
 * ARKodeSetMaxHnilWarns    | maximum number of warning messages
 *                          | issued by the solver that t+h==t on the
 *                          | next internal step. A value of -1 means
 *                          | no such messages are issued.
 *                          | [10]
 *                          |
 * ARKodeSetInitStep        | initial step size.
 *                          | [estimated by ARKODE]
 *                          |
 * ARKodeSetMinStep         | minimum absolute value of step size
 *                          | allowed.
 *                          | [0.0]
 *                          |
 * ARKodeSetMaxStep         | maximum absolute value of step size
 *                          | allowed.
 *                          | [infinity]
 *                          |
 * ARKodeSetStopTime        | the independent variable value past
 *                          | which the solution is not to proceed.
 *                          | [infinity]
 *                          |
 * ARKodeSetAdapMethod      | Method to use for time step adaptivity
 *                          | [0]
 *                          |
 * ARKodeSetAdaptivityFn    | user-provided time step adaptivity function.
 *                          | [internal]
 *                          |
 * ARKodeSetMaxErrTestFails | Maximum number of error test failures
 *                          | in attempting one step.
 *                          | [7]
 *                          |
 * ARKodeSetMaxNonlinIters  | Maximum number of nonlinear solver
 *                          | iterations at one stage solution.
 *                          | [3]
 *                          |
 * ARKodeSetMaxConvFails    | Maximum number of convergence failures
 *                          | allowed in attempting one step.
 *                          | [10]
 *                          |
 * ARKodeSetNonlinConvCoef  | Coefficient in the nonlinear
 *                          | convergence test.
 *                          | [0.1]
 *                          |
 * ARKodeSetPredictMethod   | Method to use for prediction of new-time
 *                          | solutions
 *                          | [0]
 *                          |
 * -----------------------------------------------------------------
 *                             |
 * ARKodeSetRootDirection      | Specifies the direction of zero
 *                             | crossings to be monitored
 *                             | [both directions]
 *                             |
 * ARKodeSetNoInactiveRootWarn | disable warning about possible
 *                             | g==0 at beginning of integration
 *                             | 
 * -----------------------------------------------------------------
 * Return flag:
 *   ARK_SUCCESS   if successful
 *   ARK_MEM_NULL  if the arkode memory is NULL
 *   ARK_ILL_INPUT if an argument has an illegal value
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT int ARKodeSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun, void *eh_data);
SUNDIALS_EXPORT int ARKodeSetErrFile(void *arkode_mem, FILE *errfp);
SUNDIALS_EXPORT int ARKodeSetUserData(void *arkode_mem, void *user_data);
SUNDIALS_EXPORT int ARKodeSetOrd(void *arkode_mem, int maxord);
SUNDIALS_EXPORT int ARKodeSetERK(void *arkode_mem, int truefalse);
SUNDIALS_EXPORT int ARKodeSetIRK(void *arkode_mem, int truefalse);
SUNDIALS_EXPORT int ARKodeSetERKTable(void *arkode_mem, realtype s, 
				      realtype *c, realtype *A, realtype *b, 
				      realtype *bembed, realtype *bdense);
SUNDIALS_EXPORT int ARKodeSetIRKTable(void *arkode_mem, realtype s, 
				      realtype *c, realtype *A, realtype *b, 
				      realtype *bembed, realtype *bdense);
SUNDIALS_EXPORT int ARKodeSetMaxNumSteps(void *arkode_mem, long int mxsteps);
SUNDIALS_EXPORT int ARKodeSetMaxHnilWarns(void *arkode_mem, int mxhnil);
SUNDIALS_EXPORT int ARKodeSetInitStep(void *arkode_mem, realtype hin);
SUNDIALS_EXPORT int ARKodeSetMinStep(void *arkode_mem, realtype hmin);
SUNDIALS_EXPORT int ARKodeSetMaxStep(void *arkode_mem, realtype hmax);
SUNDIALS_EXPORT int ARKodeSetStopTime(void *arkode_mem, realtype tstop);
SUNDIALS_EXPORT int ARKodeSetAdaptMethod(void *arkode_mem, int imethod, 
					 realtype *adapt_params);
SUNDIALS_EXPORT int ARKodeSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun);
SUNDIALS_EXPORT int ARKodeSetMaxErrTestFails(void *arkode_mem, int maxnef);
SUNDIALS_EXPORT int ARKodeSetMaxNonlinIters(void *arkode_mem, int maxcor);
SUNDIALS_EXPORT int ARKodeSetMaxConvFails(void *arkode_mem, int maxncf);
SUNDIALS_EXPORT int ARKodeSetNonlinConvCoef(void *arkode_mem, realtype nlscoef);
SUNDIALS_EXPORT int ARKodeSetPredictMethod(void *arkode_mem, int imethod);

SUNDIALS_EXPORT int ARKodeSetRootDirection(void *arkode_mem, int *rootdir);
SUNDIALS_EXPORT int ARKodeSetNoInactiveRootWarn(void *arkode_mem);

/* -----------------------------------------------------------------
 * Function : ARKodeInit
 * -----------------------------------------------------------------
 * ARKodeInit allocates and initializes memory for a problem to
 * to be solved by ARKODE.
 *
 * arkode_mem is pointer to ARKODE memory returned by ARKodeCreate.
 *
 * fe      is the name of the C function defining the explicit 
 *         portion of the right-hand side function in 
 *                y' = fe(t,y) + fi(t,y).
 *
 * fi      is the name of the C function defining the implicit 
 *         portion of the right-hand side function in 
 *                y' = fe(t,y) + fi(t,y).
 *
 * EStab   is the name of the C function defining the stability 
 *         restriction for fe.
 *
 * t0      is the initial value of t.
 *
 * y0      is the initial condition vector y(t0).
 *
 * Return flag:
 *  ARK_SUCCESS   if successful
 *  ARK_MEM_NULL  if the arkode memory was NULL
 *  ARK_MEM_FAIL  if a memory allocation failed
 *  ARK_ILL_INPUT if an argument has an illegal value.
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT int ARKodeInit(void *arkode_mem, ARKRhsFn fe, 
			       ARKRhsFn fi, ARKExpStabFn EStab, 
			       realtype t0, N_Vector y0);

/* -----------------------------------------------------------------
 * Function : ARKodeReInit
 * -----------------------------------------------------------------
 * ARKodeReInit re-initializes ARKode for the solution of a problem,
 * where a prior call to ARKodeInit has been made with the same
 * problem size N. ARKodeReInit performs the same input checking
 * and initializations that ARKodeInit does, but it does no memory
 * allocation, assuming that the existing internal memory is 
 * sufficient for the new problem.
 *
 * The use of ARKodeReInit requires that the method order, ord, is 
 * no larger for the new problem than for the problem specified in
 * the last call to ARKodeInit.  This condition is automatically 
 * fulfilled if the default value for ord is specified.
 *
 * All of the arguments to ARKodeReInit have names and meanings
 * identical to those of ARKodeInit.
 *
 * The return value of ARKodeReInit is equal to ARK_SUCCESS = 0 if
 * there were no errors; otherwise it is a negative int equal to:
 *   ARK_MEM_NULL     indicating arkode_mem was NULL (i.e.,
 *                    ARKodeCreate has not been called).
 *   ARK_NO_MALLOC    indicating that arkode_mem has not been
 *                    allocated (i.e., ARKodeInit has not been
 *                    called).
 *   ARK_ILL_INPUT    indicating an input argument was illegal
 *                    (including an attempt to increase maxord).
 * In case of an error return, an error message is also printed.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int ARKodeReInit(void *arkode_mem, realtype t0, N_Vector y0);

/* -----------------------------------------------------------------
 * Functions : ARKodeSStolerances
 *             ARKodeSVtolerances
 *             ARKodeWFtolerances
 * -----------------------------------------------------------------
 *
 * These functions specify the integration tolerances. One of them
 * SHOULD be called before the first call to ARKode; otherwise 
 * default values of reltol=1e-4 and abstol=1e-9 will be used, which 
 * may be entirely incorrect for a specific problem.
 *
 * ARKodeSStolerances specifies scalar relative and absolute 
 *   tolerances.
 * ARKodeSVtolerances specifies scalar relative tolerance and a 
 *   vector absolute tolerance (a potentially different absolute 
 *   tolerance for each vector component).
 * ARKodeWFtolerances specifies a user-provided function (of type 
 *   ARKEwtFn) which will be called to set the error weight vector.
 *
 * The tolerances reltol and abstol define a vector of error weights,
 * ewt, with components
 *   ewt[i] = 1/(reltol*abs(y[i]) + abstol)      (in the SS case), or
 *   ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (in the SV case).
 * This vector is used in all error and convergence tests, which
 * use a weighted RMS norm on all error-like vectors v:
 *    WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),
 * where N is the problem dimension.
 *
 * The return value of these functions is equal to ARK_SUCCESS = 0 if
 * there were no errors; otherwise it is a negative int equal to:
 *   ARK_MEM_NULL     indicating arkode_mem was NULL (i.e.,
 *                    ARKodeCreate has not been called).
 *   ARK_NO_MALLOC    indicating that arkode_mem has not been
 *                    allocated (i.e., ARKodeInit has not been
 *                    called).
 *   ARK_ILL_INPUT    indicating an input argument was illegal
 *                    (e.g. a negative tolerance)
 * In case of an error return, an error message is also printed.
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT int ARKodeSStolerances(void *arkode_mem, realtype reltol, realtype abstol);
SUNDIALS_EXPORT int ARKodeSVtolerances(void *arkode_mem, realtype reltol, N_Vector abstol);
SUNDIALS_EXPORT int ARKodeWFtolerances(void *arkode_mem, ARKEwtFn efun);

/* -----------------------------------------------------------------
 * Function : ARKodeRootInit
 * -----------------------------------------------------------------
 * ARKodeRootInit initializes a rootfinding problem to be solved
 * during the integration of the ODE system.  It must be called
 * after ARKodeCreate, and before ARKode.  The arguments are:
 *
 * arkode_mem = pointer to ARKODE memory returned by ARKodeCreate.
 *
 * nrtfn      = number of functions g_i, an integer >= 0.
 *
 * g          = name of user-supplied function, of type ARKRootFn,
 *              defining the functions g_i whose roots are sought.
 *
 * If a new problem is to be solved with a call to ARKodeReInit,
 * where the new problem has no root functions but the prior one
 * did, then call ARKodeRootInit with nrtfn = 0.
 *
 * The return value of ARKodeRootInit is ARK_SUCCESS = 0 if there were
 * no errors; otherwise it is a negative int equal to:
 *   ARK_MEM_NULL    indicating arkode_mem was NULL, or
 *   ARK_MEM_FAIL    indicating a memory allocation failed.
 *                   (including an attempt to increase maxord).
 *   ARK_ILL_INPUT   indicating nrtfn > 0 but g = NULL.
 * In case of an error return, an error message is also printed.
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT int ARKodeRootInit(void *arkode_mem, int nrtfn, ARKRootFn g);

/* -----------------------------------------------------------------
 * Function : ARKode
 * -----------------------------------------------------------------
 * ARKode integrates the ODE over an interval in t.
 *
 * ARKode may be run in one of three modes (ARK_NORMAL, ARK_ONE_STEP, 
 * ARK_FIXED_STEP), as determined by the itask argument:
 *
 * If itask is ARK_NORMAL, then the solver integrates from its
 * current internal t value to a point at or beyond tout, then
 * interpolates to t = tout and returns y(tout) in the user-
 * allocated vector yout.  This interpolation is typically less 
 * accurate than the full time step solutions produced by the 
 * solver, since the interpolating polynomial relies on the 
 * internal stage solutions, that may have reduced accuracy in 
 * comparison with the full time step solutions.  If the user 
 * wishes that this returned value have full method accuracy, they 
 * may issue a call to ARKodeSetStopTime before the call to ARKode 
 * to specify a fixed stop time to end the time step and return to 
 * the user.  Once the integrator returns at a tstop time, any 
 * future testing for tstop is disabled (and can be reenabled only 
 * though a new call to ARKodeSetStopTime).
 *
 * If itask is ARK_ONE_STEP, then the solver takes one internal 
 * time step and returns in yout the value of y at the new internal 
 * time. In this case, tout is used only during the first call to 
 * ARKode to determine the direction of integration and the rough 
 * scale of the t variable.  As with the ARK_NORMAL mode, a user 
 * may specify a specific stop time for output of this step, 
 * assuming that the requested step is smaller than the step taken 
 * by the method. 
 *
 * If itask is ARK_FIXED_STEP, then the solver will take a single 
 * step from the current time (either the t0 value specified in 
 * ARKodeInit or ARKodeReInit, or the time at which the solver 
 * last returned a value) to the time requested in tout.  In this
 * mode, no attempt is made at satisfying temporal error 
 * requirements as specified through the reltol and abstol 
 * tolerances.  It is strongly suggested that before using this 
 * mode, the user increase the allowable work limits for the 
 * nonlinear and linear solvers.
 *
 * The time reached by the solver is placed in (*tret). The
 * user is responsible for allocating the memory for this value.
 *
 * arkode_mem is the pointer to ARKODE memory returned by
 *            ARKodeCreate.
 *
 * tout  is the next time at which a computed solution is desired.
 *
 * yout  is the computed solution vector. In ARK_NORMAL mode with no
 *       errors and no roots found, yout=y(tout).
 *
 * tret  is a pointer to a real location. ARKode sets (*tret) to
 *       the time reached by the solver and returns
 *       yout=y(*tret).
 *
 * itask is ARK_NORMAL, ARK_ONE_STEP or ARK_FIXED_STEP. These three 
 * modes are described above.
 *
 * Here is a brief description of each return value:
 *
 * ARK_SUCCESS:      ARKode succeeded and no roots were found.
 *
 * ARK_ROOT_RETURN:  ARKode succeeded, and found one or more roots.
 *                   If nrtfn > 1, call ARKodeGetRootInfo to see
 *                   which g_i were found to have a root at (*tret).
 *
 * ARK_TSTOP_RETURN: ARKode succeeded and returned at tstop.
 *
 * ARK_MEM_NULL:     The arkode_mem argument was NULL.
 *
 * ARK_NO_MALLOC:    arkode_mem was not allocated.
 *
 * ARK_ILL_INPUT:    One of the inputs to ARKode is illegal. This
 *                   includes the situation when a component of the
 *                   error weight vectors becomes < 0 during
 *                   internal time-stepping.  It also includes the
 *                   situation where a root of one of the root
 *                   functions was found both at t0 and very near t0.
 *                   The ILL_INPUT flag will also be returned if the
 *                   linear solver routine ARK--- (called by the user
 *                   after calling ARKodeCreate) failed to set one of
 *                   the linear solver-related fields in arkode_mem or
 *                   if the linear solver's init routine failed. In
 *                   any case, the user should see the printed
 *                   error message for more details.
 *
 * ARK_TOO_MUCH_WORK: The solver took mxstep internal steps but
 *                   could not reach tout. The default value for
 *                   mxstep is MXSTEP_DEFAULT = 500.
 *
 * ARK_TOO_MUCH_ACC: The solver could not satisfy the accuracy
 *                   demanded by the user for some internal step.
 *
 * ARK_ERR_FAILURE:  Error test failures occurred too many times
 *                   (= MXNEF = 7) during one internal time step or
 *                   occurred with |h| = hmin.
 *
 * ARK_CONV_FAILURE: Convergence test failures occurred too many
 *                   times (= MXNCF = 10) during one internal time
 *                   step or occurred with |h| = hmin.
 *
 * ARK_LINIT_FAIL:   The linear solver's initialization function 
 *                   failed.
 *
 * ARK_LSETUP_FAIL:  The linear solver's setup routine failed in an
 *                   unrecoverable manner.
 *
 * ARK_LSOLVE_FAIL:  The linear solver's solve routine failed in an
 *                   unrecoverable manner.
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT int ARKode(void *arkode_mem, realtype tout, 
			   N_Vector yout, realtype *tret, int itask);

/* -----------------------------------------------------------------
 * Function : ARKodeGetDky
 * -----------------------------------------------------------------
 * ARKodeGetDky computes the kth derivative of the y function at
 * time t, where tn-hu <= t <= tn, tn denotes the current
 * internal time reached, and hu is the last internal step size
 * successfully used by the solver. The user may request
 * k=0, 1, ..., s-1, where s is the number of stages taken by the 
 * RK method. The derivative vector is returned in dky. This vector 
 * must be allocated by the caller. It is only legal to call this
 * function after a successful return from ARKode.
 *
 * arkode_mem is the pointer to ARKODE memory returned by
 *            ARKodeCreate.
 *
 * t   is the time at which the kth derivative of y is evaluated.
 *     The legal range for t is [tn-hu,tn] as described above.
 *
 * k   is the order of the derivative of y to be computed. The
 *     legal range for k is [0,s-1] as described above.
 *
 * dky is the output derivative vector [((d/dy)^k)y](t).
 *
 * The return value for ARKodeGetDky is one of:
 *
 *   ARK_SUCCESS:  ARKodeGetDky succeeded.
 *
 *   ARK_BAD_K:    k is not in the range 0, 1, ..., s-1.
 *
 *   ARK_BAD_T:    t is not in the interval [tn-hu,tn].
 *
 *   ARK_BAD_DKY:  The dky argument was NULL.
 *
 *   ARK_MEM_NULL: The arkode_mem argument was NULL.
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT int ARKodeGetDky(void *arkode_mem, realtype t, int k, N_Vector dky);

/* -----------------------------------------------------------------
 * Integrator optional output extraction functions
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistics related to the main integrator.
 * -----------------------------------------------------------------
 * ARKodeGetWorkSpace returns the ARKODE real and integer workspaces
 * ARKodeGetNumSteps returns the cumulative number of internal
 *                   steps taken by the solver
 * ARKodeGetNumExpSteps returns the cumulative number of stability 
 *                      limited steps taken by the solver
 * ARKodeGetNumAccSteps returns the cumulative number of accuracy 
 *                      limited steps taken by the solver
 * ARKodeGetNumConvSteps returns the cumulative number of convergence 
 *                       limited steps taken by the solver
 * ARKodeGetNumRhsEvals returns the number of calls to the user's
 *                      f function
 * ARKodeGetNumLinSolvSetups returns the number of calls made to
 *                           the linear solver's setup routine
 * ARKodeGetNumErrTestFails returns the number of local error test
 *                          failures that have occured
 * ARKodeGetActualInitStep returns the actual initial step size
 *                         used by ARKODE
 * ARKodeGetLastStep returns the step size for the last internal
 *                   step
 * ARKodeGetCurrentStep returns the step size to be attempted on
 *                      the next internal step
 * ARKodeGetCurrentTime returns the current internal time reached
 *                      by the solver
 * ARKodeGetTolScaleFactor returns a suggested factor by which the
 *                         user's tolerances should be scaled when
 *                         too much accuracy has been requested for
 *                         some internal step
 * ARKodeGetErrWeights returns the current error weight vector.
 *                     The user must allocate space for eweight.
 * ARKodeGetEstLocalErrors returns the vector of estimated local
 *                         errors. The user must allocate space
 *                         for ele.
 * ARKodeGetNumGEvals returns the number of calls to the user's
 *                    g function (for rootfinding)
 * ARKodeGetRootInfo returns the indices for which g_i was found to 
 *                   have a root. The user must allocate space for 
 *                   rootsfound. For i = 0 ... nrtfn-1, 
 *                   rootsfound[i] = 1 if g_i has a root, and = 0 if not.
 *
 * ARKodeGet* return values:
 *   ARK_SUCCESS   if succesful
 *   ARK_MEM_NULL  if the arkode memory was NULL
 *   ARK_NO_SLDET  if stability limit was not turned on
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int ARKodeGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw);
SUNDIALS_EXPORT int ARKodeGetNumSteps(void *arkode_mem, long int *nsteps);
SUNDIALS_EXPORT int ARKodeGetNumExpSteps(void *arkode_mem, long int *expsteps);
SUNDIALS_EXPORT int ARKodeGetNumAccSteps(void *arkode_mem, long int *accsteps);
SUNDIALS_EXPORT int ARKodeGetNumConvSteps(void *arkode_mem, long int *convsteps);
SUNDIALS_EXPORT int ARKodeGetNumRhsEvals(void *arkode_mem, long int *nfevals);
SUNDIALS_EXPORT int ARKodeGetNumLinSolvSetups(void *arkode_mem, long int *nlinsetups);
SUNDIALS_EXPORT int ARKodeGetNumErrTestFails(void *arkode_mem, long int *netfails);
SUNDIALS_EXPORT int ARKodeGetActualInitStep(void *arkode_mem, realtype *hinused);
SUNDIALS_EXPORT int ARKodeGetLastStep(void *arkode_mem, realtype *hlast);
SUNDIALS_EXPORT int ARKodeGetCurrentStep(void *arkode_mem, realtype *hcur);
SUNDIALS_EXPORT int ARKodeGetCurrentTime(void *arkode_mem, realtype *tcur);
SUNDIALS_EXPORT int ARKodeGetTolScaleFactor(void *arkode_mem, realtype *tolsfac);
SUNDIALS_EXPORT int ARKodeGetErrWeights(void *arkode_mem, N_Vector eweight);
SUNDIALS_EXPORT int ARKodeGetEstLocalErrors(void *arkode_mem, N_Vector ele);
SUNDIALS_EXPORT int ARKodeGetNumGEvals(void *arkode_mem, long int *ngevals);
SUNDIALS_EXPORT int ARKodeGetRootInfo(void *arkode_mem, int *rootsfound);

/* -----------------------------------------------------------------
 * As a convenience, the following functions provides the
 * optional outputs in one group.
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT int ARKodeGetIntegratorStats(void *arkode_mem, long int *nsteps,
					     long int *expsteps, long int *accsteps, 
					     long int *convsteps, long int *nfevals, 
					     long int *nlinsetups, long int *netfails,
					     realtype *hinused, realtype *hlast, 
					     realtype *hcur, realtype *tcur);

/* -----------------------------------------------------------------
 * Nonlinear solver optional output extraction functions
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistics related to the nonlinear solver.
 * -----------------------------------------------------------------
 * ARKodeGetNumNonlinSolvIters returns the number of nonlinear
 *                             solver iterations performed.
 * ARKodeGetNumNonlinSolvConvFails returns the number of nonlinear
 *                                 convergence failures.
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT int ARKodeGetNumNonlinSolvIters(void *arkode_mem, long int *nniters);
SUNDIALS_EXPORT int ARKodeGetNumNonlinSolvConvFails(void *arkode_mem, long int *nncfails);

/* -----------------------------------------------------------------
 * As a convenience, the following function provides the
 * nonlinear solver optional outputs in a group.
 * ----------------------------------------------------------------- */
  
SUNDIALS_EXPORT int ARKodeGetNonlinSolvStats(void *arkode_mem, long int *nniters,
					     long int *nncfails);

/* -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a ARKODE return flag
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT char *ARKodeGetReturnFlagName(long int flag);

/* -----------------------------------------------------------------
 * Function : ARKodeFree
 * -----------------------------------------------------------------
 * ARKodeFree frees the problem memory arkode_mem allocated by
 * ARKodeCreate and ARKodeInit. Its only argument is the pointer
 * arkode_mem returned by ARKodeCreate.
 * ----------------------------------------------------------------- */

SUNDIALS_EXPORT void ARKodeFree(void **arkode_mem);

#ifdef __cplusplus
}
#endif

#endif
