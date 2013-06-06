:tocdepth: 3

.. _Examples:

ARKode Examples
===============

ARKode comes packaged with a variety of example problems, that
exercise options including explicit, implicit and ImEx solvers,
root-finding, direct and iterative linear solvers, and the Fortran
solver interface, FARKODE.  While these examples are not an exhaustive
set of all possible usage scenarios, they are designed to show a
variety of usage scenarios, and can be used as templates for new
problems using ARKode's solvers.


**Add table differentiating between examples, so that users can
quickly ascertain which one to study/emulate**



Simple linear example (ark_analytic)
-------------------------------------

This is a very simple C example that merely shows how to use the
ARKode solver interface.

ODE system
^^^^^^^^^^^^

The problem is that of a scalar-valued initial value problem (IVP)
that is linear in the dependent variable :math:`y`, but nonlinear in
the independent variable :math:`t`:

.. math::

   \frac{dy}{dt} = \lambda y + \frac{1}{1+t^2} - \lambda \arctan(t),

where :math:`0\le t\le 10` and :math:`y(0)=0`.  The stiffness of the
problem may be tuned via the parameter :math:`\lambda`, which is
specified (along with the relative and absolute tolerances,
:math:`rtol` and :math:`atol`) in the input file
``input_analytic.txt``.  The value of :math:`\lambda` must be negative
to result in a well-posed problem; for values with magnitude larger
than 100 or so the problem becomes quite stiff.  In the provided input
file, we choose :math:`\lambda=-100` and tolerances
:math:`rtol=10^{-6}` and :math:`atol=10^{-10}`.    After each unit
time interval, the solution is output to the screen.


Numerical method
^^^^^^^^^^^^^^^^^

The example routine solves this problem using a diagonally-implicit
Runge-Kutta method.  Each stage is solved using the built-in modified
Newton iteration, but since the ODE is linear in :math:`y` these
should only require a single iteration per stage.  Internally, Newton
will use the ARKDENSE dense linear solver, which in the case of this
scalar-valued problem is just division.  The example file contains
functions to evaluate both :math:`f(t,y)` and :math:`J(t,y)=\lambda`.

Aside from the input tolerance values, this problem uses only the
default parameters for the ARKode solver.


Routines
^^^^^^^^^^^^^^^^^

We reproduce the relevant aspects of the ``main()`` routine and
auxiliary functions here for explanatory purposes (see the in-line
comments for details; error-checking has been removed for brevity):

.. code-block:: c

   /* Header files */
   #include <stdio.h>
   #include <math.h>
   #include <arkode/arkode.h>           /* prototypes for ARKODE fcts., consts. */
   #include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
   #include <arkode/arkode_dense.h>     /* prototype for ARKDense solver */
   #include <sundials/sundials_dense.h> /* definitions of DlsMat and DENSE_ELEM */
   #include <sundials/sundials_types.h> /* definition of type 'realtype' */
   
   /* User-supplied Functions Called by the Solver */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


   /* Main Program */
   int main()
   {
     /* general problem parameters */
     realtype T0 = RCONST(0.0);      /* initial time */
     realtype Tf = RCONST(10.0);     /* final time */
     realtype dTout = RCONST(1.0);   /* time between outputs */
     long int NEQ = 1;               /* number of dependent variables */
   
     /* general problem variables */
     int flag;                       /* reusable error-checking flag */
     N_Vector y = NULL;              /* empty vector for storing solution */
     void *arkode_mem = NULL;        /* empty ARKode memory structure */
   
     /* read problem parameter and tolerances from input file:
        lamda  - problem stiffness parameter
        reltol - desired relative tolerance
        abstol - desired absolute tolerance */
     double reltol_, abstol_, lamda_;
     FILE *FID;
     FID=fopen("input_analytic.txt","r");
     flag = fscanf(FID,"  lamda = %lf\n",  &lamda_);
     flag = fscanf(FID,"  reltol = %lf\n", &reltol_);
     flag = fscanf(FID,"  abstol = %lf\n", &abstol_);
     fclose(FID);
   
     /* convert the inputs to 'realtype' format */
     realtype reltol = reltol_;
     realtype abstol = abstol_;
     realtype lamda  = lamda_;
   
     /* Initial diagnostics output */
     printf("\nAnalytical ODE test problem:\n");
     printf("    lamda = %g\n",    lamda);
     printf("   reltol = %.1e\n",  reltol);
     printf("   abstol = %.1e\n\n",abstol);
   
     /* Create serial vector of length NEQ for solution */
     y = N_VNew_Serial(NEQ);
   
     /* Initialize y to 0 (to specify initial condition) */
     NV_Ith_S(y,0) = 0.0;
   
     /* Call ARKodeCreate to create the solver memory */
     arkode_mem = ARKodeCreate();
     
     /* Call ARKodeInit to initialize the integrator memory and specify the
        hand-side side function in y'=f(t,y), the inital time T0, and
        the initial dependent variable vector y.  Note: since this
	problem is fully implicit, we set f_E to NULL and f_I to f. */
     flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
   
     /* Call ARKodeSetUserData to pass lamda to user functions */
     flag = ARKodeSetUserData(arkode_mem, (void *) &lamda);
   
     /* Call ARKodeSStolerances to specify the scalar relative and absolute
        tolerances */
     flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
   
     /* Call ARKDense to specify the ARKDENSE dense linear solver */
     flag = ARKDense(arkode_mem, NEQ);
   
     /* Set the Jacobian routine to 'Jac' (user-supplied) */
     flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);
   
     /* Main time-stepping loop: calls ARKode to perform the
        integration, then prints results.  Stops when the final time
	has been reached */
     realtype t = T0;
     realtype tout = T0+dTout;
     realtype u;
     printf("        t           u\n");
     printf("   ---------------------\n");
     while (Tf - t > 1.0e-15) {
   
       /* call integrator */
       flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);

       /* access/print current solution */
       u = NV_Ith_S(y,0);
       printf("  %10.6f  %10.6f\n", t, u);

       /* check for successful solve: update time or break */
       if (flag >= 0) {
         tout += dTout;
         tout = (tout > Tf) ? Tf : tout;
       } else {
         fprintf(stderr,"Solver failure, stopping integration\n");
         break;
       }

     }
     printf("   ---------------------\n");
   
     /* Get/print some final statistics on how the solve progressed */
     long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
     flag = ARKodeGetNumSteps(arkode_mem, &nst);
     flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
     flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
     flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
     flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
     flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
     flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
     flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
     flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
   
     printf("\nFinal Solver Statistics:\n");
     printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
     printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
     printf("   Total linear solver setups = %li\n", nsetups);
     printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
     printf("   Total number of Jacobian evaluations = %li\n", nje);
     printf("   Total number of Newton iterations = %li\n", nni);
     printf("   Total number of linear solver convergence failures = %li\n", ncfn);
     printf("   Total number of error test failures = %li\n\n", netf);
   
     /* Free y vector */
     N_VDestroy_Serial(y);
   
     /* Free integrator memory */
     ARKodeFree(&arkode_mem);

     /* Return with successful completion */
     return 0;
   }
  
   
   /*-------------------------------
    * Functions called by the solver
    *-------------------------------*/
   
   /* f routine to compute the ODE RHS function f(t,y). */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
   {
     realtype *rdata = (realtype *) user_data;   /* cast user_data to realtype */
     realtype lamda = rdata[0];                  /* set shortcut for stiffness parameter */
     realtype u = NV_Ith_S(y,0);                 /* access current solution value */
   
     /* fill in the RHS function: here "NV_Ith_S" access the 0th entry
        of the vector ydot */
     NV_Ith_S(ydot,0) = lamda*u + 1.0/(1.0+t*t) - lamda*atan(t);

     /* return with success flag */
     return 0;
   }
   
   /* Jacobian routine to compute J(t,y) = df/dy. */
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
   {
     realtype *rdata = (realtype *) user_data;   /* cast user_data to realtype */
     realtype lamda = rdata[0];                  /* set shortcut for stiffness parameter */

     /* Fill in Jacobian of f: here "DENSE_ELEM" accesses the (0,0)
        entry of the Jacobian matrix J */
     DENSE_ELEM(J,0,0) = lamda;
   
     /* return with success flag */
     return 0;
   }


Solutions
^^^^^^^^^^^^

This problem is included both as a simple example, but also because it
has an analytical solution, :math:`y(t) = \arctan(t)`.  As seen in the
plots below, the computed solution tracks the analytical solution
quite well, and results in errors below those specified by the input
error tolerances.



Simple nonlinear example (ark_analytic_nonlin)
-----------------------------------------------

ODE system
^^^^^^^^^^^^

.. math::

   \frac{dy}{dt} = (t+1) e^{-y},

for the interval :math:`t \in [0.0, 10.0]`, with initial condition
:math:`y(0)=0`.  This has analytical solution :math:`y(t) =
\log\left(\frac{t^2}{2} + t + 1\right)`.  



Numerical method
^^^^^^^^^^^^^^^^^

This program solves the problem with the DIRK method,
Newton iteration with the ARKDENSE dense linear solver, and a
user-supplied Jacobian routine.
Output is printed every 1.0 units of time (10 total).
Run statistics (optional outputs) are printed at the end.


Routines
^^^^^^^^^^^^^^^^^

We reproduce the relevant aspects of the ``main()`` routine and
auxiliary functions here for explanatory purposes (see the in-line
comments for details; error-checking has been removed for brevity):

.. code-block:: c

   /* Header files */
   #include <stdio.h>
   #include <math.h>
   #include <arkode/arkode.h>
   #include <nvector/nvector_serial.h>
   #include <arkode/arkode_dense.h>
   #include <sundials/sundials_dense.h>
   #include <sundials/sundials_types.h>
   
   /* User-supplied Functions Called by the Solver */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
   

   /* Main Program */
   int main()
   {
     /* general problem parameters */
     realtype T0 = RCONST(0.0);
     realtype Tf = RCONST(10.0);
     realtype dTout = RCONST(1.0);
     long int NEQ = 1;
   
     /* general problem variables */
     int flag;
     N_Vector y = NULL;
     void *arkode_mem = NULL;
   
     /* read problem parameter and tolerances from input file:
        lamda  - problem stiffness parameter
        reltol - desired relative tolerance
        abstol - desired absolute tolerance */
     double reltol_, abstol_;
     FILE *FID;
     FID=fopen("input_analytic_nonlin.txt","r");
     fscanf(FID,"  reltol = %lf\n", &reltol_);
     fscanf(FID,"  abstol = %lf\n", &abstol_);
     fclose(FID);
   
     /* convert the inputs to 'realtype' format */
     realtype reltol = reltol_;
     realtype abstol = abstol_;
   
     /* open solver diagnostics output file for writing */
     FILE *DFID;
     DFID=fopen("diags_ark_analytic_nonlin.txt","w");
     
     /* Initial problem output */
     printf("\nAnalytical ODE test problem:\n");
     printf("   reltol = %.1e\n",  reltol);
     printf("   abstol = %.1e\n\n",abstol);
   
   
     /* Create serial vector of length NEQ for initial condition */
     y = N_VNew_Serial(NEQ);
   
     /* Initialize y to 0 */
     NV_Ith_S(y,0) = 0.0;
   
     /* Call ARKodeCreate to create the solver memory */
     arkode_mem = ARKodeCreate();
     
     /* Call ARKodeInit to initialize the integrator memory and specify the
        user's right hand side function in y'=f(t,y), the inital time T0, and
        the initial dependent variable vector y */
     flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
   
     /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
     flag = ARKodeSetDiagnostics(arkode_mem, DFID);
   
     /* Call ARKodeSStolerances to specify the scalar relative and absolute
        tolerances */
     flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
   
     /* Call ARKDense to specify the ARKDENSE dense linear solver */
     flag = ARKDense(arkode_mem, NEQ);
   
     /* Set the Jacobian routine to Jac (user-supplied) */
     flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);
   
     /* In loop, call ARKode, print results, and test for error.
        Break out of loop when the final output time has been reached */
     realtype t = T0;
     realtype tout = T0+dTout;
     realtype u;
     printf("        t           u\n");
     printf("   ---------------------\n");
     while (Tf - t > 1.0e-15) {
   
       flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
       u = NV_Ith_S(y,0);
       printf("  %10.6f  %10.6f\n", t, u);
   
       if (flag >= 0) {
         tout += dTout;
         tout = (tout > Tf) ? Tf : tout;
       }
     }
     printf("   ---------------------\n");
   
     /* Print some final statistics */
     long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
     flag = ARKodeGetNumSteps(arkode_mem, &nst);
     flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
     flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
     flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
     flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
     flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
     flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
     flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
     flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
   
     printf("\nFinal Solver Statistics:\n");
     printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
     printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
     printf("   Total linear solver setups = %li\n", nsetups);
     printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
     printf("   Total number of Jacobian evaluations = %li\n", nje);
     printf("   Total number of Newton iterations = %li\n", nni);
     printf("   Total number of linear solver convergence failures = %li\n", ncfn);
     printf("   Total number of error test failures = %li\n\n", netf);
   
     /* Free y vector */
     N_VDestroy_Serial(y);
   
     /* Free integrator memory */
     ARKodeFree(&arkode_mem);
   
     /* close solver diagnostics output file */
     fclose(DFID);
   
     return 0;
   }
   
   
   /*-------------------------------
    * Functions called by the solver
    *-------------------------------*/
   
   /* f routine to compute the ODE RHS function f(t,y). */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
   {
     NV_Ith_S(ydot,0) = (t+1.0)*exp(-NV_Ith_S(y,0));
     return 0;
   }
   
   /* Jacobian routine to compute J(t,y) = df/dy. */
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
   {
     DENSE_ELEM(J,0,0) = -(t+1.0)*exp(-NV_Ith_S(y,0));
     return 0;
   }


Solutions
^^^^^^^^^^^^



Simple linear system example (ark_analytic_sys)
------------------------------------------------

ODE system
^^^^^^^^^^^^

.. math::

   \frac{dy}{dt} = Ay

where :math:`A = V D V^{-1}`.  Here, we use

.. math::

   V = \left[\begin{array}{rrr} 1 & -1 & 1\\ -1 & 2 & 1\\ 0 * -1 *
       2\end{array}\right], \qquad
   V^{-1} = \frac14 \left[\begin{array}{rrr} 5 & 1 & -3\\ 2 & 2 & -2\\
       1 & 1 & 1 \end{array}\right], \qquad
   D = \left[\begin{array}{rrr} -1/2 & 0 & 0\\ 0 & -1/10 & 0\\ 0 & 0 &
       \lambda \end{array}\right].

where :math:`\lambda` is a large negative number. The analytical
solution to this problem is 

.. math::

   Y(t) = V e^{Dt} V^{-1} Y(0).

We evolve the problem for :math:`t` in the interval :math:`\left[0,\,
\frac{1}{20}\right]`, with initial condition :math:`Y(0) = \left[1,\,
1,\, 1\right]^T`.


Numerical method
^^^^^^^^^^^^^^^^^

The stiffness of the problem is directly proportional to the 
value of :math:`\lambda`, which is specified through an input file,
along with the desired relative and absolute tolerances.  The value of
:math:`\lambda` should be negative to result in a well-posed ODE; for
values with magnitude larger than 100 the problem becomes quite stiff.

In the example input file, we choose :math:`\lambda = -100`.
 
This program solves the problem with the DIRK method,
Newton iteration with the ARKDENSE dense linear solver, and a
user-supplied Jacobian routine.
Output is printed every 0.005 units of time (10 total).
Run statistics (optional outputs) are printed at the end.


Routines
^^^^^^^^^^^^^^^^^

We reproduce the relevant aspects of the ``main()`` routine and
auxiliary functions here for explanatory purposes (see the in-line
comments for details; error-checking has been removed for brevity):

.. code-block:: c

   /* Header files */
   #include <stdio.h>
   #include <math.h>
   #include <arkode/arkode.h>
   #include <nvector/nvector_serial.h>
   #include <arkode/arkode_dense.h>
   #include <sundials/sundials_dense.h>
   #include <sundials/sundials_types.h>
   
   /* User-supplied Functions Called by the Solver */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
   
   /* Private function to perform matrix-matrix product */
   static int dense_MM(DlsMat A, DlsMat B, DlsMat C);
   
   
   /* Main Program */
   int main()
   {
     /* general problem parameters */
     realtype T0 = RCONST(0.0);
     realtype Tf = RCONST(0.05);
     realtype dTout = RCONST(0.005);
     long int NEQ = 3;
   
     /* general problem variables */
     int flag;
     N_Vector y = NULL;
     void *arkode_mem = NULL;
   
     /* read problem parameter and tolerances from input file:
        lamda  - problem stiffness parameter
        reltol - desired relative tolerance
        abstol - desired absolute tolerance */
     double reltol_, abstol_, lamda_;
     FILE *FID;
     FID=fopen("input_analytic_sys.txt","r");
     fscanf(FID,"  lamda = %lf\n",  &lamda_);
     fscanf(FID,"  reltol = %lf\n", &reltol_);
     fscanf(FID,"  abstol = %lf\n", &abstol_);
     fclose(FID);
   
     /* convert the inputs to 'realtype' format */
     realtype reltol = reltol_;
     realtype abstol = abstol_;
     realtype lamda  = lamda_;
   
     /* open solver diagnostics output file for writing */
     FILE *DFID;
     DFID=fopen("diags_ark_analytic_sys.txt","w");
     
     /* Initial problem output */
     printf("\nAnalytical ODE test problem:\n");
     printf("    lamda = %g\n",    lamda);
     printf("   reltol = %.1e\n",  reltol);
     printf("   abstol = %.1e\n\n",abstol);
   
   
     /* Create serial vector of length NEQ for initial condition */
     y = N_VNew_Serial(NEQ);
   
     /* Initialize y to 0 */
     NV_Ith_S(y,0) = 1.0;
     NV_Ith_S(y,1) = 1.0;
     NV_Ith_S(y,2) = 1.0;
   
     /* Call ARKodeCreate to create the solver memory and specify the 
        Backward Differentiation Formula and the use of a Newton iteration */
     arkode_mem = ARKodeCreate();
     
     /* Call ARKodeInit to initialize the integrator memory and specify the
        user's right hand side function in y'=f(t,y), the inital time T0, and
        the initial dependent variable vector y */
     flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
   
     /* Call ARKodeSetUserData to pass lamda to user functions */
     flag = ARKodeSetUserData(arkode_mem, (void *) &lamda);
   
     /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
     flag = ARKodeSetDiagnostics(arkode_mem, DFID);
   
     /* Call ARKodeSStolerances to specify the scalar relative and absolute
        tolerances */
     flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
   
     /* Call ARKDense to specify the ARKDENSE dense linear solver */
     flag = ARKDense(arkode_mem, NEQ);
   
     /* Set the Jacobian routine to Jac (user-supplied) */
     flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);
   
     /* In loop, call ARKode, print results, and test for error.
        Break out of loop when the final output time has been reached */
     realtype t = T0;
     realtype tout = T0+dTout;
     realtype y0, y1, y2;
     printf("      t        y0        y1        y2\n");
     printf("   --------------------------------------\n");
     while (Tf - t > 1.0e-15) {
   
       flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
       y0 = NV_Ith_S(y,0);
       y1 = NV_Ith_S(y,1);
       y2 = NV_Ith_S(y,2);
       printf("  %8.4f  %8.5f  %8.5f  %8.5f\n", t, y0, y1, y2);
   
       if (flag >= 0) {
         tout += dTout;
         tout = (tout > Tf) ? Tf : tout;
       }
     }
     printf("   --------------------------------------\n");
   
     /* Print some final statistics */
     long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
     flag = ARKodeGetNumSteps(arkode_mem, &nst);
     flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
     flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
     flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
     flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
     flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
     flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
     flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
     flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
   
     printf("\nFinal Solver Statistics:\n");
     printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
     printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
     printf("   Total linear solver setups = %li\n", nsetups);
     printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
     printf("   Total number of Jacobian evaluations = %li\n", nje);
     printf("   Total number of Newton iterations = %li\n", nni);
     printf("   Total number of linear solver convergence failures = %li\n", ncfn);
     printf("   Total number of error test failures = %li\n\n", netf);
   
     /* Free y vector */
     N_VDestroy_Serial(y);
   
     /* Free integrator memory */
     ARKodeFree(&arkode_mem);
   
     /* close solver diagnostics output file */
     fclose(DFID);
   
     return 0;
   }
   
   
   /*-------------------------------
    * Functions called by the solver
    *-------------------------------*/
   
   /* f routine to compute the ODE RHS function f(t,y). */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
   {
     realtype *rdata = (realtype *) user_data;
     realtype lam = rdata[0];
     realtype y0 = NV_Ith_S(y,0);
     realtype y1 = NV_Ith_S(y,1);
     realtype y2 = NV_Ith_S(y,2);
     realtype yd0, yd1, yd2;
     
     /* f(t,y) = V*D*Vi*y, where 
           V = [1 -1 1; -1 2 1; 0 -1 2] 
           Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1]
           D = [-0.5 0 0; 0 -0.1 0; 0 0 lam] */
   
     /*   yd = Vi*y */
     yd0 = 0.25*(5.0*y0 + 1.0*y1 - 3.0*y2);
     yd1 = 0.25*(2.0*y0 + 2.0*y1 - 2.0*y2);
     yd2 = 0.25*(1.0*y0 + 1.0*y1 + 1.0*y2);
   
     /*   y = D*yd */
     y0  = -0.5*yd0;
     y1  = -0.1*yd1;
     y2  =  lam*yd2;
   
     /*   yd = V*y */
     yd0 =  1.0*y0 - 1.0*y1 + 1.0*y2;
     yd1 = -1.0*y0 + 2.0*y1 + 1.0*y2;
     yd2 =  0.0*y0 - 1.0*y1 + 2.0*y2;
   
     NV_Ith_S(ydot,0) = yd0;
     NV_Ith_S(ydot,1) = yd1;
     NV_Ith_S(ydot,2) = yd2;
   
     return 0;
   }
   
   /* Jacobian routine to compute J(t,y) = df/dy. */
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
   {
     realtype *rdata = (realtype *) user_data;
     realtype lam = rdata[0];
     DlsMat V  = NewDenseMat(3,3);
     DlsMat D  = NewDenseMat(3,3);
     DlsMat Vi = NewDenseMat(3,3);
   
     /* initialize temporary matrices to zero */
     DenseScale(0.0, V);
     DenseScale(0.0, D);
     DenseScale(0.0, Vi);
   
     /* J = V*D*Vi, where 
           V = [1 -1 1; -1 2 1; 0 -1 2] 
           Vi = 0.25*[5 1 -3; 2 2 -2; 1 1 1]
           D = [-0.5 0 0; 0 -0.1 0; 0 0 lam] */
     DENSE_ELEM(V,0,0) =  1.0;
     DENSE_ELEM(V,0,1) = -1.0;
     DENSE_ELEM(V,0,2) =  1.0;
     DENSE_ELEM(V,1,0) = -1.0;
     DENSE_ELEM(V,1,1) =  2.0;
     DENSE_ELEM(V,1,2) =  1.0;
     DENSE_ELEM(V,2,0) =  0.0;
     DENSE_ELEM(V,2,1) = -1.0;
     DENSE_ELEM(V,2,2) =  2.0;
   
     DENSE_ELEM(Vi,0,0) =  0.25*5.0;
     DENSE_ELEM(Vi,0,1) =  0.25*1.0;
     DENSE_ELEM(Vi,0,2) = -0.25*3.0;
     DENSE_ELEM(Vi,1,0) =  0.25*2.0;
     DENSE_ELEM(Vi,1,1) =  0.25*2.0;
     DENSE_ELEM(Vi,1,2) = -0.25*2.0;
     DENSE_ELEM(Vi,2,0) =  0.25*1.0;
     DENSE_ELEM(Vi,2,1) =  0.25*1.0;
     DENSE_ELEM(Vi,2,2) =  0.25*1.0;
   
     DENSE_ELEM(D,0,0) = -0.5;
     DENSE_ELEM(D,1,1) = -0.1;
     DENSE_ELEM(D,2,2) = lam;
   
     /* J = D*Vi */
     if (dense_MM(D,Vi,J) != 0) {
       fprintf(stderr, "matmul error\n");
       return 1;
     }
   
     /* D = V*J [= V*D*Vi] */
     if (dense_MM(V,J,D) != 0) {
       fprintf(stderr, "matmul error\n");
       return 1;
     }
   
     /* J = D [= V*D*Vi] */
     DenseCopy(D, J);
   
     return 0;
   }


Solutions
^^^^^^^^^^^^



Stiff nonlinear system example (ark_brusselator)
-------------------------------------------------

ODE system
^^^^^^^^^^^^

The following test simulates a brusselator problem from chemical 
kinetics.  This is an ODE system with 3 components, :math:`Y = [u,\,
v,\, w]`, satisfying the equations,

.. math::

   \frac{du}{dt} &= a - (w+1)u + v u^2, \\
   \frac{dv}{dt} &= w u - v u^2, \\
   \frac{dw}{dt} &= \frac{b-w}{\varepsilon} - w u.

We integrate over the interval :math:`0 \le t \le 10`, with the
initial conditions :math:`u(0) = u_0`, :math:`v(0) = v_0`, :math:`w(0) = w_0`.
After each unit time interval, the solution is output to the screen.

We have 3 different testing scenarios:

Test 1:  :math:`u_0=3.9`,  :math:`v_0=1.1`,  :math:`w_0=2.8`,
:math:`a=1.2`, :math:`b=2.5`, and :math:`\varepsilon=10^{-5}` 

Test 2:  :math:`u_0=1.2`, :math:`v_0=3.1`, :math:`w_0=3`, :math:`a=1`,
:math:`b=3.5`, and :math:`\varepsilon=5\cdot10^{-6}` 

Test 3:  :math:`u_0=3`, :math:`v_0=3`, :math:`w_0=3.5`, :math:`a=0.5`,
:math:`b=3`, and :math:`\varepsilon=5\cdot10^{-4}` 

These tests are selected within the input file (test = {1,2,3}), 
with the default set to test 2 in case the input is invalid.
Also in the input file, we allow specification of the desired 
relative and absolute tolerances.



Numerical method
^^^^^^^^^^^^^^^^^

This program solves the problem with the DIRK method, using a
Newton iteration with the ARKDENSE dense linear solver, and a
user-supplied Jacobian routine.

100 outputs are printed at equal intervals, and run statistics 
are printed at the end.


Routines
^^^^^^^^^^^^^^^^^

We reproduce the relevant aspects of the ``main()`` routine and
auxiliary functions here for explanatory purposes (see the in-line
comments for details; error-checking has been removed for brevity):

.. code-block:: c

   /* Header files */
   #include <stdio.h>
   #include <math.h>
   #include <arkode/arkode.h>
   #include <nvector/nvector_serial.h>
   #include <arkode/arkode_dense.h>
   #include <sundials/sundials_dense.h>
   #include <sundials/sundials_types.h>
   
   /* User-supplied Functions Called by the Solver */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
   
   
   /* Main Program */
   int main()
   {
     /* general problem parameters */
     realtype T0 = RCONST(0.0);
     realtype Tf = RCONST(10.0);
     realtype dTout = RCONST(1.0);
     int Nt = ceil(Tf/dTout);
     realtype a, b, ep, u0, v0, w0;
     long int NEQ = 3;
   
     /* general problem variables */
     int flag;
     N_Vector y = NULL;
     void *arkode_mem = NULL;
   
     /* read problem parameter and tolerances from input file:
        test   - test problem choice
        reltol - desired relative tolerance
        abstol - desired absolute tolerance */
     int test;
     double reltol_, abstol_;
     FILE *FID;
     FID=fopen("input_brusselator.txt","r");
     fscanf(FID,"  test = %i\n", &test);
     fscanf(FID,"  reltol = %lf\n", &reltol_);
     fscanf(FID,"  abstol = %lf\n", &abstol_);
     fclose(FID);
   
     /* convert the inputs to 'realtype' format */
     realtype reltol = reltol_;
     realtype abstol = abstol_;
   
     /* open solver diagnostics output file for writing */
     FILE *DFID;
     DFID=fopen("diags_ark_brusselator.txt","w");
     
     /* set up the test problem according to the desired input */
     if (test == 1) {
       u0 = RCONST(3.9);
       v0 = RCONST(1.1);
       w0 = RCONST(2.8);
       a  = RCONST(1.2);
       b  = RCONST(2.5);
       ep = RCONST(1.0e-5);
     } else if (test == 3) {
       u0 = RCONST(3.0);
       v0 = RCONST(3.0);
       w0 = RCONST(3.5);
       a  = RCONST(0.5);
       b  = RCONST(3.0);
       ep = RCONST(5.0e-4);
     } else {
       u0 = RCONST(1.2);
       v0 = RCONST(3.1);
       w0 = RCONST(3.0);
       a  = RCONST(1.0);
       b  = RCONST(3.5);
       ep = RCONST(5.0e-6);
     }
   
     /* set user data to contain problem-defining parameters */
     realtype rdata[3] = {a, b, ep};
   
     /* Initial problem output */
     printf("\nBrusselator ODE test problem:\n");
     printf("    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n",u0,v0,w0);
     printf("    problem parameters:  a = %g,  b = %g,  ep = %g\n",a,b,ep);
     printf("    reltol = %.1e,  abstol = %.1e\n\n",reltol,abstol);
   
   
     /* Create serial vector of length NEQ for initial condition */
     y = N_VNew_Serial(NEQ);
   
     /* Set initial conditions into y */
     NV_Ith_S(y,0) = u0;
     NV_Ith_S(y,1) = v0;
     NV_Ith_S(y,2) = w0;
   
     /* Call ARKodeCreate to create the solver memory */
     arkode_mem = ARKodeCreate();
     
     /* Call ARKodeInit to initialize the integrator memory and specify the
        user's right hand side function in y'=f(t,y), the inital time T0, and
        the initial dependent variable vector y */
     flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
   
     /* Call ARKodeSetUserData to pass rdata to user functions */
     flag = ARKodeSetUserData(arkode_mem, (void *) rdata);
   
     /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
     flag = ARKodeSetDiagnostics(arkode_mem, DFID);
   
     /* Call ARKodeSStolerances to specify the scalar relative and absolute
        tolerances */
     flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
   
     /* Call ARKDense to specify the ARKDENSE dense linear solver */
     flag = ARKDense(arkode_mem, NEQ);
   
     /* Set the Jacobian routine to Jac (user-supplied) */
     flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);
   
     /* In loop, call ARKode, print results, and test for error.
        Break out of loop when the final output time has been reached */
     realtype t = T0;
     realtype tout = T0+dTout;
     realtype u, v, w;
     printf("        t           u           v           w\n");
     printf("   -------------------------------------------\n");
     int iout;
     for (iout=0; iout<Nt; iout++) {
   
       flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
       u = NV_Ith_S(y,0);
       v = NV_Ith_S(y,1);
       w = NV_Ith_S(y,2);
       printf("  %10.6f  %10.6f  %10.6f  %10.6f\n", t, u, v, w);
   
       if (flag >= 0) {
         tout += dTout;
         tout = (tout > Tf) ? Tf : tout;
       }
     }
     printf("   -------------------------------------------\n");
   
     /* Print some final statistics */
     long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
     flag = ARKodeGetNumSteps(arkode_mem, &nst);
     flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
     flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
     flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
     flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
     flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
     flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
     flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
     flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
   
     printf("\nFinal Solver Statistics:\n");
     printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
     printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
     printf("   Total linear solver setups = %li\n", nsetups);
     printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
     printf("   Total number of Jacobian evaluations = %li\n", nje);
     printf("   Total number of Newton iterations = %li\n", nni);
     printf("   Total number of linear solver convergence failures = %li\n", ncfn);
     printf("   Total number of error test failures = %li\n\n", netf);
   
     /* Free y vector */
     N_VDestroy_Serial(y);
   
     /* Free integrator memory */
     ARKodeFree(&arkode_mem);
   
     /* close solver diagnostics output file */
     fclose(DFID);
   
     return 0;
   }
   
   
   /*-------------------------------
    * Functions called by the solver
    *-------------------------------*/
   
   /* f routine to compute the ODE RHS function f(t,y). */
   
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
   {
     realtype *rdata = (realtype *) user_data;
     realtype a  = rdata[0];
     realtype b  = rdata[1];
     realtype ep = rdata[2];
     realtype u = NV_Ith_S(y,0);
     realtype v = NV_Ith_S(y,1);
     realtype w = NV_Ith_S(y,2);
   
     /* du/dt = a - (w+1)*u + v*u^2 */
     NV_Ith_S(ydot,0) = a - (w+1.0)*u + v*u*u;
   
     /* dv/dt = w*u - v*u^2 */
     NV_Ith_S(ydot,1) = w*u - v*u*u;
   
     /* dw/dt = (b-w)/ep - w*u */
     NV_Ith_S(ydot,2) = (b-w)/ep - w*u;
   
     return 0;
   }
   
   /* Jacobian routine to compute J(t,y) = df/dy. */
   
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
   {
     realtype *rdata = (realtype *) user_data;
     realtype ep = rdata[2];
     realtype u = NV_Ith_S(y,0);
     realtype v = NV_Ith_S(y,1);
     realtype w = NV_Ith_S(y,2);
   
     /* du/dt = a - (w+1)*u + v*u^2 */
     DENSE_ELEM(J,0,0) = -(w+1.0) + 2.0*u*v;
     DENSE_ELEM(J,0,1) = u*u;
     DENSE_ELEM(J,0,2) = -u;
   
     /* dv/dt = w*u - v*u^2 */
     DENSE_ELEM(J,1,0) = w - 2.0*u*v;
     DENSE_ELEM(J,1,1) = -u*u;
     DENSE_ELEM(J,1,2) = u;
   
     /* dw/dt = (b-w)/ep - w*u */
     DENSE_ELEM(J,2,0) = -w;
     DENSE_ELEM(J,2,1) = 0.0;
     DENSE_ELEM(J,2,2) = -1.0/ep - u;
   
     return 0;
   }
   
   
Solutions
^^^^^^^^^^^^

The computed solutions will of course depend on which test is
performed:

Test 1:  Here, all three components exhibit a rapid transient change
during the first 0.2 time units, followed by a slow and smooth evolution. 

Test 2: Here, :math:`w` experiences a fast initial transient, jumping
0.5 within a few steps.  All values proceed smoothly until around
:math:`t=6.5`, when both :math:`u` and :math:`v` undergo a sharp
transition, with :math:`u` increaseing from around 0.5 to 5 and
:math:`v` decreasing from around 6 to 1 in less than 0.5 time units.
After this transition, both :math:`u` and :math:`v` continue to evolve
somewhat rapidly for another 1.4 time units, and finish off smoothly.

Test 3: Here, all components undergo very rapid initial transients
during the first 0.3 time units, and all then proceed very smoothly
for the remainder of the simulation.



Stiff nonlinear system, Fortran example (ark_bruss)
----------------------------------------------------

This is a Fortran-90 version of the same test brusselator test problem
as above.  

ODE system
^^^^^^^^^^^^

The test problem has 3 dependent variables :math:`u`, :math:`v` and
:math:`w`, that depend on the independent variable :math:`t` via the
IVP system

.. math::

   \frac{du}{dt} &= a - (w+1)u + v u^2, \\
   \frac{dv}{dt} &= w u - v u^2, \\
   \frac{dw}{dt} &= \frac{b-w}{\varepsilon} - w u.

We integrate over the interval :math:`0 \le t \le 10`, with the
initial conditions :math:`u(0) = 3.9`, :math:`v(0) = 1.1`, :math:`w(0) = 2.8`,
and parameters :math:`a=1.2`, :math:`b=2.5` and
:math:`\varepsilon=10^{-5}`.  After each unit time interval, the
solution is output to the screen.


Numerical method
^^^^^^^^^^^^^^^^^

Since this driver and utility functions are written in Fortran-90,
this example demonstrates the use of the FARKODE interface for the
ARKode solver.  For time integration, this example uses the
fourth-order additive Runge-Kutta method, where the right-hand sides
are broken up as

.. math::

   f_E(t,u,v,w) = \left(\begin{array}{c} a - (w+1)u + v u^2 \\ 
     w u - v u^2 \\ - w u  \end{array}\right), \quad\text{and}\quad 
   f_I(t,u,v,w) = \left(\begin{array}{c} 0\\0\\
     \frac{b-w}{\varepsilon}\end{array}\right). 

The implicit systems are solved using the built-in modified Newton
iteration, with the ARKDENSE dense linear solver.  Both the Jacobian
routine and right-hand side functions are supplied by functions
provided in the example file.

The only non-default solver options are the tolerances
:math:`atol=10^{-10}` and :math:`rtol=10^{-6}`, adaptivity method 2 (I
controller), a maximum of 8 Newton iterations per step, a nonlinear
solver convergence coefficient :math:`nlscoef=10^{-8}`, and a maximum
of 1000 internal time steps.



Routines
^^^^^^^^^^^^^^^^^

We reproduce the relevant aspects of the ``main()`` routine and
auxiliary functions here for explanatory purposes (see the in-line
comments for details; error-checking has been removed for brevity):

.. code-block:: fortran

   program driver
   
     ! Declarations
     implicit none
   
     ! general problem variables
     integer*8, parameter :: NEQ=3
     real*8,    parameter :: T0=0.d0, Tf=10.d0
     real*8    :: dTout, Tout, Tcur, rtol, atol, rout(6)
     integer   :: it, Nt, ier, btable2(2)
     integer*8 :: iout(22)
     real*8, dimension(NEQ) :: y
   
     ! real/integer parameters to pass through to supplied functions
     !    ipar(1) -> unused
     !    rpar(1) -> "a" parameter
     !    rpar(2) -> "b" parameter 
     !    rpar(3) -> "ep" parameter
     integer :: ipar
     real*8  :: rpar(3)
   
     ! solver parameters
     integer :: order, adapt_method, maxcor
     real*8 :: nlscoef
   
     !-----------------------
     ! set some solver parameters
     order = 4
     adapt_method = 2
     maxcor = 8
     nlscoef = 1.d-8
   
     ! time-stepping information
     dTout = (Tf-T0)/10.d0
     Nt = Tf/dTout + 0.5
   
     ! set initial conditions, problem parameters
     y(1) = 3.9d0     ! u0
     y(2) = 1.1d0     ! v0
     y(3) = 2.8d0     ! w0
     rpar(1) = 1.2    ! a
     rpar(2) = 2.5    ! b
     rpar(3) = 1.d-5  ! ep
   
     ! set tolerances according to problem specifications
     atol = 1.d-10
     rtol = 1.d-6
     
     ! initialize vector module
     call FNVInitS(4, NEQ, ier)
   
     ! initialize ARKode solver to use IMEX integrator, scalar tolerances
     call FARKMalloc(T0, y, 2, 1, rtol, atol, &
                     iout, rout, ipar, rpar, ier)
   
     ! set integrator options
     call FARKSetIin('ORDER', order, ier)
     call FARKSetIin('ADAPT_METHOD', adapt_method, ier)
     call FARKSetIin('MAX_NITERS', maxcor, ier)
     call FARKSetRin('NLCONV_COEF', nlscoef, ier)
     call FARKSetIin('MAX_NSTEPS', 1000, ier)
   
     ! specify use of dense linear solver
     call FARKDense(NEQ, ier)
     call FARKDenseSetJac(1, ier)
   
     ! loop over time outputs
     Tout = T0
     Tcur = T0
     print *, '        t           u           v           w'
     print *, '  ----------------------------------------------------'
     print '(3x,4(es12.5,1x))', Tcur, y
     do it = 1,Nt
   
        ! set next output time
        Tout = min(Tout + dTout, Tf)
   
        ! call solver
        call FARKode(Tout, Tcur, y, 1, ier)
        if (ier < 0) then
           print *, 'Error at step ',it,', FARKode return flag =',ier
           exit
        end if
   
        ! output current solution information
        print '(3x,4(es12.5,1x))', Tcur, y
   
     end do
     print *, '  ----------------------------------------------------'
   
     ! output solver statistics
     print *, '  '
     print *, 'Final Solver Statistics:'
     print '(2(A,i7),A)', '   Internal solver steps =', iout(3), &
          ' (attempted =', iout(6), ')'
     print '(2(A,i7))', '   Total RHS evals:  Fe =', iout(7), &
          ',  Fi =', iout(8)
     print '(A,i7)', '   Total linear solver setups =', iout(9)
     print '(A,i7)', '   Total RHS evals for setting up the linear system =', iout(17)
     print '(A,i7)', '   Total number of Jacobian evaluations =', iout(18)
     print '(A,i7)', '   Total number of Newton iterations =', iout(11)
     print '(A,i7)', '   Total number of nonlinear solver convergence failures =', &
          iout(12)
     print '(A,i7)', '   Total number of error test failures =', iout(10)
     print *, '  '
   
     ! output final solution
     print *, '     y(Tf) =', y
     print *, '  '
   
     ! clean up
     call FARKFree()
   
   end program driver
   !-----------------------------------------------------------------
   
   
   
   !-----------------------------------------------------------------
   ! Required subroutines for FARKODE interface
   !-----------------------------------------------------------------
   
   
   subroutine farkifun(t, y, ydot, ipar, rpar, ier)
   !-----------------------------------------------------------------
   ! Implicit portion of the right-hand side of the ODE system
   !-----------------------------------------------------------------
   
     ! Declarations
     implicit none
   
     ! Arguments
     real*8,  intent(in)  :: t, rpar(3)
     integer, intent(in)  :: ipar(1)
     integer, intent(out) :: ier
     real*8,  intent(in)  :: y(3)
     real*8,  intent(out) :: ydot(3)
   
     ! temporary variables
     real*8 :: u, v, w, a, b, ep
   
     ! set temporary values
     a  = rpar(1)
     b  = rpar(2)
     ep = rpar(3)
     u  = y(1)
     v  = y(2)
     w  = y(3)
   
     ! fill implicit RHS
     ydot(1) = 0.d0
     ydot(2) = 0.d0
     ydot(3) = (b-w)/ep
     ier = 0
     
   end subroutine farkifun
   !-----------------------------------------------------------------
   
   
   subroutine farkefun(t, y, ydot, ipar, rpar, ier)
   !-----------------------------------------------------------------
   ! Explicit portion of the right-hand side of the ODE system
   !-----------------------------------------------------------------
   
     ! Declarations
     implicit none
   
     ! Arguments
     real*8,  intent(in)  :: t, rpar(3)
     integer, intent(in)  :: ipar(1)
     integer, intent(out) :: ier
     real*8,  intent(in)  :: y(3)
     real*8,  intent(out) :: ydot(3)
   
     ! temporary variables
     real*8 :: u, v, w, a, b, ep
   
     ! set temporary values
     a  = rpar(1)
     b  = rpar(2)
     ep = rpar(3)
     u  = y(1)
     v  = y(2)
     w  = y(3)
   
     ! fill explicit RHS
     ydot(1) = a - (w+1.d0)*u + v*u*u
     ydot(2) = w*u - v*u*u
     ydot(3) = -w*u
     ier = 0
     
   end subroutine farkefun
   !-----------------------------------------------------------------
   
   
   subroutine farkdjac(neq,t,y,fy,DJac,h,ipar,rpar,wk1,wk2,wk3,ier)
   !-----------------------------------------------------------------
   ! Jacobian computation routine
   !-----------------------------------------------------------------
   
     ! Declarations
     implicit none
   
     ! Arguments
     real*8,  intent(in)  :: t, h, rpar(3)
     integer, intent(in)  :: neq, ipar(1)
     integer, intent(out) :: ier
     real*8,  intent(in), dimension(neq) :: y, fy, wk1, wk2, wk3
     real*8,  intent(out) :: DJac(neq,neq)
   
     ! temporary variables
     real*8 :: u, v, w, a, b, ep
   
     ! set temporary values
     a  = rpar(1)
     b  = rpar(2)
     ep = rpar(3)
     u  = y(1)
     v  = y(2)
     w  = y(3)
   
     ! fill implicit Jacobian
     DJac = 0.d0
     DJac(3,3) = -1.d0/ep
     ier = 0
     
     
   end subroutine farkdjac
   !-----------------------------------------------------------------
   

Solutions
^^^^^^^^^^^^

With this setup, all three solution components exhibit a rapid
transient change during the first 0.2 time units, followed by a slow
and smooth evolution, as seen in the figure below.




Stiff nonlinear system example (ark_robertson)
------------------------------------------------

ODE system
^^^^^^^^^^^^

This test simulates the Robertson problem, corresponding to the
kinetics of an autocatalytic reaction.  This is an ODE system with 3
components, :math:`Y = [u,\, v,\, w]^T`, satisfying the equations,

.. math::

   \frac{du}{dt} &= -0.04 u + 10^4 v w, \\
   \frac{dv}{dt} &= 0.04 u - 10^4 v w - 3\cdot10^7 v^2, \\
   \frac{dw}{dt} &= 3\cdot10^7 v^2.

We integrate over the interval :math:`0\le t\le 10^{11}`, with initial
conditions  :math:`Y(0) = [1,\, 0,\, 0]^T`. 


Numerical method
^^^^^^^^^^^^^^^^^

In the input file, ``input_robertson.txt``, we allow specification of
the desired relative and absolute tolerances. 
 
This program solves the problem with one of the solvers, ERK, DIRK or
ARK.  For DIRK and ARK, implicit subsystems are solved using a Newton
iteration with the ARKDENSE dense linear solver, and a user-supplied
Jacobian routine. 

100 outputs are printed at equal intervals, and run statistics are
printed at the end.


Routines
^^^^^^^^^^^^^^^^^

We reproduce the relevant aspects of the ``main()`` routine and
auxiliary functions here for explanatory purposes (see the in-line
comments for details; error-checking has been removed for brevity):

.. code-block:: c

   /* Header files */
   #include <stdio.h>
   #include <math.h>
   #include <arkode/arkode.h>
   #include <nvector/nvector_serial.h>
   #include <arkode/arkode_dense.h>
   #include <sundials/sundials_dense.h>
   #include <sundials/sundials_types.h>
   
   /* User-supplied Functions Called by the Solver */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
   
   
   /* Main Program */
   int main()
   {
     /* general problem parameters */
     realtype T0 = RCONST(0.0);
     realtype Tf = RCONST(1.e11);
     realtype dTout = (Tf-T0)/100;
     int Nt = ceil(Tf/dTout);
     realtype u0, v0, w0, h0;
     long int NEQ = 3;
   
     /* declare solver parameters */
     int flag;
   
     /* general problem variables */
     int idense;
     N_Vector y = NULL;
     void *arkode_mem = NULL;
   
     /* open solver diagnostics output file for writing */
     FILE *DFID;
     DFID=fopen("diags_ark_robertson.txt","w");
     
     /* set up the initial conditions */
     u0 = RCONST(1.0);
     v0 = RCONST(0.0);
     w0 = RCONST(0.0);
   
     /* Initial problem output */
     printf("\nRobertson ODE test problem:\n");
     printf("    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n",u0,v0,w0);
   
     /* Create serial vectors of length NEQ for initial condition */
     y = N_VNew_Serial(NEQ);
   
     /* Set initial conditions into y */
     NV_Ith_S(y,0) = u0;
     NV_Ith_S(y,1) = v0;
     NV_Ith_S(y,2) = w0;
   
     /* Call ARKodeCreate to create the solver memory */
     arkode_mem = ARKodeCreate();
     
     /* Call ARKodeInit to initialize the integrator memory and specify the
        user's right hand side function in y'=f(t,y), the inital time T0, and
        the initial dependent variable vector y */
     flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
   
     /* Set tolerances */
     realtype reltol = 1.e-4;
     realtype abstol = 1.e-8;
     h0 = 1.e-4 * reltol;
   
     /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
     flag = ARKodeSetDiagnostics(arkode_mem, DFID);
   
     /* Set custom initial step */
     flag = ARKodeSetInitStep(arkode_mem, h0);
   
     /* Increase maximum number of error test failures */
     flag = ARKodeSetMaxErrTestFails(arkode_mem, 20);
   
     /* Call ARKodeSetmaxNonlinIters to increase default for this problem*/
     flag = ARKodeSetMaxNonlinIters(arkode_mem, 8);
   
     /* Call ARKodeSetNonlinConvCoef */
     flag = ARKodeSetNonlinConvCoef(arkode_mem, 1.e-7);
   
     /* Call ARKodeSetMaxNumSteps to increase default (for testing) */
     flag = ARKodeSetMaxNumSteps(arkode_mem, 100000);
   
     /* Call ARKodeSStolerances to specify the scalar relative and absolute
        tolerances */
     flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
   
     /* Call ARKDense to specify the ARKDENSE dense linear solver */
     flag = ARKDense(arkode_mem, NEQ);
   
     /* Set the Jacobian routine to Jac (user-supplied) */
     flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);
   
     /* In loop, call ARKode, print results, and test for error.
        Break out of loop when the final output time has been reached */
     realtype t = T0;
     realtype tout = T0+dTout;
     realtype u, v, w;
     printf("        t           u           v           w\n");
     printf("   --------------------------------------------------\n");
     printf("  %10.3e  %12.5e  %12.5e  %12.5e\n", 
   	 t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
     int iout;
     for (iout=0; iout<Nt; iout++) {
   
       flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
       printf("  %10.3e  %12.5e  %12.5e  %12.5e\n", 
   	   t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
   
       if (flag >= 0) {
         tout += dTout;
         tout = (tout > Tf) ? Tf : tout;
       }
   
     }
     printf("   --------------------------------------------------\n");
   
     /* Print some final statistics */
     long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
     flag = ARKodeGetNumSteps(arkode_mem, &nst);
     flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
     flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
     flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
     flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
     flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
     flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
     flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
     flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
   
     printf("\nFinal Solver Statistics:\n");
     printf("   Internal solver steps = %li (attempted = %li)\n", 
   	 nst, nst_a);
     printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
     printf("   Total linear solver setups = %li\n", nsetups);
     printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
     printf("   Total number of Jacobian evaluations = %li\n", nje);
     printf("   Total number of Newton iterations = %li\n", nni);
     printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
     printf("   Total number of error test failures = %li\n", netf);
   
     /* Free y vector */
     N_VDestroy_Serial(y);
   
     /* Free integrator memory */
     ARKodeFree(&arkode_mem);
   
     /* close solver diagnostics output file */
     fclose(DFID);
   
     return 0;
   }
   
   
   /*-------------------------------
    * Functions called by the solver
    *-------------------------------*/
   
   /* f routine to compute the ODE RHS function f(t,y). */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
   {
     realtype u = NV_Ith_S(y,0);
     realtype v = NV_Ith_S(y,1);
     realtype w = NV_Ith_S(y,2);
   
     /* du/dt = -0.04*u + 1.e4*v*w */
     NV_Ith_S(ydot,0) = -0.04*u + 1.e4*v*w;
   
     /* dv/dt = 0.04*u - 1.e4*v*w - 3.e7*v*v */
     NV_Ith_S(ydot,1) = 0.04*u - 1.e4*v*w - 3.e7*v*v;
   
     /* dw/dt = 3.e7*v*v */
     NV_Ith_S(ydot,2) = 3.e7*v*v;
   
     return 0;
   }
   
   /* Jacobian routine to compute J(t,y) = df/dy. */
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
   {
     realtype v = NV_Ith_S(y,1);
     realtype w = NV_Ith_S(y,2);
     SetToZero(J);
   
     /* du/dt = -0.04*u + 1.e4*v*w */
     DENSE_ELEM(J,0,0) = -0.04;
     DENSE_ELEM(J,0,1) = 1.e4*w;
     DENSE_ELEM(J,0,2) = 1.e4*v;
   
     /* dv/dt = 0.04*u - 1.e4*v*w - 3.e7*v*v */
     DENSE_ELEM(J,1,0) = 0.04;
     DENSE_ELEM(J,1,1) = -1.e4*w - 6.e7*v;
     DENSE_ELEM(J,1,2) = -1.e4*v;
   
     /* dw/dt = 3.e7*v*v */
     DENSE_ELEM(J,2,1) = 6.e7*v;
   
     return 0;
   }
   

Solutions
^^^^^^^^^^^^



Stiff nonlinear system with root-finding example (ark_robertson_root)
-----------------------------------------------------------------------

ODE system
^^^^^^^^^^^^

This is the same test as in the above problem (the Robertson problem).
This is an ODE system with 3 components, :math:`Y = [u,\, v,\, w]^T`,
satisfying the equations,

.. math::

   \frac{du}{dt} &= -0.04 u + 10^4 v w, \\
   \frac{dv}{dt} &= 0.04 u - 10^4 v w - 3\cdot10^7 v^2, \\
   \frac{dw}{dt} &= 3\cdot10^7 v^2.

We integrate over the interval :math:`0\le t\le 10^{11}`, with initial
conditions  :math:`Y(0) = [1,\, 0,\, 0]^T`. 


Numerical method
^^^^^^^^^^^^^^^^^

In the input file, ``input_robertson.txt``, we allow specification of
the desired relative and absolute tolerances. 
 
This program solves the problem with one of the solvers, ERK, DIRK or
ARK.  For DIRK and ARK, implicit subsystems are solved using a Newton
iteration with the ARKDENSE dense linear solver, and a user-supplied
Jacobian routine. 

100 outputs are printed at equal intervals, and run statistics are
printed at the end.

However, unlike in the previous problem, while integrating the system,
we use the rootfinding feature of ARKode to find the times at which
either :math:`u=10^{-4}` or :math:`w=10^{-2}`.



Routines
^^^^^^^^^^^^^^^^^

We reproduce the relevant aspects of the ``main()`` routine and
auxiliary functions here for explanatory purposes (see the in-line
comments for details; error-checking has been removed for brevity):

.. code-block:: c

   /* Header files */
   #include <stdio.h>
   #include <math.h>
   #include <arkode/arkode.h>
   #include <nvector/nvector_serial.h>
   #include <arkode/arkode_dense.h>
   #include <sundials/sundials_dense.h>
   #include <sundials/sundials_types.h>
   
   /* User-supplied Functions Called by the Solver */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
   static int g(realtype t, N_Vector y, 
   	     realtype *gout, void *user_data);
   

   /* Main Program */
   int main()
   {
     /* general problem parameters */
     realtype T0 = RCONST(0.0);
     realtype T1 = RCONST(0.4);
     realtype TMult = RCONST(10.0);
     int Nt = 12;
     realtype u0, v0, w0;
     long int NEQ = 3;
     int rootsfound[2];
     long int nst, nst_a, nfe, nfi, nsetups;
     long int nje, nfeLS, nni, ncfn, netf, nge;
     
     /* declare solver parameters */
     int flag, rtflag;
   
     /* general problem variables */
     N_Vector y = NULL;
     N_Vector atols = NULL;
     void *arkode_mem = NULL;
   
     /* open solver diagnostics output file for writing */
     FILE *DFID;
     DFID=fopen("diags_ark_robertson_root.txt","w");
     
     /* set up the initial conditions */
     u0 = RCONST(1.0);
     v0 = RCONST(0.0);
     w0 = RCONST(0.0);
   
     /* Initial problem output */
     printf("\nRobertson ODE test problem (with rootfinding):\n");
     printf("    initial conditions:  u0 = %g,  v0 = %g,  w0 = %g\n",u0,v0,w0);
   
     /* Create serial vectors of length NEQ for initial condition */
     y = N_VNew_Serial(NEQ);
     atols = N_VNew_Serial(NEQ);
   
     /* Set initial conditions into y */
     NV_Ith_S(y,0) = u0;
     NV_Ith_S(y,1) = v0;
     NV_Ith_S(y,2) = w0;
   
     /* Call ARKodeCreate to create the solver memory */
     arkode_mem = ARKodeCreate();
     
     /* Call ARKodeInit to initialize the integrator memory and specify the
        user's right hand side function in y'=f(t,y), the inital time T0, and
        the initial dependent variable vector y */
     flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
   
     /* Set tolerances */
     realtype reltol = RCONST(1.0e-4);
     NV_Ith_S(atols,0) = RCONST(1.0e-8);
     NV_Ith_S(atols,1) = RCONST(1.0e-8);
     NV_Ith_S(atols,2) = RCONST(1.0e-8);
   
     /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
     flag = ARKodeSetDiagnostics(arkode_mem, DFID);
   
     /* Increase maximum number of error test failures */
     flag = ARKodeSetMaxErrTestFails(arkode_mem, 20);
   
     /* Call ARKodeSetmaxNonlinIters to increase default for this problem*/
     flag = ARKodeSetMaxNonlinIters(arkode_mem, 8);
   
     /* Call ARKodeSetNonlinConvCoef */
     flag = ARKodeSetNonlinConvCoef(arkode_mem, 1.e-7);
   
     /* Call ARKodeSetMaxNumSteps to increase default (for testing) */
     flag = ARKodeSetMaxNumSteps(arkode_mem, 100000);
   
     /* Call ARKodeSStolerances to specify the scalar relative and absolute
        tolerances */
     flag = ARKodeSVtolerances(arkode_mem, reltol, atols);
   
     /* Call ARKodeRootInit to specify the root function with 2 equations */
     flag = ARKodeRootInit(arkode_mem, 2, g);
   
     /* Call ARKDense to specify the ARKDENSE dense linear solver */
     flag = ARKDense(arkode_mem, NEQ);
   
     /* Set the Jacobian routine to Jac (user-supplied) */
     flag = ARKDlsSetDenseJacFn(arkode_mem, Jac);
   
     /* Print out root-finding expectations */
     printf("\n Roots should be found at the following times and values:\n");
     printf("   t=2.64019e-1, u=9.89965e-1, v=3.47049e-5,  w=1.00000e-2  [ 0 1]\n");
     printf("   t=2.07956e+7, u=1.00000e-4, v=3.96207e-10, w=9.99900e-1  [-1 0]\n\n");
   
     /* In loop, call ARKode, print results, and test for error.
        Break out of loop when Nt preset output times have been reached */
     realtype t = T0;
     printf("        t             u             v             w\n");
     printf("   -----------------------------------------------------\n");
     printf("  %12.5e  %12.5e  %12.5e  %12.5e\n", 
   	 t, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
     realtype tout = T1;
     int iout=0;
     while(1) {
   
       flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
       printf("  %12.5e  %12.5e  %12.5e  %12.5e\n",  t, 
   	   NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2));
   
       if (flag == ARK_ROOT_RETURN) {
         rtflag = ARKodeGetRootInfo(arkode_mem, rootsfound);
         printf("      rootsfound[] = %3d %3d\n", 
   	     rootsfound[0], rootsfound[1]);
       }
   
       if (flag >= 0) {
         iout++;
         tout *= TMult;
       }
   
       if (iout == Nt) break;
     }
     printf("   -----------------------------------------------------\n");
   
     /* Print some final statistics */
     flag = ARKodeGetNumSteps(arkode_mem, &nst);
     flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
     flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
     flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
     flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
     flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
     flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
     flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
     flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
     flag = ARKodeGetNumGEvals(arkode_mem, &nge);
   
     printf("\nFinal Solver Statistics:\n");
     printf("   Internal solver steps = %li (attempted = %li)\n", 
   	 nst, nst_a);
     printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
     printf("   Total linear solver setups = %li\n", nsetups);
     printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
     printf("   Total number of Jacobian evaluations = %li\n", nje);
     printf("   Total number of Newton iterations = %li\n", nni);
     printf("   Total root-function g evals = %li\n", nge);
     printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
     printf("   Total number of error test failures = %li\n", netf);
   
     /* Free y vector */
     N_VDestroy_Serial(y);
   
     /* Free integrator memory */
     ARKodeFree(&arkode_mem);
   
     /* close solver diagnostics output file */
     fclose(DFID);
   
     return 0;
   }
   
   
   /*-------------------------------
    * Functions called by the solver
    *-------------------------------*/
   
   /* f routine to compute the ODE RHS function f(t,y). */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
   {
     realtype u = NV_Ith_S(y,0);
     realtype v = NV_Ith_S(y,1);
     realtype w = NV_Ith_S(y,2);
   
     /* du/dt = -0.04*u + 1.e4*v*w */
     NV_Ith_S(ydot,0) = -0.04*u + 1.e4*v*w;
   
     /* dv/dt = 0.04*u - 1.e4*v*w - 3.e7*v*v */
     NV_Ith_S(ydot,1) = 0.04*u - 1.e4*v*w - 3.e7*v*v;
   
     /* dw/dt = 3.e7*v*v */
     NV_Ith_S(ydot,2) = 3.e7*v*v;
   
     return 0;
   }
   
   /* Jacobian routine to compute J(t,y) = df/dy. */
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
   {
     realtype v = NV_Ith_S(y,1);
     realtype w = NV_Ith_S(y,2);
     SetToZero(J);
   
     /* du/dt = -0.04*u + 1.e4*v*w */
     DENSE_ELEM(J,0,0) = -0.04;
     DENSE_ELEM(J,0,1) = 1.e4*w;
     DENSE_ELEM(J,0,2) = 1.e4*v;
   
     /* dv/dt = 0.04*u - 1.e4*v*w - 3.e7*v*v */
     DENSE_ELEM(J,1,0) = 0.04;
     DENSE_ELEM(J,1,1) = -1.e4*w - 6.e7*v;
     DENSE_ELEM(J,1,2) = -1.e4*v;
   
     /* dw/dt = 3.e7*v*v */
     DENSE_ELEM(J,2,1) = 6.e7*v;
   
     return 0;
   }
   
   /* g routine to compute the root-finding function g(t,y). */
   static int g(realtype t, N_Vector y, realtype *gout, void *user_data)
   {
     realtype u = NV_Ith_S(y,0);
     realtype w = NV_Ith_S(y,2);
   
     /* check for u == 1e-4 */
     gout[0] = u - RCONST(0.0001);
   
     /* check for w == 1e-2 */
     gout[1] = w - RCONST(0.01);
   
     return 0;
   }


Solutions
^^^^^^^^^^^^



Stiff PDE system example (ark_brusselator1D)
---------------------------------------------

ODE system
^^^^^^^^^^^^

This test simulates a brusselator problem from chemical kinetics, but
in PDE form.  This system has 3 components, :math:`Y = [u,\, v,\, w]^T`,  
satisfying the equations,

.. math::

   \frac{\partial u}{\partial t} &= d_u \frac{\partial^2 u}{\partial
      x^2} + a - (w+1) u + v u^2, \\
   \frac{\partial v}{\partial t} &= d_v \frac{\partial^2 v}{\partial
      x^2} + w u - v u^2, \\
   \frac{\partial w}{\partial t} &= d_w \frac{\partial^2 w}{\partial
      x^2} + \frac{b-w}{\varepsilon} - w u.

We integrate for :math:`t \in [0, 80]`, and :math:`x \in [0, 1]`, with
initial conditions 

.. math::

   u(0,x) &=  a + \frac{1}{10} \sin(\pi x),\\
   v(0,x) &= \frac{b}{a} + \frac{1}{10}\sin(\pi x),\\
   w(0,x) &=  b + \frac{1}{10}\sin(\pi x),

and with stationary boundary conditions, i.e. 

.. math::

   \frac{\partial u}{\partial t}(t,0) &= \frac{\partial u}{\partial t}(t,1) = 0,\\
   \frac{\partial v}{\partial t}(t,0) &= \frac{\partial v}{\partial t}(t,1) = 0,\\
   \frac{\partial w}{\partial t}(t,0) &= \frac{\partial w}{\partial t}(t,1) = 0.

We note that these can also be implemented as Dirichlet boundary
conditions with values identical to the initial conditions. 



Numerical method
^^^^^^^^^^^^^^^^^

The spatial derivatives are computed using second-order 
centered differences, with the data distributed over :math:`N` points 
on a uniform spatial grid.

The number of spatial points :math:`N`, the parameters :math:`a`,
:math:`b`, :math:`d_u`, :math:`d_v`, :math:`d_w` and
:math:`\varepsilon`, as well as the desired relative and absolute
solver tolerances, are provided in the input file ``input_brusselator1D.txt``.
 
This program solves the problem with a DIRK method, using a Newton
iteration with the ARKBAND banded linear solver, and a user-supplied
Jacobian routine. 

100 outputs are printed at equal intervals, and run statistics 
are printed at the end.


Routines
^^^^^^^^^^^^^^^^^

We reproduce the relevant aspects of the ``main()`` routine and
auxiliary functions here for explanatory purposes (see the in-line
comments for details; error-checking has been removed for brevity):

.. code-block:: c

   /* Header files */
   #include <stdio.h>
   #include <stdlib.h>
   #include <math.h>
   #include <arkode/arkode.h>
   #include <nvector/nvector_serial.h>
   #include <arkode/arkode_band.h>
   #include <sundials/sundials_band.h>
   #include <sundials/sundials_types.h>
   
   /* accessor macros between (x,v) location and 1D NVector array */
   #define IDX(x,v) (3*(x)+v)
   
   /* user data structure */
   typedef struct {  
     long int N;    /* number of intervals     */
     realtype dx;   /* mesh spacing            */
     realtype a;    /* constant forcing on u   */
     realtype b;    /* steady-state value of w */
     realtype du;   /* diffusion coeff for u   */
     realtype dv;   /* diffusion coeff for v   */
     realtype dw;   /* diffusion coeff for w   */
     realtype ep;   /* stiffness parameter     */
   } *UserData;
   
   /* User-supplied Functions Called by the Solver */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
   static int Jac(long int N, long int mu, long int ml,
                  realtype t, N_Vector y, N_Vector fy, 
                  DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
   
   /* Private functions  */
   static int LaplaceMatrix(realtype c, DlsMat Jac, UserData udata);
   static int ReactionJac(realtype c, N_Vector y, DlsMat Jac, UserData udata);
   
   
   /* Main Program */
   int main()
   {
     /* general problem parameters */
     realtype T0 = RCONST(0.0);
     realtype Tf = RCONST(10.0);
     int Nt = 100;
     int Nvar = 3;
     UserData udata = NULL;
     realtype *data;
     long int N, NEQ, i;
   
     /* general problem variables */
     int flag;
     N_Vector y = NULL;
     N_Vector umask = NULL;
     N_Vector vmask = NULL;
     N_Vector wmask = NULL;
     void *arkode_mem = NULL;
   
     /* allocate udata structure */
     udata = (UserData) malloc(sizeof(*udata));
   
     /* read problem parameter and tolerances from input file:
        N - number of spatial discretization points
        a - constant forcing on u
        b - steady-state value of w
        du - diffusion coefficient for u
        dv - diffusion coefficient for v
        dw - diffusion coefficient for w
        ep - stiffness parameter
        reltol - desired relative tolerance
        abstol - desired absolute tolerance */
     double a, b, du, dv, dw, ep, reltol, abstol;
     FILE *FID;
     FID=fopen("input_brusselator1D.txt","r");
     fscanf(FID,"  N = %li\n", &N);
     fscanf(FID,"  a = %lf\n", &a);
     fscanf(FID,"  b = %lf\n", &b);
     fscanf(FID,"  du = %lf\n", &du);
     fscanf(FID,"  dv = %lf\n", &dv);
     fscanf(FID,"  dw = %lf\n", &dw);
     fscanf(FID,"  ep = %lf\n", &ep);
     fscanf(FID,"  reltol = %lf\n", &reltol);
     fscanf(FID,"  abstol = %lf\n", &abstol);
     fclose(FID);
   
     /* store the inputs in the UserData structure */
     udata->N  = N;
     udata->a  = a;
     udata->b  = b;
     udata->du = du;
     udata->dv = dv;
     udata->dw = dw;
     udata->ep = ep;
   
     /* open solver diagnostics output file for writing */
     FILE *DFID;
     DFID=fopen("diags_ark_bruss1D.txt","w");
     
     /* set total allocated vector length */
     NEQ = Nvar*udata->N;
   
     /* Initial problem output */
     printf("\n1D Brusselator PDE test problem:\n");
     printf("    N = %li,  NEQ = %li\n", udata->N, NEQ);
     printf("    problem parameters:  a = %g,  b = %g,  ep = %g\n",
   	 udata->a, udata->b, udata->ep);
     printf("    diffusion coefficients:  du = %g,  dv = %g,  dw = %g\n", 
   	 udata->du, udata->dv, udata->dw);
     printf("    reltol = %.1e,  abstol = %.1e\n\n", reltol, abstol);
   
     /* Create serial vector of length NEQ for initial condition */
     y = N_VNew_Serial(NEQ);
   
     /* set spatial mesh spacing */
     udata->dx = RCONST(1.0)/(N-1);
   
     /* output mesh to disk */
     FID=fopen("bruss_mesh.txt","w");
     for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->dx*i);
     fclose(FID);
     
   
     /* Access data array for new NVector y */
     data = N_VGetArrayPointer(y);
   
     /* Set initial conditions into y */
     realtype pi = RCONST(4.0)*atan(RCONST(1.0));
     for (i=0; i<N; i++) {
       data[IDX(i,0)] =  a  + RCONST(0.1)*sin(pi*i*udata->dx);  /* u */
       data[IDX(i,1)] = b/a + RCONST(0.1)*sin(pi*i*udata->dx);  /* v */
       data[IDX(i,2)] =  b  + RCONST(0.1)*sin(pi*i*udata->dx);  /* w */
     }
   
     /* Create serial vector masks for each solution component */
     umask = N_VNew_Serial(NEQ);
     vmask = N_VNew_Serial(NEQ);
     wmask = N_VNew_Serial(NEQ);
   
     /* Set mask array values for each solution component */
     N_VConst(0.0, umask);
     data = N_VGetArrayPointer(umask);
     for (i=0; i<N; i++)  data[IDX(i,0)] = RCONST(1.0);
   
     N_VConst(0.0, vmask);
     data = N_VGetArrayPointer(vmask);
     for (i=0; i<N; i++)  data[IDX(i,1)] = RCONST(1.0);
   
     N_VConst(0.0, wmask);
     data = N_VGetArrayPointer(wmask);
     for (i=0; i<N; i++)  data[IDX(i,2)] = RCONST(1.0);
   
   
     /* Call ARKodeCreate to create the solver memory and specify the 
        Backward Differentiation Formula and the use of a Newton iteration */
     arkode_mem = ARKodeCreate();
     
     /* Call ARKodeInit to initialize the integrator memory and specify the
        user's right hand side function in y'=f(t,y), the inital time T0, and
        the initial dependent variable vector y */
     flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
     
     /* Call ARKodeSetUserData to pass rdata to user functions */
     flag = ARKodeSetUserData(arkode_mem, (void *) udata);
   
     /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
     flag = ARKodeSetDiagnostics(arkode_mem, DFID);
   
     /* Call ARKodeSStolerances to specify the scalar relative and absolute
        tolerances */
     flag = ARKodeSStolerances(arkode_mem, reltol, abstol);
   
     /* Call ARKBand to specify the ARKBAND band linear solver */
     flag = ARKBand(arkode_mem, NEQ, 4, 4);
   
     /* Set the Jacobian routine to Jac (user-supplied) */
     flag = ARKDlsSetBandJacFn(arkode_mem, Jac);
   
     /* Open output stream for results, access data arrays */
     FILE *UFID=fopen("bruss_u.txt","w");
     FILE *VFID=fopen("bruss_v.txt","w");
     FILE *WFID=fopen("bruss_w.txt","w");
     data = N_VGetArrayPointer(y);
   
     /* output initial condition to disk */
     for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[IDX(i,0)]);
     for (i=0; i<N; i++)  fprintf(VFID," %.16e", data[IDX(i,1)]);
     for (i=0; i<N; i++)  fprintf(WFID," %.16e", data[IDX(i,2)]);
     fprintf(UFID,"\n");
     fprintf(VFID,"\n");
     fprintf(WFID,"\n");
   
     /* In loop, call ARKode, print results, and test for error.
        Break out of loop when the final output time has been reached */
     realtype t  = T0;
     realtype dTout = Tf/Nt;
     realtype tout = T0+dTout;
     realtype u, v, w;
     printf("        t      ||u||_rms   ||v||_rms   ||w||_rms\n");
     printf("   ----------------------------------------------\n");
     int iout;
     for (iout=0; iout<Nt; iout++) {
   
       flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
       u = N_VWL2Norm(y,umask);
       u = sqrt(u*u/N);
       v = N_VWL2Norm(y,vmask);
       v = sqrt(v*v/N);
       w = N_VWL2Norm(y,wmask);
       w = sqrt(w*w/N);
       printf("  %10.6f  %10.6f  %10.6f  %10.6f\n", t, u, v, w);
   
       /* output results to disk */
       for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[IDX(i,0)]);
       for (i=0; i<N; i++)  fprintf(VFID," %.16e", data[IDX(i,1)]);
       for (i=0; i<N; i++)  fprintf(WFID," %.16e", data[IDX(i,2)]);
       fprintf(UFID,"\n");
       fprintf(VFID,"\n");
       fprintf(WFID,"\n");
   
       if (flag >= 0) {
         tout += dTout;
         tout = (tout > Tf) ? Tf : tout;
       }
     }
     printf("   ----------------------------------------------\n");
     fclose(UFID);
     fclose(VFID);
     fclose(WFID);
       
   
     /* Print some final statistics */
     long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
     flag = ARKodeGetNumSteps(arkode_mem, &nst);
     flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
     flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
     flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
     flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
     flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
     flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
     flag = ARKDlsGetNumJacEvals(arkode_mem, &nje);
     flag = ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
   
     printf("\nFinal Solver Statistics:\n");
     printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
     printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
     printf("   Total linear solver setups = %li\n", nsetups);
     printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
     printf("   Total number of Jacobian evaluations = %li\n", nje);
     printf("   Total number of Newton iterations = %li\n", nni);
     printf("   Total number of linear solver convergence failures = %li\n", ncfn);
     printf("   Total number of error test failures = %li\n\n", netf);
   
     /* Free vectors */
     N_VDestroy_Serial(y);
     N_VDestroy_Serial(umask);
     N_VDestroy_Serial(vmask);
     N_VDestroy_Serial(wmask);
   
     /* Free user data */
     free(udata);
   
     /* Free integrator memory */
     ARKodeFree(&arkode_mem);
   
     /* close solver diagnostics output file */
     fclose(DFID);
   
     return 0;
   }
   
   
   /*
    *-------------------------------
    * Functions called by the solver
    *-------------------------------
    */
   
   /* f routine to compute the ODE RHS function f(t,y). */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
   {
     /* clear out ydot (to be careful) */
     N_VConst(0.0, ydot);
   
     /* problem data */
     UserData udata = (UserData) user_data;
   
     /* shortcuts to number of intervals, background values */
     long int N  = udata->N;
     realtype a  = udata->a;
     realtype b  = udata->b;
     realtype ep = udata->ep;
     realtype du = udata->du;
     realtype dv = udata->dv;
     realtype dw = udata->dw;
     realtype dx = udata->dx;
   
     /* access data arrays */
     realtype *Ydata = N_VGetArrayPointer(y);
     realtype *dYdata = N_VGetArrayPointer(ydot);
   
     /* iterate over domain, computing all equations */
     realtype uconst = du/dx/dx;
     realtype vconst = dv/dx/dx;
     realtype wconst = dw/dx/dx;
     realtype u, ul, ur, v, vl, vr, w, wl, wr;
     long int i;
     for (i=1; i<N-1; i++) {
   
       /* set shortcuts */
       u = Ydata[IDX(i,0)];  ul = Ydata[IDX(i-1,0)];  ur = Ydata[IDX(i+1,0)];
       v = Ydata[IDX(i,1)];  vl = Ydata[IDX(i-1,1)];  vr = Ydata[IDX(i+1,1)];
       w = Ydata[IDX(i,2)];  wl = Ydata[IDX(i-1,2)];  wr = Ydata[IDX(i+1,2)];
   
       /* u_t = du*u_xx + a - (w+1)*u + v*u^2 */
       dYdata[IDX(i,0)] = (ul - RCONST(2.0)*u + ur)*uconst + a - (w+RCONST(1.0))*u + v*u*u;
   
       /* v_t = dv*v_xx + w*u - v*u^2 */
       dYdata[IDX(i,1)] = (vl - RCONST(2.0)*v + vr)*vconst + w*u - v*u*u;
   
       /* w_t = dw*w_xx + (b-w)/ep - w*u */
       dYdata[IDX(i,2)] = (wl - RCONST(2.0)*w + wr)*wconst + (b-w)/ep - w*u;
   
     }
   
     /* enforce stationary boundaries */
     dYdata[IDX(0,0)]   = dYdata[IDX(0,1)]   = dYdata[IDX(0,2)]   = 0.0;
     dYdata[IDX(N-1,0)] = dYdata[IDX(N-1,1)] = dYdata[IDX(N-1,2)] = 0.0;
   
     return 0;
   }
   
   
   /* Jacobian routine to compute J(t,y) = df/dy. */
   static int Jac(long int M, long int mu, long int ml,
                  realtype t, N_Vector y, N_Vector fy, 
                  DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
   {
     /* clear out Jacobian (to be careful) */
     SetToZero(J);
   
     /* problem data */
     UserData udata = (UserData) user_data;
   
     /* Fill in the Laplace matrix */
     if (LaplaceMatrix(RCONST(1.0), J, udata)) {
       printf("Jacobian calculation error in calling LaplaceMatrix!\n");
       return 1;
     }
   
     /* Add in the Jacobian of the reaction terms matrix */
     if (ReactionJac(RCONST(1.0), y, J, udata)) {
       printf("Jacobian calculation error in calling ReactionJac!\n");
       return 1;
     }
   
     return 0;
   }
   
   
   
   /*-------------------------------
    * Private helper functions
    *-------------------------------*/
   
   /* Routine to compute the stiffness matrix from (L*y), scaled by the factor c.
      We add the result into Jac and do not erase what was already there */
   static int LaplaceMatrix(realtype c, DlsMat Jac, UserData udata)
   {
     /* shortcut to number of intervals */
     long int N = udata->N;
   
     /* set shortcuts */
     long int i;
     realtype dx = udata->dx;
     
     /* iterate over intervals, filling in Jacobian entries */
     for (i=1; i<N-1; i++) {
   
       /* Jacobian of (L*y) at this node */
       BAND_ELEM(Jac,IDX(i,0),IDX(i-1,0)) += c*udata->du/dx/dx;
       BAND_ELEM(Jac,IDX(i,1),IDX(i-1,1)) += c*udata->dv/dx/dx;
       BAND_ELEM(Jac,IDX(i,2),IDX(i-1,2)) += c*udata->dw/dx/dx;
       BAND_ELEM(Jac,IDX(i,0),IDX(i,0)) += -c*RCONST(2.0)*udata->du/dx/dx;
       BAND_ELEM(Jac,IDX(i,1),IDX(i,1)) += -c*RCONST(2.0)*udata->dv/dx/dx;
       BAND_ELEM(Jac,IDX(i,2),IDX(i,2)) += -c*RCONST(2.0)*udata->dw/dx/dx;
       BAND_ELEM(Jac,IDX(i,0),IDX(i+1,0)) += c*udata->du/dx/dx;
       BAND_ELEM(Jac,IDX(i,1),IDX(i+1,1)) += c*udata->dv/dx/dx;
       BAND_ELEM(Jac,IDX(i,2),IDX(i+1,2)) += c*udata->dw/dx/dx;
     }
   
     return 0;
   }
   
   
   
   /* Routine to compute the Jacobian matrix from R(y), scaled by the factor c.
      We add the result into Jac and do not erase what was already there */
   static int ReactionJac(realtype c, N_Vector y, DlsMat Jac, UserData udata)
   {
   
     /* shortcuts to number of intervals, background values */
     long int N  = udata->N;
   
     /* access data arrays */
     realtype *Ydata = N_VGetArrayPointer(y);
   
     /* set shortcuts */
     long int i;
     realtype u, v, w;
     realtype ep = udata->ep;
     
     /* iterate over nodes, filling in Jacobian entries */
     for (i=1; i<N-1; i++) {
   
       /* set nodal value shortcuts (shifted index due to start at first interior node) */
       u = Ydata[IDX(i,0)];
       v = Ydata[IDX(i,1)];
       w = Ydata[IDX(i,2)];
   
       /* all vars wrt u */
       BAND_ELEM(Jac,IDX(i,0),IDX(i,0)) += c*(RCONST(2.0)*u*v-(w+RCONST(1.0)));
       BAND_ELEM(Jac,IDX(i,1),IDX(i,0)) += c*(w - RCONST(2.0)*u*v);
       BAND_ELEM(Jac,IDX(i,2),IDX(i,0)) += c*(-w);
   
       /* all vars wrt v */
       BAND_ELEM(Jac,IDX(i,0),IDX(i,1)) += c*(u*u);
       BAND_ELEM(Jac,IDX(i,1),IDX(i,1)) += c*(-u*u);
   
       /* all vars wrt w */
       BAND_ELEM(Jac,IDX(i,0),IDX(i,2)) += c*(-u);
       BAND_ELEM(Jac,IDX(i,1),IDX(i,2)) += c*(u);
       BAND_ELEM(Jac,IDX(i,2),IDX(i,2)) += c*(-RCONST(1.0)/ep - u);
   
     }
   
     return 0;
   }


Solutions
^^^^^^^^^^^^



PDE system example with iterative linear solver (ark_heat1D)
--------------------------------------------------------------

ODE system
^^^^^^^^^^^^

This example simulates a simple 1D heat equation,

.. math::

   \frac{\partial u}{\partial t} = k \frac{\partial^2 u}{\partial x^2} + f,

for :math:`t \in [0, 10]`, and :math:`x in [0, 1]`, with initial
condition :math:`u(0,x) = 0`, Dirichlet boundary conditions,

.. math::

   \frac{\partial u}{\partial t}(t,0) = \frac{\partial u}{\partial t}(t,1) = 0,

and a point-source heating term, 

.. math::

   f(t,x) = \begin{cases} 1 & \text{if} x=\frac12, \\
                          0 & \text{otherwise}. \end{cases}
 

Numerical method
^^^^^^^^^^^^^^^^^

The spatial derivatives are computed using second-order 
centered differences, with the data distributed over :math:`N` points
on a uniform spatial grid. 

The number of spatial points :math:`N` and the heat conductivity
parameter :math:`k`, as well as the desired relative and absolute
solver tolerances, are provided in the input file ``input_heat1D.txt``.
 
This program solves the problem with a DIRK method, utilizing a Newton
iteration with the PCG iterative linear solver, and a user-supplied
Jacobian-vector product routine.

100 outputs are printed at equal intervals, and run statistics are
printed at the end. 


Routines
^^^^^^^^^^^^^^^^^

We reproduce the relevant aspects of the ``main()`` routine and
auxiliary functions here for explanatory purposes (see the in-line
comments for details; error-checking has been removed for brevity):

.. code-block:: c

   /* Header files */
   #include <stdio.h>
   #include <stdlib.h>
   #include <math.h>
   #include <arkode/arkode.h>
   #include <nvector/nvector_serial.h>
   #include <arkode/arkode_pcg.h>
   #include <sundials/sundials_types.h>
   
   /* user data structure */
   typedef struct {  
     long int N;    /* number of intervals     */
     realtype dx;   /* mesh spacing            */
     realtype k;    /* diffusion coefficient   */
   } *UserData;
   
   /* User-supplied Functions Called by the Solver */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
   static int Jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
   	       N_Vector fy, void *user_data, N_Vector tmp);

   
   /* Main Program */
   int main() {
   
     /* general problem parameters */
     realtype T0 = RCONST(0.0);
     realtype Tf = RCONST(1.0);
     int Nt = 10;
     UserData udata = NULL;
     realtype *data;
     long int N, i;
   
     /* general problem variables */
     int flag;
     N_Vector y = NULL;
     void *arkode_mem = NULL;
   
     /* allocate udata structure */
     udata = (UserData) malloc(sizeof(*udata));
   
     /* read problem parameter and tolerances from input file:
        N - number of spatial discretization points
        k - diffusion coefficient */
     double k;
     FILE *FID;
     FID = fopen("input_heat1D.txt","r");
     flag = fscanf(FID,"  N = %li\n", &N);
     flag = fscanf(FID,"  k = %lf\n", &k);
     fclose(FID);
   
     /* store the inputs in the UserData structure */
     udata->N = N;
     udata->k = k;
   
     /* open solver diagnostics output file for writing */
     FILE *DFID;
     DFID=fopen("diags_ark_heat1D.txt","w");
     
     /* Initial problem output */
     printf("\n1D Heat PDE test problem:\n");
     printf("  N = %li\n", udata->N);
     printf("  diffusion coefficient:  k = %g\n", udata->k);
   
     /* Create serial vector of length N for initial condition */
     y = N_VNew_Serial(N);
   
     /* set spatial mesh spacing */
     udata->dx = RCONST(1.0)/(1.0*N-1.0);
   
     /* output mesh to disk */
     FID=fopen("heat_mesh.txt","w");
     for (i=0; i<N; i++)  fprintf(FID,"  %.16e\n", udata->dx*i);
     fclose(FID);
     
     /* Set initial conditions into y */
     N_VConst(0.0, y);
   
   
     /* Call ARKodeCreate to create the solver memory and specify the 
        Backward Differentiation Formula and the use of a Newton iteration */
     arkode_mem = ARKodeCreate();
   
     /* Call ARKodeInit to initialize the integrator memory and specify the
        user's right hand side function in y'=f(t,y), the inital time T0, and
        the initial dependent variable vector y */
     flag = ARKodeInit(arkode_mem, NULL, f, T0, y);
   
     /* Call init_from_file helper routine to read and set solver parameters */
     realtype rtol = 1.e-6;
     realtype atol = 1.e-10;
   
     /* Call ARKodeSetUserData to pass rdata to user functions */
     flag = ARKodeSetUserData(arkode_mem, (void *) udata);
   
     /* Call ARKodeSetDiagnostics to set diagnostics output file pointer */
     flag = ARKodeSetDiagnostics(arkode_mem, DFID);
   
     /* Call ARKodeSetMaxNumSteps to increase default (for testing) */
     flag = ARKodeSetMaxNumSteps(arkode_mem, 10000);
   
     /* Call ARKodeSStolerances to specify the scalar relative and absolute
        tolerances */
     flag = ARKodeSStolerances(arkode_mem, rtol, atol);
   
     /* Specify the linear solver */
     flag = ARKPcg(arkode_mem, 0, N);
   
     /* Set the Jacobian routine to Jac (user-supplied) */
     flag = ARKSpilsSetJacTimesVecFn(arkode_mem, Jac);
   
     /* Open output stream for results, access data arrays */
     FILE *UFID=fopen("heat.txt","w");
     data = N_VGetArrayPointer(y);
   
     /* output initial condition to disk */
     for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[i]);
     fprintf(UFID,"\n");
   
     /* In loop, call ARKode, print results, and test for error.
        Break out of loop when the final output time has been reached */
     realtype t  = T0;
     realtype dTout = Tf/Nt;
     realtype tout = T0+dTout;
     realtype u;
     printf("        t      ||u||_rms\n");
     printf("   -------------------------\n");
     u = N_VDotProd(y,y);
     u = sqrt(u/N);
     printf("  %10.6f  %10.6f\n", t, u);
     int iout;
     for (iout=0; iout<Nt; iout++) {
   
       flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);
       u = N_VDotProd(y,y);
       u = sqrt(u/N);
       printf("  %10.6f  %10.6f\n", t, u);
   
       if (flag >= 0) {
         tout += dTout;
         tout = (tout > Tf) ? Tf : tout;
       }
   
       /* output results to disk */
       for (i=0; i<N; i++)  fprintf(UFID," %.16e", data[i]);
       fprintf(UFID,"\n");
     }
     printf("   -------------------------\n");
     fclose(UFID);
       
   
     /* Print some final statistics */
     long int nst, nst_a, nfe, nfi, nsetups, nli, nJv, nlcf, nni, ncfn, netf;
     flag = ARKodeGetNumSteps(arkode_mem, &nst);
     flag = ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
     flag = ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
     flag = ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
     flag = ARKodeGetNumErrTestFails(arkode_mem, &netf);
     flag = ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
     flag = ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
     flag = ARKSpilsGetNumLinIters(arkode_mem, &nli);
     flag = ARKSpilsGetNumJtimesEvals(arkode_mem, &nJv);
     flag = ARKSpilsGetNumConvFails(arkode_mem, &nlcf);
   
     printf("\nFinal Solver Statistics:\n");
     printf("   Internal solver steps = %li (attempted = %li)\n", 
   	 nst, nst_a);
     printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
     printf("   Total linear solver setups = %li\n", nsetups);
     printf("   Total linear iterations = %li\n", nli);
     printf("   Total number of Jacobian-vector products = %li\n", nJv);
     printf("   Total number of linear solver convergence failures = %li\n", nlcf);
     printf("   Total number of Newton iterations = %li\n", nni);
     printf("   Total number of nonlinear solver convergence failures = %li\n", ncfn);
     printf("   Total number of error test failures = %li\n", netf);
   
     /* Free vectors */
     N_VDestroy_Serial(y);
   
     /* Free user data */
     free(udata);
   
     /* Free integrator memory */
     ARKodeFree(&arkode_mem);
   
     /* close solver diagnostics output file */
     fclose(DFID);
   
     return 0;
   }
   
   
   /*--------------------------------
    * Functions called by the solver
    *--------------------------------*/
   
   /* f routine to compute the ODE RHS function f(t,y). */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
   {
     /* clear out ydot (to be careful) */
     N_VConst(0.0, ydot);
   
     /* problem data */
     UserData udata = (UserData) user_data;
   
     /* shortcuts to number of intervals, background values */
     long int N  = udata->N;
     realtype k  = udata->k;
     realtype dx = udata->dx;
   
     /* access data arrays */
     realtype *Y = N_VGetArrayPointer(y);
     realtype *Ydot = N_VGetArrayPointer(ydot);
   
     /* iterate over domain, computing all equations */
     realtype c1 = k/dx/dx;
     realtype c2 = -RCONST(2.0)*k/dx/dx;
     long int i;
     long int isource = N/2;
     Ydot[0] = 0.0;                 /* left boundary condition */
     for (i=1; i<N-1; i++)
       Ydot[i] = c1*Y[i-1] + c2*Y[i] + c1*Y[i+1];
     Ydot[N-1] = 0.0;               /* right boundary condition */
     Ydot[isource] += 1.0;          /* f */
   
     return 0;
   }
   
   
   
   /* Jacobian routine to compute J(t,y) = df/dy. */
   static int Jac(N_Vector v, N_Vector Jv, realtype t, N_Vector y, 
   	       N_Vector fy, void *user_data, N_Vector tmp)
   {
     /* clear out result (to be careful) */
     N_VConst(0.0, Jv);
   
     /* shortcuts to number of intervals, background values */
     UserData udata = (UserData) user_data;
     long int N  = udata->N;
     realtype k  = udata->k;
     realtype dx = udata->dx;
   
     /* access data arrays */
     realtype *V = N_VGetArrayPointer(v);
     realtype *JV = N_VGetArrayPointer(Jv);
   
     /* iterate over domain, computing all Jacobian-vector products */
     realtype c1 = k/dx/dx;
     realtype c2 = -RCONST(2.0)*k/dx/dx;
     long int i;
     JV[0] = 0.0;
     for (i=1; i<N-1; i++)
       JV[i] = c1*V[i-1] + c2*V[i] + c1*V[i+1];
     JV[N-1] = 0.0;
   
     return 0;
   }


Solutions
^^^^^^^^^^^^


