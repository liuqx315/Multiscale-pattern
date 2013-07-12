:tocdepth: 3

.. _ark_analytic:

ark_analytic
====================================

This is a very simple C example that merely shows how to use the
ARKode solver interface.  

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
----------------

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
--------

We reproduce the relevant aspects of the ``main()`` routine and
auxiliary functions here for explanatory purposes (see the in-line
comments for details; error-checking has been removed for brevity).




Include files and function prototypes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: c

   #include <stdio.h>
   #include <math.h>
   #include <arkode/arkode.h>           /* prototypes for ARKODE fcts., consts. */
   #include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
   #include <arkode/arkode_dense.h>     /* prototype for ARKDense solver */
   #include <sundials/sundials_dense.h> /* definitions of DlsMat and DENSE_ELEM */
   #include <sundials/sundials_types.h> /* definition of type 'realtype' */
   
   /* User-supplied functions called by the solver */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);



main()
^^^^^^^

.. code-block:: c

   int main()
   {
     /* general problem parameters */
     realtype T0 = RCONST(0.0);      /* initial time */
     realtype Tf = RCONST(10.0);     /* final time */
     realtype dTout = RCONST(1.0);   /* time between outputs */
     long int NEQ = 1;               /* number of dependent vars. */
   
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
     FID = fopen("input_analytic.txt","r");
     fscanf(FID,"  lamda = %lf\n",  &lamda_);
     fscanf(FID,"  reltol = %lf\n", &reltol_);
     fscanf(FID,"  abstol = %lf\n", &abstol_);
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

     /* Initialize data structures */
     y = N_VNew_Serial(NEQ);          /* Create serial vector for solution */
     NV_Ith_S(y,0) = 0.0;             /* Specify initial condition */
     arkode_mem = ARKodeCreate();     /* Create the solver memory */
     
     /* Call ARKodeInit to initialize the integrator memory and specify the
        hand-side side function in y'=f(t,y), the inital time T0, and
        the initial dependent variable vector y.  Note: since this
	problem is fully implicit, we set f_E to NULL and f_I to f. */
     ARKodeInit(arkode_mem, NULL, f, T0, y);
   
     /* Set routines */
     ARKodeSetUserData(arkode_mem, (void *) &lamda);  /* Pass lamda to user functions */
     ARKodeSStolerances(arkode_mem, reltol, abstol);  /* Specify tolerances */

     /* Linear solver specification */
     ARKDense(arkode_mem, NEQ);                       /* Specify dense linear solver */
     ARKDlsSetDenseJacFn(arkode_mem, Jac);            /* Set Jacobian routine */
   
     /* Main time-stepping loop: calls ARKode to perform the integration, then
        prints results.  Stops when the final time has been reached */
     realtype t = T0;
     realtype tout = T0+dTout;
     printf("        t           u\n");
     printf("   ---------------------\n");
     while (Tf - t > 1.0e-15) {
   
       flag = ARKode(arkode_mem, tout, y, &t, ARK_NORMAL);      /* call integrator */
       printf("  %10.6f  %10.6f\n", t, NV_Ith_S(y,0));          /* access/print solution */
       if (flag >= 0) {                                         /* successful solve: update time */
         tout += dTout;
         tout = (tout > Tf) ? Tf : tout;
       } else {                                                 /* unsuccessful solve: break */
         fprintf(stderr,"Solver failure, stopping integration\n");
         break;
       }
     }
     printf("   ---------------------\n");
   
     /* Get/print some final statistics on how the solve progressed */
     long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
     ARKodeGetNumSteps(arkode_mem, &nst);
     ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
     ARKodeGetNumRhsEvals(arkode_mem, &nfe, &nfi);
     ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
     ARKodeGetNumErrTestFails(arkode_mem, &netf);
     ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
     ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
     ARKDlsGetNumJacEvals(arkode_mem, &nje);
     ARKDlsGetNumRhsEvals(arkode_mem, &nfeLS);
   
     printf("\nFinal Solver Statistics:\n");
     printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
     printf("   Total RHS evals:  Fe = %li,  Fi = %li\n", nfe, nfi);
     printf("   Total linear solver setups = %li\n", nsetups);
     printf("   Total RHS evals for setting up the linear system = %li\n", nfeLS);
     printf("   Total number of Jacobian evaluations = %li\n", nje);
     printf("   Total number of Newton iterations = %li\n", nni);
     printf("   Total number of linear solver convergence failures = %li\n", ncfn);
     printf("   Total number of error test failures = %li\n\n", netf);
   
     /* Clean up and return with successful completion */
     N_VDestroy_Serial(y);     /* Free y vector */
     ARKodeFree(&arkode_mem);  /* Free integrator memory */
     return 0;
   }



f() 
^^^^

.. code-block:: c

   /* f routine to compute the ODE RHS function f(t,y). */
   static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
   {
     realtype *rdata = (realtype *) user_data;   /* cast user_data to realtype */
     realtype lamda = rdata[0];                  /* set shortcut for stiffness parameter */
     realtype u = NV_Ith_S(y,0);                 /* access current solution value */
   
     /* fill in the RHS function: "NV_Ith_S" accesses the 0th entry of ydot */
     NV_Ith_S(ydot,0) = lamda*u + 1.0/(1.0+t*t) - lamda*atan(t);

     return 0;                                   /* return with success */
   }
   



Jac()
^^^^^^^

.. code-block:: c

   /* Jacobian routine to compute J(t,y) = df/dy. */
   static int Jac(long int N, realtype t,
                  N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                  N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
   {
     realtype *rdata = (realtype *) user_data;   /* cast user_data to realtype */
     realtype lamda = rdata[0];                  /* set shortcut for stiffness parameter */

     /* Fill in Jacobian of f: "DENSE_ELEM" accesses the (0,0) entry of J */
     DENSE_ELEM(J,0,0) = lamda;
   
     return 0;                                   /* return with success */
   }



Solutions
---------

This problem is included both as a simple example, but also because it
has an analytical solution, :math:`y(t) = \arctan(t)`.  As seen in the
plots below, the computed solution tracks the analytical solution
quite well (left), and results in errors below those specified by the input
error tolerances (right).

.. image:: figs/plot-ark_analytic.png
   :width: 45 %
.. image:: figs/plot-ark_analytic_error.png
   :width: 45 %
