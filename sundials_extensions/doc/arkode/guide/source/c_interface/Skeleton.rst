:tocdepth: 3


.. _CInterface.Skeleton:

A skeleton of the user's main program
============================================

The following is a skeleton of the user's main program (or calling
program) for the integration of an ODE IVP.  Some steps are
independent of the NVECTOR implementation used; where this is not
the case, usage specifications are given for the two implementations
provided with ARKode: steps marked [P] correspond to
NVECTOR_PARALLEL, while steps marked [S] correspond to
NVECTOR_SERIAL. 

1. [P] Initialize MPI 
 
   Call ``MPI_Init`` to initialize MPI if used by the user's program.

2. Set problem dimensions

   [S] Set ``N``, the problem size :math:`N`.

   [P] Set ``Nlocal``, the local vector length (the sub-vector length
   for this process); ``N``, the global vector length (the problem size
   :math:`n`, and the sum of all the values of ``Nlocal``); and the
   active set of processes. 

3. Set vector of initial values

   To set the vector ``y0`` of initial values, use the appropriate
   functions defined by the particular NVECTOR implementation.  If a
   ``realtype`` array ``ydata`` containing the initial values of :math:`y`
   already exists, then make the call: 

   [S] ``y0 = N_VMake_Serial(N, ydata);``

   [P] ``y0 = N_VMake_Parallel(comm, Nlocal, N, ydata);``

   Otherwise, make the call: 

   [S] ``y0 = N_VNew_Serial(N);``

   [P] ``y0 = N_VNew_Parallel(comm, Nlocal, N);``

   and load initial values into the structure defined by: 

   [S] ``NV_DATA_S(y0)``

   [P] ``NV_DATA_P(y0)``

   Here ``comm`` is the MPI communicator containing the set of active
   processes to be used (may be ``MPI_COMM_WORLD``). 

4. Create ARKode object

   Call ``arkode_mem = ARKodeCreate()`` to create the ARKode memory
   block. :c:func:`ARKodeCreate()` returns a pointer to the ARKode memory
   structure. See the section :ref:`CInterface.Initialization` for
   details.  

5. Initialize ARKode solver

   Call :c:func:`ARKodeInit()` to provide required problem specifications,
   allocate internal memory for ARKode, and initialize
   ARKode. :c:func:`ARKodeInit()` returns a flag, the value of which indicates
   either success or an illegal argument value. See the section
   :ref:`CInterface.Initialization` for details. 

6. Specify integration tolerances

   Call :c:func:`ARKodeSStolerances()` or :c:func:`ARKodeSVtolerances()` to
   specify either a scalar relative tolerance and scalar absolute
   tolerance, or a scalar relative tolerance and a vector of absolute
   tolerances, respectively. Alternatively, call :c:func:`ARKodeWFtolerances()`
   to specify a function which sets directly the weights used in
   evaluating WRMS vector norms. See the section
   :ref:`CInterface.Tolerances` for details. 

7. Set optional inputs 

   Call ``ARKodeSet*`` functions to change any optional inputs that
   control the behavior of ARKode from their default values. See
   the section :ref:`CInterface.OptionalInputs` for details. 

8. Attach linear solver module

   If an implicit solve is required, initialize the linear solver
   module with one of the following calls (for details see the section
   :ref:`CInterface.LinearSolvers`):  

   [S] ``ier = ARKDense(...);``

   [S] ``ier = ARKBand(...);``

   [S] ``ier = ARKLapackDense(...);`` 

   [S] ``ier = ARKLapackBand(...);``

   ``ier = ARKSpgmr(...);``

   ``ier = ARKSpbcg(...);``

   ``ier = ARKSptfqmr(...);``

   ``ier = ARKPcg(...);``

9. Set linear solver optional inputs 

   Call ``ARK*Set*`` functions from the selected linear solver module to
   change optional inputs specific to that linear solver. See the section
   :ref:`CInterface.OptionalInputs` for details. 

10. Specify rootfinding problem

    Optionally, call :c:func:`ARKodeRootInit()` to initialize a rootfinding
    problem to be solved during the integration of the ODE system. See
    the section :ref:`CInterface.RootFinding` for general details, and
    the section :ref:`CInterface.OptionalInputs` for relevant optional
    input calls. 

11. Advance solution in time

    For each point at which output is desired, call 

    ``ier = ARKode(arkode_mem, tout, yout, &tret, itask)``

    Here, :c:func:`ARKode()` requires that ``itask``
    specify the return mode. The vector ``y`` (which can be the same as
    the vector ``y0`` above) will contain :math:`y(t)`. See the section
    :ref:`CInterface.Integration` for details. 

12. Get optional outputs 

    Call ``ARK*Get*`` functions to obtain optional output. See
    the section :ref:`CInterface.OptionalInputs` for details.  

13. Deallocate memory for solution vector 

    Upon completion of the integration, deallocate memory for the
    vector ``y`` by calling the destructor function defined by the
    NVECTOR implementation:

    [S] ``N_VDestroy_Serial(y);``

    [P] ``N_VDestroy_Parallel(y);`` 

14. Free solver memory 

    Call ``ARKodeFree(&arkode_mem)`` to free the memory allocated for ARKode. 

15. [P] Finalize MPI 

    Call ``MPI_Finalize`` to terminate MPI.
