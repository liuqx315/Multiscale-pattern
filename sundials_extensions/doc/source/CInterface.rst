.. _CInterface:

Using ARKode for C and C++ Applications
=======================================

This chapter is concerned with the use of ARKode for the solution
of initial value problems (IVPs) in a C or C++ language setting.  The
following sections treat the header files and the layout of the user's
main program, and provide descriptions of the ARKode user-callable 
functions and user-supplied functions. 

The sample programs described in the companion document
[R2013]_ may also be helpful. Those codes may be used as
templates (with the removal of some lines used in testing) and are
included in the ARKode package. 

Users with applications written in Fortran should see the chapter
:ref:`FortranInterface`, that describes the Fortran/C interface
module. 

The user should be aware that not all linear solver modules are
compatible with all NVECTOR implementations.  For example,
NVECTOR_PARALLEL is not compatible with the direct dense or direct
band linear solvers since these linear solver modules need to form the
complete system Jacobian on a single processor.  Specifically, the
following ARKode modules can only be used with NVECTOR_SERIAL:
ARKDENSE, ARKBAND (using either the internal or the LAPACK
implementation) and ARKBANDPRE. Also, the preconditioner module
ARKBBDPRE can only be used with NVECTOR_PARALLEL. 

ARKode uses various constants for both input and output. These are
defined as needed in this chapter, but for convenience are also listed
separately in the section :ref:`Constants`. 


.. _CInterface.Headers:

Access to library and header files
----------------------------------

At this point, it is assumed that the installation of ARKode,
following the procedure described in the section :ref:`Installation`,
has been completed successfully. 

Regardless of where the user's application program resides, its
associated compilation and load commands must make reference to the
appropriate locations for the library and header files required by
ARKode. The relevant library files are 

- ``libdir/libsundialsarkode.lib``,
- ``libdir/libsundials_nvec*.lib`` (one or two files), 

where the file extension ``.lib`` is typically ``.so`` for shared
libraries and ``.a`` for static libraries. The relevant header files are
located in the subdirectories 

- ``incdir/include/arkode``
- ``incdir/include/sundials``
- ``incdir/include/nvector``

The directories ``libdir`` and ``incdir`` are the install library and
include directories, respectively.  For a default installation, these
are ``instdir/lib`` and ``instdir/include``, respectively, where ``instdir``
is the directory where SUNDIALS was installed (see the section
:ref:`Installation`).


.. _CInterface.DataTypes:

Data Types
----------

The ``sundials_types.h`` file contains the definition of the type
``realtype``, which is used by the SUNDIALS solvers for all
floating-point data.  The type ":index:`realtype`" can be ``float``,
``double``, or ``long double``, depending on how SUNDIALS was
installed (with the default being ``double``). The user can change the
precision of the SUNDIALS solvers floating-point arithmetic at the
configuration stage (see the section :ref:`Installation`). 

Additionally, based on the current precision, ``sundials_types.h``
defines the values :index:`BIG_REAL` to be the largest value
representable as a ``realtype``, :index:`SMALL_REAL` to be the
smallest positive value representable as a ``realtype``, and
:index:`UNIT_ROUNDOFF` to be the difference between 1.0 and the
minimum ``realtype`` greater than 1.0.  

Within SUNDIALS, real constants may be set to have the appropriate
precision by way of a macro called :index:`RCONST`.  It is this macro
that needs the ability to branch on the definition ``realtype``.  In
ANSI C, a floating-point constant with no suffix is stored as a
``double``. Placing the suffix "F" at the end of a floating point
constant makes it a ``float``, whereas using the suffix "L" makes it a
``long double``. For example,

.. code-block:: c

   #define A 1.0 
   #define B 1.0F 
   #define C 1.0L

defines ``A`` to be a ``double`` constant equal to 1.0, ``B`` to be a
``float`` constant equal to 1.0, and ``C`` to be a ``long double`` constant
equal to 1.0.  The macro call ``RCONST(1.0)`` automatically expands to
1.0 if ``realtype`` is ``double``, to ``1.0F`` if ``realtype`` is ``float``, or
to ``1.0L`` if ``realtype`` is ``long double``. SUNDIALS uses the ``RCONST``
macro internally to declare all of its floating-point constants. 

A user program which uses the type ``realtype`` and the ``RCONST`` macro
to handle floating-point constants is precision-independent except for
any calls to precision-specific standard math library functions (Our
example programs use both ``realtype`` and ``RCONST``).  Users can,
however, use the types ``double``, ``float``, or ``long double`` in their
code (assuming that this usage is consistent with the ``typedef`` for
``realtype``).  Thus, a previously existing piece of ANSI C code can use
SUNDIALS without modifying the code to use ``realtype``, so long as
the SUNDIALS libraries have been compiled using the same precision
(for details see the section :ref:`Installation`). 

SUNDIALS also defines a type ":index:`booleantype`", that can have
values ``TRUE`` and ``FALSE``, which is used for logic arguments
within the library.



Header Files
------------

The calling program must include several header files so that various
macros and data types can be used. The header file that is always
required is: 

- ``arkode.h``, the main header file for ARKode, which defines the
  several types and various constants, and includes function
  prototypes. 

Note that ``arkode.h`` includes ``sundials_types.h`` directly, which
defines the types ``realtype`` and ``booleantype`` and the
constants ``FALSE`` and ``TRUE``, so a user program does not need to
include ``sundials_types.h`` directly. 

The calling program must also include an NVECTOR implementation
header file (see the chapter :ref:`NVectors` for details).  For the two
NVECTOR implementations that are included in the ARKode package, the
corresponding header files are: 

* ``nvector_serial.h``, which defines the serial implementation
  NVECTOR_SERIAL; 
* ``nvector_parallel.h``, which defines the parallel (MPI)
  implementation, NVECTOR_PARALLEL.

Note that both these files in turn include the header file
``sundials_nvector.h`` which defines the abstract ``N_Vector`` data
type.

Finally, if the user includes a non-trivial implicit component to their
ODE system (and hence requires a Newton solver for the resulting
nonlinear systems of equations), then a linear solver module header
file will be required. The header files corresponding to the various
linear solvers availble for use with ARKode are: 

- ``arkode_dense.h``, which is used with the dense direct linear solver; 
- ``arkode_band.h``, which is used with the band direct linear solver;
- ``arkode_lapack.h``, which is used with LAPACK implementations of dense
  or band direct linear solvers; 
- ``arkode_spgmr.h``, which is used with the scaled, preconditioned GMRES
  Krylov linear solver SPGMR;
- ``arkode_spbcgs.h``, which is used with the scaled, preconditioned
  Bi-CGStab Krylov linear solver SPBCG;
- ``arkode_sptfqmr.h``, which is used with the scaled, preconditioned
  TFQMR Krylov solver SPTFQMR.

The header files for the dense and banded linear solvers (both
internal and LAPACK) include the file ``arkode_direct.h``, which defines
common functions.  This in turn includes a file (``sundials_direct.h``)
which defines the matrix type for these direct linear solvers
(``DlsMat``), as well as various functions and macros for acting on and
accessing entries of such matrices. 

The header files for the Krylov iterative solvers each include
``arkode_spils.h`` which defines common functions and which in turn
includes a header file (``sundials_iterative.h``) which enumerates the
preconditioning type and the choices for the Gram-Schmidt process (for
the SPGMR solver). 

Other headers may be needed, according to the choice of
preconditioner, etc.  For example, in the ``arkDiurnal_kry_p.c`` example
(see [R2013]_), preconditioning is done with a block-diagonal
matrix.  For this, even though the :c:func:`ARKSpgmr()` linear solver
is used, the header ``sundials_dense.h`` is included for access to the
underlying generic dense linear solver that is used for preconditioning.



.. _CInterface.Skeleton:

A skeleton of the user's main program
-------------------------------------

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



User-callable functions
-----------------------

This section describes the ARKode functions that are called by the
user to setup and then solve an IVP. Some of these are
required. However, starting with the section
:ref:`CInterface.OptionalInputs`, the functions listed involve
optional inputs/outputs or restarting, and those paragraphs may be
skipped for a casual use of ARKode. In any 
case, refer to the ssection :ref:`CInterface.Skeleton` for the correct
order of these calls. 

On an error, each user-callable function returns a negative value and
sends an error message to the error handler routine, which prints the
message on ``stderr`` by default. However, the user can set a file as
error output or can provide his own error handler function
(see the section :ref:`CInterface.OptionalInputs` for details).



.. _CInterface.Initialization:

ARKode initialization and deallocation functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: void *ARKodeCreate()

   The function ARKodeCreate creates an internal
   memory block for a problem to be solved by ARKode.

   **Arguments:**  None

   **Return value:**  If successful, a pointer to initialized problem memory
   of type ``void *``, to be passed to :c:func:`ARKodeInit()`.
   If unsuccessful, a ``NULL`` pointer, and an error
   message will be printed to ``stderr``.


.. c:function:: int ARKodeInit(void *arkode_mem, ARKRhsFn fe, ARKRhsFn fi, realtype t0, realtype y0)

   The function ARKodeInit allocates and initializes
   memory for a problem to to be solved by ARKode.

   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block
        (that was returned by :c:func:`ARKodeCreate()`)
      * `fe` -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the explicit portion of the right-hand side function in 
        :math:`\dot{y} = f_E(t,y) + f_I(t,y)` 
      * `fi` -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the implicit portion of the right-hand side function in 
        :math:`\dot{y} = f_E(t,y) + f_I(t,y)`
      * `t0` -- the initial value of :math:`t`
      * `y0` -- the initial condition vector :math:`y(t_0)`

   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL  if the ARKode memory was ``NULL``
      * ARK_MEM_FAIL  if a memory allocation failed
      * ARK_ILL_INPUT if an argument has an illegal value.


.. c:function:: void ARKodeFree(void *arkode_mem)

   The function ARKodeFree frees the problem memory
   `arkode_mem` allocated by :c:func:`ARKodeCreate()` and :c:func:`ARKodeInit()`.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
   
   **Return value:**  None



.. _CInterface.Tolerances:

ARKode tolerance specification functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These functions specify the integration tolerances. One of them
**should** be called before the first call to :c:func:`ARKode()`; otherwise
default values of ``reltol = 1e-4`` and ``abstol = 1e-9`` will be
used, which may be entirely incorrect for a specific problem.

The tolerances ``reltol`` and ``abstol`` define a vector of error
weights, ``ewt``.  In the case of :c:func:`ARKodeSStolerances()`, this vector
has components 

.. code-block:: c

   ewt[i] = 1.0/(reltol*abs(y[i]) + abstol);

whereas in the case of :c:func:`ARKodeSVtolerances()` the vector components
are given by 

.. code-block:: c

   ewt[i] = 1.0/(reltol*abs(y[i]) + abstol[i]);

This vector is used in all error and convergence tests, which use a
weighted RMS norm on all error-like vectors v:

.. math::
    \|v\|_{WRMS} = \left( \frac{1}{n} \sum_{i=1}^n (v_i*ewt_i)^2 \right)^{1/2},

where :math:`n` is the problem dimension.

Alternatively, the user may supply a custom function to supply the
``ewt`` vector, through a call to :c:func:`ARKodeWFtolerances()`.



.. c:function:: int ARKodeSStolerances(void *arkode_mem, realtype reltol, realtype abstol)

   Specifies scalar relative and absolute tolerances.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `reltol` -- scalar relative tolerance
      * `abstol` -- scalar absolute tolerance
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL  if the ARKode memory was ``NULL``
      * ARK_NO_MALLOC  if the ARKode memory was not allocated by :c:func:`ARKodeInit()`
      * ARK_ILL_INPUT if an argument has an illegal value (e.g. a
        negative tolerance).



.. c:function:: int ARKodeSVtolerances(void *arkode_mem, realtype reltol, N_Vector abstol)

   Specifies a scalar relative tolerance and a 
   vector absolute tolerance (a potentially different absolute 
   tolerance for each vector component).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `reltol` -- scalar relative tolerance
      * `abstol` -- vector containing the absolute tolerances for each
        solution component
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL  if the ARKode memory was ````NULL````
      * ARK_NO_MALLOC  if the ARKode memory was not allocated by :c:func:`ARKodeInit()`
      * ARK_ILL_INPUT if an argument has an illegal value (e.g. a
        negative tolerance).



.. c:function:: int ARKodeWFtolerances(void *arkode_mem, ARKEwtFn efun)

   Specifies a user-supplied function `efun` to compute
   the error weight vector `ewt`.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `efun` -- the name of the C function (of type :c:func:`ARKEwtFn()`)
        that implements the error weight vector computation.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL  if the ARKode memory was ``NULL``
      * ARK_NO_MALLOC  if the ARKode memory was not allocated by :c:func:`ARKodeInit()`


General advice on the choice of tolerances
""""""""""""""""""""""""""""""""""""""""""

For many users, the appropriate choices for tolerance values in reltol
and abstol are a concern. The following pieces of advice are
relevant. 

(1) The scalar relative tolerance ``reltol`` is to be set to control
    relative errors. So a value of :math:`10^{-4}` means that errors
    are controlled to .01%. We do not recommend using ``reltol`` larger
    than :math:`10^{-3}`. On the other hand, ``reltol`` should not be so small
    that it is comparable to the unit roundoff of the machine
    arithmetic (generally around :math:`10^{-15}`). 

(2) The absolute tolerances ``abstol`` (whether scalar or vector) need
    to be set to control absolute errors when any components of the
    solution vector :math:`y` may be so small that pure relative error
    control is meaningless.  For example, if :math`y_i` starts at some
    nonzero value, but in time decays to zero, then pure relative
    error control on :math:`y_i` makes no sense (and is overly costly)
    after :math:`y_i` is below some noise level. Then ``abstol`` (if
    scalar) or ``abstol[i]`` (if a vector) needs to be set to that
    noise level. If the different components have different noise
    levels, then ``abstol`` should be a vector. See the example
    ``arkRoberts_dns.c`` in the ARKode package, and the discussion
    of it in the ARKode Examples document [R2013]_. In that
    problem, the three components vary betwen 0 and 1, and have
    different noise levels; hence the ``abstol`` vector. It is
    impossible to give any general advice on ``abstol`` values,
    because the appropriate noise levels are completely
    problem-dependent. The user or modeler hopefully has some idea as
    to what those noise levels are. 

(3) Finally, it is important to pick all the tolerance values
    conservately, because they control the error committed on each
    individual time step. The final (global) errors are an
    accumulation of those per-step errors, where that accumulation
    factor is problem-dependent.  A general rule of thumb is to reduce
    the tolerances by a factor of 10 from the actual desired limits on
    errors.  So if you want .01% relative accuracy (globally), a good
    choice for ``reltol`` is :math:`10^{-5}`.  But in any case, it is
    a good idea to do a few experiments with the tolerances to see how
    the computed solution values vary as tolerances are reduced.


Advice on controlling unphysical negative values
""""""""""""""""""""""""""""""""""""""""""""""""
In many applications, some components in the true solution are always
positive or non-negative, though at times very small.  In the
numerical solution, however, small negative (hence unphysical) values
can then occur. In most cases, these values are harmless, and simply
need to be controlled, not eliminated. The following pieces of advice
are relevant. 

(1) The best way to control the size of unwanted negative computed
    values is with tighter absolute tolerances.  Again this requires
    some knowledge of the noise level of these components, which may
    or may not be different for different components. Some
    experimentation may be needed. 

(2) If output plots or tables are being generated, and it is important
    to avoid having negative numbers appear there (for the sake of
    avoiding a long explanation of them, if nothing else), then
    eliminate them, but only in the context of the output medium. Then
    the internal values carried by the solver are unaffected. Remember
    that a small negative value in :math:`y` returned by ARKode, with
    magnitude comparable to ``abstol`` or less, is equivalent to zero
    as far as the computation is concerned. 

(3) The user's right-hand side routines :math:`f_E` and :math:`f_I`
    should never change a negative value in the solution vector :math:`y`
    to a non-negative value in attempt to "solve" to this problem,
    since this can cause instability.  If the :math:`f_E` or
    :math:`f_I` routines cannot tolerate a zero or negative value
    (e.g. because there is a square root or log), then the offending
    value should be changed to zero or a tiny positive number in a
    temporary variable (not in the input :math:`y` vector) for the
    purposes of computing :math:`f_E(t, y)` or :math:`f_I(t, y)`. 

(4) Positivity and non-negativity constraints on components can be
    enforced by use of the recoverable error return feature in the
    user-supplied right-hand side function. However, because this option
    involves some additional overhead cost, it should only be exercised if
    the use of absolute tolerances to control the computed values is
    unsuccessful. 


.. _CInterface.LinearSolvers:

Linear solver specification functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As previously explained, the modified Newton iteration used in solving
implicit systems within ARKode requires the solution of linear
systems of the form 

.. math::
    A\left(y^n(m)\right) s^m = -F\left(y^n(m)\right)

where 

.. math::
    A \approx M - \gamma J, \qquad J = \frac{\partial f_I}{\partial y}.

There are five ARKode linear solvers currently available for this
task: ARKDENSE, ARKBAND, ARKSPGMR, ARKSPBCG, and ARKSPTFQMR.

The first two linear solvers are direct and derive their names from
the type of approximation used for the Jacobian :math:`J`;
ARKDENSE and ARKBAND work with dense and banded approximations
to :math:`J`, respectively. The SUNDIALS suite includes both
internal implementations of these two linear solvers and interfaces to
LAPACK implementations. Together, these linear solvers are referred to
as ARKDLS (from Direct Linear Solvers). 

The last three ARKode linear solvers, ARKSPGMR, ARKSPBCG,
and ARKSPTFQMR, are Krylov iterative solvers, which use scaled
preconditioned GMRES, scaled preconditioned Bi-CGStab, and scaled
preconditioned TFQMR, respectively. Together, they are referred to as
ARKSPILS (from Scaled Preconditioned Iterative Linear Solvers). 

With any of the Krylov methods, preconditioning can be done on the
left only, on the right only, on both the left and the right, or not
at all. For the specification of a preconditioner, see the iterative
linear solver sections in :ref:`CInterface.OptionalOutputs` and
:ref:`CInterface.UserSupplied`. 

If preconditioning is done, user-supplied functions define left and
right preconditioner matrices :math:`P_1` and :math:`P_2` (either of
which could be the identity matrix), such that the product P1P2
approximates the Newton matrix  :math:`A = M - \gamma J`. 

To specify a ARKode linear solver, after the call to
:c:func:`ARKodeCreate()` but before any calls to :c:func:`ARKode()`, the user's
program must call one of the functions
:c:func:`ARKDense()`/:c:func:`ARKLapackDense()`, :c:func:`ARKBand()`/:c:func:`ARKLapackBand()`,
:c:func:`ARKSpgmr()`, :c:func:`ARKSpbcg()`, or :c:func:`ARKSptfqmr()`, as
documented below. The first argument passed to these functions is the
ARKode memory pointer returned by :c:func:`ARKodeCreate()`. A call to one
of these functions links the main ARKode integrator to a linear
solver and allows the user to specify parameters which are specific to
a particular solver, such as the half-bandwidths in the :c:func:`ARKBand()`
case. The use of each of the linear solvers involves certain constants
and possibly some macros, that are likely to be needed in the user
code. These are available in the corresponding header file associated
with the linear solver, as specified below. 

In each case except LAPACK direct solvers, the linear solver module
used by ARKode is actually built on top of a generic linear system
solver, which may be of interest in itself.  These generic solvers,
denoted DENSE, BAND, SPGMR, SPBCG, and SPTFQMR,
are described separately in the section :ref:`LinearSolvers`.



.. c:function:: int ARKDense(void *arkode_mem, long int N)

   A call to the ARKDense function links the main
   integrator with the ARKDENSE linear solver.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `N` -- the number of components in the ODE system.
   
   **Return value:** 
       * ARKDLS_SUCCESS   if successful
       * ARKDLS_MEM_NULL  if the ARKode memory was ``NULL``
       * ARKDLS_MEM_FAIL  if there was a memory allocation failure
       * ARKDLS_ILL_INPUT if a required vector operation is missing
   
   **Notes:**  The ARKDENSE linear solver may not be compatible with the
   particular implementation of the NVECTOR module. Of the two
   nvector modules provided with SUNDIALS, only NVECTOR_SERIAL is
   compatible. 



.. c:function:: int ARKLapackDense(void *arkode_mem, int N)

   A call to the ARKLapackDense function links the main
   integrator with the ARKLAPACK linear solver dense Jacobians.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `N` -- the number of components in the ODE system.
   
   **Return value:** 
      * ARKDLS_SUCCESS   if successful
      * ARKDLS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKDLS_MEM_FAIL  if there was a memory allocation failure
      * ARKDLS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** Here `N` is restricted to be of type ``int``, because of the
   corresponding type restriction in the LAPACK solvers.



.. c:function:: int ARKBand(void *arkode_mem, long int N, long int mupper, long int mlower)

   A call to the ARKBand function links the main
   integrator with the ARKBAND linear solver.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `N` -- the number of components in the ODE system
      * `mupper` -- the upper bandwidth of the band Jacobian approximation
      * `mlower` -- is the lower bandwidth of the band Jacobian approximation.
   
   **Return value:** 
      * ARKDLS_SUCCESS   if successful
      * ARKDLS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKDLS_MEM_FAIL  if there was a memory allocation failure
      * ARKDLS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** The ARKBAND linear solver may not be compatible with the
   particular implementation of the NVECTOR module. Of the two
   NVECTOR modules provided with SUNDIALS, only
   NVECTOR_SERIAL is compatible. The half-bandwidths are to be set
   such that the nonzero locations `(i, j)` in the banded
   (approximate) Jacobian satisfy `-mlower` :math:`\le` `j-i`
   :math:`\le` `mupper`. 



.. c:function:: int ARKLapackBand(void *arkode_mem, int N, int mupper, int mlower)

   A call to the ARKLapackBand function links the main
   integrator with the ARKLAPACK linear solver using banded Jacobians.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `N` -- the number of components in the ODE system
      * `mupper` -- the upper bandwidth of the band Jacobian approximation
      * `mlower` -- is the lower bandwidth of the band Jacobian approximation.
   
   **Return value:** 
      * ARKDLS_SUCCESS   if successful
      * ARKDLS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKDLS_MEM_FAIL  if there was a memory allocation failure
      * ARKDLS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** Here, each of `N`, `mupper` and `mlower` are restricted
   to be of type ``int``, because of the corresponding type restriction
   in the LAPACK solvers.



.. c:function:: int ARKSpgmr(void *arkode_mem, int pretype, int maxl)

   A call to the ARKSpgmr function links the main
   ARKode integrator with the ARKSPGMR linear solver.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `pretype` -- the type of user preconditioning to be done.  This
        must be one of the four enumeration constants PREC_NONE,
        PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in
        ``sundials_iterative.h``. These correspond to no preconditioning,
        left preconditioning only, right preconditioning only, and both
        left and right preconditioning, respectively.
      * `maxl` -- the maximum Krylov dimension. This is an optional input
        to the ARKSPGMR solver. Pass 0 to use the default value of 5.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKSPILS_MEM_FAIL  if there was a memory allocation failure
      * ARKSPILS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** The ARKSPGMR solver uses a scaled preconditioned GMRES
   iterative method to solve the linear systems.



.. c:function:: int ARKSpbcg(void *arkode_mem, int pretype, int maxl)

   A call to the ARKSpbcg function links the main
   ARKode integrator with the ARKSPBCG linear solver.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `pretype` -- the type of user preconditioning to be done.  This
        must be one of the four enumeration constants PREC_NONE,
        PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in
        ``sundials_iterative.h``. These correspond to no preconditioning,
        left preconditioning only, right preconditioning only, and both
        left and right preconditioning, respectively.
      * `maxl` -- the maximum Krylov dimension. This is an optional input
        to the ARKSPBCG solver. Pass 0 to use the default value of 5.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKSPILS_MEM_FAIL  if there was a memory allocation failure
      * ARKSPILS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** The ARKSPBCG solver uses a scaled preconditioned Bi-CGStab 
   iterative method to solve the linear systems.
   


.. c:function:: int ARKSptfqmr(void *arkode_mem, int pretype, int maxl)

   A call to the ARKSptfqmr function links the main
   ARKode integrator with the ARKSPTFQMR linear solver.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `pretype` -- the type of user preconditioning to be done.  This
        must be one of the four enumeration constants PREC_NONE,
        PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in
        ``sundials_iterative.h``. These correspond to no preconditioning,
        left preconditioning only, right preconditioning only, and both
        left and right preconditioning, respectively.
      * `maxl` -- the maximum Krylov dimension. This is an optional input
        to the ARKSPTFMR solver. Pass 0 to use the default value of 5.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKSPILS_MEM_FAIL  if there was a memory allocation failure
      * ARKSPILS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** The ARKSPTFQMR solver uses a scaled preconditioned TFQMR
   iterative method to solve the linear systems.



.. _CInterface.RootFinding:

Rootfinding initialization function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While solving the IVP, ARKode has the capability to find the roots
of a set of user-defined functions.  To activate the root-finding
algorithm, call the following function:



.. c:function:: int ARKodeRootInit(void *arkode_mem, int nrtfn, ARKRootFn g)

   Initializes a rootfinding problem to be solved
   during the integration of the ODE system.  It must be called
   after :c:func:`ARKodeCreate()`, and before :c:func:`ARKode()`. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nrtfn` -- number of functions :math:`g_i`, an integer :math:`\ge` 0.
      * `g` -- name of user-supplied function, of type :c:func:`ARKRootFn()`,
        defining the functions :math:`g_i` whose roots are sought. 
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL  if the ARKode memory was ``NULL``
      * ARK_MEM_FAIL  if there was a memory allocation failure
      * ARK_ILL_INPUT if `nrtfn` is greater than zero but `g` = ``NULL``.
   
   **Notes:** If a new IVP is to be solved with a call to :c:func:`ARKodeReInit()`,
   where the new IVP has no rootfinding problem but the prior one did,
   then call ARKodeRootInit with `nrtfn=0`.



.. _CInterface.Integration:

ARKode solver function
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the central step in the solution process -- the call to perform
the integration of the IVP.  One of the input arguments (`itask`)
specifies one of two modes as to where ARKode is to return a
solution.  These modes are modified if the user has set a stop time
(with a call to the optional input function :c:func`ARKodeSetStopTime()`) or
has requested rootfinding. 



.. c:function:: int ARKode(void *arkode_mem, realtype tout, N_Vector yout, realtype *tret, int itask)

   Integrates the ODE over an interval in :math:`t`.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `tout` -- the next time at which a computed solution is desired
      * `yout` -- the computed solution vector
      * `tret` -- the time reached by the solver (output)
      * `itask` -- a flag indicating the job of the solver for the next
        user step. The ARK_NORMAL option causes the solver to take internal
        steps until it has reached or just passed the user-specified `tout`
        parameter. The solver then interpolates in order to return an
        approximate value of :math:`y`(`tout`). This interpolation is
        typically less accurate than the full time step solutions produced
        by the solver, since the interpolation uses a cubic Hermite
        polynomial even when the RK method is of higher order.  If the user 
        wishes that this returned value have full method accuracy, they 
        may issue a call to :c:func:`ARKodeSetStopTime()` before the call to ARKode
        to specify a fixed stop time to end the time step and return to 
        the user.  Once the integrator returns at a `tstop` time, any 
        future testing for `tstop` is disabled (and can be reenabled only 
        though a new call to :c:func:`ARKodeSetStopTime()`).  The ARK_ONE_STEP
        option tells the solver to take just one internal step and then
        return the solution at the point reached by that step. 
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_ROOT_RETURN if ARKode succeeded, and found one or more roots.
        If `nrtfn` is greater than 1, call :c:func:`ARKodeGetRootInfo()` to see
        which :math:`g_i` were found to have a root at `(*tret)`. 
      * ARK_TSTOP_RETURN if ARKode succeeded and returned at `tstop`.
      * ARK_MEM_NULL if the `arkode_mem` argument was ``NULL``.
      * ARK_NO_MALLOC if `arkode_mem` was not allocated.
      * ARK_ILL_INPUT if one of the inputs to ARKode is illegal, or
        some other input to the solver was either illegal or missing.  The
        latter category includes the following situations:  (a) The
        tolerances have not been set. (b) A component of the error weight
        vector became zero during internal time-stepping. (c) The linear
        solver initialization function (called by the user after calling
        :c:func:`ARKodeCreate()`) failed to set the linear solver-specific
        `lsolve` field in `arkode_mem`. (d) A root of one of the root
        functions was found both at a point :math:`t` and also very near
        :math:`t`. In any case, the user should see the error message for
        details.
      * ARK_TOO_MUCH_WORK if the solver took `mxstep` internal steps
        but could not reach `tout`.  The default value for `mxstep` is
        `MXSTEP_DEFAULT = 500`.
      * ARK_TOO_MUCH_ACC if the solver could not satisfy the accuracy
        demanded by the user for some internal step.
      * ARK_ERR_FAILURE if error test failures occurred either too many
        times (`ark_maxnef`) during one internal time step or occurred
        with :math:`|h| = h_{min}`. 
      * ARK_CONV_FAILURE if either convergence test failures occurred
        too many times (`ark_maxncf`) during one internal time step or
        occurred with :math:`|h| = h_{min}`. 
      * ARK_LINIT_FAIL if the linear solver's initialization function failed.
      * ARK_LSETUP_FAIL if the linear solver's setup routine failed in
        an unrecoverable manner.
      * ARK_LSOLVE_FAIL if the linear solver's solve routine failed in
        an unrecoverable manner.
   
   **Notes:** The vector `yout` can occupy the same space as the vector
   `y0` of initial conditions that was passed to :c:func:`ARKodeInit()`. 
   
   In the ARK_ONE_STEP mode, `tout` is used only on the first
   call, and only to get the direction and a rough scale of the
   independent variable. 
 
   All failure return values are negative and so testing the return
   argument for negative values will trap all ARKode failures.
   
   On any error return in which one or more internal steps were taken
   by ARKode, the returned values of `tret` and `yout`
   correspond to the farthest point reached in the integration. On all
   other error returns, `tret` and `yout` are left unchanged from
   the previous ARKode return. 




.. _CInterface.OptionalInputs:

Optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are numerous optional input parameters that control the behavior
of the ARKode solver. ARKode provides functions that can be
used to change these optional input parameters from their default
values. The following tables list all optional input functions in
ARKode which are then described in detail in the remainder of this
section, begining with those for the main ARKode solver
(:ref:`CInterface.ARKodeInputTable`), then the dense linear solver
modules (:ref:`CInterface.ARKDlsInputTable`) and finally the optional
inputs for the sparse linear solver modules
(:ref:`CInterface.ARKSpilsInputTable`).  For the most casual use of
ARKode, the reader can skip to the section
:ref:`CInterface.UserSupplied`.

We note that, on an error return, all of the optional input functions
send an error message to the error handler function.  We also note
that all error return values are negative, so a test on the return
arguments for negative values will catch all errors. 

.. _CInterface.ARKodeInputTable:


Table: Optional inputs for ARKode
"""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

===============================================  ========================================  ==============
Optional input                                   Function name                             Default
===============================================  ========================================  ==============
Set default solver parameters                    :c:func:`ARKodeSetDefaults()`             internal
Set 'optimal' adaptivity params                  :c:func:`ARKodeSetOptimalParams()`        internal
Error handler function                           :c:func:`ARKodeSetErrHandlerFn()`         internal fn
Pointer to an error file                         :c:func:`ARKodeSetErrFile()`              ``stderr``
User data                                        :c:func:`ARKodeSetUserData()`             ``NULL``
Pointer to a diagnostics file                    :c:func:`ARKodeSetDiagnostics()`          ``NULL``
Set method order                                 :c:func:`ARKodeSetOrder()`                4
Set dense output order                           :c:func:`ARKodeSetDenseOrder()`           3
Specify linearly implicit :math:`f_I`            :c:func:`ARKodeSetLinear()`               ``FALSE``
Specify nonlinearly implicit :math:`f_I`         :c:func:`ARKodeSetNonlinear()`            ``TRUE``
Specify explicit problem                         :c:func:`ARKodeSetExplicit()`             ``FALSE``
Specify implicit problem                         :c:func:`ARKodeSetImplicit()`             ``FALSE``
Specify implicit/explicit problem                :c:func:`ARKodeSetImEx()`                 ``TRUE``
Set explicit RK table                            :c:func:`ARKodeSetERKTable()`             internal
Set implicit RK table                            :c:func:`ARKodeSetIRKTable()`             internal
Set additive RK tables                           :c:func:`ARKodeSetARKTables()`            internal
Specify explicit RK table number                 :c:func:`ARKodeSetERKTableNum()`          internal
Specify implicit RK table number                 :c:func:`ARKodeSetIRKTableNum()`          internal
Specify additive RK tables number                :c:func:`ARKodeSetARKTableNum()`          internal
Maximum no. of internal steps before `tout`      :c:func:`ARKodeSetMaxNumSteps()`          500
Maximum no. of warnings for :math:`t_n+h = t_n`  :c:func:`ARKodeSetMaxNumSteps()`          10
Initial step size                                :c:func:`ARKodeSetInitStep()`             estimated
Minimum absolute step size                       :c:func:`ARKodeSetMinStep()`              0.0
Maximum absolute step size                       :c:func:`ARKodeSetMaxStep()`              :math:`\infty`
Value of :math:`t_{stop}`                        :c:func:`ARKodeSetStopTime()`             :math:`\infty`
Time step adaptivity method                      :c:func:`ARKodeSetAdaptivityMethod()`     0
Time step adaptivity function                    :c:func:`ARKodeSetAdaptivityFn()`         internal
Time step adaptivity constants                   :c:func:`ARKodeSetAdaptivityConstants()`  internal
Newton convergence constants                     :c:func:`ARKodeSetNewtonConstants()`      internal
Linear solver setup decision constants           :c:func:`ARKodeSetLSetupConstants()`      internal
Implicit predictor method                        :c:func:`ARKodeSetPredictorMethod()`      3
Explicit stability function                      :c:func:`ARKodeSetStabilityFn()`          internal
Maximum no. of error test failures               :c:func:`ARKodeSetMaxErrTestFails()`      7
Maximum no. of nonlinear iterations              :c:func:`ARKodeSetMaxNonlinIters()`       3
Maximum no. of convergence failures              :c:func:`ARKodeSetMaxConvFails()`         10
Coefficient in the nonlinear convergence test    :c:func:`ARKodeSetNonlinConvCoef()`       0.2
===============================================  ========================================  ==============




.. c:function:: int ARKodeSetDefaults(void *arkode_mem)

   Resets all optional inputs to ARKode default
   values.  
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Does not change problem-defining function pointers `fe`
   and `fi` or the `user_data` pointer.  
   
   Also leaves alone any data structures or  options related to
   root-finding (those can be reset using :c:func:`ARKodeRootInit()`).



.. c:function:: int ARKodeSetOptimalParams(void *arkode_mem)

   Sets all adaptivity and solver parameters to our 'best
   guess' values, for a given integration method (ERK, DIRK, ARK) and
   a given method order.  
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Should only be called after the method order and integration
   method have been set.



.. c:function:: int ARKodeSetErrHandlerFn(void *arkode_mem, ARKErrHandlerFn ehfun, void *eh_data)

   Specifies the optional user-defined function to be used
   in handling error messages.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `ehfun` -- name of user-supplied error handler function. 
      * `eh_data` -- pointer to user data passed to `ehfun` every time
        it is called
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Error messages indicating that the ARKode solver memory is
   ``NULL`` will always be directed to ``stderr``.



.. c:function:: int ARKodeSetErrFile(void *arkode_mem, FILE *errfp)

   Specifies a pointer to the file where all ARKode
   warning and error messages will be written if the default internal
   error handling function is used. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `errfp` -- pointer to the output file. 
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** The default value for `errfp` is ``stderr``.
    
   Passing a ``NULL`` value disables all future error message output
   (except for the case wherein the ARKode memory pointer is
   ``NULL``.  This use of the function is strongly discouraged.
   
   If used, this routine should be called before any other
   optional input functions, in order to take effect for subsequent
   error messages.



.. c:function:: int ARKodeSetUserData(void *arkode_mem, void *user_data)

   Specifies the user data block `user_data` and
   attaches it to the main ARKode memory block.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `user_data` -- pointer to the user data
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** If specified, the pointer to `user_data` is passed to all
   user-supplied functions for which it is an argument; otherwise
   ``NULL`` is passed.
   
   If `user_data` is needed in user preconditioner functions, the
   call to this function must be made *before* the call to
   specify the linear solver.



.. c:function:: int ARKodeSetDiagnostics(void *arkode_mem, FILE *diagfp)

   Specifies the file pointer for a diagnostics file where
   all ARKode step adaptivity and solver information is written.  
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `diagfp` -- pointer to the diagnostics output file
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This parameter can be ``stdout`` or ``stderr``, although the
   suggested approach is to specify a pointer to a unique file opened
   by the user and returned by ``fopen``.  If not called, or if called
   with a ``NULL`` file pointer, all diagnostics output is disabled.
   
   When run in parallel, only one process should set a non-NULL value
   for this pointer, since statistics from all processes would be
   identical.
   


.. c:function:: int ARKodeSetOrder(void *arkode_mem, int ord)

   Specifies the order of accuracy for the linear
   multistep method.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `ord` -- requested order of accuracy
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** For explicit methods, the allowed values are 2 :math:`\le`
   `ord` :math:`\le` 6.  For implicit and IMEX methods,  the allowed values are 3 :math:`\le`
   `ord` :math:`\le` 5.  An illegal input will result in the default value of 4.
   
   Since `ord` affects the memory requirements for the internal
   ARKode memory block, it cannot be increased between calls to
   :c:func:`ARKode()` unless :c:func:`ARKodeReInit()` is called.



.. c:function:: int ARKodeSetDenseOrder(void *arkode_mem, int dord)

   Specifies the order of accuracy for the polynomial
   interpolant used for dense output.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `dord` -- requested polynomial order of accuracy
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Allowed values are between 0 and ``min(q,3)``, where ``q`` is
   the order of the overall integration method.



.. c:function:: int ARKodeSetLinear(void *arkode_mem)

   Specifies that the implicit portion of the problem is linear.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Tightens the linear solver tolerances and takes only a single
   Newton iteration.



.. c:function:: int ARKodeSetNonlinear(void *arkode_mem)

   Specifies that the implicit portion of the problem is nonlinear.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This is the default behavior of ARKode, so the function
   ARKodeSetNonlinear is primarily useful to undo a previous call
   to :c:func:`ARKodeSetLinear()`. 



.. c:function:: int ARKodeSetExplicit(void *arkode_mem)

   Specifies that the implicit portion of problem is disabled,
   and to use an explicit RK method.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This is automatically deduced when the function pointer `fi`
   passed to :c:func:`ARKodeInit()` is ``NULL``, but may be set
   directly by the user if desired.



.. c:function:: int ARKodeSetImplicit(void *arkode_mem)

   Specifies that the explicit portion of problem is disabled,
   and to use a diagonally implicit RK method.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This is automatically deduced when the function pointer `fe`
   passed to :c:func:`ARKodeInit()` is ``NULL``, but may be set directly by the
   user if desired.



.. c:function:: int ARKodeSetImEx(void *arkode_mem)

   Specifies that both the implicit and explicit portions
   of problem are enabled, and to use an additive Runge Kutta method.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This is automatically deduced when neither of the function
   pointers `fe` or `fi` passed to :c:func:`ARKodeInit()` are ``NULL``, but
   may be set directly by the user if desired.



.. c:function:: int ARKodeSetERKTable(void *arkode_mem, int s, int q, int p, realtype *c, realtype *A, realtype *b, realtype *bembed)

   Specifies a customized Butcher table for the explicit portion of the system.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `s` -- number of stages in the RK method
      * `q` -- global order of accuracy for the RK method
      * `p` -- global order of accuracy for the embedded RK method
      * `c` -- array (of length `s`) of stage times for the RK method.
      * `A` -- array of coefficients defining the RK stages.  This should
        be stored as a 1D array of size `s*s`, in row-major order.
      * `b` -- array of coefficients (of length `s`) defining the time step solution.
      * `bembed` -- array of coefficients (of length `s`) defining the embedded solution.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This automatically calls :c:func:`ARKodeSetExplicit()`.
   
   No error checking is performed to ensure that either `p` or `q`
   correctly describe the coefficients that were input.
   
   Error checking is performed to ensure that `A` is strictly
   lower-triangular (i.e. that it specifies an ERK method).
   
   The embedding `bembed` is required.



.. c:function:: int ARKodeSetIRKTable(void *arkode_mem, int s, int q, int p, realtype *c, realtype *A, realtype *b, realtype *bembed)

   Specifies a customized Butcher table for the implicit portion of the system.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `s` -- number of stages in the RK method
      * `q` -- global order of accuracy for the RK method
      * `p` -- global order of accuracy for the embedded RK method
      * `c` -- array (of length `s`) of stage times for the RK method.
      * `A` -- array of coefficients defining the RK stages.  This should
        be stored as a 1D array of size `s*s`, in row-major order.
      * `b` -- array of coefficients (of length `s`) defining the time step solution.
      * `bembed` -- array of coefficients (of length `s`) defining the embedded solution.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This automatically calls :c:func:`ARKodeSetImplicit()`.
   
   No error checking is performed to ensure that either `p` or `q`
   correctly describe the coefficients that were input.
   
   Error checking is performed to ensure that `A` is 
   lower-triangular with nonzeros on at least some of the diagonal
   entries (i.e. that it specifies a DIRK method).
   
   The embedding `bembed` is required.



.. c:function:: int ARKodeSetARKTables(void *arkode_mem, int s, int q, int p, realtype *c, realtype *Ai, realtype *Ae, realtype *b, realtype *bembed)

   Specifies a customized Butcher table pair for the
   additive RK method.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `s` -- number of stages in the RK method
      * `q` -- global order of accuracy for the RK method
      * `p` -- global order of accuracy for the embedded RK method
      * `c` -- array (of length `s`) of stage times for the RK method.
      * `Ai` -- array of coefficients defining the implicit RK stages.  This should
        be stored as a 1D array of size `s*s`, in row-major order.
      * `Ae` -- array of coefficients defining the explicit RK stages.  This should
        be stored as a 1D array of size `s*s`, in row-major order.
      * `b` -- array of coefficients (of length `s`) defining the time step solution.
      * `bembed` -- array of coefficients (of length `s`) defining the embedded solution.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``   
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This automatically calls :c:func:`ARKodeSetImEx()`.
   
   No error checking is performed to ensure that either `p` or `q`
   correctly describe the coefficients that were input.
   
   Error checking is performed on both `Ai` and `Ae` to ensure
   that they specify DIRK and ERK methods, respectively.  
   
   Both RK methods must share the same `c`, `b` and `bembed` coefficients.
  
   The embedding `bembed` is required.



.. c:function:: int ARKodeSetERKTableNum(void *arkode_mem, int etable)

   Specifies to use a built-in Butcher table for the
   explicit portion of the system.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `etable` -- index of the Butcher table.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** `etable` should match an existing method in the function
   ARKodeLoadButcherTable within the file ``arkode_butcher.c``.
   Error-checking is performed to ensure that the table exists, and is
   not implicit.  
   
   This automatically calls :c:func:`ARKodeSetExplicit()`. 



.. c:function:: int ARKodeSetIRKTableNum(void *arkode_mem, int itable)

   Specifies to use a built-in Butcher table for the
   implicit portion of the system.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `itable` -- index of the Butcher table.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** `itable` should match an existing method in the function
   ARKodeLoadButcherTable within the file ``arkode_butcher.c``.
   Error-checking is performed to ensure that the table exists, and is
   not explicit.  
   
   This automatically calls :c:func:`ARKodeSetImplicit()`. 



.. c:function:: int ARKodeSetARKTableNum(void *arkode_mem, int itable, int etable)

   Specifies to use built-in Butcher tables for the ImEx system.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `itable` -- index of the DIRK Butcher table.
      * `etable` -- index of the ERK Butcher table.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Both `itable` and `etable` should match existing methods
   in the function ARKodeLoadButcherTable within the file
   ``arkode_butcher.c``. 
   
   Error-checking is performed to ensure that the tables exist.
   Subsequent error-checking is automatically performed to ensure that
   the tables' stage times and solution coefficients match.  This
   automatically calls :c:func:`ARKodeSetImEx()`. 



.. c:function:: int ARKodeSetMaxNumSteps(void *arkode_mem, long int mxsteps)

   Specifies the maximum number of steps to be taken by the
   solver in its attempt to reach the next output time.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `mxsteps` -- maximum allowed number of internal steps.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Passing `mxsteps = 0` results in ARKode using the
   default value (500).
   
   Passing `mxsteps < 0` disables the test `(not recommended)`.



.. c:function:: int ARKodeSetMaxHnilWarns(void *arkode_mem, int mxhnil)

   Specifies the maximum number of messages issued by the
   solver warning that :math:`t+h=t` on the next internal step.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `mxhnil` -- maximum allowed number of warning messages (>0).
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** The default value is 10.

   A negative value indicates that no warning messages should be issued.



.. c:function:: int ARKodeSetInitStep(void *arkode_mem, realtype hin)

   Specifies the initial time step size.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `hin` -- value of the initial step to be attempted :math:`(\ge 0)`
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Pass 0.0 to use the default value.  
   
   By default, ARKode estimates the initial step size to be the
   solution :math:`h` of the equation :math:`\left\| \frac{h^2
   \ddot{y}}{2}\right\| = 1`, where :math:`\ddot{y}` is an estimated
   value of the second derivative of the solution at `t0`.



.. c:function:: int ARKodeSetMinStep(void *arkode_mem, realtype hmin)

   Specifies the lower bound on the magnitude of the time step size.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `hmin` -- minimum absolute value of the time step size :math:`(\ge 0)`
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** The default value is 0.0.  



.. c:function:: int ARKodeSetMaxStep(void *arkode_mem, realtype hmax)

   Specifies the upper bound on the magnitude of the time step size.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `hmax` -- maximum absolute value of the time step size :math:`(\ge 0)`
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Pass `hmax = 0.0` to set the default value of :math:`\infty`.  



.. c:function:: int ARKodeSetStopTime(void *arkode_mem, realtype tstop)

   Specifies the value of the independent variable
   :math:`t` past which the solution is not to proceed.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `tstop` -- stopping time for the integrator.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** The default is that no stop time is imposed.



.. c:function:: int ARKodeSetAdaptivityMethod(void *arkode_mem, int imethod, realtype *adapt_params)

   Specifies the method (and associated parameters) used
   for time step adaptivity.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `imethod` -- accuracy-based adaptivity method choice 
        (0 :math:`\le` `imethod` :math:`\le` 5): 
        0 is PID, 1 is PI, 2 is I, 3 is explicit Gustafsson, 4 is
        implicit Gustafsson, and 5 is the ImEx Gustafsson.
      * `adapt_params[0]` -- (*cfl*) fraction of the estimated explicitly stable
        step to use (default is 0.5)
      * `adapt_params[1]` -- (*safety*) safety factor applied to accuracy-based time
        step (default is 0.96)
      * `adapt_params[2]` -- (*bias*) bias applied to error in accuracy-based time
        step estimation (default is 1.5)
      * `adapt_params[3]` -- (*growth*) maximum allowed growth factor between
        consecutive time steps (default is 20.0)
      * `adapt_params[4]` -- (*lb*) lower bound on window to leave step size fixed (default is 1.0)
      * `adapt_params[5]` -- (*ub*) upper bound on window to leave step size fixed (default is 1.5)
      * `adapt_params[6]` -- :math:`k_1` parameter within accuracy-based adaptivity algorithms.
      * `adapt_params[7]` -- :math:`k_2` parameter within accuracy-based adaptivity algorithms.
      * `adapt_params[8]` -- :math:`k_3` parameter within accuracy-based adaptivity algorithms.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Any zero-valued parameter will imply a reset to the default
   value.  
   
   Any negative parameter will be left unchanged from the previous value.


      
.. c:function:: int ARKodeSetAdaptivityFn(void *arkode_mem, ARKAdaptFn hfun, void *h_data)

   Sets a user-supplied time-step adaptivity function.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `hfun` -- name of user-supplied adaptivity function.
      * `h_data` -- pointer to user data passed to `hfun` every time
        it is called
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This function should focus on accuracy-based time step
   estimation; for stability based time steps the function
   :c:func:`ARKodeSetStabilityFn()` should be used instead.


      
.. c:function:: int ARKodeSetAdaptivityConstants(void *arkode_mem, realtype etamx1, realtype etamxf, realtype etacf, int small_nef)

   Specifies additional parameters used in time step adaptivity.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `etamx1` -- maximum allowed growth factor after the first time
        step (default is 10000.0)
      * `etamxf` -- time step reduction factor on multiple error fails (default is 0.3)
      * `etacf` -- time step reduction factor on a nonlinear solver
        convergence failure (default is 0.25)
      * `small_nef` -- bound to determine `multiple` for `etamxf` (default is 2)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This function is designed only for advanced ARKode
   usage.

   Any zero-valued parameter will imply a reset to the default value.

   Any negative parameter will be left unchanged from the previous state.



.. c:function:: int ARKodeSetNewtonConstants(void *arkode_mem, realtype crdown, realtype rdiv)

   Specifies nonlinear convergence constants.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `crdown` -- nonlinear convergence rate estimation constant (default is 0.3)
      * `rdiv` -- Tolerance on Newton correction size ratio to declare divergence (default is 2.3)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This function is designed only for advanced ARKode usage.

   Any zero-valued parameter will imply a reset to the default value.

   Any negative parameter will be left unchanged from the previous state.



.. c:function:: int ARKodeSetLSetupConstants(void *arkode_mem, realtype dgmax, int msbp)

   Specifies linear setup decision constants.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `dgmax` -- tolerance on step size ratio change before calling
        linear solver setup routine (default is 0.2)
      * `msbp` -- maximum no. of time steps between linear solver setup calls (default is 20)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This function is designed only for advanced ARKode usage.

   Any zero-valued parameter will imply a reset to the default value.

   Any negative parameter will be left unchanged from the previous state.



.. c:function:: int ARKodeSetPredictorMethod(void *arkode_mem, int method)

   Specifies the method to use for predicting implicit solutions.  
   Non-default choices are {1,2,3}, all others will use default 
   (trivial) predictor.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `method` -- method choice (0 :math:`\le` `method` :math:`\le`
        3): 0 is the trivial predictor, 1 is the dense output predictor, 2
        is the dense output predictor that decreases the polynomial degree
        for more distant RK stages, 3 is the dense output predictor to max
        order for early RK stages, and a first-order predictor for distant
        RK stages.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This function is designed only for advanced ARKode usage.



.. c:function:: int ARKodeSetStabilityFn(void *arkode_mem, ARKExpStabFn EStab, void *estab_data)

   Sets the problem-dependent function to estimate a stable
   time step size for the explicit portion of the ODE system.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `EStab` -- name of user-supplied stability function.
      * `estab_data` -- pointer to user data passed to `EStab` every time
        it is called.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This function should return an estimate of the maximum stable
   time step for the explicit portion of the IMEX system.  It is not
   required, since accuracy-based adaptivity may be sufficient at
   retaining stability, but this can be quite useful for problems
   where the IMEX splitting may retain stiff components in
   :math:`f_E(t,y)`. 



.. c:function:: int ARKodeSetMaxErrTestFails(void *arkode_mem, int maxnef)

   Specifies the maximum number of error test failures
   permitted in attempting one step.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `maxnef` -- maximum allowed number of error test failures :math:`(>0)`
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** The default value is 7.



.. c:function:: int ARKodeSetMaxNonlinIters(void *arkode_mem, int maxcor)

   Specifies the maximum number of nonlinear solver
   iterations permitted per RK stage within each time step.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `maxcor` -- maximum allowed solver iterations per stage :math:`(>0)`
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** The default value is 3.



.. c:function:: int ARKodeSetMaxConvFails(void *arkode_mem, int maxncf)

   Specifies the maximum number of nonlinear solver
   convergence failures permitted during one step.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `maxncf` -- maximum allowed nonlinear solver convergence failures
        per step :math:`(>0)`
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** The default value is 10.  Upon each convergence failure,
   ARKode will first call the Jacobian setup routine and try again;
   if a convergence failure still occurs, the time step size is reduced
   by the factor `etacf` (set within
   :c:func:`ARKodeSetAdaptivityConstants()`). 



.. c:function:: int ARKodeSetNonlinConvCoef(void *arkode_mem, realtype nlscoef)

   Specifies the safety factor used within the nonlinear
   solver convergence test.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nlscoef` -- coefficient in nonlinear solver convergence test :math:`(>0.0)`
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** The default value is 0.2.



Direct linear solvers optional input functions
""""""""""""""""""""""""""""""""""""""""""""""

.. _CInterface.ARKDlsInputTable:

Table: Optional inputs for ARKDLS
"""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

==========================  ===============================  =============
Optional input              Function name                    Default
==========================  ===============================  =============
Dense Jacobian function     :c:func:`ARKDlsSetDenseJacFn()`     ``DQ``
Band Jacobian function      :c:func:`ARKDlsSetBandJacFn()`      ``DQ``
==========================  ===============================  =============

The ARKDENSE solver needs a function to compute a dense approximation
to the Jacobian matrix :math:`J(t,y)`. This function must be of type
:c:func:`ARKDlsDenseJacFn()`. The user can supply his/her own dense Jacobian
function, or use the default internal difference quotient
approximation that comes with the ARKDENSE solver. To specify a 
user-supplied Jacobian function `djac`, ARKDENSE provides the
function :c:func:`ARKDlsSetDenseJacFn()`. The ARKDENSE solver
passes the pointer user data to the dense Jacobian function. This
allows the user to create an arbitrary structure with relevant problem
data and access it during the execution of the user-supplied Jacobian
function, without using global data in the program. The pointer user
data may be specified through :c:func:`ARKodeSetUserData()`.



.. c:function:: int ARKDlsSetDenseJacFn(void *arkode_mem, ARKDlsDenseJacFn djac)

   Specifies the dense Jacobian approximation routine to
   be used for a direct dense linear solver. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `djac` -- name of user-supplied dense Jacobian approximation function.
   
   **Return value:** 
      * ARKDLS_SUCCESS  if successful
      * ARKDLS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKDLS_LMEM_NULL if the linear solver memory was ``NULL``
   
   **Notes:** By default, ARKDENSE uses an internal difference quotient
   function.  
   
   If ``NULL`` is passed in for `djac`, this default is used.
  
   The function type :c:func:`ARKDlsDenseJacFn()` is described in the section
   :ref:`CInterface.UserSupplied`.



Similarly, the ARKBAND solver needs a function to compute a banded
approximation to the Jacobian matrix :math:`J(t,y)`. This function
must be of type :c:func:`ARKDlsBandJacFn()`. The user can supply his/her own
banded Jacobian approximation function, or use the default internal
difference quotient approximation that comes with the ARKBAND
solver. To specify a user-supplied Jacobian function `bjac`,
ARKBAND provides the function :c:func:`ARKDlsSetBandJacFn()`. The
ARKBAND solver passes the pointer user data to the banded Jacobian
approximation function.  This allows the user to create an arbitrary
structure with relevant problem data and access it during the
execution of the user-supplied Jacobian function, without using global
data in the program. The pointer user data may be specified through
:c:func:`ARKodeSetUserData()`. 



.. c:function:: int ARKDlsSetBandJacFn(void *arkode_mem, ARKDlsBandJacFn bjac)

   Specifies the band Jacobian approximation routine to be
   used for a direct band linear solver. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `bjac` -- name of user-supplied banded Jacobian approximation function.
   
   **Return value:** 
      * ARKDLS_SUCCESS  if successful
      * ARKDLS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKDLS_LMEM_NULL if the linear solver memory was ``NULL``
   
   **Notes:** By default, ARKBAND uses an internal difference quotient
   function.
   
   If ``NULL`` is passed in for `bjac`, this default is used.
   
   The function type :c:func:`ARKDlsBandJacFn()` is described in the section
   :ref:`CInterface.UserSupplied`.



Iterative linear solvers optional input functions
"""""""""""""""""""""""""""""""""""""""""""""""""""

If any preconditioning is to be done within one of the ARKSPILS
linear solvers, then the user must supply a preconditioner solve
function `psolve` and specify its name in a call to
:c:func:`ARKSpilsSetPreconditioner()`. The evaluation and preprocessing
of any Jacobian-related data needed by the user's preconditioner solve
function is done in the optional user-supplied function
`psetup`. Both of these functions are fully specified in the section
:ref:`CInterface.UserSupplied`. If used, the `psetup` function
should also be specified in the call to
:c:func:`ARKSpilsSetPreconditioner()`. The pointer user data received
through :c:func:`ARKodeSetUserData()` (or a pointer to ``NULL`` if user
data was not specified) is passed to the preconditioner `psetup` and
`psolve` functions. This allows the user to create an arbitrary
structure with relevant problem data and access it during the
execution of the user-supplied preconditioner functions without using
global data in the program. 

The ARKSPILS solvers require a function to compute an
approximation to the product between the Jacobian matrix
:math:`J(t,y)` and a vector :math:`v`. The user can supply his/her own
Jacobian-times-vector approximation function, or use the default
internal difference quotient function that comes with the ARKSPILS
solvers. A user-defined Jacobian-vector function must be of type
:c:func:`ARKSpilsJacTimesVecFn()` and can be specified through a call to
:c:func:`ARKSpilsSetJacTimesVecFn()` (see the section
:ref:`CInterface.UserSupplied` for specification details). As with the
preconditioner user-supplied functions, a pointer to the user-defined
data structure, `user_data`, specified through
:c:func:`ARKodeSetUserData()` (or a ``NULL`` pointer otherwise) is
passed to the Jacobian-times-vector function `jtimes` each time it
is called.

.. _CInterface.ARKSpilsInputTable:

Table: Optional inputs for ARKSPILS
"""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

=============================================  =====================================  ==================
Optional input                                 Function name                          Default
=============================================  =====================================  ==================
Preconditioner functions                       :c:func:`ARKSpilsSetPreconditioner()`  ``NULL``, ``NULL``
Jacobian-times-vector function                 :c:func:`ARKSpilsSetJacTimesVecFn()`   ``DQ``
Preconditioning type                           :c:func:`ARKSpilsSetPrecType()`        none
Ratio between linear and nonlinear tolerances  :c:func:`ARKSpilsSetEpsLin()`          0.05
Type of Gram-Schmidt orthogonalization `(a)`   :c:func:`ARKSpilsSetGSType()`          classical GS
Maximum Krylov subspace size `(b)`             :c:func:`ARKSpilsSetMaxl()`            5
=============================================  =====================================  ==================


`(a)` Only for ARKSPGMR

`(b)` Only for ARKSPBCG and ARMSPTFQMR



.. c:function:: int ARKSpilsSetPreconditioner(void *arkode_mem, ARKSpilsPrecSetupFn psetup, ARKSpilsPrecSolveFn psolve)

   Specifies the preconditioner setup and solve functions.  
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `psetup` -- user defined preconditioner setup function.  Pass
        ``NULL`` if no setup is to be done
      * `psolve` -- user-defined preconditioner solve function.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful.
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``.
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``.
      * ARKSPILS_ILL_INPUT if an input has an illegal value.
   
   **Notes:** The default is ``NULL`` for both arguments (i.e. no
   preconditioning).
    
   Both of the function types :c:func`ARKSpilsPrecSetupFn()` and
   c:func:`ARKSpilsPrecSolveFn()` are described in the section
   :ref:`CInterface.UserSupplied`. 



.. c:function:: int ARKSpilsSetJacTimesVecFn(void *arkode_mem, ARKSpilsJacTimesVecFn jtimes)

   Specifies the Jacobian-times-vector function. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `jtimes` -- user-defined Jacobian-vector product function.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful.
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``.
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``.
      * ARKSPILS_ILL_INPUT if an input has an illegal value.

   **Notes:** The default is to use an internal finite difference
   approximation routine.  If ``NULL`` is passed to `jtimes`, this
   default function is used.
   
   The function type :c:func:`ARKSpilsJacTimesVecFn()` is described in the
   section :ref:`CInterface.UserSupplied`.



.. c:function:: int ARKSpilsSetPrecType(void *arkode_mem, int pretype)

   Resets the type of preconditioner, `pretype`, from the value previously set.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `pretype` -- the type of preconditioning to use, must be one of
        PREC_NONE, PREC_LEFT, PREC_RIGHT or PREC_BOTH. 
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful.
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``.
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``.
      * ARKSPILS_ILL_INPUT if an input has an illegal value.
   
   **Notes:** The preconditioning type is initially set in the call to the
   linear solver's specification function (see the section
   :ref:`CInterface.LinearSolvers`).  This function call is needed
   only if `pretype` is being changed from its original value.



.. c:function:: int ARKSpilsSetEpsLin(void *arkode_mem, realtype eplifac)

   Specifies the factor by which the tolerance on the
   nonlinear iteration is multiplied to get a tolerance on the linear iteration.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `eplifac` -- linear convergence safety factor :math:`(\ge 0.0)`.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful.
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``.
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``.
      * ARKSPILS_ILL_INPUT if an input has an illegal value.
   
   **Notes:** Passing a value `eplifac` of 0.0 indicates to use the default value of 0.05.



.. c:function:: int ARKSpilsSetGSType(void *arkode_mem, int gstype)

   Specifies the type of Gram-Schmidt orthogonalization to
   be used with the ARKSPGMR linear solver. This must be one of
   the two enumeration constants MODIFIED_GS or CLASSICAL_GS
   defined in ``iterative.h``. These correspond to using modified
   Gram-Schmidt and classical Gram-Schmidt, respectively.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `gstype` -- type of Gram-Schmidt orthogonalization.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful.
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``.
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``.
      * ARKSPILS_ILL_INPUT if an input has an illegal value.
   
   **Notes:** The default value is MODIFIED_GS.
   
   This option is available only for the ARKSPGMR linear solver.



.. c:function:: int ARKSpilsSetMaxl(void *arkode_mem, int maxl)

   Resets the maximum Krylov subspace size, `maxl`, from
   the value previously set, when using the Bi-CGStab or TFQMR linear
   solver methods.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `maxl` -- maximum dimension of the Krylov subspace.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful.
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``.
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``.
      * ARKSPILS_ILL_INPUT if an input has an illegal value.
   
   **Notes:** The maximum subspace dimension is initially specified in the
   call to the linear solver specification function (see the section
   :ref:`CInterface.LinearSolvers`).  This function call is needed
   only if `maxl` is being changed from its previous value.
  
   An input value `maxl` :math:`\le 0`, gives the default value, 5.
   
   This option is available only for the ARKSPBCG and
   ARKSPTFQMR linear solvers.



Rootfinding optional input functions
"""""""""""""""""""""""""""""""""""""

The following functions can be called to set optional inputs to
control the rootfinding algorithm.

.. cssclass:: table-bordered

=============================================  =======================================  ==================
Optional input                                 Function name                            Default
=============================================  =======================================  ==================
Direction of zero-crossings to monitor         :c:func:`ARKodeSetRootDirection()`       both
Disabling inactive root warnings               :c:func:`ARKodeSetNoInactiveRootWarn()`  warning
=============================================  =======================================  ==================



.. c:function:: int ARKodeSetRootDirection(void *arkode_mem, int *rootdir)

   Specifies the direction of zero-crossings to be located
   and returned.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `rootdir` -- state array of length `nrtfn`, the number of root
        functions :math:`g_i`, as specified in the call to the function
        :c:func:`ARKodeRootInit()`. A value of 0 for ``rootdir[i]``
        indicates that crossing in either direction for :math:`g_i` should
        be reported.  A value of +1 or -1 indicates that the solver should
        report only zero-crossings where :math:`g_i` is increasing or
        decreasing, respectively.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** The default behavior is to monitor for both zero-crossing
      directions.



.. c:function:: int ARKodeSetNoInactiveRootWarn(void *arkode_mem)

   Disables issuing a warning if some root function appears
   to be identically zero at the beginning of the integration.
  
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
   
   **Notes:** ARKode will not report the initial conditions as a
   possible zero-crossing (assuming that one or more components
   :math:`g_i` are zero at the initial time).  However, if it appears
   that some :math:`g_i` is identically zero at the initial time
   (i.e., :math:`g_i` is zero at the initial time and after the first
   step), ARKode will issue a warning which can be disabled with
   this optional input function. 





.. _CInterface.InterpolatedOutput:

Interpolated output function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An optional function :c:func:`ARKodeGetDky()` is available to obtain
additional output values.  This function should only be called after a
successful return from :c:func:`ARKode()` as it provides interpolated
values either of :math:`y` or of its derivatives (up to the 3rd
derivative) interpolated to any value of :math:`t` in the last
internal step taken by :c:func:`ARKode()`. 



.. c:function:: int ARKodeGetDky(void *arkode_mem, realtype t, int k, N_Vector dky)

   Computes the `k`-th derivative of the function
   :math:`y` at the time `t`, i.e. :math:`\frac{d^(k)y}{dt^(k)}`,
   where :math:`t_n-h_n \le t \le t_n`, :math:`t_n` denotes the
   current internal time reached, and :math:`h_n` is the last internal
   step size successfully used by the solver.  The user may request
   `k` in the range 0,1,2,3.  This routine uses an interpolating
   polynomial of degree `max(dord, k)`, where `dord` is the
   argument provided to :c:func:`ARKodeSetDenseOrder()`, i.e. it will
   form a polynomial of the degree requested by the user through
   `dord`, unless higher-order derivatives are requested.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `t` -- the value of the independent variable at which the
        derivative is to be evaluated
      * `k` -- the derivative order requested
      * `dky` -- vector containing the derivative.  This vector must be
        allocated by the user.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_BAD_K if `k` is not in the range 0,1,2,3.
      * ARK_BAD_T if `t` is not in the interval :math:`[t_n-h_n, t_n]`
      * ARK_BAD_DKY if the `dky` argument was ``NULL``
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
   
   **Notes:** It is only legal to call this function after a successful
   return from :c:func:`ARKode()`.  See :c:func:`ARKodeGetCurrentTime()`
   and :c:func:`ARKodeGetLastStep()` in the next section for access to
   :math:`t_n` and :math:`h_n`, respectively.




.. _CInterface.OptionalOutputs:

Optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^

ARKode provides an extensive set of functions that can be used to
obtain solver performance information. In the tables 
:ref:`CInterface.ARKodeOutputTable`,
:ref:`CInterface.ARKodeRootOutputTable`,
:ref:`CInterface.ARKDlsOutputTable` and
:ref:`CInterface.ARKSpilsOutputTable`, we list all of the optional
output functions in ARKode, which are then described in detail in
the remainder of this section. 

Some of the optional outputs, especially the various counters, can be
very useful in determining how successful the :c:func:`ARKode()` solver
is in doing its job.  For example, the counters `nsteps`,
`nfe_evals` and `nfi_evals` provide a rough measure of the overall
cost of a given run, and can be compared among runs with differing
input options to suggest which set of options is most efficient.  The
ratio `nniters`/`nsteps` measures the performance of the modified
Newton iteration in solving the nonlinear systems at each stage;
typical values for this range from 1.1 to 1.8.  The ratio
`njevals`/`nniters` (in the case of a direct linear solver), and
the ratio `npevals`/`nniters` (in the case of an iterative linear
solver) measure the overall degree of nonlinearity in these systems,
and also the quality of the approximate Jacobian or preconditioner
being used.  Thus, for example, `njevals`/`nniters` can indicate
if a user-supplied Jacobian is inaccurate, if this ratio is larger
than for the case of the corresponding internal Jacobian.  The ratio
`nliters`/`nniters` measures the performance of the Krylov
iterative linear solver, and thus (indirectly) the quality of the
preconditioner.

Similarly, the ratio of explicit stability-limited steps to
accuracy-limited steps can measure the quality of the ImEx splitting
used (with a higher-quality splitting dominated by accuracy-limited
steps). 


Main solver optional output functions
"""""""""""""""""""""""""""""""""""""""

.. _CInterface.ARKodeOutputTable:

Table: Optional outputs for ARKode
"""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

===================================================  ============================================
Optional output                                      Function name
===================================================  ============================================
Size of ARKode real and integer workspaces           :c:func:`ARKodeGetWorkSpace()`
Cumulative number of internal steps                  :c:func:`ARKodeGetNumSteps()`
No. of explicit stability-limited steps              :c:func:`ARKodeGetNumExpSteps()`
No. of accuracy-limited steps                        :c:func:`ARKodeGetNumAccSteps()`
No. of solver convergence-limited steps              :c:func:`ARKodeGetNumConvSteps()`
No. of calls to `fe` and `fi` functions              :c:func:`ARKodeGetNumRhsEvals()`
No. of calls to linear solver setup function         :c:func:`ARKodeGetNumLinSolvSetups()`
No. of local error test failures that have occurred  :c:func:`ARKodeGetNumErrTestFails()`
Actual initial time step size used                   :c:func:`ARKodeGetActualInitStep()`
Step size used for the last successful step          :c:func:`ARKodeGetLastStep()`
Step size to be attempted on the next step           :c:func:`ARKodeGetCurrentStep()`
Current internal time reached by the solver          :c:func:`ARKodeGetCurrentTime()`
Current ERK and DIRK Butcher tables                  :c:func:`ARKodeGetCurrentButcherTables()`
Suggested factor for tolerance scaling               :c:func:`ARKodeGetTolScaleFactor()`
Error weight vector for state variables              :c:func:`ARKodeGetErrWeights()`
Estimated local truncation error vector              :c:func:`ARKodeGetEstLocalErrors()`
Single accessor to many statistics at once           :c:func:`ARKodeGetIntegratorStats()`
No. of nonlinear solver iterations                   :c:func:`ARKodeGetNumNonlinSolvIters()`
No. of nonlinear solver convergence failures         :c:func:`ARKodeGetNumNonlinSolvConvFails()`
Single accessor to all nonlinear solver statistics   :c:func:`ARKodeGetNonlinSolvStats()`
Name of constant associated with a return flag       :c:func:`ARKodeGetReturnFlagName()`
===================================================  ============================================ 




.. c:function:: int ARKodeGetWorkSpace(void *arkode_mem, long int *lenrw, long int *leniw)

   Returns the ARKode real and integer workspace sizes.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `lenrw` -- the number of `realtype` values in the ARKode workspace.
      * `leniw` -- the number of integer values in the ARKode workspace.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetNumSteps(void *arkode_mem, long int *nsteps)

   Returns the cumulative number of internal steps taken by
   the solver (so far).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nsteps` -- number of steps taken in the solver.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetNumExpSteps(void *arkode_mem, long int *expsteps)

   Returns the cumulative number of stability-limited steps
   taken by the solver (so far).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `expsteps` -- number of stability-limited steps taken in the solver.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetNumAccSteps(void *arkode_mem, long int *accsteps)

   Returns the cumulative number of accuracy-limited steps
   taken by the solver (so far).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `accsteps` -- number of accuracy-limited steps taken in the solver.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetNumConvSteps(void *arkode_mem, long int *convsteps)

   Returns the cumulative number of convergence-limited
   steps taken by the solver (so far).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `convsteps` -- number of convergence-limited steps taken in the solver.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetNumRhsEvals(void *arkode_mem, long int *nfe_evals, long int *nfi_evals)

   Returns the number of calls to the user's right-hand
   side functions, :math:`f_E` and :math:`f_I` (so far).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nfe_evals` -- number of calls to the user's :math:`f_E(t,y)` function.
      * `nfi_evals` -- number of calls to the user's :math:`f_I(t,y)` function.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``
   
   **Notes:** The `nfi_evals` value does not account for calls made to
   :math:`f_I` by a linear solver or preconditioner module.



.. c:function:: int ARKodeGetNumLinSolvSetups(void *arkode_mem, long int *nlinsetups)

   Returns the number of calls made to the linear solver's
   setup routine (so far).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nlinsetups` -- number of linear solver setup calls made
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetNumErrTestFails(void *arkode_mem, long int *netfails)

   Returns the number of local error test failures that
   have occured (so far).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `netfails` -- number of error test failures
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetActualInitStep(void *arkode_mem, realtype *hinused)

   Returns the value of the integration step size used on the first step.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `hinused` -- actual value of initial step size
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``
   
   **Notes:** Even if the value of the initial integration step was
   specified by the user through a call to
   :c:func:`ARKodeSetInitStep()`, this value may have been changed by
   ARKode to ensure that the step size fell within the prescribed
   bounds :math:`(h_{min} \le h_0 \le h_{max})`, or to satisfy the
   local error test condition, or to ensure convergence of the
   nonlinear solver.



.. c:function:: int ARKodeGetLastStep(void *arkode_mem, realtype *hlast)

   Returns the integration step size taken on the last successful internal step.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `hlast` -- step size taken on the last internal step
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetCurrentStep(void *arkode_mem, realtype *hcur)

   Returns the integration step size to be attempted on the next internal step.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `hcur` -- step size to be attempted on the next internal step
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetCurrentTime(void *arkode_mem, realtype *tcur)

   Returns the current internal time reached by the solver.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `tcur` -- current internal time reached
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetCurrentButcherTables(void *arkode_mem, int *s, int *q, int *p, realtype *Ai, realtype *Ae, realtype *c, realtype *b, realtype *bembed)

   Returns the explicit and implicit Butcher tables
   currently in use by the solver.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `s` -- number of stages in the method.
      * `q` -- global order of accuracy of the method.
      * `p` -- global order of accuracy of the embedding.
      * `Ai` -- coefficients of DIRK method.
      * `Ae` -- coefficients of ERK method.
      * `c` -- array of internal stage times.
      * `b` -- array of solution coefficients.
      * `bembed` -- array of embedding coefficients.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``
   
   **Notes:**  The user must allocate space for `Ae` and `Ai` of size
   ``ARK_S_MAX*ARK_S_MAX``, and for `c`, `b` and `bembed` of size
   ``ARK_S_MAX``. 



.. c:function:: int ARKodeGetTolScaleFactor(void *arkode_mem, realtype *tolsfac)

   Returns a suggested factor by which the user's
   tolerances should be scaled when too much accuracy has been
   requested for some internal step.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `tolsfac` -- suggested scaling factor for user-supplied tolerances.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetErrWeights(void *arkode_mem, N_Vector eweight)

   Returns the current error weight vector.  
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `eweight` -- solution error weights at the current time.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``
   
   **Notes:** The user must allocate space for `eweight`.



.. c:function:: int ARKodeGetEstLocalErrors(void *arkode_mem, N_Vector ele)

   Returns the vector of estimated local truncation errors
   for the current step.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `ele` -- vector of estimated local truncation errors.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``
   
   **Notes:**  The user must allocate space for `ele`.
   
   The values returned in `ele` are valid only if :c:func:`ARKode()`
   returned a non-negative value.
   
   The `ele` vector, together with the `eweight` vector from
   :c:func:`ARKodeGetErrWeights()`, can be used to determine how the
   various components of the system contributed to the estimated local
   error test.  Specifically, that error test uses the RMS norm of a
   vector whose components are the products of the components of these
   two vectors.  Thus, for example, if there were recent error test
   failures, the components causing the failures are those with largest
   values for the products, denoted loosely as ``eweight[i]*ele[i]``.



.. c:function:: int ARKodeGetIntegratorStats(void *arkode_mem, long int *nsteps, long int *expsteps, long int *accsteps, long int *convsteps, long int *nfe_evals, long int *nfi_evals, long int *nlinsetups, long int *netfails, realtype *hinused, realtype *hlast, realtype *hcur, realtype *tcur)

   Returns many of the most useful integrator statistics in a single call.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nsteps` -- number of steps taken in the solver.
      * `expsteps` -- number of stability-limited steps taken in the solver.
      * `accsteps` -- number of accuracy-limited steps taken in the solver.
      * `convsteps` -- number of convergence-limited steps taken in the solver.
      * `nfe_evals` -- number of calls to the user's :math:`f_E(t,y)` function.
      * `nfi_evals` -- number of calls to the user's :math:`f_I(t,y)` function.
      * `nlinsetups` -- number of linear solver setup calls made.
      * `netfails` -- number of error test failures.
      * `hinused` -- actual value of initial step size.
      * `hlast` -- step size taken on the last internal step.
      * `hcur` -- step size to be attempted on the next internal step.
      * `tcur` -- current internal time reached.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetNumNonlinSolvIters(void *arkode_mem, long int *nniters)

   Returns the number of nonlinear solver iterations
   performed (so far).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nniters` -- number of nonlinear iterations performed.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetNumNonlinSolvConvFails(void *arkode_mem, long int *nncfails)

   Returns the number of nonlinear solver convergence
   failures that have occurred (so far).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nncfails` -- number of nonlinear convergence failures
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetNonlinSolvStats(void *arkode_mem, long int *nniters, long int *nncfails)

   Returns all of the nonlinear solver statistics in a single call.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nniters` -- number of nonlinear iterations performed.
      * `nncfails` -- number of nonlinear convergence failures
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: char *ARKodeGetReturnFlagName(long int flag)

   Returns the name of the ARKode constant corresponding to `flag`.
   
   **Arguments:**
      * `flag` -- a return flag from an ARKode function.
   
   **Return value:**  
   The return value is a string containing the name of
   the corresponding constant. 



Rootfinding optional output functions
"""""""""""""""""""""""""""""""""""""""

.. _CInterface.ARKodeRootOutputTable:

Table: Optional rootfinding outputs
"""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

===================================================  ==========================================
Optional output                                      Function name
===================================================  ==========================================
Array showing roots found                            :c:func:`ARKodeGetRootInfo()`
No. of calls to user root function                   :c:func:`ARKodeGetNumGEvals()`
===================================================  ========================================== 



.. c:function:: int ARKodeGetRootInfo(void *arkode_mem, int *rootsfound)

   Returns an array showing which functions were found to
   have a root.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `rootsfound` -- array of length `nrtfn` with the indices of the
        user functions :math:`g_i` found to have a root.  For :math:`i = 0 \ldots` `nrtfn`-1, 
        ``rootsfound[i]`` is nonzero if :math:`g_i` has a root, and 0 if not.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``
   
   **Notes:** The user must allocate space for `rootsfound`. 
   
   For the components of :math:`g_i` for which a root was found, the
   sign of ``rootsfound[i]`` indicates the direction of
   zero-crossing.  A value of +1 indicates that :math:`g_i` is
   increasing, while a value of -1 indicates a decreasing :math:`g_i`.



.. c:function:: int ARKodeGetNumGEvals(void *arkode_mem, long int *ngevals)

   Returns the cumulative number of calls made to the
   user's root function :math:`g`.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `ngevals` -- number of calls made to :math:`g` so far.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``




Direct linear solvers optional output functions
""""""""""""""""""""""""""""""""""""""""""""""""

The following optional outputs are available from the ARKDLS
modules: workspace requirements, number of calls to the Jacobian
routine, number of calls to the implicit right-hand side routine for
finite-difference Jacobian approximation, and last return value from
an ARKDLS function.  Note that, where the name of an output would
otherwise conflict with the name of an optional output from the main
solver, a suffix LS (for Linear Solver) has been added here
(e.g. `lenrwLS`). 


.. _CInterface.ARKDlsOutputTable:

Table: Optional outputs for ARKDLS
""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

===================================================  ===================================
Optional output                                      Function name
===================================================  ===================================
Size of real and integer workspaces                  :c:func:`ARKDlsGetWorkSpace()`
No. of Jacobian evaluations                          :c:func:`ARKDlsGetNumJacEvals()`
No. of `fi` calls for finite diff. Jacobian evals    :c:func:`ARKDlsGetNumRhsEvals()`
Last return flag from a linear solver function       :c:func:`ARKDlsGetLastFlag()`
Name of constant associated with a return flag       :c:func:`ARKDlsGetReturnFlagName()`
===================================================  =================================== 



    
.. c:function:: int ARKDlsGetWorkSpace(void *arkode_mem, long int *lenrwLS, long int *leniwLS)

   Returns the real and integer workspace used by the
   ARKDLS linear solver (ARKDENSE or ARKBAND).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `lenrwLS` -- the number of ``realtype`` values in the ARKDLS workspace.
      * `leniwLS` -- the number of integer values in the ARKDLS workspace.
   
   **Return value:**  
      * ARKDLS_SUCCESS if successful
      * ARKDLS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKDLS_LMEM_NULL if the linear solver memory was ``NULL``
   
   **Notes:** For the ARKDENSE linear solver, in terms of the problem
   size :math:`n`, the actual size of the real workspace is
   :math:`2n^2` ``realtype`` words, and the actual size of the integer
   workspace is :math:`n` integer words. For the ARKBAND linear
   solver, in terms of :math:`n` and the Jacobian lower and upper
   half-bandwidths :math:`m_L` and :math:`m_U`, the actual size of the
   real workspace is :math:`(2m_U + 3m_L + 2)n` ``realtype`` words,
   and the actual size of the integer workspace is :math:`n` integer
   words.



.. c:function:: int ARKDlsGetNumJacEvals(void *arkode_mem, long int *njevals)

   Returns the number of calls made to the ARKDLS
   (dense or band) Jacobian approximation routine.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `njevals` -- number of calls to the Jacobian function.
   
   **Return value:**  
      * ARKDLS_SUCCESS if successful
      * ARKDLS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKDLS_LMEM_NULL if the linear solver memory was ``NULL``



.. c:function:: int ARKDlsGetNumRhsEvals(void *arkode_mem, long int *nfevalsLS)

   Returns the number of calls made to the user-supplied
   :math:`f_I` routine due to the finite difference (dense or band)
   Jacobian approximation.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nfevalsLS` -- the number of calls made to the user-supplied
        :math:`f_I` function.
   
   **Return value:**  
      * ARKDLS_SUCCESS if successful
      * ARKDLS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKDLS_LMEM_NULL if the linear solver memory was ``NULL``
   
   **Notes:** The value of `nfevalsLS` is incremented only if hte default
   internal difference quotient function is used.



.. c:function:: int ARKDlsGetLastFlag(void *arkode_mem, long int *lsflag)

   Returns the last return value from an ARKDLS routine.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `lsflag` -- the value of the last return flag from an ARKDLS function.
   
   **Return value:**  
      * ARKDLS_SUCCESS if successful
      * ARKDLS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKDLS_LMEM_NULL if the linear solver memory was ``NULL``
   
   **Notes:** If the ARKDENSE setup function failed
   (i.e. :c:func:`ARKode()` returned ARK_LSETUP_FAIL), then the
   value of `lsflag` is equal to the column index (numbered from
   one) at which a zero diagonal element was encountered during the LU
   factorization of the (dense or banded) Jacobian matrix.  For all
   other failures, `lsflag` is negative.



.. c:function:: char *ARKDlsGetReturnFlagName(long int lsflag)

   Returns the name of the ARKDLS constant
   corresponding to `lsflag`.
   
   **Arguments:**
      * `lsflag` -- a return flag from an ARKDLS function.
   
   **Return value:**  The return value is a string containing the name of
   the corresponding constant. If 1 :math:`\le` `lsflag` :math:`\le
   n` (LU factorization failed), this routine returns "NONE". 




Iterative linear solvers optional output functions
""""""""""""""""""""""""""""""""""""""""""""""""""""

The following optional outputs are available from the ARKSPILS
modules: workspace requirements, number of linear iterations, number
of linear convergence failures, number of calls to the preconditioner
setup and solve routines, number of calls to the Jacobian-vector
product routine, number of calls to the implicit right-hand side
routine for finite-difference Jacobian-vector product approximation,
and last return value from a linear solver function.  Note that, where
the name of an output would otherwise conflict with the name of an
optional output from the main solver, a suffix LS (for Linear Solver)
has been added here (e.g. `lenrwLS`). 


.. _CInterface.ARKSpilsOutputTable:

Table: Optional outputs for ARKSPILS
""""""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

===========================================================  ====================================== 
Optional output                                              Function name
===========================================================  ====================================== 
Size of real and integer workspaces                          :c:func:`ARKSpilsGetWorkSpace()`
No. of preconditioner evaluations                            :c:func:`ARKSpilsGetNumPrecEvals()`
No. of preconditioner solves                                 :c:func:`ARKSpilsGetNumPrecSolves()`
No. of linear iterations                                     :c:func:`ARKSpilsGetNumLinIters()`
No. of linear convergence failures                           :c:func:`ARKSpilsGetNumConvFails()`
No. of Jacobian-vector product evaluations                   :c:func:`ARKSpilsGetNumJtimesEvals()`
No. of `fi` calls for finite diff. Jacobian-vector evals.    :c:func:`ARKSpilsGetNumRhsEvals()`
Last return from a linear solver function                    :c:func:`ARKSpilsGetLastFlag()`
Name of constant associated with a return flag               :c:func:`ARKSpilsGetReturnFlagName()`
===========================================================  ====================================== 




.. c:function:: int ARKSpilsGetWorkSpace(void *arkode_mem, long int *lenrwLS, long int *leniwLS)

   Returns the global sizes of the ARKSPILS real and integer workspaces.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `lenrwLS` -- the number of ``realtype`` values in the ARKSPILS workspace.
      * `leniwLS` -- the number of integer values in the ARKSPILS workspace.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``
   
   **Notes:** In terms of the problem size :math:`n` and maximum Krylov subspace
   size :math:`m`, the actual size of the real workspace is roughly:
   :math:`(m+5)n+m(m+4)+1` ``realtype`` words for ARKSPGMR,
   :math:`9n` ``realtype`` words for ARKSPBCG, and :math:`11n`
   ``realtype`` words for ARKSPTFQMR.  
   
   In a parallel setting, the above values are global, summed over all
   processors.



.. c:function:: int ARKSpilsGetNumPrecEvals(void *arkode_mem, long int *npevals)

   Returns the total number of preconditioner evaluations,
   i.e. the number of calls made to `psetup` with `jok` = ``FALSE``.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `npevals` -- the current number of calls to `psetup`.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumPrecSolves(void *arkode_mem, long int *npsolves)

   Returns the number of calls made to the preconditioner
   solve function, `psolve`.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `npsolves` -- the number of calls to `psolve`.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumLinIters(void *arkode_mem, long int *nliters)

   Returns the cumulative number of linear iterations.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nliters` -- the current number of linear iterations.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumConvFails(void *arkode_mem, long int *nlcfails)

   Returns the cumulative number of linear convergence failures.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nlcfails` -- the current number of linear convergence failures.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumJtimesEvals(void *arkode_mem, long int *njvevals)

   Returns the cumulative number of calls made to the
   Jacobian-vector function, `jtimes`.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `njvevals` -- the current number of calls to `jtimes`.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``



.. c:function:: int ARKSpilsGetNumRhsEvals(void *arkode_mem, long int *nfevalsLS)

   Returns the number of calls to the user-supplied
   implicit right-hand side function :math:`f_I` for finite difference
   Jacobian-vector product approximation.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nfevalsLS` -- the number of calls to the user implicit
        right-hand side function.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``
   
   **Notes:** The value `nfevalsLS` is incremented only if the default
   ARKSpilsDQJtimes difference quotient function is used.



.. c:function:: int ARKSpilsGetLastFlag(void *arkode_mem, long int *lsflag)

   Returns the last return value from an ARKSPILS routine.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `lsflag` -- the value of the last return flag from an
        ARKSPILS function.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``
   
   **Notes:** If the ARKSPILS setup function failed (:c:func:`ARKode()`
   returned ARK_LSETUP_FAIL), then `lsflag` will be
   SPGMR_PSET_FAIL_UNREC, SPBCG_PSET_FAIL_UNREC, or
   SPTFQMR_PSET_FAIL_UNREC. 
   
   If the ARKSPGMR solve function failed (:c:func:`ARKode()`
   returned ARK_LSOLVE_FAIL), then `lsflag` contains the error
   return flag from SpgmrSolve and will be one of:
   SPGMR_MEM_NULL, indicating that the SPGMR memory is
   ``NULL``; SPGMR_ATIMES_FAIL_UNREC, indicating an unrecoverable
   failure in the :math:`J*v` function; SPGMR_PSOLVE_FAIL_UNREC,
   indicating that the preconditioner solve function `psolve` failed
   unrecoverably; SPGMR_GS_FAIL, indicating a failure in the
   Gram-Schmidt procedure; or SPGMR_QRSOL_FAIL, indicating that
   the matrix :math:`R` was found to be singular during the QR solve
   phase. 
  
   If the ARKSPBCG solve function failed (:c:func:`ARKode()`
   returned ARK_LSOLVE_FAIL), then `lsflag` contains the error
   return flag from SpbcgSolve and will be one of:
   SPBCG_MEM_NULL, indicating that the SPBCG memory is
   ``NULL``; SPBCG_ATIMES_FAIL_UNREC, indicating an unrecoverable
   failure in the :math:`J*v` function; or
   SPBCG_PSOLVE_FAIL_UNREC, indicating that the preconditioner
   solve function `psolve` failed unrecoverably. 
   
   If the ARKSPTFQMR solve function failed (:c:func:`ARKode()`
   returned ARK_LSOLVE_FAIL), then `lsflag` contains the error
   return flag from SptfqmrSolve and will be one of:
   SPTFQMR_MEM_NULL, indicating that the SPTFQMR memory is
   ``NULL``; SPTFQMR_ATIMES_FAIL_UNREC, indicating an
   unrecoverable failure in the :math:`J*v` function; or
   SPTFQMR_PSOLVE_FAIL_UNREC, indicating that the preconditioner
   solve function `psolve` failed unrecoverably.



.. c:function:: char *ARKSpilsGetReturnFlagName(long int lsflag)

   Returns the name of the ARKSPILS constant
   corresponding to `lsflag`.
   
   **Arguments:**
      * `lsflag` -- a return flag from an ARKSPILS function.
   
   **Return value:**  
   The return value is a string containing the name of
   the corresponding constant.





.. _CInterface.Reinitialization:

ARKode reinitialization function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function :c:func:`ARKodeReInit()` reinitializes the main ARKode solver
for the solution of a problem, where a prior call to
:c:func:`ARKodeInit()` been made. The new problem must have the same
size as the previous one. ARKodeReInit performs the same input
checking and initializations that :c:func:`ARKodeInit()` does, but does
no memory allocation as it assumes that the existing internal memory
is sufficient for the new problem. 

The use of ARKodeReInit requires that the number of Runge Kutta
stages, denoted by `s`, be no larger for the new problem than for
the previous problem.  This condition is automatically fulfilled if
the method order `q` and the problem type (explicit, implicit, ImEx)
are left unchanged.  If there are changes to the linear solver
specifications, make the appropriate ARK*Set* calls, as described
in the section :ref:`CInterface.LinearSolvers`.



.. c:function:: int ARKodeReInit(void *arkode_mem, ARKRhsFn fe, ARKRhsFn fi, realtype t0, N_Vector y0)

   Provides required problem specifications and
   reinitializes ARKode.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `fe` -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the explicit portion of the right-hand side function in 
        :math:`\dot{y} = f_E(t,y) + f_I(t,y)` 
      * `fi` -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the implicit portion of the right-hand side function in 
        :math:`\dot{y} = f_E(t,y) + f_I(t,y)`
      * `t0` -- the initial value of :math:`t`
      * `y0` -- the initial condition vector :math:`y(t_0)`
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL  if the ARKode memory was ``NULL``
      * ARK_MEM_FAIL  if a memory allocation failed
      * ARK_ILL_INPUT if an argument has an illegal value.
   
   **Notes:** If an error occurred, ARKodeReInit also sends an error
   message to the error handler function.




.. _CInterface.UserSupplied:

User-supplied functions
-----------------------

The user-supplied functions consist of at least one function defining
the ODE, (optionally) a function that handles error and warning
messages, (optionally) a function that provides the error weight
vector, (optionally) a function that handles adaptive time step error
control, (optionally) a function that handles explicit time step
stability, (optionally) a function that defines the root-finding
problem(s) to solve, (optionally) a function that provides
Jacobian-related information for the linear solver (if Newton
iteration is chosen), and (optionally) one or two functions that
define the preconditioner for use in any of the Krylov iterative
algorithms.


.. _CInterface.ODERHS:

ODE right-hand side
^^^^^^^^^^^^^^^^^^^

The user must supply at least one function of type :c:func:`ARKRhsFn()` to
specify the explicit and/or implicit portions of the ODE system:



.. c:function:: typedef int (*ARKRhsFn)(realtype t, N_Vector y, N_Vector ydot, void *user_data)

   These functions compute the ODE right-hand side for a given
   value of the independent variable :math:`t` and state vector :math:`y`.
   
   **Arguments:**
      * `t` -- the current value of the independent variable.
      * `y` -- the current value of teh dependent variable vector, :math:`y(t)`.
      * `ydot` -- the output vector that forms a portion of the ODE RHS :math:`f_E(t,y) + f_I(t,y)`
      * `user_data` -- the `user_data` pointer that was passed to :c:func:`ARKodeSetUserData()`.
   
   **Return value:** 
   An ARKRhsFn should return 0 if successful, a
   positive value if a recoverable error occurred (in which case
   ARKode will attempt to correct), or a negative value if it
   failed unrecoverably (in which case the integration is halted and
   ARK_RHSFUNC_FAIL is returned).
   
   **Notes:** Allocation of memory for `ydot` is handled within
   ARKode. A recoverable failure error return from the
   ARKRhsFn is typically used to flag a value of the dependent
   variable :math:`y` that is "illegal" in some way (e.g., negative
   where only a nonnegative value is physically meaningful).  If such
   a return is made, ARKode will attempt to recover (possibly
   repeating the Newton iteration, or reducing the step size) in order
   to avoid this recoverable error return.  There are some situations
   in which recovery is not possible even if the right-hand side
   function returns a recoverable error flag.  One is when this occurs
   at the very first call to the ARKRhsFn (in which case
   ARKode returns ARK_FIRST_RHSFUNC_ERR).  Another is when a
   recoverable error is reported by ARKRhsFn after the integrator
   completes a successful stage, in which case ARKode returns
   ARK_UNREC_RHSFUNC_ERR).



.. _CInterface.ErrorHandler:

Error message handler function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As an alternative to the default behavior of directing error and
warning messages to the file pointed to by `errfp` (see
:c:func:`ARKodeSetErrFile()`), the user may provide a function of type
:c:func:`ARKErrHandlerFn()` to process any such messages. 



.. c:function:: typedef void (*ARKErrHandlerFn)(int error_code, const char *module, const char *function, char *msg, void *user_data)

   This function processes error and warning messages from
   ARKode and is sub-modules.
   
   **Arguments:**
      * `error_code` -- the error code.
      * `module` -- the name of the ARKode module reporting the error.
      * `function` -- the name of the function in which the error occurred.
      * `msg` -- the error message.
      * `user_data` -- a pointer to user data, the same as the
        `eh_data` parameter that was passed to :c:func:`ARKodeSetErrHandlerFn()`.
   
   **Return value:** 
   An ARKErrHandlerFn function has no return value.
   
   **Notes:** `error_code` is negative for errors and positive
   (ARK_WARNING) for warnings.  If a function that returns a
   pointer to memory encounters an error, it sets `error_code` to
   0.




.. _CInterface.ErrorWeight:

Error weight function
^^^^^^^^^^^^^^^^^^^^^^

As an alternative to providing the relative and absolute tolerances,
the user may provide a function of type :c:func:`ARKEwtFn()` to compute a
vector `ewt` containing the weights in the WRMS norm
:math:`\|v\|_{WRMS} = \left(\frac{1}{n} \sum_{i=1}^n \left(ewt_i * v_i\right)^2
\right)^{1/2}`.  These weights will be used in place of those defined
in the section :ref:`Mathematics`.



.. c:function:: typedef int (*ARKEwtFn)(N_Vector y, N_Vector ewt, void *user_data)

   This function computes the WRMS error weights for the vector
   :math:`y`.
   
   **Arguments:**
      * `y` -- the dependent variable vector at which the
        weight vector is to be computed.
      * `ewt` -- the output vector containing the error weights.
      * `user_data` -- a pointer to user data, the same as the
        `user_data` parameter that was passed to :c:func:`ARKodeSetUserData()`.
   
   **Return value:** 
   An ARKEwtFn function must return 0 if it
   successfully set the error weights, and -1 otherwise.
   
   **Notes:** Allocation of memory for `ewt` is handled within ARKode.
   
   The error weight vector must have all components positive.  It is
   the user's responsibility to perform this test and return -1 if it
   is not satisfied.



.. _CInterface.AdaptivityFn:

Time step adaptivity function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As an alternative to using one of the built-in time step adaptivity
methods for controlling solution error, the user may provide a
function of type :c:func:`ARKAdaptFn()` to compute a target step size
:math:`h` for the next integration step.  These steps should be chosen
as the maximum value such that the error estimates remain below 1.



.. c:function:: typedef int (*ARKAdaptFn)(N_Vector y, realtype t, realtype h, realtype e1, realtype e2,  realtype e3, int q, int p, realtype *hnew, void *user_data)

   This function implements a time step adaptivity algorithm
   that chooses :math:`h` satisfying the error tolerances..
   
   **Arguments:**
      * `y` -- the current value of the dependent variable vector, :math:`y(t)`.
      * `t` -- the current value of the independent variable.
      * `h` -- the current value of the step size.
      * `e1` -- the error estimate from the current step, :math:`m`.
      * `e2` -- the error estimate from the previous step, :math:`m-1`.
      * `e3` -- the error estimate from the step :math:`m-2`.
      * `q` -- the global order of accuracy for the integration method.
      * `p` -- the global order of accuracy for the embedding.
      * `hnew` -- the output value of the next step size.
      * `user_data` -- a pointer to user data, the same as the
        `h_data` parameter that was passed to :c:func:`ARKodeSetAdaptivityFn()`.
   
   **Return value:** 
   An ARKAdaptFn function should return 0 if it
   successfuly set the next step size, and a non-zero value otherwise.




.. _CInterface.StabilityFn:

Explicit stability function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A user may supply a function to predict the maximum stable step size
for the explicit portion of the ImEx system, :math:`f_E(t,y)`.  While
the accuracy-based time step adaptivity algorithms may be sufficient
for retaining a stable solution to the ODE system, these may be
inefficient if :math:`f_E(t,y)` contains moderately stiff terms.  In
this scenario, a user may provide a function of type :c:func:`ARKExpStabFn()`
to provide this stability information to ARKode.  This function
must set the scalar step size satisfying the stability restriction for
the upcoming time step.  This value will subsequently be bounded by
the user-supplied values for the minimum and maximum allowed time
step, and the accuracy-based time step.  



.. c:function:: typedef int (*ARKExpStabFn)(N_Vector y, realtype t, realtype *hstab, void *user_data)

   This function predicts the maximum stable step size for the
   explicit portions of the ImEx ODE system.
   
   **Arguments:**
      * `y` -- the current value of the dependent variable vector, :math:`y(t)`.
      * `t` -- the current value of the independent variable
      * `hstab` -- the output value with the maximum stable step size.
      * `user_data` -- a pointer to user data, the same as the
        `estab_data` parameter that was passed to :c:func:`ARKodeSetStabilityFn()`.
   
   **Return value:** 
   An ARKExpStabFn function should return 0 if it
   successfully set the upcoming stable step size, and a non-zero
   value otherwise.
   
   **Notes:**  If this function is not supplied, or if it returns `hstab =
   0.0`, then ARKode will assume that there is no explicit
   stability restriction on the time step size.



.. _CInterface.RootfindingFn:

Rootfinding function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a rootfinding problem is to be solved during the integration of the
ODE system, the user must supply a function of type :c:func:`ARKRootFn()`.



.. c:function:: typedef int (*ARKRootFn)(realtype t, N_Vector y, realtype *gout, void *user_data)

   This function implements a vector-valued function
   :math:`g(t,y)` such that the roots of the `nrtfn` components
   :math:`g_i(t,y)` are sought.
   
   **Arguments:**
      * `t` -- the current value of the independent variable
      * `y` -- the current value of the dependent variable vector, :math:`y(t)`.
      * `gout` -- the output array, of length `nrtfn`, with components :math:`g_i(t,y)`.
      * `user_data` -- a pointer to user data, the same as the
        `user_data` parameter that was passed to :c:func:`ARKodeSetUserData()`.
   
   **Return value:** 
   An ARKRootFn function should return 0 if successful
   or a non-zero value if an error occurred (in which case the
   integration is halted and ARKode returns ARK_RTFUNC_FAIL).
   
   **Notes:** Allocation of memory for `gout` is handled within ARKode.



.. _CInterface.DenseJacobianFn:

Jacobian information (direct method with dense Jacobian)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the direct linear solver with dense treatment of the Jacobian is
used (i.e., :c:func:`ARKDense()` or :c:func:`ARKLapackDense()` is called in Step 8 of
the section :ref:`CInterface.Skeleton`), the user may provide a
function of type :c:func:`ARKDlsDenseJacFn()` to provide the Jacobian
approximation. 



.. c:function:: typedef int (*ARKDlsDenseJacFn)(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the dense Jacobian :math:`J =
   \frac{\partial f_I}{\partial y}` (or an approximation to it).
   
   **Arguments:**
      * `N` -- the size of the ODE system.
      * `t` -- the current value of the independent variable
      * `y` -- the current value of the dependent variable vector, namely
        the predicted value of :math:`y(t)`.
      * `fy` -- the current value of the vector :math:`f_I(t,y)`.
      * `Jac` -- the output dense Jacobian matrix (of type ``DlsMat``).
      * `user_data` -- a pointer to user data, the same as the
        `user_data` parameter that was passed to :c:func:`ARKodeSetUserData()`.
      * `tmp1`, `tmp2`, `tmp3` -- pointers to memory allocated to
        variables of type ``N_Vector`` which can be used by an
        ARKDlsDenseJacFn as temporary storage or work space.
   
   **Return value:** 
   An ARKDlsDenseJacFn function should return 0 if
   successful, a positive value if a recoverable error occurred (in
   which case ARKode will attempt to correct, while ARKDENSE
   sets `last_flag` to ARKDLS_JACFUNC_RECVR), or a negative
   value if it failed unrecoverably (in which case the integration is
   halted, :c:func:`ARKode()` returns ARK_LSETUP_FAIL and
   ARKDENSE sets `last_flag` to ARKDLS_JACFUNC_UNRECVR). 
   
   **Notes:** A user-supplied dense Jacobian function must load the `N` by
   `N` dense matrix `Jac` with an approximation to the Jacobian
   matrix :math:`J(t,y)` at the point :math:`(t,y)`. Only nonzero
   elements need to be loaded into `Jac` because `Jac` is set to
   the zero matrix before the call to the Jacobian function. The type
   of `Jac` is ``DlsMat``. 
   
   The accessor macros ``DENSE_ELEM`` and ``DENSE_COL`` allow the user
   to read and write dense matrix elements without making explicit
   references to the underlying representation of the ``DlsMat``
   type. ``DENSE_ELEM(J,i,j)`` references the ``(i,j)``-th element of
   the dense matrix ``J`` (for ``i``, ``j`` between 0 and
   N-1). This macro is meant for small problems for which
   efficiency of access is not a major concern. Thus, in terms of the
   indices :math:`m` and :math:`n` ranging from 1 to `N`, the
   Jacobian element :math:`J_{m,n}` can be set using the statement
   ``DENSE_ELEM(J, m-1, n-1)`` :math:`= J_{m,n}`. Alternatively,
   ``DENSE_COL(J,j)`` returns a pointer to the first element of the
   ``j``-th column of ``J`` (for ``j`` ranging from 0 to `N`-1),
   and the elements of the ``j``-th column can then be accessed using
   ordinary array indexing. Consequently, :math:`J_{m,n}` can be
   loaded using the statements ``col_n = DENSE_COL(J, n-1);
   col_n[m-1]`` :math:`= J_{m,n}`. For large problems, it is more
   efficient to use ``DENSE_COL`` than to use ``DENSE_ELEM``. Note
   that both of these macros number rows and columns starting from 0. 
   
   The ``DlsMat`` type and accessor macros ``DENSE_ELEM`` and
   ``DENSE_COL`` are documented in the section :ref:`LinearSolvers`.
   
   If the user's ARKDenseJacFn function uses difference quotient
   approximations, then it may need to access quantities not in the
   argument list. These include the current step size, the error
   weights, etc..  To obtain these, use the ARKodeGet* functions
   listed in :ref:`CInterface.ARKodeOutputTable`. The unit roundoff
   can be accessed as ``UNIT_ROUNDOFF``, which is defined in the
   header file ``sundials_types.h``.
   
   For the sake of uniformity, the argument `N` is of type ``long int``,
   even in the case that the LAPACK dense solver is to be used. 



.. _CInterface.BandJacobianFn:

Jacobian information (direct method with banded Jacobian)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the direct linear solver with banded treatment of the Jacobian is
used (i.e. :c:func:`ARKBand()` or :c:func:`ARKLapackBand()` is called in Step 8 of the
section :ref:`CInterface.Skeleton`), the user may provide a function
of type :c:func:`ARKDlsBandJacFn()` to provide the Jacobian approximation.



.. c:function:: typedef int (*ARKDlsBandJacFn)(long int N, long int mupper, long int mlower, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the banded Jacobian :math:`J =
   \frac{\partial f_I}{\partial y}` (or an approximation to it).
   
   **Arguments:**
      * `N` -- the size of the ODE system.
      * `mlower`, `mupper` -- the lower and upper half-bandwidths of
        the Jacobian.
      * `t` -- the current value of the independent variable
      * `y` -- the current value of the dependent variable vector, namely
        the predicted value of :math:`y(t)`.
      * `fy` -- the current value of the vector :math:`f_I(t,y)`.
      * `Jac` -- the output dense Jacobian matrix (of type ``DlsMat``).
      * `user_data` -- a pointer to user data, the same as the
        `user_data` parameter that was passed to :c:func:`ARKodeSetUserData()`.
      * `tmp1`, `tmp2`, `tmp3` -- pointers to memory allocated to
        variables of type ``N_Vector`` which can be used by an
        ARKDlsBandJacFn as temporary storage or work space.
   
   **Return value:** 
   An ARKDlsBandJacFn function should return 0 if
   successful, a positive value if a recoverable error occurred (in
   which case ARKode will attempt to correct, while ARKBAND
   sets `last_flag` to ARKDLS_JACFUNC_RECVR), or a negative
   value if it failed unrecoverably (in which case the integration is
   halted, :c:func:`ARKode()` returns ARK_LSETUP_FAIL and
   ARKBAND sets `last_flag` to ARKDLS_JACFUNC_UNRECVR). 
   
   **Notes:** A user-supplied banded Jacobian function must load the band
   matrix `Jac` of type ``DlsMat`` with the elements of the Jacobian
   :math:`J(t,y)` at the point :math:`(t,y)`. Only nonzero elements
   need to be loaded into `Jac` because `Jac` is initialized to
   the zero matrix before the call to the Jacobian function. 
  
   The accessor macros ``BAND_ELEM``, ``BAND_COL``, and
   ``BAND_COL_ELEM`` allow the user to read and write band matrix
   elements without making specific references to the underlying
   representation of the ``DlsMat`` type.  ``BAND_ELEM(J, i, j)``
   references the ``(i,j)``-th element of the band matrix ``J``,
   counting from 0. This macro is meant for use in small problems for
   which efficiency of access is not a major concern. Thus, in terms
   of the indices :math:`m` and :math:`n` ranging from 1 to `N` with
   :math:`(m, n)` within the band defined by `mupper` and
   `mlower`, the Jacobian element :math:`J_{m,n}` can be loaded
   using the statement ``BAND_ELEM(J, m-1, n-1)`` :math:`=
   J_{m,n}`. The elements within the band are those with `-mupper`
   :math:`\le m-n \le` `mlower`.  Alternatively, ``BAND_COL(J, j)``
   returns a pointer to the diagonal element of the ``j``-th column of
   ``J``, and if we assign this address to ``realtype *col_j``, then
   the ``i``-th element of the ``j``-th column is given by
   ``BAND_COL_ELEM(col_j, i, j)``, counting from 0. Thus, for
   :math:`(m,n)` within the band, :math:`J_{m,n}` can be loaded by
   setting ``col_n = BAND_COL(J, n-1); BAND_COL_ELEM(col_n, m-1,
   n-1)`` :math:`= J_{m,n}` . The elements of the ``j``-th column can
   also be accessed via ordinary array indexing, but this approach
   requires knowledge of the underlying storage for a band matrix of
   type ``DlsMat``. The array ``col_n`` can be indexed from
   `-mupper` to `mlower`. For large problems, it is more efficient
   to use ``BAND_COL`` and ``BAND_COL_ELEM`` than to use the
   ``BAND_ELEM`` macro. As in the dense case, these macros all number
   rows and columns starting from 0. 
   
   The ``DlsMat`` type and the accessor macros ``BAND_ELEM``,
   ``BAND_COL`` and ``BAND_COL_ELEM`` are documented in the section 
   :ref:`LinearSolvers`.

   If the user's ARKBandJacFn function uses difference quotient
   approximations, then it may need to access quantities not in the
   argument list.  These include the current step size, the error
   weights, etc.. To obtain these, use the ARKodeGet* functions
   listed in :ref:`CInterface.ARKodeOutputTable`. The unit roundoff
   can be accessed as ``UNIT_ROUNDOFF`` defined in the header file
   ``sundials_types.h``.
   
   For the sake of uniformity, the arguments `N`, `mlower`, and
   `mupper` are of type ``long int``, even in the case that the
   LAPACK band solver is to be used.  



.. _CInterface.JTimesFn:

Jacobian information (matrix-vector product)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If one of the Krylov iterative linear solvers SPGMR, SPBCG, or
SPTFQMR is selected (i.e. ARKSp* is called in step 8 of the
section :ref:`CInterface.Skeleton`), the user may provide a function
of type :c:func:`ARKSpilsJacTimesVecFn()` in the following form, to compute
matrix-vector products :math:`J*v`. If such a function is not
supplied, the default is a difference quotient approximation to these
products. 



.. c:function:: typedef int (*ARKSpilsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp)

   This function computes the product :math:`Jv =
   \left(\frac{\partial f_I}{\partial y}\right)v` (or an approximation to it).
   
   **Arguments:**
      * `v` -- the vector to multiply.
      * `Jv` -- the output vector computed.
      * `t` -- the current value of the independent variable
      * `y` -- the current value of the dependent variable vector.
      * `fy` -- the current value of the vector :math:`f_I(t,y)`.
      * `user_data` -- a pointer to user data, the same as the
        `user_data` parameter that was passed to :c:func:`ARKodeSetUserData()`.
      * `tmp` -- pointer to memory allocated to a variable of type
        ``N_Vector`` which can be used as temporary storage or work space.
   
   **Return value:** 
   The value to be returned by the Jacobian-vector product
   function should be 0 if successful. Any other return value will
   result in an unrecoverable error of the SPILS generic solver,
   in which case the integration is halted. 
   
   **Notes:** If the user's ARKSpilsJacTimesVecFn function uses
   difference quotient approximations, it may need to access
   quantities not in the argument list.  These include the current
   step size, the error weights, etc..  To obtain these, use the
   ARKodeGet* functions listed in
   :ref:`CInterface.ARKodeOutputTable`. The unit roundoff can be
   accessed as ``UNIT_ROUNDOFF`` defined in the header file
   ``sundials_types.h``. 




.. _CInterface.PrecSolveFn:

Preconditioning (linear system solution)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If one of the Krylov iterative linear solvers SPGMR, SPBCG, or
SPTFQMR is selected, and preconditioning is used, then the user
must provide a function of type :c:func:`ARKSpilsPrecSolveFn()` to solve the
linear system :math:`Pz=r`, where :math:`P` may be either a left or
right preconditioning matrix.  Here :math:`P` should approximate (at
least crudely) the Newton matrix :math:`A=M-\gamma J`, where :math:`M`
is the mass matrix (typically :math:`M=I` unless working in a
finite-element setting) and :math:`J = \frac{\partial f_I}{\partial
y}`  If preconditioning is done on both sides, the product of the two
preconditioner matrices should approximate :math:`A`. 



.. c:function:: typedef int (*ARKSpilsPrecSolveFn)(realtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z, realtype gamma, realtype delta, int lr, void *user_data, N_Vector tmp)

   This function solves the preconditioner system :math:`Pz=r`.
   
   **Arguments:**
      * `t` -- the current value of the independent variable.
      * `y` -- the current value of the dependent variable vector.
      * `fy` -- the current value of the vector :math:`f_I(t,y)`.
      * `r` -- the right-hand side vector of the linear system.
      * `z` -- the computed output solution vector 
      * `gamma` -- the scalar :math:`\gamma` appearing in the Newton
        matrix given by :math:`A=M-\gamma J`.
      * `delta` -- an input tolerance to be used if an iterative method
        is employed in the solution.  In that case, the resdual vector
        :math:`Res = r-Pz` of the system should be made to be less than `delta`
        in the weighted :math:`l_2` norm, i.e. :math:`\left(\sum_{i=1}^n
        \left(Res_i * ewt_i\right)^2 \right)^{1/2} < \delta`, where :math:`\delta =`
        `delta`.  To obtain the ``N_Vector`` `ewt`, call
        :c:func:`ARKodeGetErrWeights()`. 
      * `lr` -- an input flag indicating whether the preconditioner
        solve is to use the left preconditioner (`lr = 1`) or the right
        preconditioner (`lr = 2`).
      * `user_data` -- a pointer to user data, the same as the
        `user_data` parameter that was passed to :c:func:`ARKodeSetUserData()`.
      * `tmp` -- pointer to memory allocated to a variable of type
        ``N_Vector`` which can be used as temporary storage or work space.
   
   **Return value:** 
   The value to be returned by the preconditioner solve
   function is a flag indicating whether it was successful. This value
   should be 0 if successful, positive for a recoverable error (in
   which case the step will be retried), or negative for an
   unrecoverable error (in which case the integration is halted).  




.. _CInterface.PrecSetupFn:

Preconditioning (Jacobian data)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the user's preconditioner requires that any Jacobian-related data
be preprocessed or evaluated, then these actions need to occur within
a user-supplied function of type :c:func:`ARKSpilsPrecSetupFn()`. 



.. c:function:: typedef int (*ARKSpilsPrecSetupFn)(realtype t, N_Vector y, N_Vector fy, booleantype jok, booleantype *jcurPtr, realtype gamma, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function preprocesses and/or evaluates Jacobian-related
   data needed by the preconditioner.
   
   **Arguments:**
      * `t` -- the current value of the independent variable.
      * `y` -- the current value of the dependent variable vector.
      * `fy` -- the current value of the vector :math:`f_I(t,y)`.
      * `jok` -- is an input flag indicating whether the Jacobian-related
        data needs to be updated. The `jok` argument provides for the
        reuse of Jacobian data in the preconditioner solve function. When
        `jok` = ``FALSE``, the Jacobian-related data should be recomputed
        from scratch. When `jok` = ``TRUE`` the Jacobian data, if saved from the
        previous call to this function, can be reused (with the current
        value of `gamma`). A call with `jok` = ``TRUE`` can only occur
        after a call with `jok` = ``FALSE``. 
      * `jcurPtr` -- is a pointer to a flag which should be set to
        ``TRUE`` if Jacobian data was recomputed, or set to ``FALSE`` if
        Jacobian data was not recomputed, but saved data was still reused. 
      * `gamma` -- the scalar :math:`\gamma` appearing in the Newton
        matrix given by :math:`A=M-\gamma J`.
      * `user_data` -- a pointer to user data, the same as the
        `user_data` parameter that was passed to :c:func:`ARKodeSetUserData()`.
      * `tmp1`, `tmp2`, `tmp3` -- pointers to memory allocated to
        variables of type ``N_Vector`` which can be used as temporary
        storage or work space.
   
   **Return value:** 
   The value to be returned by the preconditioner setup
   function is a flag indicating whether it was successful. This value
   should be 0 if successful, positive for a recoverable error (in
   which case the step will be retried), or negative for an
   unrecoverable error (in which case the integration is halted). 
   
   **Notes:**  The operations performed by this function might include
   forming a crude approximate Jacobian, and performing an LU
   factorization of the resulting approximation to :math:`A = M -
   \gamma J`. 
   
   Each call to the preconditioner setup function is preceded by a
   call to the implicit :c:func:`ARKRhsFn()` user function with the same
   :math:`(t,y)` arguments.  Thus, the preconditioner setup function can
   use any auxiliary data that is computed and saved during the
   evaluation of the ODE right-hand side. 
   
   This function is not called in advance of every call to the
   preconditioner solve function, but rather is called only as often
   as needed to achieve convergence in the Newton iteration. 
   
   If the user's ARKSpilsPrecSetupFn function uses difference
   quotient approximations, it may need to access quantities not in
   the call list. These include the current step size, the error
   weights, etc. To obtain these, use the ARKodeGet* functions
   listed in :ref:`CInterface.ARKodeOutputTable`. The unit roundoff
   can be accessed as ``UNIT_ROUNDOFF`` defined in the header file
   ``sundials_types.h``. 




.. _CInterface.PreconditionerModules:

Preconditioner modules
-----------------------

The efficiency of Krylov iterative methods for the solution of linear
systems can be greatly enhanced through preconditioning. For problems
in which the user cannot define a more effective, problem-specific
preconditioner, ARKode provides a banded preconditioner in the
module ARKBANDPRE and a band-block-diagonal preconditioner module
ARKBBDPRE. 


.. _CInterface.BandPre:

A serial banded preconditioner module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This preconditioner provides a band matrix preconditioner for use with
any of the Krylov iterative linear solvers, when used in a serial
setting (i.e. with the NVECTOR_SERIAL module). It uses difference
quotients of the ODE right-hand side function :math:`f_I` to generate
a band matrix of bandwidth ``ml + mu + 1``, where the number of
super-diagonals (``mu``, the upper half-bandwidth) and sub-diagonals
(``ml``, the lower half-bandwidth) are specified by the user, and uses
this to form a preconditioner for use with the Krylov linear
solver. Although this matrix is intended to approximate the Jacobian
:math:`\frac{\partial f_I}{\partial y}`, it may be a very crude
approximation. The true Jacobian need not be banded, or its true
bandwidth may be larger than ``ml + mu + 1``, as long as the banded
approximation generated here is sufficiently accurate to speed
convergence as a preconditioner. 

In order to use the ARKBANDPRE module, the user need not define
any additional functions. Aside from the header files required for the
integration of the ODE problem (see the section
:ref:`CInterface.Headers`), to use the ARKBANDPRE module, the main
program must include the header file ``arkode_bandpre.h`` which
declares the needed function prototypes.  The following is a summary
of the usage of this module.  Steps that are unchanged from the
skeleton program presented in :ref:`CInterface.Skeleton` are
`italicized`. 

1. `Set problem dimensions`

2. `Set vector of initial values` 

3. `Create ARKode object` 

4. `Allocate internal memory` 

5. `Set optional inputs` 

6. Attach iterative linear solver, one of:

   (a) ``flag = ARKSpgmr(arkode_mem, pretype, maxl);`` 

   (b) ``flag = ARKSpbcg(arkode_mem, pretype, maxl);``

   (c) ``flag = ARKSptfqmr(arkode_mem, pretype, maxl);``

7. Initialize the ARKBANDPRE preconditioner module 

   Specify the upper and lower half-bandwidths (``mu`` and ``ml``,
   respectively) and call 

   ``flag = ARKBandPrecInit(arkode_mem, N, mu, ml);``

   to allocate memory and initialize the internal preconditioner
   data. 

8. `Set linear solver optional inputs`

    Note that the user should not overwrite the preconditioner setup
    function or solve function through calls to the ARKSpilsSet*
    optional input functions. 

9. `Advance solution in time`

10. Get optional outputs 

   Additional optional outputs associated with ARKBANDPRE are
   available by way of the two routines described below,
   :c:func:`ARKBandPrecGetWorkSpace()` and
   :c:func:`ARKBandPrecGetNumRhsEvals()`.  

11. `Deallocate memory for solution vector`

12. `Free solver memory`

The ARKBANDPRE preconditioner module is initialized and attached
by calling the following function:



.. c:function:: int ARKBandPrecInit(void *arkode_mem, long int N, long int mu, long int ml)

   Initializes the ARKBANDPRE preconditioner and
   allocates required (internal) memory for it.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `N` -- problem dimension (size of ODE system).
      * `mu` -- upper half-bandwidth of the Jacobian approximation.
      * `ml` -- lower half-bandwidth of the Jacobian approximation.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_ILL_INPUT if an input has an illegal value
      * ARKSPILS_MEM_FAIL if a memory allocation request failed

   **Notes:** The banded approximate Jacobian will have nonzero elements
   only in locations :math:`(i,j)` with `ml` :math:`\le j-i \le` `mu`.


The following two optional output functions are available for use with
the ARKBANDPRE module:



.. c:function:: int ARKBandPrecGetWorkSpace(void *arkode_mem, long int *lenrwLS, long int *leniwLS)

   Returns the sizes of the ARKBANDPRE real and integer
   workspaces.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `lenrwLS` -- the number of ``realtype`` values in the
        ARKBANDPRE workspace.
      * `leniwLS` -- the number of integer values in the  ARKBANDPRE workspace.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_PMEM_NULL if the preconditioner memory is ``NULL``
   
   **Notes:** In terms of the problem size :math:`N` and `smu` :math:`=
   \min(N-1,` `mu+ml` :math:`)`, the actual size of the real
   workspace is :math:`(2` `ml + mu + smu` :math:`+2)N` ``realtype``
   words, and the actual size of the integer workspace is :math:`N`
   integer words.
   
   The workspaces referred to here exist in addition to those given by
   the corresponding function ARKSpils*GetWorkspace.



.. c:function:: int ARKBandPrecGetNumRhsEvals(void *arkode_mem, long int *nfevalsBP)

   Returns the number of calls made to the user-supplied
   right-hand side function :math:`f_I` for constructing the
   finite-difference banded Jacobian approximation used within the
   preconditioner setup function.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nfevalsBP` -- number of calls to :math:`f_I`
   
   **Return value:**  
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_PMEM_NULL if the preconditioner memory is ``NULL``
   
   **Notes:**  The counter `nfevalsBP` is distinct from the counter
   `nfevalsLS` returned by the corresponding function
   ARKSpils*GetNumRhsEvals and also from `nfi_evals` returned by
   :c:func:`ARKodeGetNumRhsEvals()`.  The total number of right-hand
   side function evaluations is the sum of all three of these
   counters, plus the `nfe_evals` counter for :math:`f_E` calls
   returned by :c:func:`ARKodeGetNumRhsEvals()`.



.. _CInterface.BBDPre:

A parallel band-block-diagonal preconditioner module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A principal reason for using a parallel ODE solver such as ARKode
lies in the solution of partial differential equations
(PDEs). Moreover, the use of a Krylov iterative method for the
solution of many such problems is motivated by the nature of the
underlying linear system of equations that must be solved at each time
step. For many PDEs, the linear algebraic system is large, sparse and
structured.  However, if a Krylov iterative method is to be effective
in this setting, then a nontrivial preconditioner needs to be used.
Otherwise, the rate of convergence of the Krylov iterative method is
usually unacceptably slow.  Unfortunately, an effective preconditioner
tends to be problem-specific.  However, we have developed one type of
preconditioner that treats a rather broad class of PDE-based
problems. It has been successfully used with CVODE for several
realistic, large-scale problems [HT1998]_ and is included
in a software module within the ARKode package. This module works
with the parallel vector module NVECTOR_PARALLEL and is usable
with any of the Krylov iterative linear solvers. It generates a
preconditioner that is a block-diagonal matrix with each block being a
band matrix. The blocks need not have the same number of super- and
sub-diagonals and these numbers may vary from block to block. This
Band-Block-Diagonal Preconditioner module is called ARKBBDPRE. 

One way to envision these preconditioners is to think of the domain of
the computational PDE problem as being subdivided into :math:`Q`
non-overlapping subdomains. Each of these subdomains is then assigned
to one of the :math:`Q` processes to be used to solve the ODE
system. The basic idea is to isolate the preconditioning so that it is
local to each process, and also to use a (possibly cheaper)
approximate right-hand side function. This requires the definition of
a new function :math:`g(t,y)` which approximates the function
:math:`f_I(t,y)` in the definition of the ODE system, 

.. math::
   M\dot{y} = f_E(t,y) + f_I(t,y),

where :math:`f_I` corresponds to the ODE components to be treated
implicitly.  However, the user may set :math:`g = f_I`. Corresponding
to the domain decomposition, there is a decomposition of the solution
vector :math:`y` into :math:`Q` disjoint blocks :math:`y_q`, and a
decomposition of :math:`g` into blocks :math:`g_q`. The block
:math:`g_q` depends both on :math:`y_p` and on components of blocks
:math:`y_{q'}`	associated with neighboring subdomains (so-called
ghost-cell data). Let :math:`\bar{y}_q` denote :math:`y_q` augmented
with those other components on which :math:`g_q` depends. Then we have 

.. math::
   g(t,y) = \left[ g_1(t,\bar{y}_1), g_2(t,\bar{y}_2), \ldots , g_Q(t,\bar{y}_Q) \right]^T

and each of the blocks :math:`g_q(t,\bar{y}_q)` is decoupled from the
others.

The preconditioner associated with this decomposition has the form

.. math::
   P = diag[P_1, P_2, \ldots, P_Q]

where

.. math::
   P_q \approx M - \gamma J_q

and where :math:`J_q` is a difference quotient approximation to
:math:`\frac{\partial g_q}{\partial \bar{y}_q}`.  This matrix is taken
to be banded, with upper and lower half-bandwidths `mudq` and
`mldq` defined as the number of non-zero diagonals above and below
the main diagonal, respectively.  The difference quotient
approximation is computed using `mudq + mldq + 2` evaluations of
:math:`g_m`, but only a matrix of bandwidth `mukeep + mlkeep + 1` is
retained. Neither pair of parameters need be the true half-bandwidths
of the Jacobian of the local block of :math:`g`, if smaller values
provide a more efficient preconditioner. The solution of the complete
linear system 

.. math::
   Px = b

reduces to solving each of the distinct equations

.. math::
   P_q x_q = b_q, \quad q=1,\ldots,Q,

and this is done by banded LU factorization of :math:`P_q` followed by
a banded backsolve.

Similar block-diagonal preconditioners could be considered with
different treatments of the blocks :math:`P_q`.  For example,
incomplete LU factorization or an iterative method could be used
instead of banded LU factorization.

The ARKBBDPRE module calls two user-provided functions to
construct :math:`P`: a required function `gloc` (of type
:c:func:`ARKLocalFn()`) which approximates the right-hand side function
:math:`g(t,y) \approx f_I(t,y)` and which is computed locally, and an
optional function `cfn` (of type :c:func:`ARKCommFn()`) which performs all
interprocess communication necessary to evaluate the approximate
right-hand side :math:`g`. These are in addition to the user-supplied
right-hand side function :math:`f_I`. Both functions take as input the
same pointer `user_data` that is passed by the user to
:c:func:`ARKodeSetUserData()` and that was passed to the user's
function :math:`f_I`. The user is responsible for providing space
(presumably within `user_data`) for components of :math:`y` that are
communicated between processes by `cfn`, and that are then used by
`gloc`, which should not do any communication.



.. c:function:: typedef int (*ARKLocalFn)(long int Nlocal, realtype t, N_Vector y, N_Vector glocal, void *user_data)

   This `gloc` function computes :math:`g(t,y)`.  It
   loads the vector `glocal` as a function of `t` and `y`.
   
   **Arguments:**
      * `Nlocal` -- the local vector length
      * `t` -- the value of the independent variable
      * `y` -- the value of the dependent variable vector on this process
      * `glocal` -- the output vector of :math:`g(t,y)` on this process
      * `user_data` -- a pointer to user data, the same as the
        `user_data` parameter passed to :c:func:`ARKodeSetUserData()`.
   
   **Return value:**  
   An ARKLocalFn should return 0 if successful, a
   positive value if a recoverable error occurred (in which case
   ARKode will attempt to correct), or a negative value if it
   failed unrecoverably (in which case the integration is halted and
   :c:func:`ARKode()` will return ARK_LSETUP_FAIL).
   
   **Notes:**  This function must assume that all interprocess communication
   of data needed to calculate `glocal` has already been done, and that
   this data is accessible within user data. 
   
   The case where :math:`g` is mathematically identical to :math:`f_I`
   is allowed. 



.. c:function:: typedef int (*ARKCommFn)(long int Nlocal, realtype t, N_Vector y, void *user_data)

   This `cfn` function performs all interprocess
   communication necessary for the executation of the `gloc` function
   above, using the input vector `y`.
   
   **Arguments:**
      *  `Nlocal` -- the local vector length
      * `t` -- the value of the independent variable
      * `y` -- the value of the dependent variable vector on this process
      * `user_data` -- a pointer to user data, the same as the
        `user_data` parameter passed to :c:func:`ARKodeSetUserData()`.
   
   **Return value:**  
   An ARKCommFn should return 0 if successful, a
   positive value if a recoverable error occurred (in which case
   ARKode will attempt to correct), or a negative value if it
   failed unrecoverably (in which case the integration is halted and
   :c:func:`ARKode()` will return ARK_LSETUP_FAIL).
   
   **Notes:**  The `cfn` function is expected to save communicated data in
   space defined within the data structure `user_data`.
   
   Each call to the `cfn` function is preceded by a call to the
   right-hand side function :math:`f_I` with the same :math:`(t,y)`
   arguments. Thus, `cfn` can omit any communication done by
   :math:`f_I` if relevant to the evaluation of `glocal`. If all
   necessary communication was done in :math:`f_I`, then `cfn` =
   ``NULL`` can be passed in the call to :c:func:`ARKBBDPrecInit()` (see
   below).



Besides the header files required for the integration of the ODE problem (see the section
:ref:`CInterface.Headers`), to use the ARKBBDPRE module, the main
program must include the header file ``arkode_bbdpre.h`` which
declares the needed function prototypes. 

The following is a summary of the proper usage of this module. Steps
that are unchanged from the skeleton program presented in
:ref:`CInterface.Skeleton` are `italicized`.

1. `Initialize MPI`

2. `Set problem dimensions`

3. `Set vector of initial values`

4. `Create ARKode object`

5. `Allocate internal memory`

6. `Set optional inputs`

7. Attach iterative linear solver, one of:

   (a) ``flag = ARKSpgmr(arkode_mem, pretype, maxl);``

   (b) ``flag = ARKSpbcg(arkode_mem, pretype, maxl);``

   (c) ``flag = ARKSptfqmr(arkode_mem, pretype, maxl);``

8. Initialize the ARKBBDPRE preconditioner module 

   Specify the upper and lower half-bandwidths for computation
   ``mudq`` and ``mldq``, the upper and lower half-bandwidths for
   storage ``mukeep`` and ``mlkeep``, and call 

   ``flag = ARKBBDPrecInit(arkode_mem, Nlocal, mudq, mldq, mukeep, mlkeep, dqrely, gloc, cfn);``

   to allocate memory and initialize the internal preconditioner
   data. The last two arguments of :c:func:`ARKBBDPrecInit()` are the
   two user-supplied functions of type :c:func:`ARKLocalFn()` and
   :c:func:`ARKCommFn()` described above, respectivelyl. 

9. `Set the linear solver optional inputs`

   Note that the user should not overwrite the preconditioner setup
   function or solve function through calls to ARKSPILS optional
   input functions. 

10. `Advance solution in time`

11. `Get optional outputs`

    Additional optional outputs associated with ARKBBDPRE are
    available by way of the two routines described below,
    :c:func:`ARKBBDPrecGetWorkSpace()` and
    :c:func:`ARKBBDPrecGetNumGfnEvals()`. 

12. `Deallocate memory for solution vector`

13. `Free solver memory`

14. `Finalize MPI`

The user-callable functions that initialize (step 8 above) or
re-initialize the ARKBBDPRE preconditioner module are described
next.



.. c:function:: int ARKBBDPrecInit(void *arkode_mem, long int Nlocal, long int mudq, long int mldq, long int mukeep, long int mlkeep, realtype dqrely, ARKLocalFn gloc, ARKCommFn cfn)

   Initializes and allocates (internal) memory for the
   ARKBBDPRE preconditioner.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `Nlocal` -- local vector length.
      * `mudq` -- upper half-bandwidth to be used in the difference
        quotient Jacobian approximation.
      * `mldq` -- lower half-bandwidth to be used in the difference
        quotient Jacobian approximation.
      * `mukeep` -- upper half-bandwidth of the retained banded
        approximate Jacobian block.
      * `mlkeep` -- lower half-bandwidth of the retained banded
        approximate Jacobian block.
      * `dqrely` -- the relative increment in components of `y` used in
        the difference quotient approximations.  The default is `dqrely`
        = :math:`\sqrt{\text{unit roundoff}}`, which can be specified by
        passing `dqrely` = 0.0.
      * `gloc` -- the name of the C function (of type :c:func:`ARKLocalFn()`)
        which computes the approximation :math:`g(t,y) \approx f_I(t,y)`.
      * `cfn` -- the name of the C function (of type :c:func:`ARKCommFn()`) which
        performs all interprocess communication required for the
        computation of :math:`g(t,y)`.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_ILL_INPUT if an input has an illegal value
      * ARKSPILS_MEM_FAIL if a memory allocation request failed
   
   **Notes:**  If one of the half-bandwidths `mudq` or `mldq` to be used
   in the difference quotient calculation of the approximate Jacobian is
   negative or exceeds the value `Nlocal-1`, it is replaced by 0 or
   `Nlocal-1` accordingly. 
   
   The half-bandwidths `mudq` and `mldq` need not be the true
   half-bandwidths of the Jacobian of the local block of :math:`g`
   when smaller values may provide a greater efficiency. 
   
   Also, the half-bandwidths `mukeep` and `mlkeep` of the retained
   banded approximate Jacobian block may be even smaller than
   `mudq` and `mldq`, to reduce storage and computational costs
   further. 
   
   For all four half-bandwidths, the values need not be the same on
   every processor.



The ARKBBDPRE module also provides a reinitialization function to
allow solving a sequence of problems of the same size, with the same
linear solver choice, provided there is no change in `Nlocal`,
`mukeep`, or `mlkeep`. After solving one problem, and after
calling :c:func:`ARKodeReInit()` to re-initialize ARKode for a
subsequent problem, a call to :c:func:`ARKBBDPrecReInit()` can be made
to change any of the following: the half-bandwidths `mudq` and
`mldq` used in the difference-quotient Jacobian approximations, the
relative increment `dqrely`, or one of the user-supplied functions
`gloc` and `cfn`. If there is a change in any of the linear solver
inputs, an additional call to :c:func:`ARKSpgmr()`,
:c:func:`ARKSpbcg()`, or :c:func:`ARKSptfqmr()`, and/or one or more of
the corresponding ARKSpils*Set* functions, must also be made (in
the proper order).



.. c:function:: int ARKBBDPrecReInit(void *arkode_mem, long int mudq, long int mldq, realtype dqrely)

   Re-initializes the ARKBBDPRE preconditioner module.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `mudq` -- upper half-bandwidth to be used in the difference
        quotient Jacobian approximation.
      * `mldq` -- lower half-bandwidth to be used in the difference
        quotient Jacobian approximation.
      * `dqrely` -- the relative increment in components of `y` used in
        the difference quotient approximations.  The default is `dqrely`
        = :math:`\sqrt{\text{unit roundoff}}`, which can be specified by
        passing `dqrely` = 0.0.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_PMEM_NULL if the preconditioner memory is ``NULL``
   
   **Notes:**  If one of the half-bandwidths `mudq` or `mldq` is
   negative or exceeds the value `Nlocal-1`, it is replaced by 0 or
   `Nlocal-1` accordingly. 


The following two optional output functions are available for use with
the ARKBBDPRE module:



.. c:function:: int ARKBBDPrecGetWorkSpace(void *arkode_mem, long int *lenrwBBDP, long int *leniwBBDP)

   Returns the processor-local ARKBBDPRE real and
   integer workspace sizes.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `lenrwBBDP` -- the number of ``realtype`` values in the
        ARKBBDPRE workspace.
      * `leniwBBDP` -- the number of integer values in the  ARKBBDPRE workspace.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_PMEM_NULL if the preconditioner memory is ``NULL``
   
   **Notes:**  In terms of `Nlocal` and `smu = min(Nlocal-1,
   mukeep+mlkeep)`, the actual size of the real workspace is `(2
   mlkeep + mukeep + smu + 2)*Nlocal`  ``realtype`` words, and the
   actual size of the integer workspace is `Nlocal` integer
   words. These values are local to each process. 
   
   The workspaces referred to here exist in addition to those given by
   the corresponding function ARKSpils*GetWorkSpace. 



.. c:function:: int ARKBBDPrecGetNumGfnEvals(void *arkode_mem, long int *ngevalsBBDP)

   Returns the number of calls made to the user-supplied
   `gloc` function (of type :c:func:`ARKLocalFn()`) due to the finite
   difference approximation of the Jacobian blocks used within the
   preconditioner setup function. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `ngevalsBBDP` -- the number of calls made to the user-supplied
        `gloc` function. 
   
   **Return value:**  
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_PMEM_NULL if the preconditioner memory is ``NULL``
   
   
In addition to the `ngevalsBBDP` `gloc` evaluations, the costs
associated with ARKBBDPRE also include `nlinsetups` LU
factorizations, `nlinsetups` calls to `cfn`, `npsolves` banded
backsolve calls, and `nfevalsLS` right-hand side function
evaluations, where `nlinsetups` is an optional ARKode output and
`npsolves` and `nfevalsLS` are linear solver optional outputs (see
the table :ref:`CInterface.ARKSpilsOutputTable`).
