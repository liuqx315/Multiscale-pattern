..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _CInterface.UserCallable:

User-callable functions
=============================

This section describes the ARKode functions that are called by the
user to setup and then solve an IVP. Some of these are
required. However, starting with the section
:ref:`CInterface.OptionalInputs`, the functions listed involve
optional inputs/outputs or restarting, and those paragraphs may be
skipped for a casual use of ARKode. In any 
case, refer to the section :ref:`CInterface.Skeleton` for the correct
order of these calls. 

On an error, each user-callable function returns a negative value and
sends an error message to the error handler routine, which prints the
message on ``stderr`` by default. However, the user can set a file as
error output or can provide his own error handler function
(see the section :ref:`CInterface.OptionalInputs` for details).



.. _CInterface.Initialization:

ARKode initialization and deallocation functions
------------------------------------------------------

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
------------------------------------------------------

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
      * ARK_MEM_NULL  if the ARKode memory was ``NULL``
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    control is meaningless.  For example, if :math:`y_i` starts at some
    nonzero value, but in time decays to zero, then pure relative
    error control on :math:`y_i` makes no sense (and is overly costly)
    after :math:`y_i` is below some noise level. Then ``abstol`` (if
    scalar) or ``abstol[i]`` (if a vector) needs to be set to that
    noise level. If the different components have different noise
    levels, then ``abstol`` should be a vector. See the example
    ``ark_robertson.c`` in the ARKode package, and the discussion
    of it in the ARKode Examples Documentation [R2013]_. In that
    problem, the three components vary betwen 0 and 1, and have
    different noise levels; hence the ``atols`` vector therein. It is
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
-------------------------------------------

As previously explained, the modified Newton iteration used in solving
implicit systems within ARKode requires the solution of linear
systems of the form 

.. math::
    A\left(y^n(m)\right) s^m = -F\left(y^n(m)\right)

where 

.. math::
    A \approx M - \gamma J, \qquad J = \frac{\partial f_I}{\partial y}.

There are seven ARKode linear solvers currently available for this
task: ARKDENSE, ARKBAND, ARKSPGMR, ARKSPBCG, ARKSPTFQMR, ARKSPFGMR and
ARKPCG. 

The first two linear solvers are direct and derive their names from
the type of approximation used for the Jacobian :math:`J`;
ARKDENSE and ARKBAND work with dense and banded approximations
to :math:`J`, respectively. The SUNDIALS suite includes both
internal implementations of these two linear solvers and interfaces to
LAPACK implementations. Together, these linear solvers are referred to
as ARKDLS (from Direct Linear Solvers). 

The last five ARKode linear solvers, ARKSPGMR, ARKSPBCG,
ARKSPTFQMR, ARKSPFGMR and ARKPCG, are Krylov iterative solvers, which
use scaled preconditioned GMRES, scaled preconditioned Bi-CGStab,
scaled preconditioned TFQMR, scaled preconditioned flexible GMRES, and
preconditioned conjugate gradient, respectively. Together, they are
referred to as ARKSPILS (from Scaled Preconditioned Iterative Linear
Solvers). 

With any of the Krylov methods, preconditioning can be done on the
left only, on the right only, on both the left and the right, or not
at all. For the specification of a preconditioner, see the iterative
linear solver sections in :ref:`CInterface.OptionalOutputs` and
:ref:`CInterface.UserSupplied`. 

If preconditioning is done, user-supplied functions define left and
right preconditioner matrices :math:`P_1` and :math:`P_2` (either of
which could be the identity matrix), such that the product :math:`P_{1}P_{2}`
approximates the Newton matrix  :math:`A = M - \gamma J`. 

To specify a ARKode linear solver, after the call to
:c:func:`ARKodeCreate()` but before any calls to :c:func:`ARKode()`, the user's
program must call one of the functions
:c:func:`ARKDense()`/:c:func:`ARKLapackDense()`,
:c:func:`ARKBand()`/:c:func:`ARKLapackBand()`, 
:c:func:`ARKSpgmr()`, :c:func:`ARKSpbcg()`, :c:func:`ARKSptfqmr()`,
:c:func:`ARKSpfgmr()` or :c:func:`ARKPcg()` as documented below. The
first argument passed to these functions is the ARKode memory pointer
returned by :c:func:`ARKodeCreate()`. A call to one of these functions
links the main ARKode integrator to a linear solver and allows the user to
specify parameters which are specific to a particular solver, such as
the half-bandwidths in the :c:func:`ARKBand()` case. The use of each
of the linear solvers involves certain constants and possibly some
macros, that are likely to be needed in the user code. These are
available in the corresponding header file associated with the linear
solver, as specified below.  

In each case except LAPACK direct solvers, the linear solver module
used by ARKode is actually built on top of a generic linear system
solver, which may be of interest in itself.  These generic solvers,
denoted DENSE, BAND, SPGMR, SPBCG, SPTFQMR, SPFGMR and PCG,
are described separately in the section :ref:`LinearSolvers`.



.. c:function:: int ARKDense(void *arkode_mem, long int N)

   A call to the ARKDense function links the main
   ARKode integrator with the ARKDENSE linear solver.  

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
   ARKode integrator with the ARKLAPACK linear solver module.
   
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
   ARKode integrator with the ARKBAND linear solver.
   
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
   ARKode integrator with the ARKLAPACK linear solver using banded Jacobians.
   
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



.. c:function:: int ARKSpfgmr(void *arkode_mem, int pretype, int maxl)

   A call to the ARKSpfgmr function links the main ARKode integrator
   with the ARKSPFGMR linear solver. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `pretype` -- the type of user preconditioning to be done.  This
        must be one of the four enumeration constants PREC_NONE,
        PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in
        ``sundials_iterative.h``. These correspond to no preconditioning,
        left preconditioning only, right preconditioning only, and both
        left and right preconditioning, respectively.
      * `maxl` -- the maximum Krylov dimension. This is an optional input
        to the ARKSPFGMR solver. Pass 0 to use the default value of 5.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKSPILS_MEM_FAIL  if there was a memory allocation failure
      * ARKSPILS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** The ARKSPFGMR solver uses a scaled preconditioned
   flexible GMRES iterative method to solve the linear systems.



.. c:function:: int ARKPcg(void *arkode_mem, int pretype, int maxl)

   A call to the ARKPcg function links the main
   ARKode integrator with the ARKPCG linear solver.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `pretype` -- flag denoting whether to use preconditioning.  If
        set to any of the enumeration constants PREC_LEFT, PREC_RIGHT,
	or PREC_BOTH, defined in ``sundials_iterative.h``,
	preconditioning will be enabled. Due to the symmetric form of
	PCG, there is no choice between	left and right
	preconditioning.
      * `maxl` -- the maximum Krylov dimension. This is an optional input
        to the ARKPCG solver. Pass 0 to use the default value of 5.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKSPILS_MEM_FAIL  if there was a memory allocation failure
      * ARKSPILS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** The ARKPCG solver uses a preconditioned conjugate
   gradient iterative method to solve the linear systems.
   






.. _CInterface.MassMatrixSolvers:

Mass matrix solver specification functions
-------------------------------------------

As discussed in section :ref:`Mathematics.MassSolve`, if the ODE
system involves a non-identity mass matrix :math:`M\ne I`, then ARKode
must solve linear systems of the form 

.. math::
    M x = b.

The same solvers listed above in section
:ref:`CInterface.LinearSolvers` may be used for this purpose:
ARKDENSE, ARKBAND, ARKSPGMR, ARKSPBCG, ARKSPTFQMR, ARKSPFGMR and
ARKPCG.  

With the four Krylov solvers (SPGMR, SPBCG, SPTFQMR, SPFGMR and PCG),
preconditioning can be applied. For the specification of a
preconditioner, see the iterative linear solver sections in
:ref:`CInterface.OptionalOutputs` and :ref:`CInterface.UserSupplied`.
Moreover, If preconditioning is done, user-supplied functions define
left and right preconditioner matrices :math:`P_1` and :math:`P_2`
(either of which could be the identity matrix), such that the product
:math:`P_{1}P_{2}` approximates the mass matrix  :math:`M`. 

To specify a mass matrix solver, after the call to
:c:func:`ARKodeCreate()` but before any calls to :c:func:`ARKode()`,
the user's program must call one of the functions
:c:func:`ARKMassDense()`/:c:func:`ARKMassLapackDense()`,
:c:func:`ARKMassBand()`/:c:func:`ARKMassLapackBand()`, 
:c:func:`ARKMassSpgmr()`, :c:func:`ARKMassSpbcg()`,
:c:func:`ARKMassSptfqmr()`, :c:func:`ARKMassSpfgmr()` or
:c:func:`ARKMassPcg()` as documented below. The first argument passed
to these functions is the ARKode memory pointer returned by
:c:func:`ARKodeCreate()`. A call to one of these functions 
links the mass matrix solve with the specified solver module, and
allows the user to specify parameters which are specific to a
particular solver. The use of each of the linear solvers involves
certain constants and possibly some macros, that are likely to be
needed in the user code. These are available in the corresponding
header file associated with the linear solver, as specified below.



.. c:function:: int ARKMassDense(void *arkode_mem, long int N, ARKDlsDenseMassFn dmass)

   A call to the ARKMassDense function links the mass matrix solve
   with the ARKDENSE linear solver module, and specifies the dense
   mass matrix function.

   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `N` -- the number of components in the ODE system.
      * `dmass` -- name of user-supplied dense mass matrix function.
   
   **Return value:** 
       * ARKDLS_SUCCESS   if successful
       * ARKDLS_MEM_NULL  if the ARKode memory was ``NULL``
       * ARKDLS_MEM_FAIL  if there was a memory allocation failure
       * ARKDLS_ILL_INPUT if a required vector operation is missing
   
   **Notes:**  The ARKDENSE linear solver may not be compatible with the
   particular implementation of the NVECTOR module. Of the two
   nvector modules provided with SUNDIALS, only NVECTOR_SERIAL is
   compatible. 



.. c:function:: int ARKMassLapackDense(void *arkode_mem, int N, ARKDlsDenseMassFn dmass)

   A call to the ARKMassLapackDense function links the mass matrix
   solve with the ARKLAPACK linear solver module.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `N` -- the number of components in the ODE system.
      * `dmass` -- name of user-supplied dense mass matrix function.
   
   **Return value:** 
      * ARKDLS_SUCCESS   if successful
      * ARKDLS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKDLS_MEM_FAIL  if there was a memory allocation failure
      * ARKDLS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** Here `N` is restricted to be of type ``int``, because of the
   corresponding type restriction in the LAPACK solvers.



.. c:function:: int ARKMassBand(void *arkode_mem, long int N, long int mupper, long int mlower, ARKDlsBandMassFn bmass)

   A call to the ARKMassBand function links the mass matrix solve with
   the ARKBAND linear solver module.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `N` -- the number of components in the ODE system.
      * `mupper` -- the upper bandwidth of the band mass matrix.
      * `mlower` -- is the lower bandwidth of the band mass matrix.
      * `bmass` -- name of user-supplied band mass matrix function.
   
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
   mass matrix satisfy `-mlower` :math:`\le` `j-i` :math:`\le`
   `mupper`. 

   At present, if using a band solver for both the Jacobian systems
   and the mass matrix systems, it is required that the band mass
   matrix have identical band structure to the Jacobian matrix.  While
   this is typical of finite-element problems, if this is not true for
   a specific problem it can be handled by manually zero-padding the
   mass matrix.



.. c:function:: int ARKMassLapackBand(void *arkode_mem, int N, int mupper, int mlower, ARKDlsBandMassFn bmass)

   A call to the ARKMassLapackBand function links the mass matrix
   solve with the ARKLAPACK linear solver module.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `N` -- the number of components in the ODE system.
      * `mupper` -- the upper bandwidth of the band mass matrix.
      * `mlower` -- is the lower bandwidth of the band mass matrix.
      * `bmass` -- name of user-supplied band mass matrix function.
   
   **Return value:** 
      * ARKDLS_SUCCESS   if successful
      * ARKDLS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKDLS_MEM_FAIL  if there was a memory allocation failure
      * ARKDLS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** Here, each of `N`, `mupper` and `mlower` are restricted
   to be of type ``int``, because of the corresponding type restriction
   in the LAPACK solvers.

   At present, if using a band solver for both the Jacobian systems
   and the mass matrix systems, it is required that the band mass
   matrix have identical band structure to the Jacobian matrix.  While
   this is typical of finite-element problems, if this is not true for
   a specific problem it can be handled by manually zero-padding the
   mass matrix.



.. c:function:: int ARKMassSpgmr(void *arkode_mem, int pretype, int maxl, ARKSpilsMassTimesVecFn mtimes)

   A call to the ARKMassSpgmr function links the mass matrix solve
   with the ARKSPGMR linear solver module.
   
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
      * `mtimes` -- user-defined mass-matrix-vector product function.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKSPILS_MEM_FAIL  if there was a memory allocation failure
      * ARKSPILS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** The ARKSPGMR solver uses a scaled preconditioned GMRES
   iterative method to solve the linear systems.



.. c:function:: int ARKMassSpbcg(void *arkode_mem, int pretype, int maxl, ARKSpilsMassTimesVecFn mtimes)

   A call to the ARKMassSpbcg function links the mass matrix solve
   with the ARKSPBCG linear solver module.
   
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
      * `mtimes` -- user-defined mass-matrix-vector product function.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKSPILS_MEM_FAIL  if there was a memory allocation failure
      * ARKSPILS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** The ARKSPBCG solver uses a scaled preconditioned Bi-CGStab 
   iterative method to solve the linear systems.
   


.. c:function:: int ARKMassSptfqmr(void *arkode_mem, int pretype, int maxl, ARKSpilsMassTimesVecFn mtimes)

   A call to the ARKMassSptfqmr function links the mass matrix solve with
   the ARKSPTFQMR linear solver.
   
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
      * `mtimes` -- user-defined mass-matrix-vector product function.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKSPILS_MEM_FAIL  if there was a memory allocation failure
      * ARKSPILS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** The ARKSPTFQMR solver uses a scaled preconditioned TFQMR
   iterative method to solve the linear systems.



.. c:function:: int ARKMassSpfgmr(void *arkode_mem, int pretype, int maxl, ARKSpilsMassTimesVecFn mtimes)

   A call to the ARKMassSpfgmr function links the mass matrix solve with
   the ARKSPFGMR linear solver. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `pretype` -- the type of user preconditioning to be done.  This
        must be one of the four enumeration constants PREC_NONE,
        PREC_LEFT, PREC_RIGHT, or PREC_BOTH defined in
        ``sundials_iterative.h``. These correspond to no preconditioning,
        left preconditioning only, right preconditioning only, and both
        left and right preconditioning, respectively.
      * `maxl` -- the maximum Krylov dimension. This is an optional input
        to the ARKSPFGMR solver. Pass 0 to use the default value of 5.
      * `mtimes` -- user-defined mass-matrix-vector product function.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKSPILS_MEM_FAIL  if there was a memory allocation failure
      * ARKSPILS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** The ARKSPFGMR solver uses a scaled preconditioned
   flexible GMRES iterative method to solve the linear systems.



.. c:function:: int ARKMassPcg(void *arkode_mem, int pretype, int maxl, ARKSpilsMassTimesVecFn mtimes)

   A call to the ARKMassPcg function links the mass matrix solve with
   the ARKPCG linear solver.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `pretype` -- flag denoting whether to use preconditioning.  If
        set to any of the enumeration constants PREC_LEFT, PREC_RIGHT,
	or PREC_BOTH, defined in ``sundials_iterative.h``,
	preconditioning will be enabled. Due to the symmetric form of
	PCG, there is no choice between	left and right
	preconditioning.
      * `maxl` -- the maximum Krylov dimension. This is an optional input
        to the ARKPCG solver. Pass 0 to use the default value of 5.
      * `mtimes` -- user-defined mass-matrix-vector product function.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKSPILS_MEM_FAIL  if there was a memory allocation failure
      * ARKSPILS_ILL_INPUT if a required vector operation is missing
   
   **Notes:** The ARKPCG solver uses a preconditioned conjugate
   gradient iterative method to solve the linear systems.
   




.. _CInterface.RootFinding:

Rootfinding initialization function
--------------------------------------

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
-------------------------

This is the central step in the solution process -- the call to perform
the integration of the IVP.  One of the input arguments (`itask`)
specifies one of two modes as to where ARKode is to return a
solution.  These modes are modified if the user has set a stop time
(with a call to the optional input function :c:func:`ARKodeSetStopTime()`) or
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
        approximate value of :math:`y(tout)`. This interpolation is
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
      * ARK_MASSINIT_FAIL if the mass matrix solver's
	initialization function failed. 
      * ARK_MASSSETUP_FAIL if the mass matrix solver's setup routine
	failed.
      * ARK_MASSSOLVE_FAIL if the mass matrix solver's solve routine
	failed.
   
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
-------------------------

There are numerous optional input parameters that control the behavior
of the ARKode solver. ARKode provides functions that can be
used to change these optional input parameters from their default
values. The following tables list all optional input functions in
ARKode which are then described in detail in the remainder of this
section.  These optional inputs are grouped into the following
categories:

* General solver options (:ref:`CInterface.ARKodeInputTable`), 
* IVP method solver options (:ref:`CInterface.ARKodeMethodInputTable`), 
* Step adaptivity solver options (:ref:`CInterface.ARKodeAdaptivityInputTable`), 
* Implicit stage solver options (:ref:`CInterface.CInterface.ARKodeSolverInputTable`), 
* Dense linear solver options (:ref:`CInterface.ARKDlsInputs`) 
* Iterative linear solver options (:ref:`CInterface.ARKSpilsInputs`).  

For the most casual use of ARKode, the reader can skip to the section
:ref:`CInterface.UserSupplied`.

We note that, on an error return, all of the optional input functions
send an error message to the error handler function.  We also note
that all error return values are negative, so a test on the return
arguments for negative values will catch all errors. 



.. _CInterface.ARKodeInputTable:

Optional inputs for ARKode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

===============================================  ====================================  ==============
Optional input                                   Function name                         Default
===============================================  ====================================  ==============
Set dense output order                           :c:func:`ARKodeSetDenseOrder()`       3
Set default solver parameters                    :c:func:`ARKodeSetDefaults()`         internal
Pointer to a diagnostics file                    :c:func:`ARKodeSetDiagnostics()`      ``NULL``
Pointer to an error file                         :c:func:`ARKodeSetErrFile()`          ``stderr``
Error handler function                           :c:func:`ARKodeSetErrHandlerFn()`     internal fn
Initial step size                                :c:func:`ARKodeSetInitStep()`         estimated
Maximum no. of warnings for :math:`t_n+h = t_n`  :c:func:`ARKodeSetMaxHnilWarns()`     10
Maximum no. of internal steps before `tout`      :c:func:`ARKodeSetMaxNumSteps()`      500
Maximum no. of error test failures               :c:func:`ARKodeSetMaxErrTestFails()`  7
Maximum absolute step size                       :c:func:`ARKodeSetMaxStep()`          :math:`\infty`
Minimum absolute step size                       :c:func:`ARKodeSetMinStep()`          0.0
Set 'optimal' adaptivity params                  :c:func:`ARKodeSetOptimalParams()`    internal
Value of :math:`t_{stop}`                        :c:func:`ARKodeSetStopTime()`         :math:`\infty`
User data                                        :c:func:`ARKodeSetUserData()`         ``NULL``
===============================================  ====================================  ==============



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
   
   **Notes:** The default value is 10; set *mxhnil* to zero to specify
   this default.

   A negative value indicates that no warning messages should be issued.




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
   
   **Notes:** The default value is 7; set *maxnef* :math:`\le 0`
   to specify this default.



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



.. c:function:: int ARKodeSetMinStep(void *arkode_mem, realtype hmin)

   Specifies the lower bound on the magnitude of the time step size.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `hmin` -- minimum absolute value of the time step size :math:`(\ge 0)`
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Pass `hmin \le 0.0` to set the default value of 0.



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







.. _CInterface.ARKodeMethodInputTable:

Optional inputs for IVP method selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

=================================  ================================  ==============
Optional input                     Function name                     Default
=================================  ================================  ==============
Set method order                   :c:func:`ARKodeSetOrder()`        4
Specify implicit/explicit problem  :c:func:`ARKodeSetImEx()`         ``TRUE``
Specify explicit problem           :c:func:`ARKodeSetExplicit()`     ``FALSE``
Specify implicit problem           :c:func:`ARKodeSetImplicit()`     ``FALSE``
Set additive RK tables             :c:func:`ARKodeSetARKTables()`    internal
Set explicit RK table              :c:func:`ARKodeSetERKTable()`     internal
Set implicit RK table              :c:func:`ARKodeSetIRKTable()`     internal
Specify additive RK tables number  :c:func:`ARKodeSetARKTableNum()`  internal
Specify explicit RK table number   :c:func:`ARKodeSetERKTableNum()`  internal
Specify implicit RK table number   :c:func:`ARKodeSetIRKTableNum()`  internal
=================================  ================================  ==============



.. c:function:: int ARKodeSetOrder(void *arkode_mem, int ord)

   Specifies the order of accuracy for the RK method.
   
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
   
   **Notes:** Both `itable` and `etable` should match an existing
   implicit/explicit pair from the section :ref:`Mathematics.Butcher.additive`.   
   Error-checking is performed to ensure that the tables exist.
   Subsequent error-checking is automatically performed to ensure that
   the tables' stage times and solution coefficients match.  

   This automatically calls :c:func:`ARKodeSetImEx()`. 



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
   
   **Notes:** `etable` should match an existing explicit method from
   the section :ref:`Mathematics.Butcher.explicit`.
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
   
   **Notes:** `itable` should match an existing implicit method from
   the section :ref:`Mathematics.Butcher.implicit`.
   Error-checking is performed to ensure that the table exists, and is
   not explicit.  
   
   This automatically calls :c:func:`ARKodeSetImplicit()`. 







.. _CInterface.ARKodeAdaptivityInputTable:

Optional inputs for time step adaptivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

==============================================  =====================================  ========
Optional input                                  Function name                          Default
==============================================  =====================================  ========
Time step adaptivity function                   :c:func:`ARKodeSetAdaptivityFn()`      internal
Time step adaptivity method                     :c:func:`ARKodeSetAdaptivityMethod()`  0
Explicit stability safety factor                :c:func:`ARKodeSetCFLFraction()`       0.5
Time step error bias factor                     :c:func:`ARKodeSetErrorBias()`         1.5
Bounds determining no change in step size       :c:func:`ARKodeSetFixedStepBounds()`   1.0  1.5
Maximum step growth factor on error test fail   :c:func:`ARKodeSetMaxEFailGrowth()`    0.3
Maximum first step growth factor                :c:func:`ARKodeSetMaxFirstGrowth()`    10000.0
Maximum step growth factor                      :c:func:`ARKodeSetMaxGrowth()`         20.0
Maximum step growth factor on convergence fail  :c:func:`ARKodeSetMaxCFailGrowth()`    0.25
Time step safety factor                         :c:func:`ARKodeSetSafetyFactor()`      0.96
Error fails before MaxEFailGrowth takes effect  :c:func:`ARKodeSetSmallNumEFails()`    2
Explicit stability function                     :c:func:`ARKodeSetStabilityFn()`       internal
==============================================  =====================================  ========



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


      
.. c:function:: int ARKodeSetAdaptivityMethod(void *arkode_mem, int imethod, int idefault, int pq, realtype *adapt_params)

   Specifies the method (and associated parameters) used for time step adaptivity.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `imethod` -- accuracy-based adaptivity method choice 
        (0 :math:`\le` `imethod` :math:`\le` 5): 
        0 is PID, 1 is PI, 2 is I, 3 is explicit Gustafsson, 4 is
        implicit Gustafsson, and 5 is the ImEx Gustafsson.
      * `idefault` -- flag denoting whether to use default adaptivity
	parameters (1), or that they will be supplied in the
	*adapt_params* argument (0).
      * `pq` -- flag denoting whether to use the embedding order of
	accuracy `p` (0) or the method order of accuracy `q` (1)
	within the adaptivity algorithm.
      * `adapt_params[0]` -- :math:`k_1` parameter within accuracy-based adaptivity algorithms.
      * `adapt_params[1]` -- :math:`k_2` parameter within accuracy-based adaptivity algorithms.
      * `adapt_params[2]` -- :math:`k_3` parameter within accuracy-based adaptivity algorithms.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** If custom parameters are supplied, they will be checked
   for validity against published stability intervals.  If other
   parameter values are desired, it is recommended to use the
   following function, :c:func:`ARKodeSetAdaptivityFn()`.


      
.. c:function:: int ARKodeSetCFLFraction(void *arkode_mem, realtype cfl_frac)

   Specifies the fraction of the estimated explicitly stable
   step to use.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `cfl_frac` -- maximum allowed fraction of explicitly stable step (default is 0.5)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default
   value.  
   

      
.. c:function:: int ARKodeSetErrorBias(void *arkode_mem, realtype bias)

   Specifies the bias to be applied to the error estimates within
   accuracy-based adaptivity strategies.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `bias` -- bias applied to error in accuracy-based time
        step estimation (default is 1.5)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Any value below 1.0 will imply a reset to the default
   value.  


      
.. c:function:: int ARKodeSetFixedStepBounds(void *arkode_mem, realtype lb, realtype ub)

   Specifies the step growth interval in which the step size will
   remain unchanged.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `lb` -- lower bound on window to leave step size fixed (default is 1.0)
      * `ub` -- upper bound on window to leave step size fixed (default is 1.5)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Any interval *not* containing 1.0 will imply a reset to the default value.  
   

      
.. c:function:: int ARKodeSetMaxEFailGrowth(void *arkode_mem, realtype etamxf)

   Specifies the maximum step size growth factor upon multiple successive
   accuracy-based error failures in the solver.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `etamxf` -- time step reduction factor on multiple error fails (default is 0.3)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Any value outside the interval (0,1] will imply a reset to the default value.
   

      
.. c:function:: int ARKodeSetMaxFirstGrowth(void *arkode_mem, realtype etamx1)

   Specifies the maximum allowed step size change following the very
   first integration step.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `etamx1` -- maximum allowed growth factor after the first time
        step (default is 10000.0)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Any value :math:`\le 1.0` will imply a reset to the default value.



.. c:function:: int ARKodeSetMaxGrowth(void *arkode_mem, realtype mx_growth)

   Specifies the maximum growth of the step size between consecutive
   steps in the integration process.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `growth` -- maximum allowed growth factor between consecutive time steps (default is 20.0)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Any value :math:`\le 1.0` will imply a reset to the default
   value.  



.. c:function:: int ARKodeSetMaxCFailGrowth(void *arkode_mem, realtype etacf)

   Specifies the maximum step size growth factor upon a convergence
   failure on a stage solve within a step.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `etacf` -- time step reduction factor on a nonlinear solver
        convergence failure (default is 0.25)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Any value outside the interval (0,1] will imply a reset to the default value.



.. c:function:: int ARKodeSetSafetyFactor(void *arkode_mem, realtype safety)

   Specifies the safety factor to be applied to the accuracy-based
   estimated step.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `safety` -- safety factor applied to accuracy-based time step (default is 0.96)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default
   value.  


      
.. c:function:: int ARKodeSetSmallNumEFails(void *arkode_mem, int small_nef)

   Specifies the threshold for "multiple" successive error failures
   before the `etamxf` parameter from
   :c:func:`ARKodeSetMaxEFailGrowth()` is applied.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `small_nef` -- bound to determine `multiple` for `etamxf` (default is 2)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default value.



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
   
   **Notes:** This function should return an estimate of the absolute
   value of the maximum stable time step for the explicit portion of
   the IMEX system.  It is not required, since accuracy-based
   adaptivity may be sufficient at retaining stability, but this can
   be quite useful for problems where the IMEX splitting may retain
   stiff components in :math:`f_E(t,y)`. 








.. _CInterface.CInterface.ARKodeSolverInputTable:

Optional inputs for implicit stage solves
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered

=============================================  ========================================  =========
Optional input                                 Function name                             Default
=============================================  ========================================  =========
Specify use of the fixed-point stage solver    :c:func:`ARKodeSetFixedPoint()`           ``FALSE``
Specify use of the Newton stage solver         :c:func:`ARKodeSetNewton()`               ``TRUE``
Specify linearly implicit :math:`f_I`          :c:func:`ARKodeSetLinear()`               ``FALSE``
Specify nonlinearly implicit :math:`f_I`       :c:func:`ARKodeSetNonlinear()`            ``TRUE``
Implicit predictor method                      :c:func:`ARKodeSetPredictorMethod()`      3
Maximum no. of nonlinear iterations            :c:func:`ARKodeSetMaxNonlinIters()`       3
Coefficient in the nonlinear convergence test  :c:func:`ARKodeSetNonlinConvCoef()`       0.2
Nonlinear convergence rate constant            :c:func:`ARKodeSetNonlinCRDown()`         0.3
Nonlinear residual divergence ratio            :c:func:`ARKodeSetNonlinRDiv()`           2.3
Max change in step signaling new :math:`J`     :c:func:`ARKodeSetDeltaGammaMax()`        0.2
Max steps between calls to new :math:`J`       :c:func:`ARKodeSetMaxStepsBetweenLSet()`  20
Maximum no. of convergence failures            :c:func:`ARKodeSetMaxConvFails()`         10
=============================================  ========================================  =========




.. c:function:: int ARKodeSetFixedPoint(void *arkode_mem, long int fp_m)

   Specifies that the implicit portion of the problem should be solved
   using the accelerated fixed-point solver instead of the modified
   Newton iteration, and provides the maximum dimension of the
   acceleration subspace.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `fp_m` -- number of vectors to store within the Anderson
        acceleration subspace.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Since the accelerated fixed-point solver has a slower
   rate of convergence than the Newton iteration (but each iteration
   is typically much more efficient), it is recommended that the
   maximum nonlinear correction iterations be increased through a call
   to :c:func:`ARKodeSetMaxNonlinIters()`. 



.. c:function:: int ARKodeSetNewton(void *arkode_mem)

   Specifies that the implicit portion of the problem should be solved
   using the modified Newton solver.  
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This is the default behavior of ARKode, so the function
   ARKodeSetNewton is primarily useful to undo a previous call
   to :c:func:`ARKodeSetFixedPoint()`. 



.. c:function:: int ARKodeSetLinear(void *arkode_mem)

   Specifies that the implicit portion of the problem is linear.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Tightens the linear solver tolerances and takes only a single
   Newton iteration.  Only useful when used in combination with the
   modified Newton iteration (not the fixed-point solver).



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



.. c:function:: int ARKodeSetPredictorMethod(void *arkode_mem, int method)

   Specifies the method to use for predicting implicit solutions.  
   Non-default choices are {1,2,3,4}, all others will use default 
   (trivial) predictor.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `method` -- method choice (0 :math:`\le` `method` :math:`\le` 4): 

        * 0 is the trivial predictor, 

        * 1 is the dense output predictor, 

	* 2 is the dense output predictor that decreases the
	  polynomial degree for more distant RK stages, 

        * 3 is the dense output predictor to max order for early RK
	  stages, and a first-order predictor for distant RK stages, 

        * 4 is the bootstrap predictor.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** This function is designed only for advanced ARKode usage.



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
   
   **Notes:** The default value is 3; set *maxcor* :math:`\le 0`
   to specify this default.



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
   
   **Notes:** The default value is 0.1; set *nlscoef* :math:`\le 0`
   to specify this default.



.. c:function:: int ARKodeSetNonlinCRDown(void *arkode_mem, realtype crdown)

   Specifies the constant used in estimating the nonlinear convergence rate.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `crdown` -- nonlinear convergence rate estimation constant (default is 0.3)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default value.



.. c:function:: int ARKodeSetNonlinRDiv(void *arkode_mem, realtype rdiv)

   Specifies the nonlinear correction threshold beyond which the
   iteration will be declared divergent.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `rdiv` -- tolerance on nonlinear correction size ratio to
	declare divergence (default is 2.3) 
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:** Any non-positive parameter will imply a reset to the default value.



.. c:function:: int ARKodeSetDeltaGammaMax(void *arkode_mem, realtype dgmax)

   Specifies a scaled step size ratio tolerance beyond which the
   linear solver setup routine will be signaled.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `dgmax` -- tolerance on step size ratio change before calling
        linear solver setup routine (default is 0.2)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:**  Any non-positive parameter will imply a reset to the default value.



.. c:function:: int ARKodeSetMaxStepsBetweenLSet(void *arkode_mem, int msbp)

   Specifies the maximum number of steps allowed between calls to the
   linear solver setup routine.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `msbp` -- maximum no. of time steps between linear solver setup calls (default is 20)
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory is ``NULL``
      * ARK_ILL_INPUT if an argument has an illegal value
   
   **Notes:**  Any non-positive parameter will imply a reset to the default value.



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
   
   **Notes:** The default value is 10; set *maxncf* :math:`\le 0`
   to specify this default.  Upon each convergence failure,
   ARKode will first call the Jacobian setup routine and try again;
   if a convergence failure still occurs, the time step size is reduced
   by the factor `etacf` (set within
   :c:func:`ARKodeSetMaxCFailGrowth()`).





.. _CInterface.ARKDlsInputs:


Direct linear solvers optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Table: Optional inputs for ARKDLS
"""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

==========================  ================================  =============
Optional input              Function name                     Default
==========================  ================================  =============
Dense Jacobian function     :c:func:`ARKDlsSetDenseJacFn()`   ``DQ``
Band Jacobian function      :c:func:`ARKDlsSetBandJacFn()`    ``DQ``
Dense mass matrix function  :c:func:`ARKDlsSetDenseMassFn()`  none
Band mass matrix function   :c:func:`ARKDlsSetBandMassFn()`   none
==========================  ================================  =============

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

If the ODE system involves a non-identity mass matrix, :math:`M\ne I`,
the ARKDENSE solver needs a function to compute a dense approximation
to the mass matrix :math:`M(t)`. This function must be of type
:c:func:`ARKDlsDenseMassFn()`. The user must supply his/her own dense
Jacobian function since there is no default value.  To specify this
user-supplied Jacobian function `dmass`, ARKDENSE provides the
function :c:func:`ARKDlsSetDenseMassFn()`. Alternatively, since the mass
matrix systems must also be solved separately from the implicit stage
systems, if a user chooses the ARKDENSE solver for these mass
matrix-specific solves, then the mass matrix function will already be
supplied to :c:func:`ARKMassDense()` or
:c:func:`ARKMassLapackDense()`, and does not need to be specified
again.  We note that the ARKDENSE solver passes the pointer user data
to the dense mass matrix function. This allows the user to create an
arbitrary structure with relevant problem data and access it during
the execution of the user-supplied mass matrix function, without using 
global data in the program. The pointer user data may be specified
through :c:func:`ARKodeSetUserData()`.



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



.. c:function:: int ARKDlsSetDenseMassFn(void *arkode_mem, ARKDlsDenseMassFn dmass)

   Specifies the dense mass matrix approximation routine to
   be used for a direct dense linear solver. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `dmass` -- name of user-supplied dense mass matrix approximation function.
   
   **Return value:** 
      * ARKDLS_SUCCESS  if successful
      * ARKDLS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKDLS_MASSMEM_NULL if the mass matrix solver memory was ``NULL``
   
   **Notes:** This routine must be called after the mass matrix solver
   has been initialized through a call to one of
   :c:func:`ARKMassDense()`, :c:func:`ARKMassLapackDense()`, 
   :c:func:`ARKMassBand()`, :c:func:`ARKMassLapackBand()`,
   :c:func:`ARKMassSpgmr()`, :c:func:`ARKMassSpbcg()`,
   :c:func:`ARKMassSptfqmr()`, :c:func:`ARKMassSpfgmr()` or
   :c:func:`ARKMassPcg()`. 
   
   The function type :c:func:`ARKDlsDenseMassFn()` is described in the section
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

Moreover, if the ODE system involves a non-identity mass matrix,
:math:`M\ne I`, the ARKBAND solver needs a function to compute a banded
approximation to the mass matrix :math:`M(t)`. This function
must be of type :c:func:`ARKDlsBandMassFn()`. The user must supply his/her own
banded mass matrix approximation function, since there is no default.
To specify this user-supplied mass matrix function `bmass`,
ARKBAND provides the function :c:func:`ARKDlsSetBandMassFn()`. 
Alternatively, since the mass matrix systems must also be solved
separately from the implicit stage systems, if a user chooses the
ARKBAND solver for these mass matrix-specific solves, then the mass
matrix function will already be supplied to :c:func:`ARKMassBand()` or
:c:func:`ARKMassLapackBand()`, and does not need to be specified
again.  We note that the ARKBAND solver passes the pointer user data
to the banded mass matrix approximation function.  This allows the
user to create an arbitrary structure with relevant problem data and
access it during the execution of the user-supplied mass matrix
function, without using global data in the program. The pointer user
data may be specified through :c:func:`ARKodeSetUserData()`. 



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


.. c:function:: int ARKDlsSetBandMassFn(void *arkode_mem, ARKDlsBandMassFn bmass)

   Specifies the band mass matrix approximation routine to be
   used for a direct band linear solver. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `bmass` -- name of user-supplied banded mass matrix approximation function.
   
   **Return value:** 
      * ARKDLS_SUCCESS  if successful
      * ARKDLS_MEM_NULL  if the ARKode memory was ``NULL``
      * ARKDLS_MASSMEM_NULL if the mass matrix solver memory was ``NULL``
   
   **Notes:** This routine must be called after the mass matrix solver
   has been initialized through a call to one of
   :c:func:`ARKMassDense()`, :c:func:`ARKMassLapackDense()`, 
   :c:func:`ARKMassBand()`, :c:func:`ARKMassLapackBand()`,
   :c:func:`ARKMassSpgmr()`, :c:func:`ARKMassSpbcg()`,
   :c:func:`ARKMassSptfqmr()`, :c:func:`ARKMassSpfgmr()` or
   :c:func:`ARKMassPcg()`. 
   
   The function type :c:func:`ARKDlsBandMassFn()` is described in the section
   :ref:`CInterface.UserSupplied`.



.. _CInterface.ARKSpilsInputs:

Iterative linear solvers optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Similarly, if a problem involves a non-identity mass matrix,
:math:`M\ne I`, then the ARKSPILS solvers require a function to
compute an approximation to the product between the mass matrix
:math:`M(t)` and a vector :math:`v`.  This function must be
user-supplied, and must be of type :c:func:`ARKSpilsMassTimesVecFn()`,
and can be specified through a call to 
:c:func:`ARKSpilsSetMassTimesVecFn()`.  Alternatively, since the mass
matrix systems must also be solved separately from the implicit stage
systems, if a user chooses an ARKSPILS solver for these mass matrix
systems then the mass matrix vector product function must be supplied
when setting up that linear solver module in any of the functions
:c:func:`ARKMassSpgmr()`, :c:func:`ARKMassSpbcg()`,
:c:func:`ARKMassSptfqmr()`, :c:func:`ARKMassPcg()` or
:c:func:`ARKMassSpfgmr()`. 


Table: Optional inputs for ARKSPILS
"""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

==============================================  =========================================  ==================
Optional input                                  Function name                              Default
==============================================  =========================================  ==================
Jacobian-times-vector function                  :c:func:`ARKSpilsSetJacTimesVecFn()`       ``DQ``
Mass matrix-times-vector function               :c:func:`ARKSpilsSetMassTimesVecFn()`      none
Maximum Krylov subspace size `(b)`              :c:func:`ARKSpilsSetMaxl()`                5
Max mass matrix Krylov subspace size `(b)`      :c:func:`ARKSpilsSetMassMaxl()`            5
Type of Gram-Schmidt orthogonalization `(a)`    :c:func:`ARKSpilsSetGSType()`              classical GS
Type of mass matrix Gram-Schmidt orthog. `(a)`  :c:func:`ARKSpilsSetMassGSType()`          classical GS
Preconditioning type                            :c:func:`ARKSpilsSetPrecType()`            none
Mass matrix preconditioning type                :c:func:`ARKSpilsSetMassPrecType()`        none
Preconditioner functions                        :c:func:`ARKSpilsSetPreconditioner()`      ``NULL``, ``NULL``
Mass matrix preconditioner functions            :c:func:`ARKSpilsSetMassPreconditioner()`  ``NULL``, ``NULL``
Ratio between linear and nonlinear tolerances   :c:func:`ARKSpilsSetEpsLin()`              0.05
Ratio between mass matrix and nonlinear tols    :c:func:`ARKSpilsSetMassEpsLin()`          0.05
==============================================  =========================================  ==================


`(a)` Only for ARKSPGMR and ARKSPFGMR

`(b)` Only for ARKSPBCG, ARMSPTFQMR and ARKPCG



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



.. c:function:: int ARKSpilsSetMassTimesVecFn(void *arkode_mem, ARKSpilsMassTimesVecFn mtimes)

   Specifies the mass matrix-times-vector function. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `mtimes` -- user-defined mass matrix-vector product function.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful.
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``.
      * ARKSPILS_MASSMEM_NULL if the mass matrix solver memory was ``NULL``.
      * ARKSPILS_ILL_INPUT if an input has an illegal value.

   **Notes:** This function must be called *after* the mass matrix
   solver has been initialized, through a call to one of
   :c:func:`ARKMassDense()`, :c:func:`ARKMassLapackDense()`,
   :c:func:`ARKMassBand()`, :c:func:`ARKMassLapackBand()`,
   :c:func:`ARKMassSpgmr()`, :c:func:`ARKMassSpbcg()`,
   :c:func:`ARKMassSptfqmr()`, :c:func:`ARKMassSpfgmr()` or
   :c:func:`ARKMassPcg()`. 
   
   The function type :c:func:`ARKSpilsMassTimesVecFn()` is described
   in the section :ref:`CInterface.UserSupplied`.



.. c:function:: int ARKSpilsSetMaxl(void *arkode_mem, int maxl)

   Resets the maximum Krylov subspace size, `maxl`, from
   the value previously set, when using the Bi-CGStab, TFQMR or PCG linear
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
   
   This option is available only for the ARKSPBCG, ARKSPTFQMR and
   ARKPCG linear solvers. 



.. c:function:: int ARKSpilsSetMassMaxl(void *arkode_mem, int maxl)

   Resets the maximum mass matrix Krylov subspace size, `maxl`, from
   the value previously set, when using the Bi-CGStab, TFQMR or PCG linear
   solver methods.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `maxl` -- maximum dimension of the mass matrix Krylov subspace.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful.
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``.
      * ARKSPILS_MASSMEM_NULL if the mass matrix solver memory was ``NULL``.
      * ARKSPILS_ILL_INPUT if an input has an illegal value.
   
   **Notes:** The maximum subspace dimension is initially specified in the
   call to the linear solver specification function (see the section
   :ref:`CInterface.LinearSolvers`).  This function call is needed
   only if `maxl` is being changed from its previous value.
  
   An input value `maxl` :math:`\le 0`, gives the default value, 5.
   
   This option is available only for the ARKSPBCG, ARKSPTFQMR and
   ARKPCG linear solvers. 



.. c:function:: int ARKSpilsSetGSType(void *arkode_mem, int gstype)

   Specifies the type of Gram-Schmidt orthogonalization to
   be used with the ARKSPGMR or ARKSPFGMR linear solvers. This must be
   one of the two enumeration constants MODIFIED_GS or CLASSICAL_GS
   defined in ``sundials_iterative.h``. These correspond to using modified
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
   
   This option is available only for the ARKSPGMR and ARKSPFGMR linear
   solvers. 



.. c:function:: int ARKSpilsSetMassGSType(void *arkode_mem, int gstype)

   Specifies the type of Gram-Schmidt orthogonalization to
   be used with the ARKSPGMR or ARKSPFGMR linear mass matrix
   solvers. This must be one of the two enumeration constants
   MODIFIED_GS or CLASSICAL_GS defined in ``sundials_iterative.h``. These
   correspond to using modified Gram-Schmidt and classical
   Gram-Schmidt, respectively. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `gstype` -- type of Gram-Schmidt orthogonalization.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful.
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``.
      * ARKSPILS_MASSMEM_NULL if the mass matrix solver memory was ``NULL``.
      * ARKSPILS_ILL_INPUT if an input has an illegal value.
   
   **Notes:** The default value is MODIFIED_GS.
   
   This option is available only for the ARKSPGMR and ARKSPFGMR linear
   solvers. 



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



.. c:function:: int ARKSpilsSetMassPrecType(void *arkode_mem, int pretype)

   Resets the type of mass matrix preconditioner, `pretype`, from the value previously set.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `pretype` -- the type of preconditioning to use, must be one of
        PREC_NONE, PREC_LEFT, PREC_RIGHT or PREC_BOTH. 
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful.
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``.
      * ARKSPILS_MASSMEM_NULL if the mass matrix solver memory was ``NULL``.
      * ARKSPILS_ILL_INPUT if an input has an illegal value.
   
   **Notes:** The preconditioning type is initially set in the call to
   the mass matrix solver's specification function (see the section
   :ref:`CInterface.LinearSolvers`).  This function call is needed
   only if `pretype` is being changed from its original value.



.. c:function:: int ARKSpilsSetPreconditioner(void *arkode_mem, ARKSpilsPrecSetupFn psetup, ARKSpilsPrecSolveFn psolve)

   Specifies the preconditioner setup and solve functions.  
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `psetup` -- user defined preconditioner setup function.  Pass
        ``NULL`` if no setup is to be done.
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



.. c:function:: int ARKSpilsSetMassPreconditioner(void *arkode_mem, ARKSpilsMassPrecSetupFn psetup, ARKSpilsMassPrecSolveFn psolve)

   Specifies the mass matrix preconditioner setup and solve functions.  
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `psetup` -- user defined preconditioner setup function.  Pass
        ``NULL`` if no setup is to be done.
      * `psolve` -- user-defined preconditioner solve function.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful.
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``.
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``.
      * ARKSPILS_ILL_INPUT if an input has an illegal value.
   
   **Notes:** The default is ``NULL`` for both arguments (i.e. no
   preconditioning).
    
   Both of the function types :c:func`ARKSpilsMassPrecSetupFn()` and
   c:func:`ARKSpilsMassPrecSolveFn()` are described in the section
   :ref:`CInterface.UserSupplied`. 



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



.. c:function:: int ARKSpilsSetMassEpsLin(void *arkode_mem, realtype eplifac)

   Specifies the factor by which the tolerance on the
   nonlinear iteration is multiplied to get a tolerance on the mass
   matrix linear iteration.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `eplifac` -- linear convergence safety factor :math:`(\ge 0.0)`.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if successful.
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``.
      * ARKSPILS_MASSMEM_NULL if the mass matrix solver memory was ``NULL``.
      * ARKSPILS_ILL_INPUT if an input has an illegal value.
   
   **Notes:** Passing a value `eplifac` of 0.0 indicates to use the default value of 0.05.





Rootfinding optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
--------------------------------

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
------------------------------

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
ratio `nniters`/`nsteps` measures the performance of the nonlinear
iteration in solving the nonlinear systems at each stage; 
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
No. of attempted steps                               :c:func:`ARKodeGetNumStepAttempts()`
No. of calls to `fe` and `fi` functions              :c:func:`ARKodeGetNumRhsEvals()`
No. of calls to linear solver setup function         :c:func:`ARKodeGetNumLinSolvSetups()`
No. of calls to mass matrix solver                   :c:func:`ARKodeGetNumMassSolves()`
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



.. c:function:: int ARKodeGetNumStepAttempts(void *arkode_mem, long int *step_attempts)

   Returns the cumulative number of steps attempted by the solver (so far).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `step_attempts` -- number of steps attempted by solver.
   
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
      * `nlinsetups` -- number of linear solver setup calls made.
   
   **Return value:**  
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL if the ARKode memory was ``NULL``



.. c:function:: int ARKodeGetNumMassSolves(void *arkode_mem, long int *nMassSolves)

   Returns the number of calls made to the mass matrix solver (so far).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nMassSolves` -- number of mass matrix solves made.
   
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



.. c:function:: int ARKodeGetIntegratorStats(void *arkode_mem, long int *nsteps, long int *expsteps, long int *accsteps, long int *step_attempts, long int *nfe_evals, long int *nfi_evals, long int *nlinsetups, long int *netfails, realtype *hinused, realtype *hlast, realtype *hcur, realtype *tcur)

   Returns many of the most useful integrator statistics in a single call.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nsteps` -- number of steps taken in the solver.
      * `expsteps` -- number of stability-limited steps taken in the solver.
      * `accsteps` -- number of accuracy-limited steps taken in the solver.
      * `step_attempts` -- number of steps attempted by the solver.
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional outputs are available from the ARKDLS
modules: workspace requirements, number of calls to the Jacobian
routine, number of calls to the mass matrix routine, number of calls
to the implicit right-hand side routine for finite-difference Jacobian
approximation, and last return value from an ARKDLS function.  Note
that, where the name of an output would otherwise conflict with the
name of an optional output from the main solver, a suffix LS (for
Linear Solver) or MLS (for Mass Linear Solver) has been added here
(e.g. `lenrwLS`).  


.. _CInterface.ARKDlsOutputTable:

Table: Optional outputs for ARKDLS
""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

===================================================  ===================================
Optional output                                      Function name
===================================================  ===================================
Size of real and integer workspaces                  :c:func:`ARKDlsGetWorkSpace()`
Size of real and integer workspaces                  :c:func:`ARKDlsGetMassWorkSpace()`
No. of Jacobian evaluations                          :c:func:`ARKDlsGetNumJacEvals()`
No. of mass matrix evaluations                       :c:func:`ARKDlsGetNumMassEvals()`
No. of `fi` calls for finite diff. Jacobian evals    :c:func:`ARKDlsGetNumRhsEvals()`
Last return flag from a linear solver function       :c:func:`ARKDlsGetLastFlag()`
Last return flag from a mass matrix solver function  :c:func:`ARKDlsGetLastMassFlag()`
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



.. c:function:: int ARKDlsGetMassWorkSpace(void *arkode_mem, long int *lenrwMLS, long int *leniwMLS)

   Returns the real and integer workspace used by the
   ARKDLS mass matrix linear solver (ARKDENSE or ARKBAND).
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `lenrwMLS` -- the number of ``realtype`` values in the ARKDLS workspace.
      * `leniwMLS` -- the number of integer values in the ARKDLS workspace.
   
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



.. c:function:: int ARKDlsGetNumMassEvals(void *arkode_mem, long int *nmevals)

   Returns the number of calls made to the ARKDLS
   (dense or band) mass matrix construction routine.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nmevals` -- number of calls to the mass matrix function.
   
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



.. c:function:: int ARKDlsGetLastMassFlag(void *arkode_mem, long int *mlsflag)

   Returns the last return value from an ARKDLS mass matrix solve routine.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `mlsflag` -- the value of the last return flag from an ARKDLS
	mass matrix solver function.
   
   **Return value:**  
      * ARKDLS_SUCCESS if successful
      * ARKDLS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKDLS_LMEM_NULL if the linear solver memory was ``NULL``
   
   **Notes:** If the ARKDENSE setup function failed
   (i.e. :c:func:`ARKode()` returned ARK_LSETUP_FAIL), then the
   value of `lsflag` is equal to the column index (numbered from
   one) at which a zero diagonal element was encountered during the LU
   factorization of the (dense or banded) mass matrix.  For all
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional outputs are available from the ARKSPILS
modules: workspace requirements, number of linear iterations, number
of linear convergence failures, number of calls to the preconditioner
setup and solve routines, number of calls to the Jacobian-vector
product routine, number of calls to the mass-matrix-vector
product routine, number of calls to the implicit right-hand side
routine for finite-difference Jacobian-vector product approximation,
and last return value from a linear solver function.  Note that, where
the name of an output would otherwise conflict with the name of an
optional output from the main solver, a suffix LS (for Linear Solver)
or MLS (for Mass Linear Solver) has been added here (e.g. `lenrwLS`). 


.. _CInterface.ARKSpilsOutputTable:

Table: Optional outputs for ARKSPILS
""""""""""""""""""""""""""""""""""""""""""""

.. cssclass:: table-bordered

===========================================================  ======================================== 
Optional output                                              Function name
===========================================================  ========================================
Size of real and integer workspaces                          :c:func:`ARKSpilsGetWorkSpace()`
Size of real and integer mass matrix solve workspaces        :c:func:`ARKSpilsGetMassWorkSpace()`
No. of preconditioner evaluations                            :c:func:`ARKSpilsGetNumPrecEvals()`
No. of mass matrix preconditioner evaluations                :c:func:`ARKSpilsGetNumMassPrecEvals()`
No. of preconditioner solves                                 :c:func:`ARKSpilsGetNumPrecSolves()`
No. of mass matrix preconditioner solves                     :c:func:`ARKSpilsGetNumMassPrecSolves()`
No. of linear iterations                                     :c:func:`ARKSpilsGetNumLinIters()`
No. of mass matrix linear iterations                         :c:func:`ARKSpilsGetNumMassIters()`
No. of linear convergence failures                           :c:func:`ARKSpilsGetNumConvFails()`
No. of mass matrix solver convergence failures               :c:func:`ARKSpilsGetNumMassConvFails()`
No. of Jacobian-vector product evaluations                   :c:func:`ARKSpilsGetNumJtimesEvals()`
No. of mass-matrix-vector product evaluations                :c:func:`ARKSpilsGetNumMtimesEvals()`
No. of `fi` calls for finite diff. Jacobian-vector evals.    :c:func:`ARKSpilsGetNumRhsEvals()`
Last return from a linear solver function                    :c:func:`ARKSpilsGetLastFlag()`
Last return from a mass matrix solver function               :c:func:`ARKSpilsGetLastMassFlag()`
Name of constant associated with a return flag               :c:func:`ARKSpilsGetReturnFlagName()`
===========================================================  ========================================




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
   :math:`9n` ``realtype`` words for ARKSPBCG, :math:`11n`
   ``realtype`` words for ARKSPTFQMR, :math:`(2m+4)n+m(m+4)+1`
   ``realtype`` words for ARKSPFGMR, and :math:`4n+1`
   ``realtype`` words for ARKPCG.
   
   In a parallel setting, the above values are global, summed over all
   processors.



.. c:function:: int ARKSpilsGetMassWorkSpace(void *arkode_mem, long int *lenrwMLS, long int *leniwMLS)

   Returns the global sizes of the ARKSPILS real and integer workspaces.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `lenrwMLS` -- the number of ``realtype`` values in the ARKSPILS workspace.
      * `leniwMLS` -- the number of integer values in the ARKSPILS workspace.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``
   
   **Notes:** In terms of the problem size :math:`n` and maximum Krylov subspace
   size :math:`m`, the actual size of the real workspace is roughly:
   :math:`(m+5)n+m(m+4)+1` ``realtype`` words for ARKSPGMR,
   :math:`9n` ``realtype`` words for ARKSPBCG, :math:`11n`
   ``realtype`` words for ARKSPTFQMR, :math:`(2m+4)n+m(m+4)+1`
   ``realtype`` words for ARKSPFGMR, and :math:`4n+1`
   ``realtype`` words for ARKPCG.
   
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



.. c:function:: int ARKSpilsGetNumMassPrecEvals(void *arkode_mem, long int *nmpevals)

   Returns the total number of mass matrix preconditioner evaluations,
   i.e. the number of calls made to `psetup`.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nmpevals` -- the current number of calls to `psetup`.
   
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



.. c:function:: int ARKSpilsGetNumMassPrecSolves(void *arkode_mem, long int *nmpsolves)

   Returns the number of calls made to the mass matrix preconditioner
   solve function, `psolve`.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nmpsolves` -- the number of calls to `psolve`.
   
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



.. c:function:: int ARKSpilsGetNumMassIters(void *arkode_mem, long int *nmiters)

   Returns the cumulative number of mass matrix solver iterations.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nmiters` -- the current number of mass matrix solver linear iterations.
   
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



.. c:function:: int ARKSpilsGetNumMassConvFails(void *arkode_mem, long int *nmcfails)

   Returns the cumulative number of mass matrix solver convergence failures.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nmcfails` -- the current number of mass matrix solver convergence failures.
   
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



.. c:function:: int ARKSpilsGetNumMtimesEvals(void *arkode_mem, long int *nmvevals)

   Returns the cumulative number of calls made to the
   mass-matrix-vector product function, `mtimes`.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nmvevals` -- the current number of calls to `mtimes`.
   
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
   SPGMR_PSET_FAIL_UNREC, SPBCG_PSET_FAIL_UNREC,
   SPTFQMR_PSET_FAIL_UNREC, SPFGMR_PSET_FAIL_UNREC, or
   PCG_PSET_FAIL_UNREC.  
   
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

   If the ARKSPFGMR solve function failed (:c:func:`ARKode()`
   returned ARK_LSOLVE_FAIL), then `lsflag` contains the error
   return flag from SpfgmrSolve and will be one of:
   SPFGMR_MEM_NULL, indicating that the SPFGMR memory is
   ``NULL``; SPFGMR_ATIMES_FAIL_UNREC, indicating an unrecoverable
   failure in the :math:`J*v` function; SPFGMR_PSOLVE_FAIL_UNREC,
   indicating that the preconditioner solve function `psolve` failed
   unrecoverably; SPFGMR_GS_FAIL, indicating a failure in the
   Gram-Schmidt procedure; or SPFGMR_QRSOL_FAIL, indicating that
   the matrix :math:`R` was found to be singular during the QR solve
   phase. 
  
   If the ARKPCG solve function failed (:c:func:`ARKode()`
   returned ARK_LSOLVE_FAIL), then `lsflag` contains the error
   return flag from PcgSolve and will be one of:
   PCG_MEM_NULL, indicating that the PCG memory is
   ``NULL``; PCG_ATIMES_FAIL_UNREC, indicating an unrecoverable
   failure in the :math:`J*v` function; or
   PCG_PSOLVE_FAIL_UNREC, indicating that the preconditioner
   solve function `psolve` failed unrecoverably. 
   



.. c:function:: int ARKSpilsGetLastMassFlag(void *arkode_mem, long int *msflag)

   Returns the last return value from an ARKSPILS mass matrix solver routine.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `msflag` -- the value of the last return flag from an
        ARKSPILS mass matrix solver function.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if successful
      * ARKSPILS_MEM_NULL if the ARKode memory was ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory was ``NULL``
   
   **Notes:** The values of `msflag` for each of the various solvers
   will match those described above for the function
   :c:func:`ARKSpilsGetLastFlag()`.   



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
-------------------------------------

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
   
   **Notes:** If an error occurred, :c:func:`ARKodeReInit()` also
   sends an error message to the error handler function.




.. _CInterface.Resizing:

ARKode system resize function
-------------------------------------

For simulations involving changes to the number of equations and
unknowns in the ODE system (e.g. when using spatially-adaptive finite
elements), the ARKode integrator may be "resized" between integration
steps, through calls to the :c:func:`ARKodeResize()` function.
This function modifies ARKode's internal memory structures to
use the new problem size, without destruction of the temporal
adaptivity heuristics.  It is assumed that the dynamical time scales
before and after the vector resize will be comparable, so that all
time-stepping heuristics prior to calling :c:func:`ARKodeResize()`
remain valid after the call.  If instead the dynamics should be
re-calibrated, the ARKode memory structure should be deleted with a
call to :c:func:`ARKodeFree()`, and re-created with calls to
:c:func:`ARKodeCreate()` and :c:func:`ARKodeInit()`. 

To aid in the vector-resize operation, the user can supply a vector
resize function that will take as input a vector with the previous
size, and transform it in-place to return a corresponding vector of
the new size.  If this function (of type :c:func:`ARKVecResizeFn()`)
is not supplied (i.e. is set to ``NULL``), then all existing vectors
internal to ARKode will be destroyed and re-cloned from the input
vector. 

In the case that the dynamical time scale should be modified slightly
from the previous time scale, an input `hscale` is allowed, that will
re-scale the upcoming time step by the specified factor.  If a value
:math:`hscale \le 0` is specified, the default of 1.0 will be used.



.. c:function:: int ARKodeResize(void *arkode_mem, N_Vector ynew, realtype hscale, realtype t0, ARKVecResizeFn resize, void *resize_data)

   Re-initializes ARKode with a different state vector but with
   comparable dynamical time scale.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `ynew` -- the newly-sized solution vector, holding the current
	dependent variable values :math:`y(t_0)`. 
      * `hscale` -- the desired scaling factor for the dynamical time
	scale (i.e. the next step will be of size :math:`h*hscale`)
      * `t0` -- the current value of the independent variable
	:math:`t_0` (this must be consistent with `ynew`.
      * `resize` -- the user-supplied vector resize function (of type
	:c:func:`ARKVecResizeFn()`.
      * `resize_data` -- the user-supplied data structure to be passed
	to `resize` when modifying internal ARKode vectors.
   
   **Return value:** 
      * ARK_SUCCESS if successful
      * ARK_MEM_NULL  if the ARKode memory was ``NULL``
      * ARK_NO_MALLOC if `arkode_mem` was not allocated.
      * ARK_ILL_INPUT if an argument has an illegal value.
   
   **Notes:** If an error occurred, :c:func:`ARKodeResize()` also sends an error
   message to the error handler function.


Resizing the linear solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using any of the built-in linear solver modules, the linear
solver memory structures must also be resized.  At present, none of
these include a solver-specific 'resize' function, so the linear
solver memory must be destroyed and re-allocated **following** each
call to :c:func:`ARKodeResize()`.  For each of the built-in ARKDLS and
ARKSPILS linear solvers, the specification call itself
(e.g. :c:func:`ARKDense()` or :c:func:`ARKSpgmr()`) will internally
destroy the solver-specific memory prior to re-allocation.   

If any user-supplied routines are provided to aid the linear solver
(e.g. Jacobian construction, Jacobian-vector product,
mass-matrix-vector product, preconditioning), then the corresponding
"Set" routines must be called again **following** the solver
re-specification.


Resizing the absolute tolerance array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If using array-valued absolute tolerances, the absolute tolerance
vector will be invalid after the call to :c:func:`ARKodeResize()`, so
the new absolute tolerance vector should be re-set **following** each
call to :c:func:`ARKodeResize()` through a new call to
:c:func:`ARKodeSVtolerances()`.   

If scalar-valued tolerances or a tolerance function was specified
through either :c:func:`ARKodeSStolerances()` or
:c:func:`ARKodeWFtolerances()`, then these will remain valid so no
further action is necessary. 


.. note:: For an example of :c:func:`ARKodeResize()` usage, see the
	  supplied serial C example problem, ``ark_heat1D_adapt.c``.
