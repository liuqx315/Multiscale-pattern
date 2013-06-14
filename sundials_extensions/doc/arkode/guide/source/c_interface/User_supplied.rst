:tocdepth: 3



.. _CInterface.UserSupplied:

User-supplied functions
=============================

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
-----------------------------

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
--------------------------------------

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
--------------------------------------

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
--------------------------------------

As an alternative to using one of the built-in time step adaptivity
methods for controlling solution error, the user may provide a
function of type :c:func:`ARKAdaptFn()` to compute a target step size
:math:`h` for the next integration step.  These steps should be chosen
as the maximum value such that the error estimates remain below 1.



.. c:function:: typedef int (*ARKAdaptFn)(N_Vector y, realtype t, realtype h1, realtype h2, realtype h3, realtype e1, realtype e2,  realtype e3, int q, int p, realtype *hnew, void *user_data)

   This function implements a time step adaptivity algorithm
   that chooses :math:`h` satisfying the error tolerances..
   
   **Arguments:**
      * `y` -- the current value of the dependent variable vector, :math:`y(t)`.
      * `t` -- the current value of the independent variable.
      * `h1` -- the current step size, :math:`t_m - t_{m-1}`.
      * `h2` -- the previous step size, :math:`t_{m-1} - t_{m-2}`.
      * `h3` -- the step size :math:`t_{m-2}-t_{m-3}`.
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
--------------------------------------

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
      * `hstab` -- the output value with the absolute value of the
 	maximum stable step size. 
      * `user_data` -- a pointer to user data, the same as the
        `estab_data` parameter that was passed to :c:func:`ARKodeSetStabilityFn()`.
   
   **Return value:** 
   An ARKExpStabFn function should return 0 if it
   successfully set the upcoming stable step size, and a non-zero
   value otherwise.
   
   **Notes:**  If this function is not supplied, or if it returns
   `hstab \le 0.0`, then ARKode will assume that there is no explicit
   stability restriction on the time step size.



.. _CInterface.RootfindingFn:

Rootfinding function
--------------------------------------

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
--------------------------------------------------------------

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
--------------------------------------------------------------

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
--------------------------------------------------------------

If one of the Krylov iterative linear solvers SPGMR, SPBCG, 
SPTFQMR, or PCG is selected (i.e. ARKSp* is called in step 8 of the
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

Preconditioning (linear system solution)^
--------------------------------------------------------------

If one of the Krylov iterative linear solvers SPGMR, SPBCG,
SPTFQMR, or PCG is selected, and preconditioning is used, then the user
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
--------------------------------------------------------------

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
