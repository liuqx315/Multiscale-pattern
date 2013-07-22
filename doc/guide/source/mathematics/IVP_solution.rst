:tocdepth: 3


.. _Mathematics.IVP:

IVP solution
=================

The methods used in ARKode are variable-step 
:index:`additive Runge-Kutta methods` (ARK), based on formulas of the
form 

.. math::
   M z_i &= M y_{n-1} + h_n \sum_{j=0}^{i-1} A^E_{i,j} f_E(t_{n-1} + c_j h_n, z_j) 
                 + h_n \sum_{j=0}^{i}   A^I_{i,j} f_I(t_{n-1} + c_j h_n, z_j),
   \quad i=1,\ldots,s, \\
   M y_n &= M y_{n-1} + h_n \sum_{i=0}^{s} b_i \left(f_E(t_{n-1} + c_i h_n, z_i) 
                 + f_I(t_{n-1} + c_i h_n, z_i)\right), \\
   M \tilde{y}_n &= M y_{n-1} + h_n \sum_{i=0}^{s} \tilde{b}_i 
       \left(f_E(t_{n-1} + c_i h_n, z_i) + f_I(t_{n-1} + c_i h_n, z_i)\right).
   :label: ARK

Here the :math:`y_n` are computed approximations to :math:`y(t_n)`,
:math:`\tilde{y}_n` are lower-order embedded solutions (used in error
estimation), and :math:`h_n \equiv t_n - t_{n-1}` is the step size.
The coefficients :math:`A^E \in \Re^{s\times s}`, :math:`A^I \in
\Re^{s\times s}`, :math:`b \in \Re^{s}` and :math:`c \in \Re^{s}` 
correspond with the explicit and implicit Butcher tables (that must
share :math:`b` and :math:`c`) in an ARK pair.  The user of ARKode
must choose appropriately between one of three classes of methods:
*multi-rate*, *nonstiff* and *stiff*.  All of the Butcher tables
encoding the coefficients :math:`c`, :math:`A`, :math:`b` and
:math:`\tilde{b}` are further described in the section
:ref:`Mathematics.Butcher`. 

For multi-rate problems, a user must provide both of the functions
:math:`f_E` and :math:`f_I`.  On such problems, ARKode implements the
ARK methods proposed in [KC2003]_, allowing for methods
having order :math:`q = \{3,4,5\}`.

For nonstiff problems, a user may specify that :math:`f_I = 0`, or in
other words that the equation :eq:`IVP` reduces down to the non-split
IVP 

.. math::
   M\dot{y} = f_E(t,y), \qquad y(t_0) = y_0.
   :label: IVP_explicit

In this scenario, the Butcher table :math:`A^I=0` in :eq:`ARK`, and
the ARK methods reduce to classical :index:`explicit Runge-Kutta methods` 
(ERK).  For these classes of methods, ARKode allows orders of accuracy
:math:`q = \{2,3,4,5,6\}`.  These default to the Heun-Euler,
Bogacki-Shampine, Zonneveld, Cash-Karp and Verner methods,
respectively.

Finally, for stiff problems the user may specify that :math:`f_E = 0`,
so the equation :eq:`IVP` reduces to the non-split IVP 

.. math::
   M\dot{y} = f_I(t,y), \qquad y(t_0) = y_0,
   :label: IVP_implicit

the Butcher table :math:`A^E=0` in :eq:`ARK`, and the ARK methods
reduce to classical :index:`diagonally-implicit Runge-Kutta methods` 
(DIRK).  For these classes of methods, ARKode allows orders of
accuracy :math:`q = \{3,4,5\}`, that default to the SDIRK 2(1),
ARK3(2)4L[2]SA, SDIRK 5(4) and ARK5(4)8L[2]SA methods, respectively.

For both the DIRK and ARK methods corresponding to :eq:`IVP` and
:eq:`IVP_implicit`, a nonlinear system

.. math::
   G(z_i) \equiv M z_i - h_n A^I_{i,i} f_I(t_{n-1} + c_i h_n, z_i) - a_i = 0
   :label: Residual

must be solved for each stage :math:`z_i, i=1,\ldots,s`, where 

.. math::
   a_i \equiv M y_{n-1} + h_n \sum_{j=0}^{i-1} \left[
      A^E_{i,j} f_E(t_{n-1} + c_j h_n, z_j) +
      A^I_{i,j} f_I(t_{n-1} + c_j h_n, z_j) \right]
   
for the ARK methods, or 

.. math::
   a_i \equiv M y_{n-1} + h_n \sum_{j=0}^{i-1} 
      A^I_{i,j} f_I(t_{n-1} + c_j h_n, z_j)
   
for the DIRK methods.  For these nonlinear systems, ARKode can use a
type of :index:`Newton iteration`,

.. math::
   z_i^{(m+1)} = z_i^{(m)} + \delta^{(m+1)},
   :label: Newton_iteration

where :math:`m` is the Newton iteration index, and the :index:`Newton
update` :math:`\delta^{(m+1)}` in turn requires the solution of the
linear :index:`Newton system` 

.. math::
   A\left(z_i^{(m)}\right) \delta^{(m+1)} = -G\left(z_i^{(m)}\right), 
   :label: Newton_system

in which

.. math::
   A \approx M - \gamma J, \quad J = \frac{\partial f_I}{\partial y},
   \quad\text{and}\quad \gamma = h_n A^I_{i,i}.
   :label: NewtonMatrix

Alternately, ARKode may solver for each stage :math:`z_i,
i=1,\ldots,s` using an :index:`Anderson-accelerated fixed point
iteration`

.. math::
   z_i^{(m+1)} = g(z_i^{(m)}).
   :label: AAFP_iteration

For details on how this iteration is performed, we refer the reader to
the reference [WN2011]_.

Unlike Newton-based methods, these fixed-point solvers do not require
the solution of any linear systems.  As a result, they may converge
somewhat more slowly than Newton-like methods, but each iteration can
be much less expensive.  Moreover, these solvers do allow for user
specification of the Anderson-accelerates subspace size, :math:`m_k`.
While the required amount of solver memory grows proportionately to
:math:`m_k n` where :math:`n` is the size of the ODE system, larger
values of :math:`m_k` may result in faster convergence.  In our
experience, this improvement may be significant even for "small"
values, e.g. :math:`1\le m_k\le 5`, and that convergence may even
deteriorate for larger values of :math:`m_k`.

While ARKode uses the Newton iteration as its default solver due to
its increased robustness on very stiff problems, it is highly
recommended that users also consider the fixed-point solver for their
when attempting a new problem.


It is well-known that both the efficiency and robustness of any
nonlinear solver algorithm intimately depends on the choice of a good
initial guess.  In ARKode, the initial guess for either nonlinear
solution method is a predicted value :math:`z_i^{(0)}` that is
computed explicitly from the previously-computed data
(e.g. :math:`y_{n-2}`, :math:`y_{n-1}`, and :math:`z_j` where
:math:`j<i`).  For further information on the predictor algorithms
implemented in ARKode, see the section :ref:`Mathematics.Predictors`.

When a Newton-based method is chosen for solving each nonlinear
system, a linear system of equations must be solved at each nonlinear
iteration.  For this solve ARKode provides several choices, including
the option of a user-supplied linear solver module.  The linear solver
modules distributed with SUNDIALS are organized into two families: a
*direct* family comprising direct linear solvers for dense or banded
matrices, and a *spils* family comprising scaled, preconditioned,
iterative (Krylov) linear solvers.  The methods offered through these
modules are as follows:

* dense direct solvers, using either an internal implementation or a
  BLAS/LAPACK implementation (serial version only),
* band direct solvers, using either an internal implementation or a
  BLAS/LAPACK implementation (serial version only),
* SPGMR, a scaled, preconditioned GMRES (Generalized Minimal Residual
  method) solver without restarts,
* SPBCG, a scaled, preconditioned Bi-CGStab (Bi-Conjugate Gradient
  Stable method) solver,
* SPTFQMR, a scaled, preconditioned TFQMR (Transpose-free
  Quasi-Minimal Residual method) solver, or
* SPFGMR, a scaled, preconditioned Flexible GMRES (Generalized Minimal
  Residual method) solver without restarts,
* PCG, a preconditioned conjugate gradient solver for symmetric linear
  systems.

For large stiff systems where direct methods are infeasible, the
combination of an implicit Runge-Kutta integrator and a preconditioned
Krylov method (SPGMR, SPBCG, SPTFQMR, SPFGMR or PCG) can yield a
powerful tool because it combines established methods for stiff
integration, nonlinear solver iteration, and Krylov (linear) iteration
with a problem-specific treatment of the dominant sources of
stiffness, in the form of a user-supplied preconditioner matrix
[BH1989]_.  We note that the direct linear solvers
provided by SUNDIALS (dense and band) can only be used with the serial
vector representations.

In the process of controlling errors at various levels (time
integration, nonlinear solution, linear solution), ARKode uses a
:index:`weighted root-mean-square norm`, denoted
:math:`\|\cdot\|_{WRMS}`, for all error-like quantities,

.. math::
   \|v\|_{WRMS} = \left( \frac{1}{N} \sum_{i=1}^N \left(v_i\,
   w_i\right)^2\right)^{1/2}. 
   :label: WRMS_NORM

The multiplicative :index:`error weight vector`  :math:`w` is based
on the current solution and on the relative and absolute tolerances
input by the user, namely

.. math::
   w_i = \frac{1}{RTOL\cdot |y_i| + ATOL_i}.
   :label: EWT

Since :math:`1/w_i` represents a tolerance in the component
:math:`y_i`, a vector whose WRMS norm is 1 is regarded as "small."
For brevity, we will typically drop the subscript WRMS on norms in the
remainder of this section.

In the case of a direct solver (dense or band), ARKode utilizes a
variant on the Newton method called a *modified Newton iteration*. In
such methods, the matrix :math:`A` is held fixed, i.e. each Newton
iteration is computed from the modified equation

.. math::
   \tilde{A}\left(z_i^{(m)}\right) \delta^{(m+1)} = -G\left(z_i^{(m)}\right), 
   :label: modified_Newton_system

in which

.. math::
   \tilde{A} \approx M - \tilde{\gamma} \tilde{J}, \quad \tilde{J} =
   \frac{\partial f_I}{\partial y}(\tilde y), \quad\text{and}\quad
   \tilde{\gamma} = \tilde{h} A^I_{i,i}. 
   :label: modified_NewtonMatrix

Here, the solution :math:`\tilde{y}` and step size :math:`\tilde{h}`
upon which the modified Jacobian rely, are merely values of the
solution and step size from a previous iteration.  In other words, the
matrix :math:`\tilde{A}` is only computed rarely, and reused for
repeated stage solves.  

In the case that a Krylov linear solver is chosen, ARKode utilizes a
Newton method variant called an *Inexact Newton iteration*.  Here, the
matrix :math:`A` is not itself required since the algorithms only
require the product of this matrix with a given vector.  Additionally,
each Newton system :eq:`Newton_system` is not solved completely, since
the Krylov solvers are iterative (hence the "inexact" in the name).
Resultingly. for these linear solvers :math:`A` is applied in a
matrix-free manner, with matrix-vector products :math:`Jv` obtained by
either difference quotients or a user-supplied routine.  As with the
direct linear solver scenario, in which :math:`A` was reused between
solves, ARKode's inexact Newton iteration also recomputes the
preconditioner matrix :math:`P` as infrequently as possible to balance
the high costs of matrix operations against other costs.  

Specifically, we update the Newton matrix :math:`\tilde{A}` or
preconditioner matrix :math:`P` only in the following circumstances:

* when starting the problem,
* when more than 20 steps have been taken since the last update (this
  value may be changed via the ``msbp`` argument to
  :c:func:`ARKodeSetLSetupConstants()`), 
* when the value :math:`\bar{\gamma}` of :math:`\gamma` at the last
  update satisfies :math:`\left|\gamma/\bar{\gamma} - 1\right| > 0.2`
  (this tolerance may be changed via the ``dgmax`` argument to 
  :c:func:`ARKodeSetLSetupConstants()`), 
* when a non-fatal convergence failure just occurred, or
* when an error test failure just occurred.

When an update is forced due to a convergence failure, an update of
:math:`\tilde{A}` or :math:`P` may or may not involve a reevaluation of
:math:`J` (in :math:`\tilde{A}`) or of Jacobian data (in :math:`P`),
depending on whether errors in the Jacobian were the likely cause of the
failure.  More generally, the decision is made to reevaluate :math:`J`
(or instruct the user to reevaluate Jacobian data in :math:`P`) when:

* starting the problem,
* more than 50 steps have been taken since the last evaluation,
* a convergence failure occurred with an outdated matrix, and the
  value :math:`\bar{\gamma}` of :math:`\gamma` at the last update
  satisfies :math:`\left|\gamma/\bar{\gamma} - 1\right| > 0.2`,
* a convergence failure occurred that forced a step size reduction.



The stopping test for all of ARKode's nonlinear solvers is related to
the subsequent local error test, with the goal of keeping the
nonlinear iteration errors from interfering with local error control.
As described below, the final computed value of each stage solution
:math:`z_i^{(m)}` must satisfy a local error test
:math:`\|z_i^{(m)} - z_i^{(0)}\| \le \epsilon`.  Letting
:math:`z_i` denote the true solution to the nonlinear problem
:eq:`Residual`, we want to ensure that the iteration error
:math:`z_i - z_i^{(m)}` is small relative to :math:`\epsilon`,
specifically that it is less than :math:`0.2\epsilon` (the safety
factor 0.2 may be changed by the user via the
:c:func:`ARKodeSetNonlinConvCoef()` function).  For this, we also
estimate the linear convergence rate :math:`R_i` of the nonlinear
iteration as follows.  We first initialize :math:`R_i` to 1, and reset
:math:`R_i=1` when either :math:`A` or :math:`P` are updated.  After
computing a nonlinear correction :math:`\delta^{(m)} = z_i^{(m)} -
z_i^{(m-1)}`, we update :math:`R_i` if :math:`m>1` as

.. math:: 
   R_i \leftarrow \max\{ 0.3 R_i, \left\|\delta^{(m)}\right\| / \left\|\delta^{(m-1)}\right\| \}.

where the factor 0.3 is user-modifiable as the ``crdown`` input to the
the function :c:func:`ARKodeSetNonlinConstants()`.  Denoting the
combined time step solution from the true stage solutions :math:`z_i`
as :math:`y_n`, and the combined time step solution from the computed
stage solutions :math:`z_i^{(m)}` as :math:`\tilde{y}_n` we use the
estimate 

.. math::
   \left\| y_n - \tilde{y}_n \right\| \approx 
   \max_i \left\| z_i^{(m+1)} - z_i^{(m)} \right\| \approx
   \max_i R_i \left\| z_i^{(m)} - z_i^{(m-1)} \right\| =
   \max_i R_i \left\| \delta^{(m)} \right\|.

Therefore the convergence (stopping) test for the nonlinear iteration
for each stage is 

.. math::
   R_i \left\|\delta^{(m)} \right\| < 0.2\epsilon,

where the factor 0.2 is user-modifiable as the ``nlscoef`` input to the
the function :c:func:`ARKodeSetNonlinConvCoef()`.  We allow at most 3
nonlinear iterations (this may be modified through the function
:c:func:`ARKodeSetMaxNonlinIters()`).  We also declare the 
nonlinear iteration to be divergent if any of the ratios
:math:`\|\delta^{(m)}\| / \|\delta^{(m-1)}\| > 2.3` with :math:`m>1`
(the value 2.3 may be modified as the ``rdiv`` input to the function 
:c:func:`ARKodeSetNonlinRDiv()`).  If convergence fails in the fixed
point iteration, or in the Newton iteration with :math:`J` or
:math:`A` current, we must then reduce the step size by a factor of
0.25 (modifiable via the ``etacf`` input to the
:c:func:`ARKodeSetMaxCFailGrowth()` function).  The integration 
is halted after 10 convergence failures (modifiable via the
:c:func:`ARKodeSetMaxConvFails()` function).

When a Krylov method is used to solve the linear systems
:eq:`Newton_system`, its errors must also be controlled; this error
control also uses the local error test constant.  To this end, we
approximate the linear iteration error in the solution vector
:math:`\delta^{(m)}` using the preconditioned residual vector.  In an
attempt to ensure that the linear iteration errors do not interfere
with the nonlinear solution error and local time integration error
controls, we require that the norm of the preconditioned residual be
less than :math:`0.05\cdot(0.2\epsilon)`.  Here 0.2 is the same value
as that used above for the nonlinear error control; the value 0.05 is
not currently modifiable by the user.

With the direct and band solvers for the linear systems
:eq:`Newton_system`, the Jacobian may be supplied by a user routine,
or approximated by finite-differences.  In the case of differencing,
we use the standard approximation

.. math::
   A_{i,j}(t,y) = (M\,e_j)_i - \gamma 
   \frac{f_{I,i}(t,y+\sigma_j e_j) - f_{I,i}(t,y)}{\sigma_j}.

Here :math:`e_j` is the jth unit vector, and the increments
:math:`\sigma_j` are given by 

.. math::
   \sigma_j = \max\left\{ \sqrt{U}\, |y_j|, \sigma_0/w_j \right\},

where :math:`U` is the unit roundoff, :math:`\sigma_0` is a
dimensionless value, and :math:`w_j` is the error weight defined in
:eq:`EWT`.  In the dense case, this approach requires :math:`N`
evaluations of :math:`f_I`, one for each column of :math:`J`.  In the
band case, the columns of :math:`J` are computed in groups, using the
Curtis-Powell-Reid algorithm, with the number of :math:`f_I`
evaluations equal to the bandwidth.

As will be further discussed in the section
:ref:`Mathematics.Preconditioning`, in the case of a Krylov method,
preconditioning may be applied on the left, right, or on both sides of
:math:`A`, with user-supplied routines for the preconditioner setup
and solve operations.  Optionally, a user may supply a routine to
compute the required matrix-vector products :math:`Jv`.
If a routine for :math:`Jv` is not supplied, these products will be
computed with directional differencing using the formula

.. math::
   Jv = \frac{f_I(t,y+\sigma_j v) - f_I(t,y)}{\sigma_j},

where the increment :math:`\sigma = 1/\|v\|` to ensure that 
:math:`\|\sigma v\| = 1`.

In the following four sub-sections (:ref:`Mathematics.Preconditioning`,
:ref:`Mathematics.Predictors`, :ref:`Mathematics.Adaptivity` and
:ref:`Mathematics.Stability`), we provide details on optional
user-supplied information that can be used to better control the
behavior of ARKode.  In these sections, we also discuss the algorithms
currently provided by ARKode.  Finally, in the last sub-section of this
chapter, :ref:`Mathematics.Rootfinding`, we discuss the algorithms
providing root-finding capabilities within ARKode.
