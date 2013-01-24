.. _Mathematics:

Mathematical Considerations
===========================

ARKode solves ODE initial value problems (IVPs) in real :math:`N`
-space, which we write in the abstract form

.. math::
   M\dot{y} = f_E(t,y) + f_I(t,y), \qquad y(t_0) = y_0,
   :label: IVP

where :math:`y \in \Re^N`, and :math:`M` is a user-specified
nonsingular operator from :math:`\Re^N \to \Re^N`. Here we use
:math:`\dot{y}` to denote :math:`dy/dt`, and the two right-hand side
components may be described as:

* :math:`f_E(t,y)` contains the "slow" time scale components of the
  system, that should be integrated using explicit methods.

* :math:`f_I(t,y)` contains the "fast" time scale components of the
  system, that should be integrated using implicit methods.

While we use :math:`t` to denote the independent variable, and usually
this is time, it certainly need not be.  ARKode may be used to solve
stiff, nonstiff and multi-rate problems.  Roughly speaking, stiffness
is characterized by the presence of at least one rapidly damped mode,
whose time constant is small compared to the time scale of the
solution itself.



.. _Mathematics.IVP:

IVP solution
---------------

The methods used in ARKode are variable-step 
:index:`additive Runge-Kutta methods` (ARK), based on formulas of the
form 

.. math::
   z_i &= y_{n-1} + h_n \sum_{j=0}^{i-1} A^E_{i,j} f_E(t_n + c_j h_n, z_j) 
                 + h_n \sum_{j=0}^{i}   A^I_{i,j} f_I(t_n + c_j h_n, z_j),
   \quad i=1,\ldots,s, \\
   y_n &= y_{n-1} + h_n \sum_{i=0}^{s} b_i \left(f_E(t_n + c_i h_n, z_i) 
                 + f_I(t_n + c_i h_n, z_i)\right), \\
   \tilde{y}_n &= y_{n-1} + h_n \sum_{i=0}^{s} \tilde{b}_i 
       \left(f_E(t_n + c_i h_n, z_i) + f_I(t_n + c_i h_n, z_i)\right).
   :label: ARK

Here the :math:`y_n` are computed approximations to :math:`y(t_n)`,
:math:`\tilde{y}_n` are lower-order embedded solutions (used in error
estimation), and :math:`h_n \equiv t_n - t_{n-1}` is the step size.
The coefficients :math:`A^E \in \Re^{s\times s}`, :math:`A^I \in
\Re^{s\times s}`, :math:`b \in \Re^{s}` and :math:`c \in \Re^{s}` 
correspond with the explicit and implicit Butcher tables (that must
share :math:`b` and :math:`c`) in an ARK pair.  The user of ARKode
must choose appropriately between one of three classes of methods:
*multi-rate*, *nonstiff* and *stiff*.

For multi-rate problems, a user must provide both of the functions
:math:`f_E` and :math:`f_I`.  On such problems, ARKode implements the
ARK methods proposed in [KennedyCarpenter2003]_, allowing for methods
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
accuracy :math:`q = \{3,4,5\}`, that default to the Billington, SDIRK
5(4) and Kvaerno(7,4,5) methods, respectively.

For both the DIRK and ARK methods corresponding to :eq:`IVP` and
:eq:`IVP_implicit`, a nonlinear system

.. math::
   G(z_i) \equiv z_i - h_n A^I_{i,i} f_I(t_n + c_i h_n, z_i) - a_i = 0
   :label: Residual

must be solved for each stage :math:`z_i, i=1,\ldots,s`, where 

.. math::
   a_i \equiv y_{n-1} + h_n \sum_{j=0}^{i-1} \left[
      A^E_{i,j} f_E(t_n + c_j h_n, z_j) +
      A^I_{i,j} f_I(t_n + c_j h_n, z_j) \right]
   
for the ARK methods, or 

.. math::
   a_i \equiv y_{n-1} + h_n \sum_{j=0}^{i-1} 
      A^I_{i,j} f_I(t_n + c_j h_n, z_j)
   
for the DIRK methods.  For these nonlinear systems, ARKode uses a
type of :index:`Newton iteration`, 

.. math::
   z_i^{(m+1)} = z_i^{(m)} + \delta^{(m+1)},
   :label: Newton_iteration

where :math:`m` is the Newton iteration index.  Here, the 
update :math:`\delta^{(m+1)}` in turn requires the solution of linear 
:index:`Newton systems`

.. math::
   A\left(z_i^{(m)}\right) \delta^{(m+1)} = -G\left(z_i^{(m)}\right), 
   :label: Newton_system

where

.. math::
   A \approx M - \gamma J, \quad J = \frac{\partial f_I}{\partial y},
   \quad\text{and}\quad \gamma = h_n A^I_{i,i}.
   :label: NewtonMatrix

The initial guess for the iteration is a predicted value
:math:`z_i^{(0)}` that is computed explicitly from the
previously-computed data (e.g. :math:`y_{n-2}`, :math:`y_{n-1}`,
and :math:`z_j` where :math:`j<i`).  For further information on the
predictor algorithms implemented in ARKode, see the section
:ref:`Mathematics.Predictors`.

For the solution of the linear systems within the Newton
iteration, ARKode provides several choices, including the option of a
user-supplied linear solver module.  The linear solver modules
distributed with SUNDIALS are organized into two families: a *direct*
family comprising direct linear solvers for dense or banded matrices,
and a *spils* family comprising scaled, preconditioned, iterative
(Krylov) linear solvers.  The methods offered through these modules
are as follows:

* dense direct solvers, using either an internal implementation or a
  BLAS/LAPACK implementation (serial version only),
* band direct solvers, using either an internal implementation or a
  BLAS/LAPACK implementation (serial version only),
* SPGMR, a scaled, preconditioned GMRES (Generalized Minimal Residual
  method) solver without restarts,
* SPBCG, a scaled, preconditioned Bi-CGStab (Bi-Conjugate Gradient
  Stable method) solver, or
* SPTFQMR, a scaled, preconditioned TFQMR (Transpose-free
  Quasi-Minimal Residual method) solver.

For large stiff systems where direct methods are infeasible, the
combination of an implicit Runge-Kutta integrator and a preconditioned
Krylov method (SPGMR, SPBCG or SPTFQMR) can yield a powerful tool
because it combines established methods for stiff integration,
nonlinear solver iteration, and Krylov (linear) iteration with a
problem-specific treatment of the dominant sources of stiffness, in
the form of a user-supplied preconditioner matrix
[BrownHindmarsh1989]_.  We note that the direct linear solvers
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

In the case of a direct solver (dense or band), the iteration is a
modified Newton iteration, in that the matrix :math:`A` is fixed
throughout the nonlinear iterations for a given stage :math:`z_i`.
However, for any of the Krylov methods, it is an Inexact Newton
iteration, in which :math:`A` is applied in a matrix-free manner, with
matrix-vector products :math:`Jv` obtained by either difference
quotients or a user-supplied routine.  The matrix :math:`A` (direct
cases) or a preconditioner matrix :math:`P` (Krylov cases) is obtained
as infrequently as possible to balance the high costs of matrix
operations against other costs.  Specifically, this matrix update
occurs when:

* starting the problem,
* more than 20 steps have been taken since the last update (this may
  be changed via the ``msbp`` argument to
  :c:func:`ARKodeSetLSetupConstants()`), 
* the value :math:`\bar{\gamma}` of :math:`\gamma` at the last update
  satisfies :math:`\left|\gamma/\bar{\gamma} - 1\right| > 0.2` (this
  tolerance may be changed via the ``dgmax`` argument to 
  :c:func:`ARKodeSetLSetupConstants()`), 
* a non-fatal convergence failure just occurred, or
* an error test failure just occurred.

When an update is forced due to a convergence failure, an update of
:math:`A` or :math:`P` may or may not involve a reevaluation of
:math:`J` (in :math:`A`) or of Jacobian data (in :math:`P`), depending
on whether errors in the Jacobian were the likely cause of the
failure.  More generally, the decision is made to reevaluate :math:`J`
(or instruct the user to reevaluate Jacobian data in :math:`P`) when:

* starting the problem,
* more than 50 steps have been taken since the last evaluation,
* a convergence failure occurred with an outdated matrix, and the
  value :math:`\bar{\gamma}` of :math:`\gamma` at the last update
  satisfies :math:`\left|\gamma/\bar{\gamma} - 1\right| > 0.2`,
* a convergence failure occurred that forced a step size reduction.



The stopping test for the Newton iteration is related to the
subsequent local error test, with the goal of keeping the nonlinear
iteration errors from interfering with local error control.  As
described below, the final computed value of each stage solution
:math:`z_i^{(m)}` will have to satisfy a local error test
:math:`\|z_i^{(m)} - z_i^{(0)}\| \le \epsilon`.  Letting
:math:`z_i` denote the true solution to the nonlinear problem
:eq:`Residual`, we want to ensure that the iteration error
:math:`z_i - z_i^{(m)}` is small relative to :math:`\epsilon`,
specifically that it is less than :math:`0.2\epsilon` (the safety
factor 0.2 may be changed by the user via the
:c:func:`ARKodeSetNonlinConvCoef()` function).  For this, we also
estimate the linear convergence rate :math:`R_i` of the modified Newton
iteration as follows.  We first initialize :math:`R_i` to 1, and reset
:math:`R_i=1` when either :math:`A` or :math:`P` are updated.  After
computing a Newton correction :math:`\delta^{(m)} = z_i^{(m)} -
z_i^{(m-1)}`, we update :math:`R_i` if :math:`m>1` as

.. math:: 
   R_i \leftarrow \max\{ 0.3 R_i, \left\|\delta^{(m)}\right\| / \left\|\delta^{(m-1)}\right\| \}.

where the factor 0.3 is user-modifiable as the ``crdown`` input to the
the function :c:func:`ARKodeSetNewtonConstants()`.  Denoting the
combined time step solution from the true stage solutions :math:`z_i`
as :math:`y_n`, and the combined time step solution from the computed
stage solutions :math:`z_i^{(m)}` as :math:`\tilde{y}_n` we use the
estimate 

.. math::
   \left\| y_n - \tilde{y}_n \right\| \approx 
   \max_i \left\| z_i^{(m+1)} - z_i^{(m)} \right\| \approx
   \max_i R_i \left\| z_i^{(m)} - z_i^{(m-1)} \right\| =
   \max_i R_i \left\| \delta^{(m)} \right\|.

Therefore the convergence (stopping) test for the modified Newton
iteration for each stage is

.. math::
   R_i \left\|\delta^{(m)} \right\| < 0.2\epsilon.

We allow at most 3 Newton iterations (this may be modified through the
function :c:func:`ARKodeSetMaxNonlinIters()`).  We also declare the
Newton iteration to be divergent if any of the ratios
:math:`\|\delta^{(m)}\| / \|\delta^{(m-1)}\| > 2.3` with :math:`m>1`
(the value 2.3 may be modified as the ``rdiv`` input to the function 
:c:func:`ARKodeSetNewtonConstants()`).  If convergence fails with
:math:`J` or :math:`A` current, we must then reduce the step size by a
factor of 0.25 (modifiable via the ``etacf`` input to the
:c:func:`ARKodeSetAdaptivityConstants()` function).  The integration
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
   J_{i,j}(t,y) = \frac{f_{I,i}(t,y+\sigma_j e_j) - f_{I,i}(t,y)}{\sigma_j}.

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

In the following four sections (:ref:`Mathematics.Preconditioning`,
:ref:`Mathematics.Predictors`, :ref:`Mathematics.Adaptivity` and
:ref:`Mathematics.Stability`), we provide details on optional
user-supplied information that can be used to better control the
behavior of ARKode.  In these sections, we also discuss the algorithms
currently provided by ARKode.  Finally, in the final section of this
chapter, :ref:`Mathematics.Rootfinding`, we discuss the algorithms
providing root-finding capabilities within ARKode.




.. _Mathematics.Preconditioning:

Preconditioning
------------------

When using a Newton method to solve the nonlinear system
:eq:`Residual`, ARKode makes repeated use of a linear solver to solve
linear systems of the form :math:`Ax = b`, where :math:`x` is a
correction vector and :math:`b` is a residual vector.  If this linear
system solve is done with one of the scaled preconditioned iterative
linear solvers, these solvers are rarely efficient if used without
preconditioning. A system :math:`Ax=b` can be preconditioned as one of:

.. math::
   (P^{-1}A)x = P^{-1}b & \qquad\text{[left preconditioning]}, \\
   (AP^{-1})Px = b  & \qquad\text{[right preconditioning]}, \\
   (P_L^{-1} A P_R^{-1}) P_R x = P_L^{-1} & \qquad\text{[left and right
   preconditioning]}.

The Krylov method is then applied to a system with the
matrix :math:`P^{-1}A`, :math:`AP^{-1}`, or :math:`P_L^{-1} A P_R^{-1}`,
instead of :math:`A`.  In order to improve the convergence of the
Krylov iteration, the preconditioner matrix :math:`P`, or the product
:math:`P_L P_R` in the third case, should in some sense approximate
the system matrix :math:`A`.  Yet at the same time, in order to be
cost-effective the matrix :math:`P` (or matrices :math:`P_L` and
:math:`P_R`) should be reasonably efficient to evaluate and
solve.  Finding an optimal point in this tradeoff between rapid
convergence and low cost can be quite challenging.  Good choices are
often problem-dependent (for example, see [BrownHindmarsh1989]_ for an
extensive study of preconditioners for reaction-transport systems). 

The ARKode solver allow for preconditioning either side, or on both
sides, although for non-symmetric matrices :math:`A` we know of few
situations where preconditioning on both sides is superior to
preconditioning on one side only (with the product :math:`P = P_L P_R`).
Moreover, for a given preconditioner matrix, the merits of left
vs. right preconditioning are unclear in general, and the user should
experiment with both choices.  Performance will differ between these
choices because the inverse of the left preconditioner is included in
the linear system residual whose norm is being tested in the Krylov
algorithm.  As a rule, however, if the preconditioner is the product
of two matrices, we recommend that preconditioning be done either on
the left only or the right only, rather than using one factor on each
side. 

Typical preconditioners used with ARKode are based on approximations
to the system Jacobian, :math:`J = \partial f_I / \partial y`.  Since
the Newton iteration matrix involved is :math:`A = M - \gamma J`, any
approximation :math:`\bar{J}` to :math:`J` yields a matrix that is of
potential use as a preconditioner, namely :math:`P = M - \gamma
\bar{J}`. Because the Krylov iteration occurs within a Newton
iteration and further also within a time integration, and since each
of these iterations has its own test for convergence, the
preconditioner may use a very crude approximation, as long as it
captures the dominant numerical feature(s) of the system.  We have
found that the combination of a preconditioner with the Newton-Krylov
iteration, using even a relatively poor approximation to the Jacobian,
can be surprisingly superior to using the same matrix without Krylov
acceleration (i.e., a modified Newton iteration), as well as to using
the Newton-Krylov method with no preconditioning.



.. _Mathematics.Predictors:

Implicit predictors
----------------------

As mentioned in the previous section, :ref:`Mathematics.IVP`, for
problems with implicit components, ARKode will employ a prediction
algorithm for constructing the initial guesses for each Runge-Kutta
stage, :math:`z_i^{(0)}`.  As is well-known with Newton-like methods,
the selection of an initial guess can have dramatic effects on both
the speed and robustness of the nonlinear solve, enabling the
difference between divergence and quadratic convergence of the
iteration.  To this end, ARKode implements a variety of prediction
algorithms that may be selected by the user.  In each case, the stages
guesses :math:`z_i^{(0)}` are constructed explicitly using
readily-available information, including the previous step solutions
:math:`y_{n-1}` and :math:`y_{n-2}`, as well as any previous stage
solutions :math:`z_j, \quad j<i`.  In all cases, prediction is
performed through construction of an interpolating polynomial through
existing data, which is then evaluated at the subsequent stage times
to hopefully provide a reasonable prediction of the future solution
value.  Specifically, for all of the Runge-Kutta methods implemented
in ARKode (and the vast majority in general), each stage solution
satisfies

.. math::
   z_i \approx y(t_{n-1} + c_i h_n),

so by constructing an interpolating polynomial :math:`p_q(t)` through
a set of existing data, the initial guess at stage solutions may be
approximated as 

.. math::
   z_i^{(0)} = p_q(t_{n-1} + c_i h_n).

Denoting :math:`[a,b]` as the interval containing the data used to
construct :math:`p_q(t)`, it is typically the case that 
:math:`t_{n-1} + c_i h_n > b`.  The dangers of using a polynomial
interpolant to extrapolate values outside the interpolation interval
are well-known, with higher-order polynomials and predictions further
outside the interval causing the greatest potential inaccuracies.  

Each prediction algorithm therefore constructs a different type of
interpolant :math:`p_q(t)`, as described below.



.. _Mathematics.Predictors.Trivial:

Trivial predictor
^^^^^^^^^^^^^^^^^^^^

The so-called "trivial predictor" is given by the formula

.. math::
   z_i^{(0)} = y_{n-1}.

While this piecewise-constant interpolant is clearly not an ideal
candidate for problems with time-varying solutions, it is often the
most robust approach for problems with constraints whose violation may
cause illegal solution values (e.g. a negative temperature).


.. _Mathematics.Predictors.Max:

Maximum order predictor
^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the opposite end of the spectrum, ARKode will construct an
interpolant :math:`p_q(t)` of polynomial order up to :math:`q=3`.
Here, the function :math:`p_q(t)` is identical to the one used for
interpolation of output solution values between time steps, i.e. for
:math:`y(t)` values where :math:`t_{n-1} < t < t_n`.  The order of
this polynomial, :math:`q`, may be specified by the user with the
function :c:func:`ARKodeSetDenseOrder()`.

The interpolants generated are either of Lagrange or Hermite form, and
use the data :math:`\left\{ y_{n-2}, f_{n-2}, y_{n-1}, f_{n-1}
\right\}`, where by :math:`f_{k}` we mean :math:`f(t_k,y_k)`.  Defining
a scaled and shifted "time" variables for the interval
:math:`[t_{n-2}, t_{n-1}]` as

.. math::
   \tau(t) = (t-t_n)/h_{n-1},

we may denote the predicted stage times in the time interval
:math:`[t_{n-1}, t_{n}]` as 

.. math::
   \tau_i = c_i \frac{h_n}{h_{n-1}}.

We then construct the interpolants :math:`p(t)` as follows:

* :math:`q=0`: this chooses the constant function

  .. math::
     p_0(\tau) = \frac{y_{n-2} + y_{n-1}}{2}.

* :math:`q=1`: this chooses the linear Lagrange interpolant

  .. math::
     p_1(\tau) = -\tau\, y_{n-2} + (1+\tau)\, y_{n-1}.

* :math:`q=2`: this chooses the quadratic Hermite interpolant

  .. math::
     p_2(\tau) =  \tau^2\,y_{n-2} + (1-\tau^2)\,y_{n-1} + h(\tau+\tau^2)\,f_{n-1}.

* :math:`q=3`: this chooses the cubic Hermite interpolant

  .. math::
     p_3(\tau) =  (3\tau^2 + 2\tau^3)\,y_{n-2} +
     (1-3\tau^2-2\tau^3)\,y_{n-1} + h(\tau^2+\tau^3)\,f_{n-2} +
     h(\tau+2\tau^2+\tau^3)\,f_{n-1}. 

These higher-order predictors may be useful when using lower-order
methods in which :math:`h_n` is not too large.  We further note that
although higher-order interpolants are possible, these are not
implemented due to the greater chance of error in predicting late
stage solutions.



.. _Mathematics.Predictors.Decreasing:

Variable order predictor
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This predictor attempts to use higher-order interpolations
:math:`p_q(t)` for predicting earlier stages in the subsequent time
interval, and lower-order interpolants for later stages.  It uses the
same formulas as described above, but chooses :math:`q` adaptively
based on the stage index :math:`i`, under the (rather tenuous)
assumption that the stage times are increasing, i.e. :math:`c_j < c_k`
for :math:`j<k`:

.. math::
   q = \max\{ q_{max} - i,\; 1 \}.



.. _Mathematics.Predictors.Cutoff:

Cutoff order predictor
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This predictor follows a similar idea as the previous algorithm, but
monitors the actual stage times to determine the polynomial
interpolant to use for prediction:

.. math::
   q = \begin{cases}
      q_{max}, & \text{if}\quad \tau < \tfrac12,\\
      1, & \text{otherwise}.
   \end{cases}





.. _Mathematics.Adaptivity:

Time step adaptivity
----------------------

A critical part of ARKode, making it an IVP "solver" rather than just
an integrator, is its adaptive control of local truncation error
(LTE).  At every step, the local error is estimated and required to
satisfy tolerance conditions.  If this local error test fails, then
the step is redone with a reduced step size.  All of the Runge-Kutta
methods implemented within ARKode admit an embedded solution
:math:`\tilde{y}_n`, as shown in equation :eq:`ARK`.  Generally, these
embedded solutions attain a slightly lower order of accuracy than the
computed solution :math:`y_n`.  Denoting these orders of accuracy as
:math:`p` and :math:`q`, where :math:`p` corresponds to the embedding
and :math:`q` corresponds to the method, for the majority of embedded
methods :math:`p = q-1`, but in all methods :math:`p<q`.  These values
of :math:`p` and :math:`q` correspond to the global order of accuracy
for the method, hence each admit local errors satisfying 

.. math::
   \| y_1 - y(t_1) \| = C h_1^{q+1} + \mathcal O(h_1^{q+2}), \\
   \| \tilde{y}_1 - y(t_1) \| = D h_1^{p+1} + \mathcal O(h_1^{p+2}),
   :label: AsymptoticErrors

where :math:`C` and :math:`D` are constants independent of :math:`h`,
and where we have assumed that :math:`y_0 = y(t_0)`.  Combining these
estimates, we have

.. math::
   \| y_1 - \tilde{y}_1 \| = \| y_1 - y(t_1) - \tilde{y}_1 + y(t_1) \| 
   \le \| y_1 - y(t_1) \| + \| \tilde{y}_1 - y(t_1) \| 
   \le D h_1^{p+1} + \mathcal O(h_1^{p+2}).

We therefore use this difference norm as an estimate for the local
error,

.. math::
   \text{LTE} = \beta \left(y_1 - \tilde{y}_1\right) = 
   \beta h_1 \sum_{i=0}^{s} \left(b_i - \tilde{b}_i\right) 
   \left(f_E(t_n + c_i h_n, z_i) + f_I(t_n + c_i h_n, z_i)\right).

Here, :math:`\beta>0` is an error *bias* to help account for the error
constant :math:`D`; the default value of this is :math:`\beta = 1.5`,
and may be modified by the user through the function
:c:func:`ARKodeSetAdaptivityMethod()`.  

With this LTE estimate, the local error test is simply :math:`\|LTE\|
< 1`, where we remind the user that this is actually the WRMS norm
defined in equation :eq:`WRMS_NORM` that includes the user-specified
relative and absolute tolerances.  If this error test passes, the step
is considered successful, and the LTE is subsequently used to estimate
the next step size, as will be described below in the section
:ref:`Mathematics.Adaptivity.ErrorControl`.  If the error test fails,
the step is rejected and a new step size :math:`h'` is then computed
using the error control algorithms described in
:ref:`Mathematics.Adaptivity.ErrorControl`.  A new attempt at the step
is made, and the error test is repeated.  If it fails multiple times
(as specified through the `small_nef` input to
:c:func:`ARKodeSetAdaptivityConstants()`, which defaults to 2), then
:math:`h'/h` is limited above to 0.3 (this is modifiable via the
etamxf` argument to :c:func:`ARKodeSetAdaptivityConstants()`), and
limited below to 0.1 after an additional step failure.  After
seven error test failures (modifiable via the function
:c:func:`ARKodeSetMaxErrTestFails()`), ARKode returns to the user
with a give-up message.

We define the step size ratio between a prospective step :math:`h'`
and a completed step :math:`h` as :math:`\eta`, i.e.

.. math::
   \eta = h' / h.

This is bounded above by :math:`\eta_{max}` to ensure that step size
adjustments are not overly aggressive.  This value is modified
according to the step and history,

.. math::
   \eta_{max} = \begin{cases}
     \text{etamx1}, & \quad\text{on the first step (default is 10000)}, \\
     \text{growth}, & \quad\text{on general steps (default is 20)}, \\
     1, & \quad\text{if the previous step had an error test failure}.
   \end{cases}

Here, the values of *etamx1* and *growth* may be modified by the user
in the functions :c:func:`ARKodeSetAdaptivityConstants()` and
:c:func:`ARKodeSetAdaptivityMethod()`, respectively.

For some problems it may be preferrable to avoid small step size
adjustments.  This can be especially true for problems that construct
and factor the Newton Jacobian matrix :math:`A` from equation
:eq:`NewtonMatrix`, where this construction is computationally
expensive, and where Newton convergence can be seriously hindered
through use of a somewhat incorrect :math:`A`.  In these scenarios,
the step is not changed when :math:`\eta \in [\eta_L, \eta_U]`.  The
default values for these parameters are :math:`\eta_L = 1` and
:math:`\eta_U = 1.5`, though these are modifiable through the function
:c:func:`ARKodeSetAdaptivityMethod()`.

The user may supply external bounds on the step sizes within ARKode,
through defining the values :math:`h_{min}` and :math:`h_{max}` with
the functions :c:func:`ARKodeSetMinStep()` and
:c:func:`ARKodeSetMaxStep()`, respectively.  These default to
:math:`h_{min}=0` and :math:`h_{max}=\infty`.  

Normally, ARKode takes steps until a user-defined output value
:math:`t = t_{out}` is overtaken, and then it computes
:math:`y(t_{out})` by interpolation (using the same dense output
routines described in the section
:ref:`Mathematics.Predictors.Max`). However, a "one step" mode option 
is available, where control returns to the calling program after each
step. There are also options to force ARKode not to integrate past a
given stopping point :math:`t = t_{stop}`, through the function
:c:func:`ARKodeSetStopTime()`.  



.. _Mathematics.Adaptivity.ErrorControl:

Asymptotic error control
^^^^^^^^^^^^^^^^^^^^^^^^^^^

As mentioned above, ARKode adapts the step size in order to attain
local errors within desired tolerances of the true solution.  These
adaptivity algorithms estimate the prospective step size :math:`h'`
based on the asymptotic local error estimates :eq:`AsymptoticErrors`.
We define the values :math:`\varepsilon_n`, :math:`\varepsilon_{n-1}`
and :math:`\varepsilon_{n-2}` as

.. math::
   \varepsilon_k &\ \equiv \ \|\text{LTE}_k\| 
      \ = \ \beta \|y_n - \tilde{y}_n\|,

corresponding to the local error estimates for three consecutive
steps, :math:`t_{n-3} \to t_{n-2} \to t_{n-1} \to t_n`.  With these
estimates, ARKode implements a variety of error control algorithms, as
specified in the subsections below. 


.. _Mathematics.Adaptivity.ErrorControl.PID:

PID controller
""""""""""""""""""

This is the default time adaptivity controller used by ARKode.  It
derives from those found in [KennedyCarpenter2003]_, [Soderlind1998]_,
[Soderlind2003]_ and  [Soderlind2006]_.  It uses all three of the
local error estimates :math:`\varepsilon_n`, :math:`\varepsilon_{n-1}`
and :math:`\varepsilon_{n-2}` in determination of a prospective step
size,

.. math::
   h' = h_n \varepsilon_n^{-k_1/p} \varepsilon_{n-1}^{k_2/p} 
        \varepsilon_{n-2}^{-k_3/p},

where the constants :math:`k_1`, :math:`k_2` and :math:`k_3` default
to 0.58, 0.21 and 0.1, respectively, though each may be changed via a
call to the function :c:func:`ARKodeSetAdaptivityMethod()`.  In this
estimate, a floor of :math:`\varepsilon > 10^{-10}` is enforced to
avoid division-by-zero errors.  These local error history values are
all initialized to 1.0 upon program initialization, to accomodate the
few initial time steps of a calculation where some of these error
estimates are undefined.



.. _Mathematics.Adaptivity.ErrorControl.PI:

PI controller
""""""""""""""""""

Like with the previous method, the PI controller derives from those
found in [KennedyCarpenter2003]_, [Soderlind1998]_,  [Soderlind2003]_
and  [Soderlind2006]_, but it differs in that it only uses the two
most recent step sizes in its adaptivity algorithm,

.. math::
   h' = h_n \varepsilon_n^{-k_1/p} \varepsilon_{n-1}^{k_2/p}.

Here, the default values of :math:`k_1` and :math:`k_2` default
to 0.8 and 0.31, respectively, though they may be changed via a
call to the function :c:func:`ARKodeSetAdaptivityMethod()`.  As with
the previous controller, at initialization :math:`k_1 = k_2 = 1.0` and
the floor of :math:`10^{-10}` is enforced on the local error
estimates.  



.. _Mathematics.Adaptivity.ErrorControl.I:

I controller
""""""""""""""""""

The so-called I controller is the standard time adaptivity control
algorithm in use by most available ODE solvers.  It bases the
prospective time step estimate entirely off of the current local error
estimate, 

.. math::
   h' = h_n \varepsilon_n^{-k_1/p}.

By default, :math:`k_1=1`, but that may be overridden by the user with
the function :c:func:`ARKodeSetAdaptivityMethod()`.



.. _Mathematics.Adaptivity.ErrorControl.eGus:

Explicit Gustafsson controller
""""""""""""""""""""""""""""""""""

This step adaptivity algorithm was proposed in [Gustafsson1991]_, and
is primarily useful in combination with explicit Runge-Kutta methods.
Using the notation of our earlier controllers, it has the form

.. math::
   h' = \begin{cases}
      h_1 \varepsilon_1^{-1/p}, &\quad\text{on the first step}, \\
      h_n \varepsilon_n^{-k_1/p} 
        \left(\varepsilon_n/\varepsilon_{n-1}\right)^{k_2/p}, &
      \quad\text{on subsequent steps}.
   \end{cases}
   :label: expGus

The default values of :math:`k_1` and :math:`k_2` are 0.4 and 0.33,
respectively, which may be changed via the function
:c:func:`ARKodeSetAdaptivityMethod()`.


.. _Mathematics.Adaptivity.ErrorControl.iGus:

Implicit Gustafsson controller
""""""""""""""""""""""""""""""""""

A version of the above controller suitable for implicit Runge-Kutta
methods was introduced in [Gustafsson1994]_, and has the form

.. math::
   h' = \begin{cases}
      h_1 \varepsilon_1^{-1/p}, &\quad\text{on the first step}, \\
      h_n \left(h_n / h_{n-1}\right) \varepsilon_n^{-k_1/p} 
        \left(\varepsilon_n/\varepsilon_{n-1}\right)^{-k_2/p}, &
      \quad\text{on subsequent steps}.
   \end{cases}
   :label: impGus

The algorithm parameters default to :math:`k_1 = 0.98` and 
:math:`k_2 = 0.95`, but may be modified by the user with
:c:func:`ARKodeSetAdaptivityMethod()`. 


.. _Mathematics.Adaptivity.ErrorControl.ieGus:

IMEX Gustafsson controller
"""""""""""""""""""""""""""""

An IMEX version of these two preceding controllers is available in
ARKode.  This approach computes the estimates :math:`h'_1` arising from
equation :eq:`expGus` and the estimate :math:`h'_2` arising from
equation :eq:`impGus`, and selects :math:`h' = \min\left\{h'_1,
h'_2\right\}`.  Here, equation :eq:`expGus` uses :math:`k_1` and
:math:`k_2` with default values of 0.4 and 0.25, while equation
:eq:`impGus` sets both parameters to the input :math:`k_3` that
defaults to 0.95.  All three of these parameters may be modified with
the function :c:func:`ARKodeSetAdaptivityMethod()`. 


.. _Mathematics.Adaptivity.ErrorControl.User:

User-supplied controller
""""""""""""""""""""""""""""

Finally, ARKode allows the user to define their own time step
adaptivity function,

.. math::
   h' = f(y, t, h_n, \varepsilon_n, \varepsilon_{n-1}, \varepsilon_{n-2}, q, p),

via a call to :c:func:`ARKodeSetAdaptivityFn()`.


.. _Mathematics.Stability:

Explicit stability
----------------------

For problems that involve an explicit component in :math:`f_E(t,y)`,
explicit and additive Runge-Kutta methods may benefit from addition
user-supplied information regarding the explicit stability region.
All of the methods in ARKode utilize step adaptivity based on
estimates of the local error.  It is often the case that such local
error control will automatically adapt the steps such that the method
remains stable (since unstable steps will typically exceed the error
control tolerances).  However, for problems in which :math:`f_E(t,y)`
includes some moderately stiff components, and especially for
higher-order integration methods, it is quite likely that a
significant number of attempted steps will exceed the error
tolerances.  In these scenarios, a stability-based time step
controller may also be useful.

Since the explicit stability region for any method is highly
problem-dependent, as it results from the eigenvalues of the
linearized operator :math:`\frac{\partial f_E}{\partial y}`,
information on the maximum stable step size is not computed internally
within ARKode.  However, for many applications such information is
readily available.  For example, in an advection-diffusion calculation
:math:`f_I` may contain the stiff diffusive components and
:math:`f_E` may contain the comparably nonstiff advection terms.  In
this scenario, an explicitly stable step :math:`h_{exp}` would be
predicted as one satisfying the Courant-Friedrichs-Lewy (CFL)
stability condition,

.. math::
   h_{exp} < \frac{\Delta x}{\lambda}

where :math:`\Delta x` is the spatial mesh size and :math:`\lambda` is
the fastest advective wave speed.

In the case that a user has supplied a routine to predict these
explicitly stable step sizes, :math:`h_{exp}` is compared against that
resulting from the local error adaptivity, :math:`h_{acc}`, and the
step used by ARKode will satisfy 

.. math::
   h = \min\{c\, h_{exp},\, h_{acc}\},

where the explicit stability step factor :math:`c` may be modified
throug the function :c:func:`ARKodeSetAdaptivityMethod()`, and has a
default value of :math:`1/2`.



.. _Mathematics.Rootfinding:

Rootfinding
--------------

The ARKode solver has been augmented to include a rootfinding
feature. This means that, while integrating the IVP :eq:`IVP`, ARKode
can also find the roots of a set of user-defined functions
:math:`g_i(t,y)` that depend on :math:`t` and the solution vector
:math:`y = y(t)`. The number of these root functions is arbitrary, and
if more than one :math:`g_i` is found to have a root in any given
interval, the various root locations are found and reported in the
order that they occur on the :math:`t` axis, in the direction of
integration. 

Generally, this rootfinding feature finds only roots of odd
multiplicity, corresponding to changes in sign of :math:`g_i(t,
y(t))`, denoted :math:`g_i(t)` for short. If a user root function has
a root of even multiplicity (no sign change), it will probably be
missed by ARKode. If such a root is desired, the user should
reformulate the root function so that it changes sign at the desired
root. 

The basic scheme used is to check for sign changes of any
:math:`g_i(t)` over each time step taken, and then (when a sign change
is found) to home in on the root (or roots) with a modified secant
method [HiebertShampine1980]_.  In addition, each time :math:`g` is
computed, ARKode checks to see if :math:`g_i(t) = 0` exactly, and if
so it reports this as a root. However, if an exact zero of any
:math:`g_i` is found at a point :math:`t`, ARKode computes
:math:`g(t+\delta)` for a small increment :math:`\delta`, slightly
further in the direction of integration, and if any
:math:`g_i(t+\delta) = 0` also, ARKode stops and reports an
error. This way, each time ARKode takes a time step, it is guaranteed
that the values of all :math:`g_i` are nonzero at some past value of
:math:`t`, beyond which a search for roots is to be done. 

At any given time in the course of the time-stepping, after suitable
checking and adjusting has been done, ARKode has an interval
:math:`(t_{lo}, t_{hi}]` in which roots of the :math:`g_i(t)` are to
be sought, such that :math:`t_{hi}` is further ahead in the direction
of integration, and all :math:`g_i(t_{lo}) \ne 0`. The endpoint
:math:`t_{hi}` is either :math:`t_n`, the end of the time step last
taken, or the next requested output time :math:`t_{out}` if this comes 
sooner. The endpoint :math:`t_{lo}` is either :math:`t_{n-1}`, or the
last output time :math:`t_{out}` (if this occurred within the last
step), or the last root location (if a root was just located within
this step), possibly adjusted slightly toward :math:`t_n` if an exact 
zero was found. The algorithm checks :math:`g(t_{hi})` for zeros, and
it checks for sign changes in :math:`(t_{lo}, t_{hi})`. If no sign
changes are found, then either a root is reported (if some
:math:`g_i(t_{hi}) = 0`) or we proceed to the next time interval
(starting at :math:`t_{hi}`). If one or more sign changes were found,
then a loop is entered to locate the root to within a rather tight
tolerance, given by 

.. math::
   \tau = 100\, U\, (|t_n| + |h|)\qquad (\text{where}\; U = \text{unit roundoff}).

Whenever sign changes are seen in two or more root functions, the one
deemed most likely to have its root occur first is the one with the
largest value of 
:math:`\left|g_i(t_{hi})\right| / \left| g_i(t_{hi}) - g_i(t_{lo})\right|`, 
corresponding to the closest to :math:`t_{lo}` of the secant method
values. At each pass through the loop, a new value :math:`t_{mid}` is
set, strictly within the search interval, and the values of
:math:`g_i(t_{mid})` are checked. Then either :math:`t_{lo}` or
:math:`t_{hi}` is reset to :math:`t_{mid}` according to which
subinterval is found to have the sign change. If there is none in
:math:`(t_{lo}, t_{mid})` but some :math:`g_i(t_{mid}) = 0`, then that
root is reported. The loop continues until :math:`\left|t_{hi} -
t_{lo} \right| < \tau`, and then the reported root location is
:math:`t_{hi}`.  In the loop to locate the root of :math:`g_i(t)`, the
formula for :math:`t_{mid}` is 

.. math::
   t_{mid} = t_{hi} - 
   \frac{g_i(t_{hi}) (t_{hi} - t_{lo})}{g_i(t_{hi}) - \alpha g_i(t_{lo})} ,

where :math:`\alpha` is a weight parameter. On the first two passes
through the loop, :math:`\alpha` is set to 1, making :math:`t_{mid}`
the secant method value. Thereafter, :math:`\alpha` is reset according
to the side of the subinterval (low vs high, i.e. toward
:math:`t_{lo}` vs toward :math:`t_{hi}`) in which the sign change was
found in the previous two passes. If the two sides were opposite,
:math:`\alpha` is set to 1. If the two sides were the same, :math:`\alpha` 
is halved (if on the low side) or doubled (if on the high side). The
value of :math:`t_{mid}` is closer to :math:`t_{lo}` when
:math:`\alpha < 1` and closer to :math:`t_{hi}` when :math:`\alpha > 1`. 
If the above value of :math:`t_{mid}` is within :math:`\tau /2` of
:math:`t_{lo}` or :math:`t_{hi}`, it is adjusted inward, such that its
fractional distance from the endpoint (relative to the interval size)
is between 0.1 and 0.5 (with 0.5 being the midpoint), and the actual
distance from the endpoint is at least :math:`\tau/2`. 







