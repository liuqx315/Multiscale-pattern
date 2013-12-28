..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3

.. _Mathematics:

===========================
Mathematical Considerations
===========================

ARKode solves ODE initial value problems (IVPs) in :math:`\mathbb{R}^N`.
These problems should be posed in explicit form, as

.. math::
   M\dot{y} = f_E(t,y) + f_I(t,y), \qquad y(t_0) = y_0.
   :label: IVP

Here, :math:`t` is the independent variable (e.g. time), and the
dependent variables are given by :math:`y \in \mathbb{R}^N`, where we
use the notation :math:`\dot{y}` to denote :math:`\frac{dy}{dt}`.

:math:`M` is a user-specified nonsingular operator from
:math:`\mathbb{R}^N \to \mathbb{R}^N`.  This operator may depend on
:math:`t` but is currently assumed to be independent of :math:`y`.
For standard systems of ordinary differential equations and for
problems arising from the spatial semi-discretization of partial
differential equations using finite difference or finite volume
methods, :math:`M` is typically the identity matrix, :math:`I`.  For
PDEs using a finite-element spatial semi-discretization :math:`M` is
typically a well-conditioned mass matrix.  

The two right-hand side functions may be described as:  

* :math:`f_E(t,y)` contains the "slow" time scale components of the
  system.  This will be integrated using explicit methods.

* :math:`f_I(t,y)` contains the "fast" time scale components of the
  system.  This will be integrated using implicit methods.

ARKode may be used to solve stiff, nonstiff and multi-rate problems.
Roughly speaking, stiffness is characterized by the presence of at
least one rapidly damped mode, whose time constant is small compared
to the time scale of the solution itself.  In the implicit/explicit
(ImEx) splitting above, these stiff components should be included in
the right-hand side function :math:`f_I(t,y)`.

In the sub-sections that follow, we elaborate on the numerical
methods that comprise the ARKode solvers.  We first discuss the
general :ref:`formulation of additive Runge-Kutta methods
<Mathematics.ARK>`, including the resulting implicit systems that must
be solved at each stage.  We then discuss the solver strategies that
ARKode uses in solving these systems: :ref:`nonlinear solvers
<Mathematics.Nonlinear>`, :ref:`linear solvers <Mathematics.Linear>`
and :ref:`preconditioners <Mathematics.Preconditioning>`.  We then
describe our approaches for :ref:`error control <Mathematics.Error>`
within the iterative nonlinear and linear solvers, including
discussion on our choice of norms used within ARKode for measuring
errors within various components of the solver.  We then discuss
specific enhancements available in ARKode, including an array of
:ref:`prediction algorithms <Mathematics.Predictors>` for the solution
at each stage, :ref:`adaptive error controllers
<Mathematics.Adaptivity>`, :ref:`mass-matrix handling
<Mathematics.MassSolve>`, and :ref:`rootfinding capabilities
<Mathematics.Rootfinding>`.





.. _Mathematics.ARK:

Additive Runge-Kutta methods
===============================

The methods used in ARKode are variable-step, embedded, 
:index:`additive Runge-Kutta methods` (ARK), based on formulas of the
form 

.. math::
   M z_i &= M y_{n-1} + h_n \sum_{j=0}^{i-1} A^E_{i,j} f_E(t_{n,j}, z_j) 
                 + h_n \sum_{j=0}^{i}   A^I_{i,j} f_I(t_{n,j}, z_j),
   \quad i=1,\ldots,s, \\
   M y_n &= M y_{n-1} + h_n \sum_{i=0}^{s} b_i \left(f_E(t_{n,i}, z_i) 
                 + f_I(t_{n,i}, z_i)\right), \\
   M \tilde{y}_n &= M y_{n-1} + h_n \sum_{i=0}^{s} \tilde{b}_i 
       \left(f_E(t_{n,i}, z_i) + f_I(t_{n,i}, z_i)\right).
   :label: ARK

Here the :math:`y_n` are computed approximations to :math:`y(t_n)`,
:math:`\tilde{y}_n` are lower-order embedded solutions (used in error
estimation), and :math:`h_n \equiv t_n - t_{n-1}` is the step size.
The internal stage times are abbreviated using the notation
:math:`t_{n,j} = t_{n-1} + c_j h_n`.  The ARK method is primarily
defined through the coefficients :math:`A^E \in \mathbb{R}^{s\times s}`, 
:math:`A^I \in \mathbb{R}^{s\times s}`, :math:`b \in \mathbb{R}^{s}` and 
:math:`c \in \mathbb{R}^{s}`, that correspond with the explicit and
implicit Butcher tables.  We note that ARKode enforces the constraint
that these tables must share :math:`b` and :math:`c` between the
explicit and implicit methods in an ARK pair.  

The user of ARKode must choose appropriately between one of three
classes of methods: *multi-rate*, *nonstiff* and *stiff*.  All of
ARKode's available Butcher tables encoding the coefficients :math:`c`,
:math:`A^E`, :math:`A^I`, :math:`b` and :math:`\tilde{b}` are further
described in the :ref:`Butcher`.  

For multi-rate problems, a user should provide both of the functions
:math:`f_E` and :math:`f_I` that define the IVP system.  For such
problems, ARKode currently implements the ARK methods proposed in
[KC2003]_, allowing for methods having order :math:`q = \{3,4,5\}`.
The tables for these methods are given in the section
:ref:`Butcher.additive`.

For nonstiff problems, a user may specify that :math:`f_I = 0`,
i.e. the equation :eq:`IVP` reduces to the non-split IVP 

.. math::
   M\dot{y} = f_E(t,y), \qquad y(t_0) = y_0.
   :label: IVP_explicit

In this scenario, the Butcher table :math:`A^I=0` in :eq:`ARK`, and
the ARK methods reduce to classical :index:`explicit Runge-Kutta
methods` (ERK).  For these classes of methods, ARKode allows orders of
accuracy :math:`q = \{2,3,4,5,6\}`, with embeddings of orders :math:`p
= \{1,2,3,4,5\}`.  These default to the :ref:`Butcher.Heun_Euler`,
:ref:`Butcher.Bogacki_Shampine`, :ref:`Butcher.Zonneveld`,
:ref:`Butcher.Cash-Karp` and :ref:`Butcher.Verner-6-5` methods,
respectively. 

Finally, for stiff problems the user may specify that :math:`f_E = 0`,
so the equation :eq:`IVP` reduces to the non-split IVP 

.. math::
   M\dot{y} = f_I(t,y), \qquad y(t_0) = y_0.
   :label: IVP_implicit

Similarly to ERK methods, in this scenario the Butcher table
:math:`A^E=0` in :eq:`ARK`, and the ARK methods reduce to classical
:index:`diagonally-implicit Runge-Kutta methods` (DIRK).  For these
classes of methods, ARKode allows orders of accuracy :math:`q =
\{2,3,4,5\}`, with embeddings of orders :math:`p = \{1,2,3,4\}`.
These default to the :ref:`Butcher.SDIRK-2-1`,
:ref:`Butcher.ARK_4_2_3_I`, :ref:`Butcher.SDIRK-5-4` and 
:ref:`Butcher.ARK_8_4_5_I` methods, respectively. 




.. _Mathematics.Nonlinear:

Nonlinear solver methods
===============================


For both the DIRK and ARK methods corresponding to :eq:`IVP` and
:eq:`IVP_implicit`, a nonlinear system

.. math::
   G(z_i) \equiv M z_i - h_n A^I_{i,i} f_I(t_{n,i}, z_i) - a_i = 0
   :label: Residual

must be solved for each stage :math:`z_i, i=1,\ldots,s`, where we have
the data 

.. math::
   a_i \equiv M y_{n-1} + h_n \sum_{j=0}^{i-1} \left[
      A^E_{i,j} f_E(t_{n,j}, z_j) +
      A^I_{i,j} f_I(t_{n,j}, z_j) \right]
   
for the ARK methods, or 

.. math::
   a_i \equiv M y_{n-1} + h_n \sum_{j=0}^{i-1} 
      A^I_{i,j} f_I(t_{n,j}, z_j)
   
for the DIRK methods.  For these nonlinear systems, ARKode allows a
choice of solution strategy. 

The default solver choice is a variant of :index:`Newton's method`,

.. math::
   z_i^{(m+1)} = z_i^{(m)} + \delta^{(m+1)},
   :label: Newton_iteration

where :math:`m` is the Newton iteration index, and the :index:`Newton
update` :math:`\delta^{(m+1)}` in turn requires the solution of the
linear :index:`Newton system` 

.. math::
   {\mathcal A}\left(z_i^{(m)}\right) \delta^{(m+1)} = -G\left(z_i^{(m)}\right), 
   :label: Newton_system

in which

.. math::
   {\mathcal A} \approx M - \gamma J, \quad J = \frac{\partial f_I}{\partial y},
   \quad\text{and}\quad \gamma = h_n A^I_{i,i}.
   :label: NewtonMatrix

As an alternate to Newton's method, ARKode may solve for each stage
:math:`z_i, i=1,\ldots,s` using an :index:`Anderson-accelerated fixed
point iteration`

.. math::
   z_i^{(m+1)} = g(z_i^{(m)}), \quad m=0,1,\ldots
   :label: AAFP_iteration

Unlike with Newton's method, this method *does not* require the
solution of a linear system at each iteration, instead opting for
solution of a low-dimensional least-squares solution to construct the
nonlinear update.  For details on how this iteration is performed, we
refer the reader to the reference [WN2011]_.

The optimal solver (Newton vs fixed-point) is highly
problem-dependent.  Since fixed-point solvers do not require
the solution of any linear systems, each iteration may be
significantly less costly than their Newton counterparts.  However,
this can come at the cost of slower convergence (or even divergence)
in comparison with Newton-like methods.  However, these fixed-point
solvers do allow for user specification of the Anderson-accelerated
subspace size, :math:`m_k`.  While the required amount of solver
memory grows proportionately to :math:`m_k N`, larger values of
:math:`m_k` may result in faster convergence.  In our experience, this
improvement may be significant even for "small" values,
e.g. :math:`1\le m_k\le 5`, and that convergence may not improve (or
even deteriorate) for larger values of :math:`m_k`.

While ARKode uses the Newton iteration as its default solver due to
its increased robustness on very stiff problems, it is highly
recommended that users also consider the fixed-point solver for their
when attempting a new problem.




For either the Newton or fixed-point solvers, it is well-known that
both the efficiency and robustness of the algorithm intimately depends
on the choice of a good initial guess.  In ARKode, the initial guess
for either nonlinear solution method is a predicted value
:math:`z_i^{(0)}` that is computed explicitly from the
previously-computed data (e.g. :math:`y_{n-2}`, :math:`y_{n-1}`, and
:math:`z_j` where :math:`j<i`).  Additional information on the
specific predictor algorithms implemented in ARKode is provided in the
following section, :ref:`Mathematics.Predictors`.



.. _Mathematics.Linear:

Linear solver methods
===============================

When a Newton-based method is chosen for solving each nonlinear
system, a linear system of equations must be solved at each nonlinear
iteration.  For this solve ARKode provides several choices, including
the option of a user-supplied linear solver module.  The linear solver
modules distributed with ARKode are organized into two families: a
*direct* family comprising direct linear solvers for dense or banded
matrices, and a *spils* family comprising scaled, preconditioned,
iterative (Krylov) linear solvers.  The methods offered through these
modules are as follows:

* dense direct solvers, using either an internal SUNDIALS
  implementation or a BLAS/LAPACK implementation (serial version
  only), 
* band direct solvers, using either an internal SUNDIALS
  implementation or a BLAS/LAPACK implementation (serial version
  only), 
* SPGMR, a scaled, preconditioned GMRES (Generalized Minimal Residual)
  solver without restarts, 
* SPBCG, a scaled, preconditioned Bi-CGStab (Bi-Conjugate Gradient
  Stable) solver,
* SPTFQMR, a scaled, preconditioned TFQMR (Transpose-free
  Quasi-Minimal Residual) solver,
* SPFGMR, a scaled, preconditioned Flexible GMRES (Generalized Minimal
  Residual) solver without restarts, or
* PCG, a preconditioned conjugate gradient solver for symmetric linear
  systems.

For large stiff systems where direct methods are infeasible, the
combination of an implicit integrator and a preconditioned
Krylov method (SPGMR, SPBCG, SPTFQMR, SPFGMR or PCG) can yield a
powerful tool because it combines established methods for stiff
integration, nonlinear solver iteration, and Krylov (linear) iteration
with a problem-specific treatment of the dominant sources of
stiffness, in the form of a user-supplied preconditioner matrix
[BH1989]_.  We note that the direct linear solvers
provided by SUNDIALS (dense and band), as well as the direct linear
solvers accessible through LAPACK, can only be used with the serial
vector representations.


.. index:: modified Newton iteration

In the case that a direct linear solver is used (dense or band),
ARKode utilizes a *modified Newton iteration*. In such methods, the
matrix :math:`{\mathcal A}` is held fixed for multiple Newton
iterations.  More precisely, each Newton iteration is computed from
the modified equation 

.. math::
   \tilde{\mathcal A}\left(z_i^{(m)}\right) \delta^{(m+1)} = -G\left(z_i^{(m)}\right), 
   :label: modified_Newton_system

in which

.. math::
   \tilde{\mathcal A} \approx M - \tilde{\gamma} \tilde{J}, \quad \tilde{J} =
   \frac{\partial f_I}{\partial y}(\tilde y), \quad\text{and}\quad
   \tilde{\gamma} = \tilde{h} A^I_{i,i}. 
   :label: modified_NewtonMatrix

Here, the solution :math:`\tilde{y}` and step size :math:`\tilde{h}`
upon which the modified Jacobian rely, are merely values of the
solution and step size from a previous iteration.  In other words, the
matrix :math:`\tilde{\mathcal A}` is only computed rarely, and reused for
repeated stage solves.  

When using the direct and band solvers for the linear systems
:eq:`modified_Newton_system`, the Jacobian may be supplied by a user
routine or approximated by finite-differences.  In the case of
differencing, we use the standard approximation

.. math::
   J_{i,j}(t,y) = \frac{f_{I,i}(t,y+\sigma_j e_j) - f_{I,i}(t,y)}{\sigma_j},

where :math:`e_j` is the jth unit vector, and the increments
:math:`\sigma_j` are given by 

.. math::
   \sigma_j = \max\left\{ \sqrt{U}\, |y_j|, \frac{\sigma_0}{w_j} \right\}.

Here :math:`U` is the unit roundoff, :math:`\sigma_0` is a
dimensionless value, and :math:`w_j` is the error weight defined in
:eq:`EWT`.  In the dense case, this approach requires :math:`N`
evaluations of :math:`f_I`, one for each column of :math:`J`.  In the
band case, the columns of :math:`J` are computed in groups, using the
Curtis-Powell-Reid algorithm, with the number of :math:`f_I`
evaluations equal to the bandwidth.




.. index:: inexact Newton iteration

In the case that an iterative linear solver is chosen, ARKode utilizes a
Newton method variant called an *Inexact Newton iteration*.  Here, the
matrix :math:`{\mathcal A}` is not itself constructed since the
algorithms only require the product of this matrix with a given
vector.  Additionally, each Newton system :eq:`Newton_system` is not
solved completely, since these linear solvers are iterative (hence the
"inexact" in the name). Resultingly. for these linear solvers
:math:`{\mathcal A}` is applied in a matrix-free manner, 

.. math::

   {\mathcal A}v = Mv - \gamma Jv.

The matrix-vector products :math:`Jv` are obtained by either calling
an optional user-supplied routine, or through directional differencing
using the formula 

.. math::
   Jv = \frac{f_I(t,y+\sigma v) - f_I(t,y)}{\sigma},

where the increment :math:`\sigma = 1/\|v\|` to ensure that 
:math:`\|\sigma v\| = 1`.


As with the modified Newton method that reused :math:`{\mathcal A}` between solves,
ARKode's inexact Newton iteration also recomputes the preconditioner
matrix :math:`P` as infrequently as possible to balance the high costs
of matrix construction and factorization against the reduced
convergence rate that may result from a stale preconditioner. 

More specifically, in both of the Newton-based solvers, we update the
Newton matrix :math:`\tilde{\mathcal A}` or preconditioner matrix :math:`P`
only in the following circumstances: 

* when starting the problem,
* when more than 20 steps have been taken since the last update (this
  value may be changed via the *msbp* argument to
  :c:func:`ARKodeSetMaxStepsBetweenLSet()`) or the *LSETUP_MSBP*
  argument to :f:func:`FARKSETIIN()`, 
* when the value :math:`\bar{\gamma}` of :math:`\gamma` at the last
  update satisfies :math:`\left|\gamma/\bar{\gamma} - 1\right| > 0.2`
  (this tolerance may be changed via the *dgmax* argument to 
  :c:func:`ARKodeSetDeltaGammaMax()`) or the *LSETUP_DGMAX*
  argument to :f:func:`FARKSETRIN()`, 
* when a non-fatal convergence failure just occurred, or
* when an error test failure just occurred.

When an update is forced due to a convergence failure, an update of
:math:`\tilde{\mathcal A}` or :math:`P` may or may not involve a reevaluation of
:math:`J` (in :math:`\tilde{\mathcal A}`) or of Jacobian data (in :math:`P`),
depending on whether errors in the Jacobian were the likely cause of the
failure.  More generally, the decision is made to reevaluate :math:`J`
(or instruct the user to reevaluate Jacobian data in :math:`P`) when:

* starting the problem,
* more than 50 steps have been taken since the last evaluation,
* a convergence failure occurred with an outdated matrix, and the
  value :math:`\bar{\gamma}` of :math:`\gamma` at the last update
  satisfies :math:`\left|\gamma/\bar{\gamma} - 1\right| > 0.2`,
* a convergence failure occurred that forced a step size reduction.


As will be further discussed in the section
:ref:`Mathematics.Preconditioning`, in the case of a Krylov method, 
preconditioning may be applied on the left, right, or on both sides of
:math:`{\mathcal A}`, with user-supplied routines for the preconditioner setup
and solve operations.  




.. _Mathematics.Error:

Iteration Error Control
========================


.. _Mathematics.Error.Norm:

Choice of norm
----------------

In the process of controlling errors at various levels (time
integration, nonlinear solution, linear solution), ARKode uses a
:index:`weighted root-mean-square norm`, denoted
:math:`\|\cdot\|_\text{WRMS}`, for all error-like quantities,

.. math::
   \|v\|_\text{WRMS} = \left( \frac{1}{N} \sum_{i=1}^N \left(v_i\,
   w_i\right)^2\right)^{1/2}. 
   :label: WRMS_NORM

The power of this choice of norm arises in the specification of the
weighting vector :math:`w`, that combines the units of the problem
with the user-supplied measure of "acceptable" error.  To this end,
ARKode constructs and :index:`error weight vector` using the
most-recent step solution and the relative and absolute tolerances
input by the user, namely

.. math::
   w_i = \frac{1}{RTOL\cdot |y_i| + ATOL_i}.
   :label: EWT

Since :math:`1/w_i` represents a tolerance in the component
:math:`y_i`, a vector whose WRMS norm is 1 is regarded as "small." For
brevity, we will typically drop the subscript WRMS on norms in the
remainder of this section. 

Additionally, for problems involving a non-identity mass matrix,
:math:`M\ne I`, the units of equation :eq:`IVP` may differ from the
units of the solution :math:`y`.  In this case, ARKode may also
construct a :index:`residual weight vector`,

.. math::
   w_i = \frac{1}{RTOL\cdot |My_i| + ATOL'_i},
   :label: RWT

where the user may specify a separate absolute residual tolerance
value or array, :math:`ATOL'_i`.  The choice of weighting vector used
in any given norm is determined by the quantity being measured: values
having solution units use :eq:`EWT`, whereas values having equation
units use :eq:`RWT`.  Obviously, for problems with :math:`M=I`, the
weighting vectors are identical. 




.. _Mathematics.Error.Nonlinear:

Nonlinear iteration error control
-----------------------------------

The stopping test for all of ARKode's nonlinear solvers is related to
the subsequent local error test, with the goal of keeping the
nonlinear iteration errors from interfering with local error control.
Denoting the final computed value of each stage solution as
:math:`z_i^{(m)}`, and the true stage solution solving :eq:`Residual`
as :math:`z_i`, we want to ensure that the iteration error
:math:`z_i - z_i^{(m)}` is "small" (recall that a norm less than 1 is
already considered "small").

To this end, we first estimate the linear convergence rate :math:`R_i`
of the nonlinear iteration.  We initialize :math:`R_i=1`, and reset it 
to this value whenever :math:`\tilde{\mathcal A}` or :math:`P` are
updated.  After computing a nonlinear correction :math:`\delta^{(m)} =
z_i^{(m)} - z_i^{(m-1)}`, if :math:`m>1` we update :math:`R_i` as

.. math:: 
   R_i \leftarrow \max\{ 0.3 R_i, \left\|\delta^{(m)}\right\| / \left\|\delta^{(m-1)}\right\| \}.

where the factor 0.3 is user-modifiable as the *crdown* input to the
the function :c:func:`ARKodeSetNonlinCRDown()` or the *NEWT_CRDOWN*
argument to :f:func:`FARKSETRIN()`.  

Denoting the true time step solution as :math:`y_n`, and the computed
time step solution (computed using the stage solutions
:math:`z_i^{(m)}`) as :math:`\tilde{y}_n`, we use the estimate 

.. math::
   \left\| y_n - \tilde{y}_n \right\| \approx 
   \max_i \left\| z_i^{(m+1)} - z_i^{(m)} \right\| \approx
   \max_i R_i \left\| z_i^{(m)} - z_i^{(m-1)} \right\| =
   \max_i R_i \left\| \delta^{(m)} \right\|.

Therefore our convergence (stopping) test for the nonlinear iteration
for each stage is 

.. math::
   R_i \left\|\delta^{(m)} \right\| < \epsilon,

where the factor :math:`\epsilon` has default value 0.1, and is
user-modifiable as the *nlscoef* input to the the function
:c:func:`ARKodeSetNonlinConvCoef()` or the *NLCONV_COEF* input to the
function :f:func:`FARKSETRIN()`.  We allow at most 3 nonlinear
iterations (modifiable through :c:func:`ARKodeSetMaxNonlinIters()`, or
as the *MAX_NSTEPS* argument to :f:func:`FARKSETIIN()`).  We also
declare the nonlinear iteration to be divergent if any of the ratios
:math:`\|\delta^{(m)}\| / \|\delta^{(m-1)}\| > 2.3` with :math:`m>1`
(the value 2.3 may be modified as the *rdiv* input to 
:c:func:`ARKodeSetNonlinRDiv()` or the *NEWT_RDIV* input to
:f:func:`FARKSETRIN()`).  If convergence fails in the fixed 
point iteration, or in the Newton iteration with :math:`J` or
:math:`{\mathcal A}` current, we must then reduce the step size by a
factor of 0.25 (modifiable via the *etacf* input to the
:c:func:`ARKodeSetMaxCFailGrowth()` function or the *ADAPT_ETACF*
input to :f:func:`FARKSETRIN()`).  The integration is halted after 10
convergence failures (modifiable via the
:c:func:`ARKodeSetMaxConvFails()` function or the *MAX_CONVFAIL*
argument to :f:func:`FARKSETIIN()`).



.. _Mathematics.Error.Linear:

Linear iteration error control
-----------------------------------

When a Krylov method is used to solve the linear systems
:eq:`Newton_system`, its errors must also be controlled.  To this end,
we approximate the linear iteration error in the solution vector
:math:`\delta^{(m)}` using the preconditioned residual vector,
e.g. :math:`r = P{\mathcal A}\delta^{(m)} + PG` for the case of left
preconditioning (the role of the preconditioner is further elaborated
on in the next section).  In an attempt to ensure that the linear
iteration errors do not interfere with the nonlinear solution error
and local time integration error controls, we require that the norm of
the preconditioned linear residual satisfies

.. math::
  
   \|r\| \le 0.05\epsilon.

Here :math:`\epsilon` is the same value as that used above for the
nonlinear error control.  The value 0.05 may be modified by the user
through the :c:func:`ARKSpilsSetEpsLin()` function; it cannot
currently be modified from Fortran applications.




.. _Mathematics.Preconditioning:

Preconditioning
===================

When using an inexact Newton method to solve the nonlinear system
:eq:`Residual`, ARKode makes repeated use of a linear solver to solve
linear systems of the form :math:`{\mathcal A}x = b`, where :math:`x` is a
correction vector and :math:`b` is a residual vector.  If this linear
system solve is done with one of the scaled preconditioned iterative
linear solvers, the efficiency of such solvers may benefit
tremendously from preconditioning. A system :math:`{\mathcal A}x=b` can be
preconditioned as one of: 

.. math::
   (P^{-1}{\mathcal A})x = P^{-1}b & \qquad\text{[left preconditioning]}, \\
   ({\mathcal A}P^{-1})Px = b  & \qquad\text{[right preconditioning]}, \\
   (P_L^{-1} {\mathcal A} P_R^{-1}) P_R x = P_L^{-1}b & \qquad\text{[left and right
   preconditioning]}.

The Krylov method is then applied to a system with the
matrix :math:`P^{-1}{\mathcal A}`, :math:`{\mathcal A}P^{-1}`, or
:math:`P_L^{-1} {\mathcal A} P_R^{-1}`, instead of :math:`{\mathcal
A}`.  In order to improve the convergence of the Krylov iteration, the
preconditioner matrix :math:`P`, or the product :math:`P_L P_R` in the
third case, should in some sense approximate the system matrix
:math:`{\mathcal A}`.  Yet at the same time, in order to be
cost-effective the matrix :math:`P` (or matrices :math:`P_L` and
:math:`P_R`) should be reasonably efficient to evaluate and solve.
Finding an optimal point in this tradeoff between rapid 
convergence and low cost can be quite challenging.  Good choices are
often problem-dependent (for example, see [BH1989]_ for an
extensive study of preconditioners for reaction-transport systems). 

The ARKode solver allows for preconditioning either side, or on both
sides, although for non-symmetric matrices :math:`{\mathcal A}` we
know of few situations where preconditioning on both sides is superior
to preconditioning on one side only (with the product :math:`P = P_L
P_R`).  Moreover, for a given preconditioner matrix, the merits of left 
vs. right preconditioning are unclear in general, and the user should
experiment with both choices.  Performance will differ between these
choices because the inverse of the left preconditioner is included in
the linear system residual whose norm is being tested in the Krylov
algorithm.  As a rule, however, if the preconditioner is the product
of two matrices, we recommend that preconditioning be done either on
the left only or the right only, rather than using one factor on each
side.  An exception to this rule is the PCG solver, that itself
assumes a symmetric matrix :math:`{\mathcal A}`, since the PCG
algorithm in fact applies the single preconditioner matrix :math:`P`
in both left/right fashion as :math:`P^{-1/2} {\mathcal A} P^{-1/2}`.

Typical preconditioners used with ARKode are based on approximations
to the system Jacobian, :math:`J = \partial f_I / \partial y`.  Since
the Newton iteration matrix involved is :math:`{\mathcal A} = M - \gamma J`, any
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
===========================

For problems with implicit components, ARKode will employ a prediction 
algorithm for constructing the initial guesses for each Runge-Kutta
stage, :math:`z_i^{(0)}`.  As is well-known with nonlinear solvers,
the selection of a good initial guess can have dramatic effects on both
the speed and robustness of the nonlinear solve, enabling the
difference between rapid quadratic convergence versus divergence of
the iteration.  To this end, ARKode implements a variety of prediction
algorithms that may be selected by the user.  In each case, the stage
guesses :math:`z_i^{(0)}` are constructed explicitly using
readily-available information, including the previous step solutions
:math:`y_{n-1}` and :math:`y_{n-2}`, as well as any previous stage
solutions :math:`z_j, \quad j<i`.  In all cases, prediction is
performed by constructing an interpolating polynomial through
existing data, which is then evaluated at the subsequent stage times
to provide an inexpensive but (hopefully) reasonable prediction of the
subsequent solution value.  Specifically, for all of the Runge-Kutta
methods implemented in ARKode (and the vast majority in general), each
stage solution satisfies

.. math::
   z_i \approx y(t_{n,i}),

so by constructing an interpolating polynomial :math:`p_q(t)` through
a set of existing data, the initial guess at stage solutions may be
approximated as 

.. math::
   z_i^{(0)} = p_q(t_{n,i}).

Denoting :math:`[a,b]` as the interval containing the data used to
construct :math:`p_q(t)`, and assuming forward integration from
:math:`a\to b`, it is typically the case that :math:`t_{n,j} > b`.
The dangers of using a polynomial interpolant to extrapolate values
outside the interpolation interval are well-known, with higher-order
polynomials and predictions further outside the interval resulting in
the greatest potential inaccuracies.

Each prediction algorithm therefore constructs a different type of
interpolant :math:`p_q(t)`, as described below.



.. _Mathematics.Predictors.Trivial:

Trivial predictor
--------------------

The so-called "trivial predictor" is given by the formula

.. math::

   p_0(\tau) = y_{n-1}.

While this piecewise-constant interpolant is clearly not a highly
accurate candidate for problems with time-varying solutions, it is
often the most robust approach for either highly stiff problems, or
problems with implicit constraints whose violation may cause illegal
solution values (e.g. a negative density or temperature). 


.. _Mathematics.Predictors.Max:

Maximum order predictor
---------------------------

At the opposite end of the spectrum, ARKode can construct an
interpolant :math:`p_q(t)` of polynomial order up to :math:`q=3`.
Here, the function :math:`p_q(t)` is identical to the one used for
interpolation of output solution values between time steps, i.e. for
":index:`dense output`" of :math:`y(t)` for :math:`t_{n-1} < t < t_n`.
The order of this polynomial, :math:`q`, may be specified by the user
with the function :c:func:`ARKodeSetDenseOrder()` or with the
*DENSE_ORDER* argument to :f:func:`FARKSETIIN()`.

The interpolants generated are either of Lagrange or Hermite form, and
use the data :math:`\left\{ y_{n-2}, f_{n-2}, y_{n-1}, f_{n-1}
\right\}`, where we use :math:`f_{k}` to denote :math:`M^{-1}
\left(f_E(t_k,y_k) + f_I(t_k,y_k)\right)`.  Defining a scaled and
shifted "time" variable :math:`\tau` for the interval :math:`[t_{n-2},
t_{n-1}]` as 

.. math::

   \tau(t) = (t-t_n)/h_{n-1},

we may denote the predicted stage times in the subsequent time
interval :math:`[t_{n-1}, t_{n}]` as 

.. math::

   \tau_i = c_i \frac{h_n}{h_{n-1}}.

We then construct the interpolants :math:`p(t)` as follows:

* :math:`q=0`: constant interpolant

  .. math::

     p_0(\tau) = \frac{y_{n-2} + y_{n-1}}{2}.

* :math:`q=1`: linear Lagrange interpolant

  .. math::

     p_1(\tau) = -\tau\, y_{n-2} + (1+\tau)\, y_{n-1}.

* :math:`q=2`: quadratic Hermite interpolant

  .. math::

     p_2(\tau) =  \tau^2\,y_{n-2} + (1-\tau^2)\,y_{n-1} + h(\tau+\tau^2)\,f_{n-1}.

* :math:`q=3`: cubic Hermite interpolant

  .. math::

     p_3(\tau) =  (3\tau^2 + 2\tau^3)\,y_{n-2} +
     (1-3\tau^2-2\tau^3)\,y_{n-1} + h(\tau^2+\tau^3)\,f_{n-2} +
     h(\tau+2\tau^2+\tau^3)\,f_{n-1}. 

These higher-order predictors may be useful when using lower-order
methods in which :math:`h_n` is not too large.  We further note that
although interpolants of order :math:`> 3` are possible, these are not
implemented due to their increased computing and storage costs, along
with their diminishing returns due to increased extrapolation error.



.. _Mathematics.Predictors.Decreasing:

Variable order predictor
---------------------------

This predictor attempts to use higher-order interpolations
:math:`p_q(t)` for predicting earlier stages in the subsequent time
interval, and lower-order interpolants for later stages.  It uses the
same formulas as described above, but chooses :math:`q` adaptively
based on the stage index :math:`i`, under the (rather tenuous)
assumption that the stage times are increasing, i.e. :math:`c_j < c_k`
for :math:`j<k`:

.. math::
   q = \max\{ q_\text{max} - i,\; 1 \}.



.. _Mathematics.Predictors.Cutoff:

Cutoff order predictor
---------------------------

This predictor follows a similar idea as the previous algorithm, but
monitors the actual stage times to determine the polynomial
interpolant to use for prediction:

.. math::
   q = \begin{cases}
      q_\text{max}, & \text{if}\quad \tau < \tfrac12,\\
      1, & \text{otherwise}.
   \end{cases}



.. _Mathematics.Predictors.Bootstrap:

Bootstrap predictor
---------------------------

This predictor does not use any information from within the preceding
step, instead using information only within the current step
:math:`[t_{n-1},t_n]` (including :math:`y_{n-1}` and
:math:`f_{n-1}`).  Instead, this approach uses the right-hand side
from a previously computed stage solution in the same step,
:math:`f(t_{n-1}+c_j h,z_j)` to construct a quadratic Hermite
interpolant for the prediction.  If we define the constants
:math:`\tilde{h} = c_j h` and :math:`\tau = c_i h`, the predictor is
given by 

.. math::
 
   z_i^{(0)} = y_{n-1} + \left(\tau - \frac{\tau^2}{2\tilde{h}}\right)
      f(t_{n-1},y_{n-1}) + \frac{\tau^2}{2\tilde{h}} f(t_{n-1}+c_j h,z_j).

For stages in which :math:`c_j=0` for all previous stages 
:math:`j = 0,\ldots,i-1`, and for the first stage of any time step
:math:`(i=0)`, this method reduces to using the trivial predictor  
:math:`z_i^{(0)} = y_{n-1}`.  For stages having multiple precdding
nonzero :math:`c_j`, we choose the stage having largest :math:`c_j`
value, to minimize the amount of extrapolation induced through the
prediction.





.. _Mathematics.Adaptivity:

Time step adaptivity
=======================

A critical component of ARKode, making it an IVP "solver" rather than
just an integrator, is its adaptive control of local truncation error.
At every step, we estimate the local error, and ensure that it
satisfies tolerance conditions.  If this local error test fails, then
the step is recomputed with a reduced step size.  To this end, every
Runge-Kutta method packaged within ARKode admit an embedded solution
:math:`\tilde{y}_n`, as shown in equation :eq:`ARK`. Generally, these
embedded solutions attain a lower order of accuracy than the computed
solution :math:`y_n`.  Denoting these orders of accuracy as :math:`p`
and :math:`q`, where :math:`p` corresponds to the embedding and
:math:`q` corresponds to the method, for the majority of embedded 
methods :math:`p = q-1`.  These values of :math:`p` and :math:`q`
correspond to the global order of accuracy for the method and
embedding, hence each admit local errors satisfying [HW1993]_

.. math::
   \| y_n - y(t_n) \| = C h_n^{q+1} + \mathcal O(h_n^{q+2}), \\
   \| \tilde{y}_n - y(t_n) \| = D h_n^{p+1} + \mathcal O(h_n^{p+2}),
   :label: AsymptoticErrors

where :math:`C` and :math:`D` are constants independent of :math:`h`,
and where we have assumed exact initial conditions for the step,
:math:`y_{n-1} = y(t_{n-1})`. Combining these estimates, we have

.. math::
   \| y_n - \tilde{y}_n \| = \| y_n - y(t_n) - \tilde{y}_n + y(t_n) \| 
   \le \| y_n - y(t_n) \| + \| \tilde{y}_n - y(t_n) \| 
   \le D h_n^{p+1} + \mathcal O(h_n^{p+2}).

We therefore use this difference norm as an estimate for the local
truncation error at the step :math:`n`,

.. math::
   T_n = \beta \left(y_n - \tilde{y}_n\right) = 
   \beta h_n M^{-1} \sum_{i=0}^{s} \left(b_i - \tilde{b}_i\right) 
   \left(f_E(t_{n-1} + c_i h_n, z_i) + f_I(t_{n-1} + c_i h_n, z_i)\right).
   :label: LTE

Here, :math:`\beta>0` is an error *bias* to help account for the error
constant :math:`D`; the default value of this is :math:`\beta = 1.5`,
and may be modified by the user through the function
:c:func:`ARKodeSetErrorBias()` or through the input *ADAPT_BIAS* to
:f:func:`FARKSETRIN()`.

With this LTE estimate, the local error test is simply :math:`\|T_n\|
< 1`, where we remind that this norm includes the user-specified
relative and absolute tolerances.  If this error test passes, the step
is considered successful, and the estimate is subsequently used to
estimate the next step size, as will be described below in the section
:ref:`Mathematics.Adaptivity.ErrorControl`.  If the error test fails,
the step is rejected and a new step size :math:`h'` is then computed
using the error control algorithms described in
:ref:`Mathematics.Adaptivity.ErrorControl`.  A new attempt at the step
is made, and the error test is repeated.  If it fails multiple times
(as specified through the *small_nef* input to
:c:func:`ARKodeSetSmallNumEFails()` or the *ADAPT_SMALL_NEF* argument
to :f:func:`FARKSETIIN()`, which defaults to 2), then
:math:`h'/h` is limited above to 0.3 (this is modifiable via the
*etamxf* argument to :c:func:`ARKodeSetMaxEFailGrowth()` or the
*ADAPT_ETAMXF* argument to :f:func:`FARKSETRIN()`), and
limited below to 0.1 after an additional step failure.  After
seven error test failures (modifiable via the function
:c:func:`ARKodeSetMaxErrTestFails()` or the *MAX_ERRFAIL* argument to
:f:func:`FARKSETIIN()`), ARKode returns to the user with a give-up
message. 

We define the step size ratio between a prospective step :math:`h'`
and a completed step :math:`h` as :math:`\eta`, i.e.

.. math::
   \eta = h' / h.

This is bounded above by :math:`\eta_\text{max}` to ensure that step size
adjustments are not overly aggressive.  This value is modified
according to the step and history,

.. math::
   \eta_\text{max} = \begin{cases}
     \text{etamx1}, & \quad\text{on the first step (default is 10000)}, \\
     \text{growth}, & \quad\text{on general steps (default is 20)}, \\
     1, & \quad\text{if the previous step had an error test failure}.
   \end{cases}

Here, the values of *etamx1* and *growth* may be modified by the user
in the functions :c:func:`ARKodeSetMaxFirstGrowth()` and
:c:func:`ARKodeSetMaxGrowth()`, respectively, or through the inputs
*ADAPT_ETAMX1* and *ADAPT_GROWTH* to the function
:f:func:`FARKSETRIN()`. 

A flowchart detailing how the time steps are modified at each
iteration to ensure solver convergence and successful steps is given
in the figure below.  Here, all norms correspond to the WRMS norm, and
the error adaptivity function **arkAdapt** is supplied by one of the
error control algorithms discussed in the subsections below. 

.. _adaptivity_figure:

.. figure:: figs/time_adaptivity.png
   :scale: 40 %
   :align: center


For some problems it may be preferrable to avoid small step size
adjustments.  This can be especially true for problems that construct
and factor the Newton Jacobian matrix :math:`{\mathcal A}` from
equation :eq:`NewtonMatrix` for either a direct solve, or as a
preconditioner for an iterative solve, where this construction is
computationally expensive, and where Newton convergence can be
seriously hindered through use of a somewhat incorrect
:math:`{\mathcal A}`.  In these scenarios, the step is not changed
when :math:`\eta \in [\eta_L, \eta_U]`.  The default values for these
parameters are :math:`\eta_L = 1` and :math:`\eta_U = 1.5`, though
these are modifiable through the function
:c:func:`ARKodeSetFixedStepBounds()` or through the input
*ADAPT_BOUNDS* to the function :f:func:`FARKSETRIN()`.

The user may supply external bounds on the step sizes within ARKode,
through defining the values :math:`h_\text{min}` and :math:`h_\text{max}` with
the functions :c:func:`ARKodeSetMinStep()` and
:c:func:`ARKodeSetMaxStep()`, or through the inputs *MIN_STEP* and
*MAX_STEP* to the function :f:func:`FARKSETRIN()`, respectively. 
These default to :math:`h_\text{min}=0` and :math:`h_\text{max}=\infty`.  

Normally, ARKode takes steps until a user-defined output value
:math:`t = t_\text{out}` is overtaken, and then it computes
:math:`y(t_\text{out})` by interpolation (using the same dense output
routines described in the section
:ref:`Mathematics.Predictors.Max`). However, a "one step" mode option 
is available, where control returns to the calling program after each
step. There are also options to force ARKode not to integrate past a
given stopping point :math:`t = t_\text{stop}`, through the function
:c:func:`ARKodeSetStopTime()` or through the input *STOP_TIME* to
:f:func:`FARKSETRIN()`. 



.. _Mathematics.Adaptivity.ErrorControl:

Asymptotic error control
---------------------------

As mentioned above, ARKode adapts the step size in order to attain
local errors within desired tolerances of the true solution.  These
adaptivity algorithms estimate the prospective step size :math:`h'`
based on the asymptotic local error estimates :eq:`AsymptoticErrors`.
We define the values :math:`\varepsilon_n`, :math:`\varepsilon_{n-1}`
and :math:`\varepsilon_{n-2}` as

.. math::
   \varepsilon_k &\ \equiv \ \|T_k\| 
      \ = \ \beta \|y_n - \tilde{y}_n\|,

corresponding to the local error estimates for three consecutive
steps, :math:`t_{n-3} \to t_{n-2} \to t_{n-1} \to t_n`.  These local
error history values are all initialized to 1.0 upon program
initialization, to accomodate the few initial time steps of a
calculation where some of these error estimates are undefined.  With
these estimates, ARKode implements a variety of error control
algorithms, as specified in the subsections below.


.. _Mathematics.Adaptivity.ErrorControl.PID:

PID controller
^^^^^^^^^^^^^^^^^^

This is the default time adaptivity controller used by ARKode.  It
derives from those found in [KC2003]_, [S1998]_, [S2003]_ and
[S2006]_.  It uses all three of the local error estimates
:math:`\varepsilon_n`, :math:`\varepsilon_{n-1}` and
:math:`\varepsilon_{n-2}` in determination of a prospective step size,

.. math::
   h' \;=\; h_n\; \varepsilon_n^{-k_1/p}\; \varepsilon_{n-1}^{k_2/p}\; 
        \varepsilon_{n-2}^{-k_3/p},

where the constants :math:`k_1`, :math:`k_2` and :math:`k_3` default
to 0.58, 0.21 and 0.1, respectively, though each may be changed via a
call to the C/C++ function :c:func:`ARKodeSetAdaptivityMethod()`, or
to the Fortran function :f:func:`FARKSETADAPTIVITYMETHOD()`.  In this
estimate, a floor of :math:`\varepsilon > 10^{-10}` is enforced to
avoid division-by-zero errors.  



.. _Mathematics.Adaptivity.ErrorControl.PI:

PI controller
^^^^^^^^^^^^^^^^^

Like with the previous method, the PI controller derives from those
found in [KC2003]_, [S1998]_, [S2003]_ and [S2006]_, but it differs in
that it only uses the two most recent step sizes in its adaptivity
algorithm, 

.. math::
   h' \;=\; h_n\; \varepsilon_n^{-k_1/p}\; \varepsilon_{n-1}^{k_2/p}.

Here, the default values of :math:`k_1` and :math:`k_2` default
to 0.8 and 0.31, respectively, though they may be changed via a
call to :c:func:`ARKodeSetAdaptivityMethod()` or
:f:func:`FARKSETADAPTIVITYMETHOD()`.  As with the previous controller,
at initialization :math:`k_1 = k_2 = 1.0` and the floor of
:math:`10^{-10}` is enforced on the local error estimates.  



.. _Mathematics.Adaptivity.ErrorControl.I:

I controller
^^^^^^^^^^^^^^^^

The so-called I controller is the standard time adaptivity control
algorithm in use by most available ODE solvers.  It bases the
prospective time step estimate entirely off of the current local error
estimate, 

.. math::
   h' \;=\; h_n\; \varepsilon_n^{-k_1/p}.

By default, :math:`k_1=1`, but that may be overridden by the user with
the function :c:func:`ARKodeSetAdaptivityMethod()` or the function
:f:func:`FARKSETADAPTIVITYMETHOD()`. 




.. _Mathematics.Adaptivity.ErrorControl.eGus:

Explicit Gustafsson controller
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step adaptivity algorithm was proposed in [G1991]_, and
is primarily useful in combination with explicit Runge-Kutta methods.
Using the notation of our earlier controllers, it has the form

.. math::
   h' \;=\; \begin{cases}
      h_1\; \varepsilon_1^{-1/p}, &\quad\text{on the first step}, \\
      h_n\; \varepsilon_n^{-k_1/p}\; 
        \left(\varepsilon_n/\varepsilon_{n-1}\right)^{k_2/p}, &
      \quad\text{on subsequent steps}.
   \end{cases}
   :label: expGus

The default values of :math:`k_1` and :math:`k_2` are 0.367 and 0.268,
respectively, which may be changed bhy calling either 
:c:func:`ARKodeSetAdaptivityMethod()` or :f:func:`FARKSETADAPTIVITYMETHOD()`.




.. _Mathematics.Adaptivity.ErrorControl.iGus:

Implicit Gustafsson controller
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A version of the above controller suitable for implicit Runge-Kutta
methods was introduced in [G1994]_, and has the form

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
:c:func:`ARKodeSetAdaptivityMethod()` or :f:func:`FARKSETADAPTIVITYMETHOD()`. 




.. _Mathematics.Adaptivity.ErrorControl.ieGus:

ImEx Gustafsson controller
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An ImEx version of these two preceding controllers is available in
ARKode.  This approach computes the estimates :math:`h'_1` arising from
equation :eq:`expGus` and the estimate :math:`h'_2` arising from
equation :eq:`impGus`, and selects 

.. math::
   h' = \frac{h}{|h|}\min\left\{|h'_1|, |h'_2|\right\}.  

Here, equation :eq:`expGus` uses :math:`k_1` and
:math:`k_2` with default values of 0.367 and 0.268, while equation
:eq:`impGus` sets both parameters to the input :math:`k_3` that
defaults to 0.95.  All three of these parameters may be modified with
the C/C++ function :c:func:`ARKodeSetAdaptivityMethod()` or the
Fortran function :f:func:`FARKSETADAPTIVITYMETHOD()`. 



.. _Mathematics.Adaptivity.ErrorControl.User:

User-supplied controller
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, ARKode allows the user to define their own time step
adaptivity function,

.. math::
   h' = H(y, t, h_n, h_{n-1}, h_{n-2}, \varepsilon_n, \varepsilon_{n-1}, \varepsilon_{n-2}, q, p),

via a call to the C/C++ routine :c:func:`ARKodeSetAdaptivityFn()` or
the Fortran routine :f:func:`FARKADAPTSET()`.





.. _Mathematics.Stability:

Explicit stability
======================

For problems that involve a nonzero explicit component,
:math:`f_E(t,y) \ne 0`, explicit and ImEx Runge-Kutta methods may 
benefit from addition user-supplied information regarding the explicit
stability region.  All ARKode adaptivity methods utilize estimates of
the local error.  It is often the case that such local error control
will be sufficient for method stability, since unstable steps will
typically exceed the error control tolerances.  However, for problems
in which :math:`f_E(t,y)` includes even moderately stiff components,
and especially for higher-order integration methods, it may occur that
a significant number of attempted steps will exceed the error
tolerances.  While these steps will automatically be recomputed, such
trial-and-error may be costlier than desired.  In these scenarios, a
stability-based time step controller may also be useful.

Since the explicit stability region for any method depends on the
problem under consideration, as it results from the eigenvalues of the
linearized operator :math:`\frac{\partial f_E}{\partial y}`,
information on the maximum stable step size is not computed internally
within ARKode.  However, for many problems such information is
readily available.  For example, in an advection-diffusion calculation,
:math:`f_I` may contain the stiff diffusive components and
:math:`f_E` may contain the comparably nonstiff advection terms.  In
this scenario, an explicitly stable step :math:`h_\text{exp}` would be
predicted as one satisfying the Courant-Friedrichs-Lewy (CFL)
stability condition,

.. math::
   |h_\text{exp}| < \frac{\Delta x}{|\lambda|}

where :math:`\Delta x` is the spatial mesh size and :math:`\lambda` is
the fastest advective wave speed.

In these scenarios, a user may supply a routine to predict this
maximum explicitly stable step size, :math:`|h_\text{exp}|`, by calling the
C/C++ function :c:func:`ARKodeSetStabilityFn()` or the Fortran
function :f:func:`FARKEXPSTABSET()`.  If a value for
:math:`|h_\text{exp}|` is supplied, it is compared against the value
resulting from the local error controller, :math:`|h_\text{acc}|`, and
the step used by ARKode will satisfy  

.. math::
   h' = \frac{h}{|h|}\min\{c\, |h_\text{exp}|,\, |h_\text{acc}|\}.

Here the explicit stability step factor (often called the "CFL
factor") :math:`c>0` may be modified through the function
:c:func:`ARKodeSetCFLFraction()` or through the input *ADAPT_CFL* to
the function :f:func:`FARKSETRIN()`, and has a default value of
:math:`1/2`. 




.. _Mathematics.MassSolve:

Mass matrix solver
=======================

Within the algorithms described above, there are three locations where a
linear solve of the form

.. math::
   M x = b

is required: (a) in constructing the time-evolved solution
:math:`y_n`, (b) in estimating the local temporal truncation error,
and (c) in constructing predictors for the implicit solver iteration
(see section :ref:`Mathematics.Predictors.Max`).  Specifically, to
construct the time-evolved solution :math:`y_n` from equation
:eq:`ARK` we must solve

.. math::
   &M y_n \ = \ M y_{n-1} + h_n \sum_{i=0}^{s} b_i \left(f_E(t_{n,i}, z_i) 
                 + f_I(t_{n,i}, z_i)\right), \\
   \Leftrightarrow \qquad & \\
   &M (y_n -y_{n-1}) \ = \ h_n \sum_{i=0}^{s} b_i \left(f_E(t_{n,i}, z_i) 
                 + f_I(t_{n,i}, z_i)\right), \\
   \Leftrightarrow \qquad & \\
   &M \nu \ = \ h_n \sum_{i=0}^{s} b_i \left(f_E(t_{n,i}, z_i) 
                 + f_I(t_{n,i}, z_i)\right),

for the update :math:`\nu = y_n - y_{n-1}`.  Similarly, in computing
the local temporal error estimate :math:`T_n` from equation :eq:`LTE`
we must solve systems of the form 

.. math::
   M\, T_n = h \sum_{i=0}^{s} \left(b_i - \tilde{b}_i\right) 
   \left(f_E(t_{n,i}, z_i) + f_I(t_{n,i}, z_i)\right).

Lastly, in constructing dense output and implicit predictors of order
2 or higher (as in the section :ref:`Mathematics.Predictors.Max` above),
we must compute the derivative information :math:`f_k` from the equation 

.. math::
   M f_k = f_E(t_k, y_k) + f_I(t_k, y_k).

Of course, for problems in which :math:`M=I` these solves are not
required; however for problems with non-identity :math:`M`, ARKode may
use either an iterative linear solver or a dense linear solver, in the
same manner as described in the section :ref:`Mathematics.Linear` for solving
the linear Newton systems.  We note that at present, the matrix
:math:`M` may depend on time :math:`t` but must be independent of the
solution :math:`y`, since we assume that each of the above systems are
linear.

At present, for DIRK and ARK problems using a dense or band solver for
the Newton nonlinear iterations, the type of linear solver (dense or
band) for the Newton systems :math:`{\mathcal A}\delta = -G` must
match the type of linear solver used for these mass-matrix systems,
since :math:`M` is included inside :math:`{\mathcal A}`.  When direct
methods (dense and band) are employed, the user must supply a routine
to compute :math:`M` in either dense or band form to match the
structure of :math:`{\mathcal A}`, using either the routine
:c:func:`ARKDlsDenseMassFn()` or :c:func:`ARKDlsBandMassFn()`.  When
iterative methods are used, a routine must be supplied to perform the
mass-matrix-vector product, :math:`Mv`, through a call to the routine
:c:func:`ARKSpilsMassTimesVecFn()`.  As with iterative solvers for the
Newton systems, preconditioning may be applied to aid in solution of
the mass matrix systems :math:`Mx=b`.

We further note that non-identity mass matrices, :math:`M\ne I`, are
only supported by the C and C++ ARKode interfaces, although Fortran
support is planned for the near future.




.. _Mathematics.Rootfinding:

Rootfinding
===============

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
method [HS1980]_.  In addition, each time :math:`g` is
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
:math:`(t_\text{lo}, t_\text{hi}]` in which roots of the :math:`g_i(t)` are to
be sought, such that :math:`t_\text{hi}` is further ahead in the direction
of integration, and all :math:`g_i(t_\text{lo}) \ne 0`. The endpoint
:math:`t_\text{hi}` is either :math:`t_n`, the end of the time step last
taken, or the next requested output time :math:`t_\text{out}` if this comes 
sooner. The endpoint :math:`t_\text{lo}` is either :math:`t_{n-1}`, or the
last output time :math:`t_\text{out}` (if this occurred within the last
step), or the last root location (if a root was just located within
this step), possibly adjusted slightly toward :math:`t_n` if an exact 
zero was found. The algorithm checks :math:`g(t_\text{hi})` for zeros, and
it checks for sign changes in :math:`(t_\text{lo}, t_\text{hi})`. If no sign
changes are found, then either a root is reported (if some
:math:`g_i(t_\text{hi}) = 0`) or we proceed to the next time interval
(starting at :math:`t_\text{hi}`). If one or more sign changes were found,
then a loop is entered to locate the root to within a rather tight
tolerance, given by 

.. math::
   \tau = 100\, U\, (|t_n| + |h|)\qquad (\text{where}\; U = \text{unit roundoff}).

Whenever sign changes are seen in two or more root functions, the one
deemed most likely to have its root occur first is the one with the
largest value of 
:math:`\left|g_i(t_\text{hi})\right| / \left| g_i(t_\text{hi}) - g_i(t_\text{lo})\right|`, 
corresponding to the closest to :math:`t_\text{lo}` of the secant method
values. At each pass through the loop, a new value :math:`t_\text{mid}` is
set, strictly within the search interval, and the values of
:math:`g_i(t_\text{mid})` are checked. Then either :math:`t_\text{lo}` or
:math:`t_\text{hi}` is reset to :math:`t_\text{mid}` according to which
subinterval is found to have the sign change. If there is none in
:math:`(t_\text{lo}, t_\text{mid})` but some :math:`g_i(t_\text{mid}) = 0`, then that
root is reported. The loop continues until :math:`\left|t_\text{hi} -
t_\text{lo} \right| < \tau`, and then the reported root location is
:math:`t_\text{hi}`.  In the loop to locate the root of :math:`g_i(t)`, the
formula for :math:`t_\text{mid}` is 

.. math::
   t_\text{mid} = t_\text{hi} - 
   \frac{g_i(t_\text{hi}) (t_\text{hi} - t_\text{lo})}{g_i(t_\text{hi}) - \alpha g_i(t_\text{lo})} ,

where :math:`\alpha` is a weight parameter. On the first two passes
through the loop, :math:`\alpha` is set to 1, making :math:`t_\text{mid}`
the secant method value. Thereafter, :math:`\alpha` is reset according
to the side of the subinterval (low vs high, i.e. toward
:math:`t_\text{lo}` vs toward :math:`t_\text{hi}`) in which the sign change was
found in the previous two passes. If the two sides were opposite,
:math:`\alpha` is set to 1. If the two sides were the same, :math:`\alpha` 
is halved (if on the low side) or doubled (if on the high side). The
value of :math:`t_\text{mid}` is closer to :math:`t_\text{lo}` when
:math:`\alpha < 1` and closer to :math:`t_\text{hi}` when :math:`\alpha > 1`. 
If the above value of :math:`t_\text{mid}` is within :math:`\tau /2` of
:math:`t_\text{lo}` or :math:`t_\text{hi}`, it is adjusted inward, such that its
fractional distance from the endpoint (relative to the interval size)
is between 0.1 and 0.5 (with 0.5 being the midpoint), and the actual
distance from the endpoint is at least :math:`\tau/2`. 

Finally, we note that when running in parallel, the ARKode rootfinding
module assumes that the entire set of root defining functions
:math:`g_i(t,y)` is replicated on every MPI task.  Since in these
cases the vector :math:`y` is distributed across tasks, it is the
user's responsibility to perform any necessary inter-task
communication to ensure that :math:`g_i(t,y)` is identical on each task.
