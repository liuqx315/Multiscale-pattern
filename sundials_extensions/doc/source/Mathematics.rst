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
   z_i^{(m+1)} = z_i^{(m)} + s^{(m)},
   :label: Newton_iteration

where :math:`m` is the Newton iteration index.  Here, the 
update :math:`s^{(m)}` in turn requires the solution of linear 
:index:`Newton systems`

.. math::
   A\left(z_i^{(m)}\right) s^{(m)} = -G\left(z_i^{(m)}\right), 
   :label: Newton_system

where

.. math::
   A \approx M - \gamma J, \quad J = \frac{\partial f_I}{\partial y},
   \quad\text{and}\quad \gamma = h_n A^I_{i,i}.

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
* more than MSBP steps have been taken since the last update,
* the value :math:`\bar{\gamma}` of :math:`\gamma` at the last update
  satisfies :math:`|\gamma/\bar{\gamma} - 1| >` DGMAX,
* a non-fatal convergence failure just occurred, or
* an error test failure just occurred,

where the parameters :index:`MSBP` and :index:`DGMAX` are described
further in the section :ref:`CInterface.OptionalInputs`.  When an
update is forced due to a convergence failure, an update of :math:`A`
or :math:`P` may or may not involve a reevaluation of :math:`J` (in
:math:`A`) or of Jacobian data (in :math:`P`), depending on whether
errors in the Jacobian were the likely cause of the failure.  More
generally, the decision is made to reevaluate :math:`J` (or instruct
the user to reevaluate Jacobian data in :math:`P`) when:

* starting the problem,
* more than 50 steps have been taken since the last evaluation,
* a convergence failure occurred with an outdated matrix, and the
  value :math:`\bar{\gamma}` of :math:`\gamma` at the last update
  satisfies :math:`|\gamma/\bar{\gamma} - 1| > 0.2`,
* a convergence failure occurred that forced a step size reduction.



[continue with discussion of the Newton stopping criteria, akin to
page 7 of the CVODE manual]




.. _Mathematics.Preconditioning:

Preconditioning
------------------




.. _Mathematics.Stability:

Explicit stability
----------------------




.. _Mathematics.Adaptivity:

Time step adaptivity
----------------------




.. _Mathematics.Predictors:

Implicit predictors
----------------------




.. _Rootfinding:

Rootfinding
--------------








