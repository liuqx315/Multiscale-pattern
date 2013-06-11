:tocdepth: 3



.. _Mathematics.Preconditioning:

Preconditioning
===================

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
often problem-dependent (for example, see [BH1989]_ for an
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

