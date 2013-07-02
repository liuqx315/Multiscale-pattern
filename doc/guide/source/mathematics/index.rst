:tocdepth: 3

.. _Mathematics:

===========================
Mathematical Considerations
===========================

ARKode solves ODE initial value problems (IVPs) in real :math:`N`
-space, which we write in the abstract form

.. math::
   M\dot{y} = f_E(t,y) + f_I(t,y), \qquad y(t_0) = y_0.
   :label: IVP

Here, :math:`t` is the independent variable (e.g. time), and the
dependent variables are given by :math:`y \in \Re^N`.  We use the
notation :math:`\dot{y}` to denote :math:`dy/dt`.  :math:`M` is a
user-specified nonsingular linear operator from :math:`\Re^N \to
\Re^N`.  For standard systems of ordinary differential equations and
for problems arising from the spatial semi-discretization of partial
differential equations using finite difference or finite volume
methods, :math:`M` is typically the identity matrix :math:`I`;  for
PDEs using a finite-element spatial semi-discretization :math:`M` is
typically a well-conditioned mass matrix. The two right-hand side
functions may be described as: 

* :math:`f_E(t,y)` contains the "slow" time scale components of the
  system, that should be integrated using explicit methods.

* :math:`f_I(t,y)` contains the "fast" time scale components of the
  system, that should be integrated using implicit methods.

ARKode may be used to solve stiff, nonstiff and multi-rate problems.
Roughly speaking, stiffness is characterized by the presence of at
least one rapidly damped mode, whose time constant is small compared
to the time scale of the solution itself.  In the implicit/explicit
(IMEX) splitting above, these stiff components should be included in
the right-hand side function :math:`f_I(t,y)`.

In the sub-sections listed below, we elaborate on the numerical
methods that comprise the ARKode solvers:

.. toctree::
   :maxdepth: 1

   IVP_solution
   Preconditioning
   Predictors
   Adaptivity
   Stability
   Mass_matrix
   Rootfinding


