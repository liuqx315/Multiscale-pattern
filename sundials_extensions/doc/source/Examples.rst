:tocdepth: 3

.. _Examples:

ARKode Examples
===============

ARKode comes packaged with a variety of example problems, that
exercise options including explicit, implicit and ImEx solvers,
root-finding, direct and iterative linear solvers, and the Fortran
solver interface, FARKODE.  While these examples are not an exhaustive
set of all possible usage scenarios, they are designed to show a
variety of usage scenarios, and can be used as templates for new
problems using ARKode's solvers.



Simple linear example (ark_analytic)
-------------------------------------

This is a very simple C example that merely shows how to use the
ARKode solver interface.

ODE system
^^^^^^^^^^^^

The problem is that of a scalar-valued initial value problem (IVP)
that is linear in the dependent variable :math:`y`, but nonlinear in
the independent variable :math:`t`:

.. math::

   \frac{dy}{dt} = \lambda y + \frac{1}{1+t^2} - \lambda \arctan(t),

where :math:`0\le t\le 10` and :math:`y(0)=0`.  The stiffness of the
problem may be tuned via the parameter :math:`\lambda`, which is
specified (along with the relative and absolute tolerances,
:math:`rtol` and :math:`atol`) in the input file
``input_analytic.txt``.  The value of :math:`\lambda` must be negative
to result in a well-posed problem; for values with magnitude larger
than 100 or so the problem becomes quite stiff.  In the provided input
file, we choose :math:`\lambda=-100` and tolerances
:math:`rtol=10^{-6}` and :math:`atol=10^{-10}`.    After each unit
time interval, the solution is output to the screen.


Numerical method
^^^^^^^^^^^^^^^^^

The example routine solves this problem using a diagonally-implicit
Runge-Kutta method.  Each stage is solved using the built-in modified
Newton iteration, but since the ODE is linear in :math:`y` these
should only require a single iteration per stage.  Internally, Newton
will use the ARKDENSE dense linear solver, which in the case of this
scalar-valued problem is just division.  The example file contains
functions to evaluate both :math:`f(t,y)` and :math:`J(t,y)=\lambda`.

Aside from the input tolerance values, this problem uses only the
default parameters for the ARKode solver.


Solutions
^^^^^^^^^^^^

This problem is included both as a simple example, but also because it
has an analytical solution, :math:`y(t) = \arctan(t)`.  As seen in the
plots below, the computed solution tracks the analytical solution
quite well, and results in errors below those specified by the input
error tolerances.



Simple nonlinear example (ark_analytic_nonlin)
-----------------------------------------------

ODE system
^^^^^^^^^^^^

Numerical method
^^^^^^^^^^^^^^^^^

Solutions
^^^^^^^^^^^^



Simple linear system example (ark_analytic_sys)
------------------------------------------------

ODE system
^^^^^^^^^^^^

Numerical method
^^^^^^^^^^^^^^^^^

Solutions
^^^^^^^^^^^^



Stiff nonlinear system example (ark_brusselator)
-------------------------------------------------

ODE system
^^^^^^^^^^^^

Numerical method
^^^^^^^^^^^^^^^^^

Solutions
^^^^^^^^^^^^



Stiff nonlinear system, Fortran example (ark_bruss)
----------------------------------------------------

This is a Fortran-90 version of the same test brusselator test problem
as above.  

ODE system
^^^^^^^^^^^^

The test problem has 3 dependent variables :math:`u`, :math:`v` and
:math:`w`, that depend on the independent variable :math:`t` via the
IVP system

.. math::

   \frac{du}{dt} &= a - (w+1)u + v u^2, \\
   \frac{dv}{dt} &= w u - v u^2, \\
   \frac{dw}{dt} &= \frac{b-w}{\varepsilon} - w u.

We integrate over the interval :math:`0 \le t \le 10`, with the
initial conditions :math:`u(0) = 3.9`, :math:`v(0) = 1.1`, :math:`w(0) = 2.8`,
and parameters :math:`a=1.2`, :math:`b=2.5` and
:math:`\varepsilon=10^{-5}`.  After each unit time interval, the
solution is output to the screen.


Numerical method
^^^^^^^^^^^^^^^^^

Since this driver and utility functions are written in Fortran-90,
this example demonstrates the use of the FARKODE interface for the
ARKode solver.  For time integration, this example uses the
fourth-order additive Runge-Kutta method, where the right-hand sides
are broken up as

.. math::

   f_E(t,u,v,w) = \left(\begin{array}{c} a - (w+1)u + v u^2 \\ 
     w u - v u^2 \\ - w u  \end{array}\right), \quad\text{and}\quad 
   f_I(t,u,v,w) = \left(\begin{array}{c} 0\\0\\ \frac{b-w}{\varepsilon}\end{array}\right).

The implicit systems are solved using the built-in modified Newton
iteration, with the ARKDENSE dense linear solver.  Both the Jacobian
routine and right-hand side functions are supplied by functions
provided in the example file.

The only non-default solver options are the tolerances
:math:`atol=10^{-10}` and :math:`rtol=10^{-6}`, adaptivity method 2 (I
controller), a maximum of 8 Newton iterations per step, a nonlinear
solver convergence coefficient :math:`nlscoef=10^{-8}`, and a maximum
of 1000 internal time steps.



Solutions
^^^^^^^^^^^^

With this setup, all three solution components exhibit a rapid
transient change during the first 0.2 time units, followed by a slow
and smooth evolution, as seen in the figure below.




Stiff nonlinear system example (ark_robertson)
------------------------------------------------

ODE system
^^^^^^^^^^^^

Numerical method
^^^^^^^^^^^^^^^^^

Solutions
^^^^^^^^^^^^



Stiff nonlinear system with root-finding example (ark_robertson_root)
-----------------------------------------------------------------------

ODE system
^^^^^^^^^^^^

Numerical method
^^^^^^^^^^^^^^^^^

Solutions
^^^^^^^^^^^^



Stiff PDE system example (ark_brusselator1D)
---------------------------------------------

ODE system
^^^^^^^^^^^^

Numerical method
^^^^^^^^^^^^^^^^^

Solutions
^^^^^^^^^^^^



PDE system example with iterative linear solver (ark_heat1D)
--------------------------------------------------------------

ODE system
^^^^^^^^^^^^

Numerical method
^^^^^^^^^^^^^^^^^

Solutions
^^^^^^^^^^^^


