:tocdepth: 3


.. _Mathematics.Stability:

Explicit stability
======================

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
   |h_{exp}| < \frac{\Delta x}{|\lambda|}

where :math:`\Delta x` is the spatial mesh size and :math:`\lambda` is
the fastest advective wave speed.

In the case that a user has supplied a routine to predict these
explicitly stable step sizes (by calling the function
:c:func:`ARKodeSetStabilityFn()`), the value :math:`|h_{exp}|` is
compared against that resulting from the local error adaptivity,
:math:`|h_{acc}|`, and the step used by ARKode will satisfy 

.. math::
   h' = \frac{h}{|h|}\min\{c\, |h_{exp}|,\, |h_{acc}|\},

where the explicit stability step factor :math:`c>0` may be modified
through the function :c:func:`ARKodeSetAdaptivityMethod()`, and has a
default value of :math:`1/2`.
