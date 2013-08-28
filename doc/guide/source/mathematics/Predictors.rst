:tocdepth: 3


.. _Mathematics.Predictors:

Implicit predictors
===========================

As mentioned in the previous section, :ref:`Mathematics.IVP`, for
problems with implicit components, ARKode will employ a prediction
algorithm for constructing the initial guesses for each Runge-Kutta
stage, :math:`z_i^{(0)}`.  As is well-known with nonlinear solvers,
the selection of an initial guess can have dramatic effects on both
the speed and robustness of the nonlinear solve, enabling the
difference between divergence and quadratic convergence of the
iteration.  To this end, ARKode implements a variety of prediction
algorithms that may be selected by the user.  In each case, the stages
guesses :math:`z_i^{(0)}` are constructed explicitly using
readily-available information, including the previous step solutions
:math:`y_{n-1}` and :math:`y_{n-2}`, as well as any previous stage
solutions :math:`z_j, \quad j<i`.  In all cases, prediction is
performed by constructing an interpolating polynomial through
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
--------------------

The so-called "trivial predictor" is given by the formula

.. math::
   z_i^{(0)} = y_{n-1}.

While this piecewise-constant interpolant is clearly not an ideal
candidate for problems with time-varying solutions, it is often the
most robust approach for problems with constraints whose violation may
cause illegal solution values (e.g. a negative temperature).


.. _Mathematics.Predictors.Max:

Maximum order predictor
---------------------------

At the opposite end of the spectrum, ARKode will construct an
interpolant :math:`p_q(t)` of polynomial order up to :math:`q=3`.
Here, the function :math:`p_q(t)` is identical to the one used for
interpolation of output solution values between time steps, i.e. for
:math:`y(t)` values where :math:`t_{n-1} < t < t_n`.  The order of
this polynomial, :math:`q`, may be specified by the user with the
function :c:func:`ARKodeSetDenseOrder()`.

The interpolants generated are either of Lagrange or Hermite form, and
use the data :math:`\left\{ y_{n-2}, f_{n-2}, y_{n-1}, f_{n-1}
\right\}`, where by :math:`f_{k}` we mean 
:math:`M^{-1} \left(f_E(t_k,y_k) + f_I(t_k,y_k)\right)`.  Defining a
scaled and shifted "time" variables for the interval :math:`[t_{n-2},
t_{n-1}]` as

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
---------------------------

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
---------------------------

This predictor follows a similar idea as the previous algorithm, but
monitors the actual stage times to determine the polynomial
interpolant to use for prediction:

.. math::
   q = \begin{cases}
      q_{max}, & \text{if}\quad \tau < \tfrac12,\\
      1, & \text{otherwise}.
   \end{cases}



.. _Mathematics.Predictors.Bootstrap:

Bootstrap predictor
---------------------------

This predictor does not use any information from the preceding step,
aside from the previous solution :math:`y_{n-1}` and right-hand side
:math:`f(t_{n-1},y_{n-1})`.  However, unlike the trivial predictor,
in computing the predictor :math:`z_i^{(0)}` this approach will use
the right-hand side from a previously computed stage solution
:math:`f(t_{n-1}+c_j h,z_j)` to construct a quadratic Hermite
interpolant for the prediction.  For stages in which :math:`c_j=0` for
all previous stages :math:`j`, and for the first stage of any time step,
this method reduces to using the trivial predictor 
:math:`z_i^{(0)} = y_{n-1}`.  For stages in which multiple
:math:`c_j\ne 0`, :math:`j` will be chosen as the stage with the
largest :math:`c_j` value.
