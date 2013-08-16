:tocdepth: 3


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

