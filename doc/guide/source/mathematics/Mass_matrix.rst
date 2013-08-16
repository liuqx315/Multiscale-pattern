:tocdepth: 3


.. _Mathematics.MassSolve:

Mass matrix solver
=======================

Within the approach described above, there are three locations where a
linear solve of the form

.. math::
   M x = b

is required: in constructing the time-evolved solution :math:`y_n`, in
estimating the local temporal truncation error, and in constructing
predictors for the implicit solver iteration (see section
:ref:`Mathematics.Predictors.Max`).  Specifically, to construct the
time-evolved solution :math:`y_n` we must solve

.. math::
   &M y_n \ = \ M y_{n-1} + h_n \sum_{i=0}^{s} b_i \left(f_E(t_{n-1} + c_i h_n, z_i) 
                 + f_I(t_{n-1} + c_i h_n, z_i)\right), \\
   \Leftrightarrow \qquad & \\
   &M (y_n -y_{n-1}) \ = \ h_n \sum_{i=0}^{s} b_i \left(f_E(t_{n-1} + c_i h_n, z_i) 
                 + f_I(t_{n-1} + c_i h_n, z_i)\right), \\
   \Leftrightarrow \qquad & \\
   &M \nu \ = \ h_n \sum_{i=0}^{s} b_i \left(f_E(t_{n-1} + c_i h_n, z_i) 
                 + f_I(t_{n-1} + c_i h_n, z_i)\right),

for the update :math:`\nu = y_n - y_{n-1}`.  Similarly, in computing
the local temporal error estimate :math:`T_n` we must solve systems of
the form 

.. math::
   M\, T_n = h \sum_{i=0}^{s} \left(b_i - \tilde{b}_i\right) 
   \left(f_E(t_{n-1} + c_i h_n, z_i) + f_I(t_{n-1} + c_i h_n, z_i)\right).

Lastly, in constructing implicit predictors of order 2 or higher, we
must compute the derivative information :math:`f_k` from the equation

.. math::
   M f_k = f_E(t_k, y_k) + f_I(t_k, y_k).

Of course, for problems in which :math:`M=I` these solves are not
required; however for problems with non-identity :math:`M`, ARKode may
use either an iterative linear solver or a dense linear solver, in the
same manner as described in section :ref:`Mathematics.IVP` for solving
the linear Newton systems.  We note that at present, the matrix
:math:`M` may depend on time :math:`t` but must be independent of the
solution :math:`y`, since we assume that each of the above systems are
linear.

**[Continue with documentation here]**

