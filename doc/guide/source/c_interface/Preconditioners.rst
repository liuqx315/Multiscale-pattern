..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Copyright (c) 2013, Southern Methodist University.
   All rights reserved.
   For details, see the LICENSE file.
   ----------------------------------------------------------------

:tocdepth: 3


.. _CInterface.PreconditionerModules:

Preconditioner modules
============================

The efficiency of Krylov iterative methods for the solution of linear
systems can be greatly enhanced through preconditioning. For problems
in which the user cannot define a more effective, problem-specific
preconditioner, ARKode provides a banded preconditioner in the
module ARKBANDPRE and a band-block-diagonal preconditioner module
ARKBBDPRE. 


.. _CInterface.BandPre:

A serial banded preconditioner module
-------------------------------------------

This preconditioner provides a band matrix preconditioner for use with
any of the Krylov iterative linear solvers, when used in a serial
setting (i.e. with the NVECTOR_SERIAL module). It uses difference
quotients of the ODE right-hand side function :math:`f_I` to generate
a band matrix of bandwidth ``ml + mu + 1``, where the number of
super-diagonals (``mu``, the upper half-bandwidth) and sub-diagonals
(``ml``, the lower half-bandwidth) are specified by the user, and uses
this to form a preconditioner for use with the Krylov linear
solver. Although this matrix is intended to approximate the Jacobian
:math:`\frac{\partial f_I}{\partial y}`, it may be a very crude
approximation. The true Jacobian need not be banded, or its true
bandwidth may be larger than ``ml + mu + 1``, as long as the banded
approximation generated here is sufficiently accurate to speed
convergence as a preconditioner. 

In order to use the ARKBANDPRE module, the user need not define
any additional functions. Aside from the header files required for the
integration of the ODE problem (see the section
:ref:`CInterface.Headers`), to use the ARKBANDPRE module, the main
program must include the header file ``arkode_bandpre.h`` which
declares the needed function prototypes.  The following is a summary
of the usage of this module.  Steps that are unchanged from the
skeleton program presented in :ref:`CInterface.Skeleton` are
`italicized`. 

1. `Set problem dimensions`

2. `Set vector of initial values` 

3. `Create ARKode object` 

4. `Allocate internal memory` 

5. `Set optional inputs` 

6. Attach iterative linear solver, one of:

   (a) ``flag = ARKSpgmr(arkode_mem, pretype, maxl);`` 

   (b) ``flag = ARKSpbcg(arkode_mem, pretype, maxl);``

   (c) ``flag = ARKSptfqmr(arkode_mem, pretype, maxl);``

   (d) ``flag = ARKSpfgmr(arkode_mem, pretype, maxl);`` 

   (e) ``flag = ARKPcg(arkode_mem, pretype, maxl);``

7. Initialize the ARKBANDPRE preconditioner module 

   Specify the upper and lower half-bandwidths (``mu`` and ``ml``,
   respectively) and call 

   ``flag = ARKBandPrecInit(arkode_mem, N, mu, ml);``

   to allocate memory and initialize the internal preconditioner
   data. 

8. `Set linear solver optional inputs`

    Note that the user should not overwrite the preconditioner setup
    function or solve function through calls to the ARKSpilsSet*
    optional input functions. 

9. `Advance solution in time`

10. Get optional outputs 

   Additional optional outputs associated with ARKBANDPRE are
   available by way of the two routines described below,
   :c:func:`ARKBandPrecGetWorkSpace()` and
   :c:func:`ARKBandPrecGetNumRhsEvals()`.  

11. `Deallocate memory for solution vector`

12. `Free solver memory`

The ARKBANDPRE preconditioner module is initialized and attached
by calling the following function:



.. c:function:: int ARKBandPrecInit(void *arkode_mem, long int N, long int mu, long int ml)

   Initializes the ARKBANDPRE preconditioner and
   allocates required (internal) memory for it.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `N` -- problem dimension (size of ODE system).
      * `mu` -- upper half-bandwidth of the Jacobian approximation.
      * `ml` -- lower half-bandwidth of the Jacobian approximation.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_ILL_INPUT if an input has an illegal value
      * ARKSPILS_MEM_FAIL if a memory allocation request failed

   **Notes:** The banded approximate Jacobian will have nonzero elements
   only in locations :math:`(i,j)` with `ml` :math:`\le j-i \le` `mu`.


The following two optional output functions are available for use with
the ARKBANDPRE module:



.. c:function:: int ARKBandPrecGetWorkSpace(void *arkode_mem, long int *lenrwLS, long int *leniwLS)

   Returns the sizes of the ARKBANDPRE real and integer
   workspaces.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `lenrwLS` -- the number of ``realtype`` values in the
        ARKBANDPRE workspace.
      * `leniwLS` -- the number of integer values in the  ARKBANDPRE workspace.
   
   **Return value:** 
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_PMEM_NULL if the preconditioner memory is ``NULL``
   
   **Notes:** In terms of the problem size :math:`N` and `smu` :math:`=
   \min(N-1,` `mu+ml` :math:`)`, the actual size of the real
   workspace is :math:`(2` `ml + mu + smu` :math:`+2)N` ``realtype``
   words, and the actual size of the integer workspace is :math:`N`
   integer words.
   
   The workspaces referred to here exist in addition to those given by
   the corresponding function ARKSpils*GetWorkspace.



.. c:function:: int ARKBandPrecGetNumRhsEvals(void *arkode_mem, long int *nfevalsBP)

   Returns the number of calls made to the user-supplied
   right-hand side function :math:`f_I` for constructing the
   finite-difference banded Jacobian approximation used within the
   preconditioner setup function.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `nfevalsBP` -- number of calls to :math:`f_I`
   
   **Return value:**  
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_PMEM_NULL if the preconditioner memory is ``NULL``
   
   **Notes:**  The counter `nfevalsBP` is distinct from the counter
   `nfevalsLS` returned by the corresponding function
   ARKSpils*GetNumRhsEvals and also from `nfi_evals` returned by
   :c:func:`ARKodeGetNumRhsEvals()`.  The total number of right-hand
   side function evaluations is the sum of all three of these
   counters, plus the `nfe_evals` counter for :math:`f_E` calls
   returned by :c:func:`ARKodeGetNumRhsEvals()`.



.. _CInterface.BBDPre:

A parallel band-block-diagonal preconditioner module
---------------------------------------------------------

A principal reason for using a parallel ODE solver such as ARKode
lies in the solution of partial differential equations
(PDEs). Moreover, the use of a Krylov iterative method for the
solution of many such problems is motivated by the nature of the
underlying linear system of equations that must be solved at each time
step. For many PDEs, the linear algebraic system is large, sparse and
structured.  However, if a Krylov iterative method is to be effective
in this setting, then a nontrivial preconditioner needs to be used.
Otherwise, the rate of convergence of the Krylov iterative method is
usually unacceptably slow.  Unfortunately, an effective preconditioner
tends to be problem-specific.  However, we have developed one type of
preconditioner that treats a rather broad class of PDE-based
problems. It has been successfully used with CVODE for several
realistic, large-scale problems [HT1998]_ and is included
in a software module within the ARKode package. This module works
with the parallel vector module NVECTOR_PARALLEL and is usable
with any of the Krylov iterative linear solvers. It generates a
preconditioner that is a block-diagonal matrix with each block being a
band matrix. The blocks need not have the same number of super- and
sub-diagonals and these numbers may vary from block to block. This
Band-Block-Diagonal Preconditioner module is called ARKBBDPRE. 

One way to envision these preconditioners is to think of the domain of
the computational PDE problem as being subdivided into :math:`Q`
non-overlapping subdomains. Each of these subdomains is then assigned
to one of the :math:`Q` processes to be used to solve the ODE
system. The basic idea is to isolate the preconditioning so that it is
local to each process, and also to use a (possibly cheaper)
approximate right-hand side function. This requires the definition of
a new function :math:`g(t,y)` which approximates the function
:math:`f_I(t,y)` in the definition of the ODE system, 

.. math::
   M\dot{y} = f_E(t,y) + f_I(t,y),

where :math:`f_I` corresponds to the ODE components to be treated
implicitly.  However, the user may set :math:`g = f_I`. Corresponding
to the domain decomposition, there is a decomposition of the solution
vector :math:`y` into :math:`Q` disjoint blocks :math:`y_q`, and a
decomposition of :math:`g` into blocks :math:`g_q`. The block
:math:`g_q` depends both on :math:`y_p` and on components of blocks
:math:`y_{q'}`	associated with neighboring subdomains (so-called
ghost-cell data). Let :math:`\bar{y}_q` denote :math:`y_q` augmented
with those other components on which :math:`g_q` depends. Then we have 

.. math::
   g(t,y) = \left[ g_1(t,\bar{y}_1), g_2(t,\bar{y}_2), \ldots , g_Q(t,\bar{y}_Q) \right]^T

and each of the blocks :math:`g_q(t,\bar{y}_q)` is decoupled from the
others.

The preconditioner associated with this decomposition has the form

.. math::
   P = diag[P_1, P_2, \ldots, P_Q]

where

.. math::
   P_q \approx M - \gamma J_q

and where :math:`J_q` is a difference quotient approximation to
:math:`\frac{\partial g_q}{\partial \bar{y}_q}`.  This matrix is taken
to be banded, with upper and lower half-bandwidths `mudq` and
`mldq` defined as the number of non-zero diagonals above and below
the main diagonal, respectively.  The difference quotient
approximation is computed using `mudq + mldq + 2` evaluations of
:math:`g_m`, but only a matrix of bandwidth `mukeep + mlkeep + 1` is
retained. Neither pair of parameters need be the true half-bandwidths
of the Jacobian of the local block of :math:`g`, if smaller values
provide a more efficient preconditioner. The solution of the complete
linear system 

.. math::
   Px = b

reduces to solving each of the distinct equations

.. math::
   P_q x_q = b_q, \quad q=1,\ldots,Q,

and this is done by banded LU factorization of :math:`P_q` followed by
a banded backsolve.

Similar block-diagonal preconditioners could be considered with
different treatments of the blocks :math:`P_q`.  For example,
incomplete LU factorization or an iterative method could be used
instead of banded LU factorization.

The ARKBBDPRE module calls two user-provided functions to
construct :math:`P`: a required function `gloc` (of type
:c:func:`ARKLocalFn()`) which approximates the right-hand side function
:math:`g(t,y) \approx f_I(t,y)` and which is computed locally, and an
optional function `cfn` (of type :c:func:`ARKCommFn()`) which performs all
interprocess communication necessary to evaluate the approximate
right-hand side :math:`g`. These are in addition to the user-supplied
right-hand side function :math:`f_I`. Both functions take as input the
same pointer `user_data` that is passed by the user to
:c:func:`ARKodeSetUserData()` and that was passed to the user's
function :math:`f_I`. The user is responsible for providing space
(presumably within `user_data`) for components of :math:`y` that are
communicated between processes by `cfn`, and that are then used by
`gloc`, which should not do any communication.



.. c:function:: typedef int (*ARKLocalFn)(long int Nlocal, realtype t, N_Vector y, N_Vector glocal, void *user_data)

   This `gloc` function computes :math:`g(t,y)`.  It
   loads the vector `glocal` as a function of `t` and `y`.
   
   **Arguments:**
      * `Nlocal` -- the local vector length
      * `t` -- the value of the independent variable
      * `y` -- the value of the dependent variable vector on this process
      * `glocal` -- the output vector of :math:`g(t,y)` on this process
      * `user_data` -- a pointer to user data, the same as the
        `user_data` parameter passed to :c:func:`ARKodeSetUserData()`.
   
   **Return value:**  
   An ARKLocalFn should return 0 if successful, a
   positive value if a recoverable error occurred (in which case
   ARKode will attempt to correct), or a negative value if it
   failed unrecoverably (in which case the integration is halted and
   :c:func:`ARKode()` will return ARK_LSETUP_FAIL).
   
   **Notes:**  This function must assume that all interprocess communication
   of data needed to calculate `glocal` has already been done, and that
   this data is accessible within user data. 
   
   The case where :math:`g` is mathematically identical to :math:`f_I`
   is allowed. 



.. c:function:: typedef int (*ARKCommFn)(long int Nlocal, realtype t, N_Vector y, void *user_data)

   This `cfn` function performs all interprocess
   communication necessary for the executation of the `gloc` function
   above, using the input vector `y`.
   
   **Arguments:**
      *  `Nlocal` -- the local vector length
      * `t` -- the value of the independent variable
      * `y` -- the value of the dependent variable vector on this process
      * `user_data` -- a pointer to user data, the same as the
        `user_data` parameter passed to :c:func:`ARKodeSetUserData()`.
   
   **Return value:**  
   An ARKCommFn should return 0 if successful, a
   positive value if a recoverable error occurred (in which case
   ARKode will attempt to correct), or a negative value if it
   failed unrecoverably (in which case the integration is halted and
   :c:func:`ARKode()` will return ARK_LSETUP_FAIL).
   
   **Notes:**  The `cfn` function is expected to save communicated data in
   space defined within the data structure `user_data`.
   
   Each call to the `cfn` function is preceded by a call to the
   right-hand side function :math:`f_I` with the same :math:`(t,y)`
   arguments. Thus, `cfn` can omit any communication done by
   :math:`f_I` if relevant to the evaluation of `glocal`. If all
   necessary communication was done in :math:`f_I`, then `cfn` =
   ``NULL`` can be passed in the call to :c:func:`ARKBBDPrecInit()` (see
   below).



Besides the header files required for the integration of the ODE problem (see the section
:ref:`CInterface.Headers`), to use the ARKBBDPRE module, the main
program must include the header file ``arkode_bbdpre.h`` which
declares the needed function prototypes. 

The following is a summary of the proper usage of this module. Steps
that are unchanged from the skeleton program presented in
:ref:`CInterface.Skeleton` are `italicized`.

1. `Initialize MPI`

2. `Set problem dimensions`

3. `Set vector of initial values`

4. `Create ARKode object`

5. `Allocate internal memory`

6. `Set optional inputs`

7. Attach iterative linear solver, one of:

   (a) ``flag = ARKSpgmr(arkode_mem, pretype, maxl);``

   (b) ``flag = ARKSpbcg(arkode_mem, pretype, maxl);``

   (c) ``flag = ARKSptfqmr(arkode_mem, pretype, maxl);``

   (d) ``flag = ARKSpfgmr(arkode_mem, pretype, maxl);``

   (e) ``flag = ARKPcg(arkode_mem, pretype, maxl);``

8. Initialize the ARKBBDPRE preconditioner module 

   Specify the upper and lower half-bandwidths for computation
   ``mudq`` and ``mldq``, the upper and lower half-bandwidths for
   storage ``mukeep`` and ``mlkeep``, and call 

   ``flag = ARKBBDPrecInit(arkode_mem, Nlocal, mudq, mldq, mukeep, mlkeep, dqrely, gloc, cfn);``

   to allocate memory and initialize the internal preconditioner
   data. The last two arguments of :c:func:`ARKBBDPrecInit()` are the
   two user-supplied functions of type :c:func:`ARKLocalFn()` and
   :c:func:`ARKCommFn()` described above, respectivelyl. 

9. `Set the linear solver optional inputs`

   Note that the user should not overwrite the preconditioner setup
   function or solve function through calls to ARKSPILS optional
   input functions. 

10. `Advance solution in time`

11. `Get optional outputs`

    Additional optional outputs associated with ARKBBDPRE are
    available by way of the two routines described below,
    :c:func:`ARKBBDPrecGetWorkSpace()` and
    :c:func:`ARKBBDPrecGetNumGfnEvals()`. 

12. `Deallocate memory for solution vector`

13. `Free solver memory`

14. `Finalize MPI`

The user-callable functions that initialize (step 8 above) or
re-initialize the ARKBBDPRE preconditioner module are described
next.



.. c:function:: int ARKBBDPrecInit(void *arkode_mem, long int Nlocal, long int mudq, long int mldq, long int mukeep, long int mlkeep, realtype dqrely, ARKLocalFn gloc, ARKCommFn cfn)

   Initializes and allocates (internal) memory for the
   ARKBBDPRE preconditioner.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `Nlocal` -- local vector length.
      * `mudq` -- upper half-bandwidth to be used in the difference
        quotient Jacobian approximation.
      * `mldq` -- lower half-bandwidth to be used in the difference
        quotient Jacobian approximation.
      * `mukeep` -- upper half-bandwidth of the retained banded
        approximate Jacobian block.
      * `mlkeep` -- lower half-bandwidth of the retained banded
        approximate Jacobian block.
      * `dqrely` -- the relative increment in components of `y` used in
        the difference quotient approximations.  The default is `dqrely`
        = :math:`\sqrt{\text{unit roundoff}}`, which can be specified by
        passing `dqrely` = 0.0.
      * `gloc` -- the name of the C function (of type :c:func:`ARKLocalFn()`)
        which computes the approximation :math:`g(t,y) \approx f_I(t,y)`.
      * `cfn` -- the name of the C function (of type :c:func:`ARKCommFn()`) which
        performs all interprocess communication required for the
        computation of :math:`g(t,y)`.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_ILL_INPUT if an input has an illegal value
      * ARKSPILS_MEM_FAIL if a memory allocation request failed
   
   **Notes:**  If one of the half-bandwidths `mudq` or `mldq` to be used
   in the difference quotient calculation of the approximate Jacobian is
   negative or exceeds the value `Nlocal-1`, it is replaced by 0 or
   `Nlocal-1` accordingly. 
   
   The half-bandwidths `mudq` and `mldq` need not be the true
   half-bandwidths of the Jacobian of the local block of :math:`g`
   when smaller values may provide a greater efficiency. 
   
   Also, the half-bandwidths `mukeep` and `mlkeep` of the retained
   banded approximate Jacobian block may be even smaller than
   `mudq` and `mldq`, to reduce storage and computational costs
   further. 
   
   For all four half-bandwidths, the values need not be the same on
   every processor.



The ARKBBDPRE module also provides a reinitialization function to
allow solving a sequence of problems of the same size, with the same
linear solver choice, provided there is no change in `Nlocal`,
`mukeep`, or `mlkeep`. After solving one problem, and after
calling :c:func:`ARKodeReInit()` to re-initialize ARKode for a
subsequent problem, a call to :c:func:`ARKBBDPrecReInit()` can be made
to change any of the following: the half-bandwidths `mudq` and
`mldq` used in the difference-quotient Jacobian approximations, the
relative increment `dqrely`, or one of the user-supplied functions
`gloc` and `cfn`. If there is a change in any of the linear solver
inputs, an additional call to :c:func:`ARKSpgmr()`,
:c:func:`ARKSpbcg()`, :c:func:`ARKSptfqmr()`, :c:func:`ARKSpfgmr()`,
or :c:func:`ARKPcg()`, and/or one or more of the corresponding
ARKSpils*Set* functions, must also be made (in the proper order).



.. c:function:: int ARKBBDPrecReInit(void *arkode_mem, long int mudq, long int mldq, realtype dqrely)

   Re-initializes the ARKBBDPRE preconditioner module.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `mudq` -- upper half-bandwidth to be used in the difference
        quotient Jacobian approximation.
      * `mldq` -- lower half-bandwidth to be used in the difference
        quotient Jacobian approximation.
      * `dqrely` -- the relative increment in components of `y` used in
        the difference quotient approximations.  The default is `dqrely`
        = :math:`\sqrt{\text{unit roundoff}}`, which can be specified by
        passing `dqrely` = 0.0.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_PMEM_NULL if the preconditioner memory is ``NULL``
   
   **Notes:**  If one of the half-bandwidths `mudq` or `mldq` is
   negative or exceeds the value `Nlocal-1`, it is replaced by 0 or
   `Nlocal-1` accordingly. 


The following two optional output functions are available for use with
the ARKBBDPRE module:



.. c:function:: int ARKBBDPrecGetWorkSpace(void *arkode_mem, long int *lenrwBBDP, long int *leniwBBDP)

   Returns the processor-local ARKBBDPRE real and
   integer workspace sizes.
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `lenrwBBDP` -- the number of ``realtype`` values in the
        ARKBBDPRE workspace.
      * `leniwBBDP` -- the number of integer values in the  ARKBBDPRE workspace.
   
   **Return value:**  
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_PMEM_NULL if the preconditioner memory is ``NULL``
   
   **Notes:**  In terms of `Nlocal` and `smu = min(Nlocal-1,
   mukeep+mlkeep)`, the actual size of the real workspace is `(2
   mlkeep + mukeep + smu + 2)*Nlocal`  ``realtype`` words, and the
   actual size of the integer workspace is `Nlocal` integer
   words. These values are local to each process. 
   
   The workspaces referred to here exist in addition to those given by
   the corresponding function ARKSpils*GetWorkSpace. 



.. c:function:: int ARKBBDPrecGetNumGfnEvals(void *arkode_mem, long int *ngevalsBBDP)

   Returns the number of calls made to the user-supplied
   `gloc` function (of type :c:func:`ARKLocalFn()`) due to the finite
   difference approximation of the Jacobian blocks used within the
   preconditioner setup function. 
   
   **Arguments:**
      * `arkode_mem` -- pointer to the ARKode memory block.
      * `ngevalsBBDP` -- the number of calls made to the user-supplied
        `gloc` function. 
   
   **Return value:**  
      * ARKSPILS_SUCCESS if no errors occurred
      * ARKSPILS_MEM_NULL if the integrator memory is ``NULL``
      * ARKSPILS_LMEM_NULL if the linear solver memory is ``NULL``
      * ARKSPILS_PMEM_NULL if the preconditioner memory is ``NULL``
   
   
In addition to the `ngevalsBBDP` `gloc` evaluations, the costs
associated with ARKBBDPRE also include `nlinsetups` LU
factorizations, `nlinsetups` calls to `cfn`, `npsolves` banded
backsolve calls, and `nfevalsLS` right-hand side function
evaluations, where `nlinsetups` is an optional ARKode output and
`npsolves` and `nfevalsLS` are linear solver optional outputs (see
the table :ref:`CInterface.ARKSpilsOutputTable`).
