/*---------------------------------------------------------------
  $Revision: 1.0 $
  $Date: $
 ---------------------------------------------------------------- 
  Programmer(s): Daniel R. Reynolds @ SMU
 ----------------------------------------------------------------
  This is the header file for FARKODE, the Fortran interface to
  the ARKODE package.                                            
 --------------------------------------------------------------*/

/*===============================================================
                FARKODE Interface Package

 The FARKODE Interface Package is a package of C functions which 
 support the use of the ARKODE solver, for the solution of ODE 
 systems 
         dy/dt = fe(t,y) + fi(t,y), 
 in a mixed Fortran/C setting.  While ARKODE is written in C, it
 is assumed here that the user's calling program and user-supplied
 problem-defining routines are written in Fortran.  This package 
 provides the necessary interface to ARKODE any acceptable NVECTOR
 implementation.
 
 A summary of the user-callable functions, with the corresponding 
 ARKODE functions, are as follows:
 
   Fortran                    ARKODE
   ---------------------      --------------------------------
   FNVINITS                   N_VNew_Serial
   FNVINITP                   N_VNew_Parallel

   FARKMALLOC                 ARKodeCreate, ARKodeSetUserData, 
                                 and ARKodeInit
   FARKREINIT                 ARKReInit
   FARKSETIIN                 ARKodeSet* (integer arguments)
   FARKSETRIN                 ARKodeSet* (real arguments)
   FARKEWTSET                 ARKodeWFtolerances

   FARKDENSE                  ARKDense
   FARKDENSESETJAC            ARKDenseSetJacFn

   FARKBAND                   ARKBand
   FARKBANDSETJAC             ARKBandSetJacFn

   FARKLAPACKDENSE            ARKLapackDense
   FARKLAPACKDENSESETJAC      ARKLapackSetJacFn

   FARKLAPACKBAND             ARKLapackBand
   FARKLAPACKBANDSETJAC       ARKLapackSetJacFn

   FARKSPGMR                  ARKSpgmr and ARKSpilsSet*
   FARKSPGMRREINIT            ARKSpilsSet*

   FARKSPBCG                  ARKSpbcg and ARKSpilsSet*
   FARKSPBCGREINIT            ARKSpilsSet*

   FARKSPTFQMR                ARKSptfqmr and ARKSpilsSet*
   FARKSPTFQMRREINIT          ARKSpilsSet*

   FARKSPILSSETJAC            ARKSpilsSetJacTimesVecFn
   FARKSPILSSETPREC           ARKSpilsSetPreconditioner
 
   FARKODE                    ARKode, ARKodeGet*, and ARK*Get*
   FARKDKY                    ARKodeGetDky

   FARKGETERRWEIGHTS          ARKodeGetErrWeights
   FARKGETESTLOCALERR         ARKodeGetEstLocalErrors

   FARKFREE                   ARKodeFree
   ---------------------      --------------------------------

 
 The user-supplied functions, each listed with the corresponding interface
 function which calls it (and its type within ARKODE), are as follows:

   Fortran:         Interface Fcn:           ARKODE Type:
   ----------       ------------------       ---------------------
   FARKIFUN         FARKfi                   ARKRhsFn
   FARKEFUN         FARKfe                   ARKRhsFn
   FARKDJAC         FARKDenseJac             ARKDenseJacFn
   FARKBJAC         FARKBandJac              ARKBandJacFn
   FARKLDJAC        FARKLapackDenseJac       ARKLapackJacFn
   FARKLBJAC        FARKLapackBandJac        ARKLapackJacFn
   FARKPSOL         FARKPSol                 ARKSpilsPrecSolveFn
   FARKPSET         FARKPSet                 ARKSpilsPrecSetupFn
   FARKJTIMES       FARKJtimes               ARKSpilsJacTimesVecFn
   FARKEWT          FARKEwtSet               ARKEwtFn
   ----------       ------------------       ---------------------

 In contrast to the case of direct use of ARKODE, and of most Fortran ODE
 solvers, the names of all user-supplied routines here are fixed, in
 order to maximize portability for the resulting mixed-language program.
 
 Important note on portability:  In this package, the names of the interface 
 functions, and the names of the Fortran user routines called by them, appear 
 as dummy names which are mapped to actual values by a series of definitions, 
 in this and other header files.
 
 =============================================================================
 
                  Usage of the FARKODE Interface Package
 
 The usage of FARKODE requires calls to five or more interface
 functions, depending on the method options selected, and one or more
 user-supplied routines which define the problem to be solved.  These
 function calls and user routines are summarized separately below.
 
 Some details are omitted, and the user is referred to the ARKODE user 
 documentation for more complete information.  Information on the arguments 
 of any given user-callable interface routine, or of a given user-supplied 
 function called by an interface function, can be found in the documentation 
 on the corresponding function in the ARKODE package.
 
 The number labels on the instructions below end with s for instructions
 that are specific to use with the N_VSerial package; similarly those that 
 end with p are specific to use with the N_VParallel package.

 -----------------------------------------------------------------------------

                               Data Types

 Throughout this documentation, we will refer to data types according to their
 usage in SUNDIALS.  The equivalent types to these may vary, depending on your
 computer architecture and on how SUNDIALS was compiled.  A Fortran user 
 should take care that all arguments passed through this Fortran/C interface 
 are declared of the appropriate type.
 
 Integers: SUNDIALS uses both 'int' and 'long int' types:
   int      -- equivalent to an INTEGER or INTEGER*4 in Fortran
   long int -- this will depend on the computer architecture:
                 32-bit -- equivalent to an INTEGER or INTEGER*4 in Fortran
                 64-bit -- equivalent to an INTEGER*8 in Fortran
	      
 Real numbers:  At compilation, SUNDIALS allows the configuration option 
 '--with-precision', that accepts values of 'single', 'double' or 'extended' 
 (the default is 'double').  This choice dictates the size of a SUNDIALS 
 'realtype' variable.  The corresponding Fortran types for these 'realtype' 
 sizes are:
   single   -- equivalent to a REAL or REAL*4 in Fortran
   double   -- equivalent to a DOUBLE PRECISION or REAL*8 in Fortran
   extended -- equivalent to a REAL*16 in Fortran

 -----------------------------------------------------------------------------

 (1) User-supplied right-hand side routines: FARKIFUN and FARKEFUN
     The user must in all cases supply at least one of the following Fortran 
     routines:

       SUBROUTINE FARKIFUN(T, Y, YDOT, IPAR, RPAR, IER)

     Sets the YDOT array to fi(T,Y), the implicit portion of the right-hand 
     side of the ODE system, as function of time T and the state variable 
     array Y.

       SUBROUTINE FARKEFUN(T, Y, YDOT, IPAR, RPAR, IER)

     Sets the YDOT array to fe(t,y), the explicit portion of the right-hand 
     side of the ODE system, as function of time T and the state variable 
     array Y.  

     The arguments are:
       Y    -- array containing state variables [realtype, input]
       YDOT -- array containing state derivatives [realtype, output]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful, 
                 >0 if a recoverable error occurred,
                 <0 if an unrecoverable error ocurred.
 
 (2s) Optional user-supplied dense Jacobian approximation routine: FARKDJAC

     As an option when using the DENSE linear solver, the user may supply a
     routine that computes a dense approximation of the system Jacobian 
     J = dfi(t,y)/dy.  If supplied, it must have the following form:

       SUBROUTINE FARKDJAC(NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, 
      &                    WK3, IER)

     Typically this routine will use only NEQ, T, Y, and DJAC. It must compute
     the Jacobian and store it column-wise in DJAC.

     The arguments are:
       NEQ  -- number of rows in the matrix [long int, input]
       T    -- current time [realtype, input]
       Y    -- array containing state variables [realtype, input]
       FY   -- array containing state derivatives [realtype, input]
       DJAC -- 2D array containing the jacobian entries [realtype of size
               (NEQ,NEQ), output]
       H    -- current step size [realtype, input]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       WK*  -- array containing temporary workspace of same size as Y 
               [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful, 
                 >0 if a recoverable error occurred,
                 <0 if an unrecoverable error ocurred.
 
 (3s) Optional user-supplied band Jacobian approximation routine: FARKBJAC

     As an option when using the BAND linear solver, the user may supply a
     routine that computes a band approximation of the system Jacobian 
     J = dfi(t,y)/dy. If supplied, it must have the following form:

       SUBROUTINE FARKBJAC(NEQ, MU, ML, MDIM, T, Y, FY, BJAC, H,
      &                    IPAR, RPAR, WK1, WK2, WK3, IER)

     Typically this routine will use only NEQ, MU, ML, T, Y, and BJAC. It 
     must load the MDIM by N array BJAC with the Jacobian matrix at the
     current (t,y) in band form.  Store in BJAC(k,j) the Jacobian element 
     J(i,j)  with k = i - j + MU + 1 (k = 1 ... ML+MU+1) and j = 1 ... N.

     The arguments are:
       NEQ  -- number of rows in the matrix [long int, input]
       MU   -- upper half-bandwidth of the matrix [long int, input]
       ML   -- lower half-bandwidth of the matrix [long int, input]
       MDIM -- leading dimension of BJAC array [long int, input]
       T    -- current time [realtype, input]
       Y    -- array containing state variables [realtype, input]
       FY   -- array containing state derivatives [realtype, input]
       BJAC -- 2D array containing the jacobian entries [realtype of size
               (MDIM,NEQ), output]
       H    -- current step size [realtype, input]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       WK*  -- array containing temporary workspace of same size as Y 
               [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful, 
                 >0 if a recoverable error occurred,
                 <0 if an unrecoverable error ocurred.
 

 (4s) Optional user-supplied Lapack dense Jacobian routine: FARKLDJAC

     See the description for FARKDJAC. NOTE: the dense Jacobian matrix
     is NOT set to zero before calling the user's FARKLDJAC.

 (5s) Optional user-supplied Lapack band Jacobian routine: FARKLBJAC

     See the description for FARKBJAC. NOTE: the band Jacobian matrix
     is NOT set to zero before calling the user's FARKLBJAC.

 (6) Optional user-supplied Jacobian-vector product routine: FARKJTIMES

     As an option when using the SP* linear solver, the user may supply a
     routine that computes the product of the system Jacobian 
     J = dfi(t,y)/dy and a given vector v.  If supplied, it must have the 
     following form:

       SUBROUTINE FARKJTIMES(V, FJV, T, Y, FY, H, IPAR, RPAR, WORK, IER)

     Typically this routine will use only NEQ, T, Y, V, and FJV.  It must
     compute the product vector J*v where the vector V, and store the product
     in FJV.  

     The arguments are:
       V    -- array containing vector to multiply [realtype, input]
       FJV  -- array containing product vector [realtype, output]
       T    -- current time [realtype, input]
       Y    -- array containing state variables [realtype, input]
       FY   -- array containing state derivatives [realtype, input]
       H    -- current step size [realtype, input]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       WORK -- array containing temporary workspace of same size as Y 
               [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful, 
                  nonzero if an error.
 
 (7) Optional user-supplied error weight vector routine: FARKEWT
 
     As an option to providing the relative and absolute tolerances, the user
     may supply a routine that computes the weights used in the WRMS norms.
     If supplied, it must have the following form:

       SUBROUTINE FARKEWT(Y, EWT, IPAR, RPAR, IER)

     It must store the error weights in EWT, given the current solution Y. 

     The arguments are:
       Y    -- array containing state variables [realtype, input]
       EWT  -- array containing the error weight vector [realtype, output]
       IPAR -- array containing integer user data that was passed to
               FARKMALLOC [long int, input]
       RPAR -- array containing real user data that was passed to
               FARKMALLOC [realtype, input]
       IER  -- return flag [int, output]:
                  0 if successful, 
                  nonzero if an error.

 -----------------------------------------------------------------------------

 (8) Initialization:  FNVINITS / FNVINITP, FARKMALLOC, FARKREINIT
 
 (8.1s) To initialize the serial vector specification, the user must make the
     following call:

        CALL FNVINITS(4, NEQ, IER)

     where the first argument is an int containing the ARKODE solver ID (4). 
     The other arguments are:
        NEQ = size of vectors [long int, input]
	IER = return completion flag [int, output]:
	          0 = success, 
		 -1 = failure.
 
 (8.1p) To initialize the parallel machine environment, the user must make 
     the following call:

        CALL FNVINITP(COMM, 4, NLOCAL, NGLOBAL, IER)

     The arguments are:
        COMM = the MPI communicator [int, input]
        NLOCAL = local size of vectors on this processor [long int, input]
        NGLOBAL = the system size, and the global size of vectors (the sum 
           of all values of NLOCAL) [long int, input]
        IER = return completion flag [int, ouptut]. 
                  0 = success, 
                 -1 = failure.

     Note: If MPI was initialized by the user, the communicator must be set to
     MPI_COMM_WORLD.  If not, this routine initializes MPI and sets the
     communicator equal to MPI_COMM_WORLD.
 
 (8.2) To set various problem and solution parameters and allocate
     internal memory, make the following call:

       CALL FARKMALLOC(T0, Y0, IMEX, IATOL, RTOL, ATOL, IOUT, ROUT, 
      &                IPAR, RPAR, IER)

     The arguments are:
        T0 = initial value of t [realtype, input]
	Y0 = array of initial conditions [realtype, input]
	IMEX = flag denoting basic integration method [int, input]: 
                  0 = implicit, 
                  1 = explicit, 
                  2 = imex
        IATOL = type for absolute tolerance ATOL [int, input]: 
                  1 = scalar, 
                  2 = array,
                  3 = user-supplied function; the user must supply a routine
		  FARKEWT to compute the error weight vector.
        RTOL = scalar relative tolerance [realtype, input]
	ATOL = scalar or array absolute tolerance [realtype, input]
	IOUT = array of length 22 for integer optional outputs
	   [long int, output]
	ROUT = array of length 6 for real optional outputs [realtype, output]
	IPAR = array of user integer data [long int, input/output]
	RPAR = array with user real data [realtype, input/output]
	IER  = return completion flag [int, output]:
                  0 = SUCCESS,
                 -1 = failure (see printed message for failure details).

     The user data arrays IPAR and RPAR are passed unmodified to all 
     subsequent calls to user-provided routines. Modifications to either array
     inside a user-provided routine will be propagated. Using these two 
     arrays, the user can dispense with COMMON blocks to pass data betwen 
     user-provided routines. 
 
     The optional outputs are:
           LENRW   = IOUT( 1) from ARKodeGetWorkSpace
           LENIW   = IOUT( 2) from ARKodeGetWorkSpace
           NST     = IOUT( 3) from ARKodeGetNumSteps
           NST_STB = IOUT( 4) from ARKodeGetNumExpSteps
           NST_ACC = IOUT( 5) from ARKodeGetNumAccSteps
           NST_CNV = IOUT( 6) from ARKodeGetNumConvSteps
           NFE     = IOUT( 7) from ARKodeGetNumRhsEvals (fe(t,y) calls)
           NFI     = IOUT( 8) from ARKodeGetNumRhsEvals (fi(t,y) calls)
           NSETUPS = IOUT( 9) from ARKodeGetNumLinSolvSetups
           NETF    = IOUT(10) from ARKodeGetNumErrTestFails
           NNI     = IOUT(11) from ARKodeGetNumNonlinSolvIters
           NCFN    = IOUT(12) from ARKodeGetNumNonlinSolvConvFails
           NGE     = IOUT(13) from ARKodeGetNumGEvals

           H0U     = ROUT( 1) from ARKodeGetActualInitStep
           HU      = ROUT( 2) from ARKodeGetLastStep
           HCUR    = ROUT( 3) from ARKodeGetCurrentStep
           TCUR    = ROUT( 4) from ARKodeGetCurrentTime
           TOLSF   = ROUT( 5) from ARKodeGetTolScaleFactor
           UROUND  = ROUT( 6) from UNIT_ROUNDOFF
     See the ARKODE manual for details. 

     If the user program includes the FARKEWT routine for the evaluation of 
     the error weights, the following call must be made

       CALL FARKEWTSET(FLAG, IER)

     with the int argument FLAG = 1 to specify that FARKEWT is provided.
     The int return flag IER is 0 if successful, and nonzero otherwise.

 (8.3) To re-initialize the ARKODE solver for the solution of a new problem
     of the same size as one already solved, make the following call:

       CALL FARKREINIT(T0, Y0, IMEX, IATOL, RTOL, ATOL, IER)

     The arguments have the same names and meanings as those of FARKMALLOC.  
     FARKREINIT performs the same initializations as FARKMALLOC, but does no 
     memory allocation, using instead the existing internal memory created by 
     the previous FARKMALLOC call.  The subsequent call to specify the linear 
     system solution method may or may not be needed; see paragraph (9) below.
 
 (8.4) To set various integer optional inputs, make the folowing call:

       CALL FARKSETIIN(KEY, VALUE, IER)

     to set the integer value VAL to the optional input specified by the 
     quoted character string KEY. KEY must be one of the following: 
     ORDER, DENSE_ORDER, LINEAR, NONLINEAR, EXPLICIT, IMPLICIT, IMEX, 
     IRK_TABLE_NUM, ERK_TABLE_NUM, MAX_NSTEPS, HNIL_WARNS, PREDICT_METHOD, 
     ARK_TABLE_NUM (pass in an int array of length 2, implicit method first), 
     MAX_ERRFAIL, MAX_NITERS or MAX_CONVFAIL.
     The int return flag IER is 0 if successful, and nonzero otherwise.

     To set various real optional inputs, make the following call:

       CALL FARKSETRIN(KEY, VALUE, IER)

     to set the realtype value VAL to the optional input specified by the 
     quoted character string KEY.  KEY must one of the following: 
     INIT_STEP, MAX_STEP, MIN_STEP, STOP_TIME or NLCONV_COEF.

     To reset all optional inputs to their default values, make the following 
     call:

       CALL FARKSETDEFAULTS(IER)

     To set a custom explicit Runge-Kutta table, make the following call:

       CALL FARKSETERKTABLE(S, Q, P, C, A, B, B2, IER)

     The arguments are:
       S = the number of stages in the table [int, input]
       Q = the global order of accuracy of the method [int, input]
       P = the global order of accuracy of the embedding [int, input]
       C = array of length S containing the stage times [realtype, input]
       A = array of length S*S containing the ERK coefficients (stored in 
           row-major, "C", order) [realtype, input]
       B = array of length S containing the solution coefficients
           [realtype, input]
       B2 = array of length S containing the embedding coefficients
           [realtype, input]

     To set a custom diagonally-implicit Runge-Kutta table, make the following
     call:

       CALL FARKSETIRKTABLE(S, Q, P, C, A, B, B2, IER)

     The arguments are:
       S = the number of stages in the table [int, input]
       Q = the global order of accuracy of the method [int, input]
       P = the global order of accuracy of the embedding [int, input]
       C = array of length S containing the stage times [realtype, input]
       A = array of length S*S containing the DIRK coefficients (stored in 
           row-major, "C", order) [realtype, input]
       B = array of length S containing the solution coefficients
           [realtype, input]
       B2 = array of length S containing the embedding coefficients
           [realtype, input]

     To set a custom additive Runge-Kutta table, make the following call:

       CALL FARKSETARKTABLES(S, Q, P, C, AI, AE, B, B2, IER)

     The arguments are:
       S = the number of stages in the table [int, input]
       Q = the global order of accuracy of the method [int, input]
       P = the global order of accuracy of the embedding [int, input]
       C = array of length S containing the stage times [realtype, input]
       AI = array of length S*S containing the DIRK coefficients (stored in 
           row-major, "C", order) [realtype, input]
       AE = array of length S*S containing the ERK coefficients (stored in 
           row-major, "C", order) [realtype, input]
       B = array of length S containing the solution coefficients
           [realtype, input]
       B2 = array of length S containing the embedding coefficients
           [realtype, input]

     To set the time step adaptivity algorithm (and it's associated 
     parameters -- ARKodeSetAdaptivityMethod), make the following call:

       CALL FARKSETADAPTIVITYMETHOD(METHOD, PARAMS, IER)

     The arguments are:
       METHOD = integer between 0 and 5 to specify the method [int, input]
       PARAMS = array of length 9 containing the adaptivity parameters 
           [realtype, input]

     To set the time step change heuristics (see 
     ARKodeSetAdaptivityConstants), make the following call:

       CALL FARKSETADAPTIVITYCONSTANTS(ETAMX1, ETAMXF, ETACF, SMALLNEF, IER)

     The arguments are:
       ETAMX1 = maximum step change for the first step [realtype, input]
       ETAMXF = step change on an error failure [realtype, input]
       ETACF = step change on a convergence failure [realtype, input]
       SMALLNEF = number of error failures before ETAMXF is enforced 
           [int, input]

     To set the Newton convergence heuristics (see ARKodeSetNewtonConstants), 
     make the following call:

       CALL FARKSETNEWTONCONSTANTS(CRDOWN, RDIV, IER)

     The arguments are:
       CRDOWN = convergence rate estimation constant [realtype, input]
       RDIV = divergence bound [realtype, input]

     To set the heuristics governing when to rebuild the Jacobian or
     preconditioner (see ARKodeSetLSetupConstants), make the following call: 

       CALL FARKSETLSETUPCONSTANTS(DGMAX, MSBP, IER)

     The arguments are:
       DGMAX = maximum allowable gamma ratio [realtype, input]
       MSBP = maximum number of time steps between lsetup calls [int, input]

 -----------------------------------------------------------------------------

 (9) Specification of linear system solution method.

     In the case of using either an implicit or ImEx method, the solution of 
     each Runge-Kutta stage may involve the solution of linear systems related 
     to the Jacobian J = dfi(t,y)/dy of the implicit portion of the ODE 
     system. ARKode presently includes seven choices for the treatment of
     these systems, and the user of FARKODE must call a routine with a
     specific name to make the desired choice. 
 
 (9.1s) DENSE treatment of the linear system.

     The user must make the call

       CALL FARKDENSE(NEQ, IER)

     The arguments are:
        NEQ = the problem size [long int; input]
	IER = error return flag [int, output]: 
	         0 = success, 
		 negative = error.
 
     If the user program includes the FARKDJAC routine for the evaluation of 
     the dense approximation to the Jacobian, the following call must be made

       CALL FARKDENSESETJAC(FLAG, IER)

     with the int FLAG=1 to specify that FARKDJAC is provided (FLAG=0 
     specifies to use the internal finite differences approximation to the
     Jacobian).  The int return flag IER=0 if successful, and nonzero 
     otherwise.
 
     Optional outputs specific to the DENSE case are:
        LENRWLS = IOUT(14) from ARKDlsGetWorkSpace (realtype space)
        LENIWLS = IOUT(15) from ARKDlsGetWorkSpace (integer space)
        LSTF    = IOUT(16) from ARKDlsGetLastFlag
        NFELS   = IOUT(17) from ARKDlsGetNumRhsEvals
        NJED    = IOUT(18) from ARKDlsGetNumJacEvals
     See the ARKODE manual for descriptions.
 
 (9.2s) BAND treatment of the linear system

     The user must make the call

       CALL FARKBAND(NEQ, MU, ML, IER)

     The arguments are:
        NEQ = problem size [long int, input]
	MU = upper bandwidth [long int, input]
	ML = lower bandwidth [long int, input]
        IER = return flag [int, output]; 0 if successful, nonzero otherwise.
 
     If the user program includes the FARKBJAC routine for the evaluation of 
     the band approximation to the Jacobian, the following call must be made

       CALL FARKBANDSETJAC(FLAG, IER)

     with the int FLAG=1 to specify that FARKBJAC is provided (FLAG=0 
     specifies to use the internal finite differences approximation to the 
     Jacobian).  The int return flag IER=0 if successful, nonzero otherwise.
 
     Optional outputs specific to the BAND case are:
        LENRWLS = IOUT(14) from ARKDlsGetWorkSpace (realtype space)
        LENIWLS = IOUT(15) from ARKDlsGetWorkSpace (integer space)
        LSTF    = IOUT(16) from ARKDlsGetLastFlag
        NFELS   = IOUT(17) from ARKDlsGetNumRhsEvals
        NJED    = IOUT(18) from ARKDlsGetNumJacEvals
     See the ARKODE manual for descriptions.

 (9.3s) LAPACK dense treatment of the linear system

     The user must make the call

       CALL FARKLAPACKDENSE(NEQ, IER)

     The arguments match those for FARKDENSE, except that NEQ is now a 
     normal int (and not a long int).

     The user may optionally call

       CALL FARKLAPACKDENSESETJAC(FLAG, IER)
       
     with the int FLAG=1 if the user provides the function FARKLDJAC. 

     The optional outputs when using FARKLAPACKDENSE match those from 
     FARKDENSE.

 (9.4s) LAPACK band treatment of the linear system

     The user must make the call

       CALL FARKLAPACKBAND(NEQ, MUPPER, MLOWER, IER)

     The arguments match those for FARKBAND, except that now all arguments 
     have type 'int'.

     The user may optionally call

       CALL FARKLAPACKBANDSETJAC(FLAG, IER)

     with the int FLAG=1 if the user provides the function FARKLBJAC. 

     The optional outputs when using FARKLAPACKBAND match those from FARKBAND.

 (9.5) SPGMR treatment of the linear systems.

     For the Scaled Preconditioned GMRES solution of the linear systems,
     the user must make the following call:

       CALL FARKSPGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)              

     The arguments are:
        IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
	IGSTYPE = Gram-schmidt process type [int, input]: 
              1 = modified G-S
              2 = classical G-S.
	MAXL = maximum Krylov subspace dimension [int; input]; 
	      0 = default
	DELT = linear convergence tolerance factor [realtype, input]; 
	      0.0 = default.
	IER = error return flag [int, output]: 
	       0 = success; 
	      <0 = an error occured
 
     Optional outputs specific to the SPGMR case are:
        LENRWLS = IOUT(14) from ARKSpilsGetWorkSpace
        LENIWLS = IOUT(15) from ARKSpilsGetWorkSpace
        LSTF    = IOUT(16) from ARKSpilsGetLastFlag
        NFELS   = IOUT(17) from ARKSpilsGetNumRhsEvals
        NJTV    = IOUT(18) from ARKSpilsGetNumJtimesEvals
        NPE     = IOUT(19) from ARKSpilsGetNumPrecEvals
        NPS     = IOUT(20) from ARKSpilsGetNumPrecSolves
        NLI     = IOUT(21) from ARKSpilsGetNumLinIters
        NCFL    = IOUT(22) from ARKSpilsGetNumConvFails
     See the ARKODE manual for descriptions.
 
     If a sequence of problems of the same size is being solved using the
     SPGMR linear solver, then following the call to FARKREINIT, a call to 
     the FARKSPGMRREINIT routine is needed if any of IPRETYPE, IGSTYPE, DELT 
     is being changed.  In that case, call FARKSPGMRREINIT as follows:

       CALL FARKSPGMRREINIT(IPRETYPE, IGSTYPE, DELT, IER)

     The arguments have the same meanings as for FARKSPGMR.  If MAXL is being
     changed, then the user should call FARKSPGMR instead.
 
 (9.6) SPBCG treatment of the linear systems.

     For the Scaled Preconditioned Bi-CGSTAB solution of the linear systems,
     the user must make the following call:

       CALL FARKSPBCG(IPRETYPE, MAXL, DELT, IER)              

     The arguments are:
       IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
       MAXL = maximum Krylov subspace dimension [int, input]; 0 = default.
       DELT = linear convergence tolerance factor [realtype, input]; 
              0.0 = default.
       IER = error return flag [int, output]: 
              0 = success; 
	     <0 = an error occured
 
     Optional outputs specific to the SPBCG case are:
        LENRWLS = IOUT(14) from ARKSpilsGetWorkSpace
        LENIWLS = IOUT(15) from ARKSpilsGetWorkSpace
        LSTF    = IOUT(16) from ARKSpilsGetLastFlag
        NFELS   = IOUT(17) from ARKSpilsGetNumRhsEvals
        NJTV    = IOUT(18) from ARKSpilsGetNumJtimesEvals
        NPE     = IOUT(19) from ARKSpilsGetNumPrecEvals
        NPS     = IOUT(20) from ARKSpilsGetNumPrecSolves
        NLI     = IOUT(21) from ARKSpilsGetNumLinIters
        NCFL    = IOUT(22) from ARKSpilsGetNumConvFails
     See the ARKODE manual for descriptions.
 
     If a sequence of problems of the same size is being solved using the
     SPBCG linear solver, then following the call to FARKREINIT, a call to the
     FARKSPBCGREINIT routine is needed if any of its arguments is being 
     changed.  The call is:

       CALL FARKSPBCGREINIT(IPRETYPE, MAXL, DELT, IER)              

     The arguments have the same meanings as for FARKSPBCG.

 (9.7) SPTFQMR treatment of the linear systems.

     For the Scaled Preconditioned TFQMR solution of the linear systems, the
     user must make the following call:

       CALL FARKSPTFQMR(IPRETYPE, MAXL, DELT, IER)              

     The arguments are:
       IPRETYPE = preconditioner type [int, input]: 
              0 = none 
              1 = left only
              2 = right only
              3 = both sides
       MAXL = maximum Krylov subspace dimension [int, input]; 0 = default.
       DELT = linear convergence tolerance factor [realtype, input]
	      0.0 = default.
       IER = error return flag [int, output]: 
              0 = success; 
	     <0 = an error occured
 
     Optional outputs specific to the SPTFQMR case are:
        LENRWLS = IOUT(14) from ARKSpilsGetWorkSpace
        LENIWLS = IOUT(15) from ARKSpilsGetWorkSpace
        LSTF    = IOUT(16) from ARKSpilsGetLastFlag
        NFELS   = IOUT(17) from ARKSpilsGetNumRhsEvals
        NJTV    = IOUT(18) from ARKSpilsGetNumJtimesEvals
        NPE     = IOUT(19) from ARKSpilsGetNumPrecEvals
        NPS     = IOUT(20) from ARKSpilsGetNumPrecSolves
        NLI     = IOUT(21) from ARKSpilsGetNumLinIters
        NCFL    = IOUT(22) from ARKSpilsGetNumConvFails
     See the ARKODE manual for descriptions.

     If a sequence of problems of the same size is being solved using the
     SPTFQMR linear solver, then following the call to FARKREINIT, a call to 
     the FARKSPTFQMRREINIT routine is needed if any of its arguments is
     being changed.  The call is:

       CALL FARKSPTFQMRREINIT(IPRETYPE, MAXL, DELT, IER)              

     The arguments have the same meanings as for FARKSPTFQMR.

 (9.8) Usage of user-supplied routines for the Krylov solvers

     If the user program includes the FARKJTIMES routine for the evaluation of
     the Jacobian vector product, the following call must be made

       CALL FARKSPILSSETJAC(FLAG, IER)

     with the int FLAG=1 to specify that FARKJTIMES is provided (FLAG=0 
     specifies to use and internal finite difference approximation to this 
     product).  The int return flag IER=0 if successful, and nonzero otherwise.
 
     Usage of the user-supplied routines FARKPSOL and FARKPSET for solution of
     the preconditioner linear system requires the following call:

       CALL FARKSPILSSETPREC(FLAG, IER)

     with the int FLAG=1. The return flag IER=0 if successful, nonzero 
     otherwise.

     The user-supplied routine FARKPSET must have the form:

       SUBROUTINE FARKPSET(T,Y,FY,JOK,JCUR,GAMMA,H,IPAR,RPAR,V1,V2,V3,IER)

     This routine must set up the preconditioner P to be used in the 
     subsequent call to FARKPSOL.  The preconditioner (or the product of the 
     left and right preconditioners if using both) should be an approximation 
     to the matrix  I - GAMMA*J  (I = identity, J = Jacobian),

     The arguments are:
       T = current time [realtype, input]
       Y = current state variable array [realtype, input]
       FY = current state variable derivative array [realtype, input]
       JOK = flag indicating whether Jacobian-related data needs to be 
           recomputed [int, input]:
                  0 = recompute, 
		  1 = reuse with the current value of GAMMA
       JCUR = return flag to denote if Jacobian data was recomputed
           [realtype, output], 1=yes, 0=no
       GAMMA = Jacobian scaling factor [realtype, input]
       H = current time step [realtype, input]
       IPAR = array of user integer data [long int, input/output]
       RPAR = array with user real data [realtype, input/output]
       V* -- array containing temporary workspace of same size as Y 
               [realtype, input]
       IER  = return completion flag [int, output]:
                  0 = SUCCESS,
		 >0 = recoverable failure
                 <0 = non-recoverable failure

     The user-supplied routine FARKPSOL must have the form:

       SUBROUTINE FARKPSOL(T,Y,FY,R,Z,GAMMA,DELTA,LR,IPAR,RPAR,VT,IER)

       DIMENSION Y(*), FY(*), VT(*), R(*), Z(*), IPAR(*), RPAR(*)

     Typically this routine will use only NEQ, T, Y, GAMMA, R, LR, and Z.  It
     must solve the preconditioner linear system Pz = r.  The preconditioner
     (or the product of the left and right preconditioners if both are 
     nontrivial) should be an approximation to the matrix  I - GAMMA*J  
     (I = identity, J = Jacobian).

     The arguments are:
       T = current time [realtype, input]
       Y = current state variable array [realtype, input]
       FY = current state variable derivative array [realtype, input]
       R = right-hand side array [realtype, input]
       Z = solution array [realtype, output]
       GAMMA = Jacobian scaling factor [realtype, input]
       DELTA = desired residual tolerance [realtype, input]
       LR = flag denoting to solve the right or left preconditioner system
                  1 = left preconditioner
		  2 = right preconditioner
       IPAR = array of user integer data [long int, input/output]
       RPAR = array with user real data [realtype, input/output]
       VT -- array containing temporary workspace of same size as Y 
               [realtype, input]
       IER  = return completion flag [int, output]:
                  0 = SUCCESS,
		 >0 = recoverable failure
                 <0 = non-recoverable failure

 -----------------------------------------------------------------------------

 (10) The integrator: FARKODE

     Carrying out the integration is accomplished by making calls as follows:

       CALL FARKODE(TOUT, T, Y, ITASK, IER)

     The arguments are:
       TOUT = next value of t at which a solution is desired [realtype, input]
       T = value of t reached by the solver [realtype, output]
       Y = array containing state variables on output (realtype, output]
       ITASK = task indicator [int, input]:
                  1 = normal mode (overshoot TOUT and interpolate)
		  2 = one-step mode (return after each internal step taken)
		  3 = normal tstop mode (like 1, but integration never 
		      proceeds past TSTOP, which must be specified through a 
		      call to FARKSETRIN using the key 'STOP_TIME')
		  4 = one step tstop (like 2, but integration never goes 
		      past TSTOP)
       IER = completion flag [int, output]: 
                  0 = success, 
		  1 = tstop return, 
		  2 = root return, 
                  values -1 ... -10 are failure modes (see ARKODE manual).
     The current values of the optional outputs are immediately available in
     the IOUT and ROUT arrays.
 
 -----------------------------------------------------------------------------

 (11) Computing solution derivatives: FARKDKY

     To obtain a derivative of the solution, of order up to the method order,
     make the following call:

       CALL FARKDKY(T, K, DKY, IER)

     The arguments are:
       T = time at which solution derivative is desired, within the interval
           [TCUR-HU,TCUR], [realtype, input].
       K = derivative order (0 .le. K .le. QU) [int, input]
       DKY = array containing computed K-th derivative of y [realtype, output]
       IER = return flag [int, output]: 0=success, <0 = illegal argument.
 
 -----------------------------------------------------------------------------

 (12) Get the current weight vector: FARKGETERRWEIGHTS

     To obtain the current weight vector, make the following call:

       CALL FARKGETERRWEIGHTS(EWT, IER)

     The arguments are:
       EWT = array containing the error weight vector [realtype, output]
       IER = return flag [int, output]: 0=success, nonzero if an error.
 
 -----------------------------------------------------------------------------

 (13) Get an estimate of the local error: FARKGETESTLOCALERR

     To obtain the current error estimate vector, make the following call:

       CALL FARKGETESTLOCALERR(ELE, IER)

     The arguments are:
       ELE = array with the estimated local error vector [realtype, output]
       IER = return flag [int, output]: 0=success, nonzero if an error.
 
 -----------------------------------------------------------------------------

 (14) Memory freeing: FARKFREE 

     To free the internal memory created by the calls to FARKMALLOC and
     FNVINITS or FNVINITP, make the call

       CALL FARKFREE()
 
===============================================================*/

#ifndef _FARKODE_H
#define _FARKODE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* header files  */
#include <arkode/arkode.h>
#include <sundials/sundials_direct.h>  /* definition of type DlsMat   */
#include <sundials/sundials_nvector.h> /* definition of type N_Vector */
#include <sundials/sundials_types.h>   /* definition of type realtype */

/* Definitions of interface function names */
#if defined(SUNDIALS_F77_FUNC)

#define FARK_MALLOC              SUNDIALS_F77_FUNC(farkmalloc,              FARKMALLOC)
#define FARK_REINIT              SUNDIALS_F77_FUNC(farkreinit,              FARKREINIT)
#define FARK_SETDEFAULTS         SUNDIALS_F77_FUNC(farksetdefaults,         FARKSETDEFAULTS)
#define FARK_SETIIN              SUNDIALS_F77_FUNC(farksetiin,              FARKSETIIN)
#define FARK_SETRIN              SUNDIALS_F77_FUNC(farksetrin,              FARKSETRIN)
#define FARK_SETERKTABLE         SUNDIALS_F77_FUNC(farkseterktable,         FARKSETERKTABLE)
#define FARK_SETIRKTABLE         SUNDIALS_F77_FUNC(farksetirktable,         FARKSETIRKTABLE)
#define FARK_SETARKTABLES        SUNDIALS_F77_FUNC(farksetarktables,        FARKSETARKTABLES)
#define FARK_EWTSET              SUNDIALS_F77_FUNC(farkewtset,              FARKEWTSET)
#define FARK_DENSE               SUNDIALS_F77_FUNC(farkdense,               FARKDENSE)
#define FARK_DENSESETJAC         SUNDIALS_F77_FUNC(farkdensesetjac,         FARKDENSESETJAC)
#define FARK_BAND                SUNDIALS_F77_FUNC(farkband,                FARKBAND)
#define FARK_BANDSETJAC          SUNDIALS_F77_FUNC(farkbandsetjac,          FARKBANDSETJAC)
#define FARK_LAPACKDENSE         SUNDIALS_F77_FUNC(farklapackdense,         FARKLAPACKDENSE)
#define FARK_LAPACKDENSESETJAC   SUNDIALS_F77_FUNC(farklapackdensesetjac,   FARKLAPACKDENSESETJAC)
#define FARK_LAPACKBAND          SUNDIALS_F77_FUNC(farklapackband,          FARKLAPACKBAND)
#define FARK_LAPACKBANDSETJAC    SUNDIALS_F77_FUNC(farklapackbandsetjac,    FARKLAPACKBANDSETJAC)
#define FARK_SPTFQMR             SUNDIALS_F77_FUNC(farksptfqmr,             FARKSPTFQMR)
#define FARK_SPTFQMRREINIT       SUNDIALS_F77_FUNC(farksptfqmrreinit,       FARKSPTFQMRREINIT)
#define FARK_SPBCG               SUNDIALS_F77_FUNC(farkspbcg,               FARKSPBCG)
#define FARK_SPBCGREINIT         SUNDIALS_F77_FUNC(farkspbcgreinit,         FARKSPBCGREINIT)
#define FARK_SPGMR               SUNDIALS_F77_FUNC(farkspgmr,               FARKSPGMR)
#define FARK_SPGMRREINIT         SUNDIALS_F77_FUNC(farkspgmrreinit,         FARKSPGMRREINIT)
#define FARK_SPILSSETJAC         SUNDIALS_F77_FUNC(farkspilssetjac,         FARKSPILSSETJAC)
#define FARK_SPILSSETPREC        SUNDIALS_F77_FUNC(farkspilssetprec,        FARKSPILSSETPREC)
#define FARK_ARKODE              SUNDIALS_F77_FUNC(farkode,                 FARKODE)
#define FARK_DKY                 SUNDIALS_F77_FUNC(farkdky,                 FARKDKY)
#define FARK_FREE                SUNDIALS_F77_FUNC(farkfree,                FARKFREE)
#define FARK_IMP_FUN             SUNDIALS_F77_FUNC(farkifun,                FARKIFUN)
#define FARK_EXP_FUN             SUNDIALS_F77_FUNC(farkefun,                FARKEFUN)
#define FARK_DJAC                SUNDIALS_F77_FUNC(farkdjac,                FARKDJAC)
#define FARK_LDJAC               SUNDIALS_F77_FUNC(farkldjac,               FARKLDJAC)
#define FARK_BJAC                SUNDIALS_F77_FUNC(farkbjac,                FARKBJAC)
#define FARK_PSOL                SUNDIALS_F77_FUNC(farkpsol,                FARKPSOL)
#define FARK_PSET                SUNDIALS_F77_FUNC(farkpset,                FARKPSET)
#define FARK_JTIMES              SUNDIALS_F77_FUNC(farkjtimes,              FARKJTIMES)
#define FARK_EWT                 SUNDIALS_F77_FUNC(farkewt,                 FARKEWT)
#define FARK_GETERRWEIGHTS       SUNDIALS_F77_FUNC(farkgeterrweights,       FARKGETERRWEIGHTS)
#define FARK_GETESTLOCALERR      SUNDIALS_F77_FUNC(farkgetestlocalerr,      FARKGETESTLOCALERR)
#define FARK_WRITEPARAMETERS     SUNDIALS_F77_FUNC(farkwriteparameters,     FARKWRITEPARAMETERS)

#else

#define FARK_MALLOC              farkmalloc_
#define FARK_REINIT              farkreinit_
#define FARK_SETDEFAULTS         farksetdefaults_
#define FARK_SETIIN              farksetiin_
#define FARK_SETRIN              farksetrin_
#define FARK_SETERKTABLE         farkseterktable_
#define FARK_SETIRKTABLE         farksetirktable_
#define FARK_SETARKTABLES        farksetarktables_
#define FARK_EWTSET              farkewtset_
#define FARK_DENSE               farkdense_
#define FARK_DENSESETJAC         farkdensesetjac_
#define FARK_BAND                farkband_
#define FARK_BANDSETJAC          farkbandsetjac_
#define FARK_LAPACKDENSE         farklapackdense_
#define FARK_LAPACKDENSESETJAC   farklapackdensesetjac_
#define FARK_LAPACKBAND          farklapackband_
#define FARK_LAPACKBANDSETJAC    farklapackbandsetjac_
#define FARK_SPTFQMR             farksptfqmr_
#define FARK_SPTFQMRREINIT       farksptfqmrreinit_
#define FARK_SPBCG               farkspbcg_
#define FARK_SPBCGREINIT         farkspbcgreinit_
#define FARK_SPGMR               farkspgmr_
#define FARK_SPGMRREINIT         farkspgmrreinit_
#define FARK_SPILSSETJAC         farkspilssetjac_
#define FARK_SPILSSETPREC        farkspilssetprec_
#define FARK_ARKODE              farkode_
#define FARK_DKY                 farkdky_
#define FARK_FREE                farkfree_
#define FARK_IMP_FUN             farkifun_
#define FARK_EXP_FUN             farkefun_
#define FARK_DJAC                farkdjac_
#define FARK_LDJAC               farkldjac_
#define FARK_BJAC                farkbjac_
#define FARK_PSOL                farkpsol_
#define FARK_PSET                farkpset_
#define FARK_JTIMES              farkjtimes_
#define FARK_EWT                 farkewt_
#define FARK_GETERRWEIGHTS       farkgeterrweights_
#define FARK_GETESTLOCALERR      farkgetestlocalerr_
#define FARK_WRITEPARAMETERS     farkwriteparameters_

#endif

  /* Type for user data */
  typedef struct {
    realtype *rpar;
    long int *ipar;
  } *FARKUserData;
  
  /* Prototypes of exported functions */
  void FARK_MALLOC(realtype *t0, realtype *y0, int *imex, 
		   int *iatol, realtype *rtol, realtype *atol, 
		   long int *iout, realtype *rout, 
		   long int *ipar, realtype *rpar, int *ier);

  void FARK_REINIT(realtype *t0, realtype *y0, int *imex,
		   int *iatol, realtype *rtol, realtype *atol,
		   int *ier);

  void FARK_SETIIN(char key_name[], long int *ival, int *ier);
  void FARK_SETRIN(char key_name[], realtype *rval, int *ier);
  void FARK_SETDEFAULTS(int *ier);
  void FARK_SETERKTABLE(int *s, int *q, int *p, realtype *c, realtype *A, 
			realtype *b, realtype *b2, int *ier);
  void FARK_SETIRKTABLE(int *s, int *q, int *p, realtype *c, 
			realtype *A, realtype *b, realtype *b2, int *ier);
  void FARK_SETARKTABLES(int *s, int *q, int *p, realtype *c, realtype *Ai, 
			 realtype *Ae, realtype *b, realtype *b2, int *ier);

  void FARK_EWTSET(int *flag, int *ier);

  void FARK_DENSE(long int *neq, int *ier);
  void FARK_DENSESETJAC(int *flag, int *ier);

  void FARK_BAND(long int *neq, long int *mupper, long int *mlower, int *ier);
  void FARK_BANDSETJAC(int *flag, int *ier);

  void FARK_LAPACKDENSE(int *neq, int *ier);
  void FARK_LAPACKDENSESETJAC(int *flag, int *ier);
  void FARK_LAPACKBAND(int *neq, int *mupper, int *mlower, int *ier);
  void FARK_LAPACKBANDSETJAC(int *flag, int *ier);

  void FARK_SPGMR(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier);
  void FARK_SPGMRREINIT(int *pretype, int *gstype, realtype *delt, int *ier);

  void FARK_SPBCG(int *pretype, int *maxl, realtype *delt, int *ier);
  void FARK_SPBCGREINIT(int *pretype, int *maxl, realtype *delt, int *ier);

  void FARK_SPTFQMR(int *pretype, int *maxl, realtype *delt, int *ier);
  void FARK_SPTFQMRREINIT(int *pretype, int *maxl, realtype *delt, int *ier);

  void FARK_SPILSSETJAC(int *flag, int *ier);
  void FARK_SPILSSETPREC(int *flag, int *ier);
  
  void FARK_ARKODE(realtype *tout, realtype *t, realtype *y, int *itask, int *ier);

  void FARK_DKY(realtype *t, int *k, realtype *dky, int *ier);

  void FARK_GETERRWEIGHTS(realtype *eweight, int *ier);
  void FARK_GETESTLOCALERR(realtype *ele, int *ier);

  void FARK_FREE(void);

  void FARK_WRITEPARAMETERS(int *ier);

  /* Prototypes: Functions Called by the ARKODE Solver */
  int FARKfi(realtype t, N_Vector y, N_Vector ydot, void *user_data);
  int FARKfe(realtype t, N_Vector y, N_Vector ydot, void *user_data);
  
  int FARKDenseJac(long int N, realtype t, 
		   N_Vector y, N_Vector fy, 
		   DlsMat J, void *user_data,
		   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKBandJac(long int N, long int mupper, long int mlower,
		  realtype t, N_Vector y, N_Vector fy,
		  DlsMat J, void *user_data,
		  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKLapackDenseJac(long int N, realtype t,
			 N_Vector y, N_Vector fy, 
			 DlsMat Jac, void *user_data,
			 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  int FARKLapackBandJac(long int N, long int mupper, long int mlower,
			realtype t, N_Vector y, N_Vector fy, 
			DlsMat Jac, void *user_data,
			N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  int FARKPSet(realtype tn, N_Vector y,N_Vector fy, booleantype jok,
	       booleantype *jcurPtr, realtype gamma, void *user_data,
	       N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FARKPSol(realtype tn, N_Vector y, N_Vector fy, 
	       N_Vector r, N_Vector z,
	       realtype gamma, realtype delta,
	       int lr, void *user_data, N_Vector vtemp);
  
  int FARKJtimes(N_Vector v, N_Vector Jv, realtype t, 
		 N_Vector y, N_Vector fy,
		 void *user_data, N_Vector work);
  
  int FARKEwtSet(N_Vector y, N_Vector ewt, void *user_data);

  /* Declarations for global variables shared amongst various routines */
  extern N_Vector F2C_ARKODE_vec;   /* defined in FNVECTOR module */

  extern void *ARK_arkodemem;       /* defined in farkode.c */
  extern long int *ARK_iout;        /* defined in farkode.c */
  extern realtype *ARK_rout;        /* defined in farkode.c */
  extern int ARK_nrtfn;             /* defined in farkode.c */
  extern int ARK_ls;                /* defined in farkode.c */

  /* Linear solver IDs */
  enum { ARK_LS_DENSE       = 1, 
	 ARK_LS_BAND        = 2, 
         ARK_LS_LAPACKDENSE = 3, 
	 ARK_LS_LAPACKBAND  = 4,
	 ARK_LS_SPGMR       = 5, 
	 ARK_LS_SPBCG       = 6, 
	 ARK_LS_SPTFQMR     = 7 };

#ifdef __cplusplus
}
#endif

#endif

/*===============================================================
   EOF
===============================================================*/
