!-----------------------------------------------------------------
! $Revision: $
! $Date: $
!-----------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ SMU
!-----------------------------------------------------------------
! Example problem:
! 
! This program solves the Fortran ODE test problem defined in the 
! file rober.f, using the FARKODE interface for the ARKode ODE 
! solver module.
! 
! Based on the inputs in the file fsolve_params.txt, this program
! will switch between ERK, DIRK and ARK solvers.  For DIRK and ARK
! methods, the implicit systems are solved with a modified Newton
! iteration with the ARKDENSE dense linear solver.  The Jacobian 
! routine and right-hand side routines come from the file 
! user-supplied Jacobian routine.
!
! Output is printed 10 times throughout the defined time interval.
! Run statistics (optional outputs) are printed at the end.
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Main driver program
!-----------------------------------------------------------------
program driver

  ! Declarations
  implicit none

  ! general problem variables
  real*8    :: T0, Tf, dTout, Tout, Tcur, rtol, rout(6), t(2), h0, apar(3)
  integer   :: it, Nt, ier, neqn, ndisc, mlmas, mumas, mljac, mujac, idef
  integer   :: ind(1)
  integer*8 :: NEQ, MU, ML, iout(22), btable2(2)
  real*8, allocatable :: y(:), ytrue(:), rtols(:), atols(:)
  logical   :: numjac, nummas, consis, tolvec, denseout
  character :: fullnm*17, problm*5, type*3

  ! real/integer parameters to pass through to supplied functions
  !    ipar(1) -> problem size
  integer :: ipar(1)
  real*8  :: rpar(1)

  ! solver parameters
  integer*8 :: order, dense_order, imex, btable, adapt_method, &
       small_nef, predictor, msbp, maxcor, pq, fixedpt, m_aa
  real*8 :: cflfac, safety, bias, growth, hfixed_lb, hfixed_ub, &
       k1, k2, k3, etamx1, etamxf, etacf, crdown, rdiv, dgmax, nlscoef

  namelist /inputs/ order, dense_order, imex, btable, adapt_method, &
       cflfac, safety, bias, growth, hfixed_lb, hfixed_ub, pq, k1,  &
       k2, k3, etamx1, etamxf, etacf, small_nef, crdown, rdiv,      &
       dgmax, predictor, msbp, fixedpt, m_aa, maxcor, nlscoef

  !-----------------------
  ! read solver inputs
  open(100,file='fsolve_params.txt',form='formatted')
  read(100,inputs)
  close(100)
  denseout = .true.
  if (dense_order == -1)  denseout = .false.

  ! quit for ImEx or Explicit, since this requires fully implicit solve
  if (imex /= 0) then
     print *, 'ark_rober.f90 error: problem requires fully implicit solver'
     stop
  end if

  ! call problem setup routine, store relevant details
  call prob(fullnm, problm, type, neqn, ndisc, t, numjac, &
            mljac, mujac, nummas, mlmas, mumas, ind)
  print *,'full problem name: ',fullnm
  print *,'short problem name: ',problm
  print *,'problem type: ',type
  T0 = t(1)
  Tf = t(2)
  dTout = (Tf-T0)/100.d0
  Nt = Tf/dTout + 0.5
  NEQ = neqn
  MU = mujac
  ML = mljac
  ipar(1) = NEQ
  allocate(y(neqn), ytrue(neqn), rtols(neqn), atols(neqn))

  ! set initial conditions according to problem specifications
  call init(neqn, T0, y, ytrue, consis)

  ! set tolerances according to problem specifications
  rtols = 0.d0
  atols = 0.d0
  call settolerances(neqn, rtols, atols, tolvec)
!!$  if (.not. tolvec)  atols = atols(1)
!!$  rtol = rtols(1)
  rtol  = 1.d-4
  atols = 1.d-4*rtol
  h0    = 1.d-4*rtol
  
  ! initialize vector module
  call FNVInitS(4, NEQ, ier)

  ! initialize ARKode solver
  call FARKMalloc(T0, y, imex, 2, rtol, atols, &
                  iout, rout, ipar, rpar, ier)

  ! set diagnostics output file
  call FARKSetDiagnostics("diags_ark_rober.txt", 19, ier);
  if (ier /= 0) then
     print *, 'FARKSetDiagnostics error = ',ier
     stop
  endif
  
  ! set optional inputs
  if (order /= 0) then
     call FARKSetIin('ORDER', order, ier)
  else if (btable /= -1) then
     if (imex == 0) then
        call FARKSetIin('IRK_TABLE_NUM', btable, ier)
     else if (imex == 1) then
        call FARKSetIin('ERK_TABLE_NUM', btable, ier)
     else 
        if (btable == 2)  btable2 = (/ btable, 15_8 /)
        if (btable == 4)  btable2 = (/ btable, 20_8 /)
        if (btable == 9)  btable2 = (/ btable, 22_8 /)
        call FARKSetIin('ARK_TABLE_NUM', btable2, ier)
     end if
  end if
  call FARKSetIin('DENSE_ORDER', dense_order, ier)
  if (imex == 0) then
     call FARKSetIin('IMPLICIT', 1, ier)
  else if (imex == 1) then
     call FARKSetIin('EXPLICIT', 1, ier)
  else
     call FARKSetIin('IMEX', 1, ier)
  end if
  apar = (/ k1, k2, k3 /)
  if (sum(abs(apar)) > 0.d0) then
     idef = 0
  else
     idef = 1
  endif
  call FARKSetAdaptivityMethod(adapt_method, idef, pq, apar, ier)
  call FARKSetRin('ADAPT_CFL',       cflfac, ier)
  call FARKSetRin('ADAPT_SAFETY',    safety, ier)
  call FARKSetRin('ADAPT_BIAS',      bias, ier)
  call FARKSetRin('ADAPT_GROWTH',    growth, ier)
  apar = (/ hfixed_lb, hfixed_ub, 0.d0 /)
  call FARKSetRin('ADAPT_BOUNDS',    apar, ier)
  call FARKSetRin('ADAPT_ETAMX1',    etamx1, ier)
  call FARKSetRin('ADAPT_ETAMXF',    etamxf, ier)
  call FARKSetRin('ADAPT_ETACF',     etacf, ier)
  call FARKSetIin('ADAPT_SMALL_NEF', small_nef, ier)
  call FARKSetRin('NONLIN_CRDOWN',   crdown, ier)
  call FARKSetRin('NONLIN_RDIV',     rdiv, ier)
  call FARKSetRin('LSETUP_DGMAX',    dgmax, ier)
  call FARKSetIin('LSETUP_MSBP',     msbp, ier)
  call FARKSetIin('PREDICT_METHOD',  predictor, ier)
  if (fixedpt > 0)  &
       call FARKSetIin('FIXEDPOINT', m_aa, ier)

  ! the following seem required to run this problem successfully
  call FARKSetIin('MAX_NITERS',      5, ier)
  call FARKSetRin('NLCONV_COEF',     1.d-6, ier)
!  call FARKSetIin('MAX_NSTEPS',      1000, ier)
  call FARKSetRin('INIT_STEP',       h0, ier)
  call FARKSetIin('MAX_ERRFAIL',     20, ier)
  
  ! output solver parameters to screen
  call FARKWriteParameters(ier)

  ! specify use of dense linear solver
  call FARKDense(NEQ, ier)
  call FARKDenseSetJac(1, ier)

  ! loop over time outputs
  Tout = T0
  Tcur = T0
  print *, '     t            y1           y2           y3'
  print *, '  ----------------------------------------------------'
  print '(3x,4(es12.5,1x))', Tcur, y
  do it = 1,Nt

     ! set next output time
     Tout = min(Tout + dTout, Tf)

     ! check whether we're interpolating output or not
     if (.not. denseout)  call FARKSetRin('STOP_TIME', Tout, ier)

     ! call solver
     call FARKode(Tout, Tcur, y, 1, ier)
     if (ier < 0) then
        write(0,*) 'Solver failure, stopping integration'
        stop
     end if

     ! output current solution information
     print '(3x,4(es12.5,1x))', Tcur, y

  end do
  print *, '  ----------------------------------------------------'

  ! output solver statistics
  print *, '  '
  print *, 'Final Solver Statistics:'
  print *, '   Internal solver steps = ', iout(3), &
       ' (attempted = ', iout(6), ')'
  print *, '   Total RHS evals:  Fe = ', iout(7), &
       '  Fi = ', iout(8)
  print *, '   Total linear solver setups = ', iout(9)
  print *, '   Total RHS evals for setting up the linear system = ', iout(17)
  print *, '   Total number of Jacobian evaluations = ', iout(18)
  print *, '   Total number of nonlinear iterations = ', iout(11)
  print *, '   Total number of nonlinear solver convergence failures = ', &
       iout(12)
  print *, '   Total number of error test failures = ', iout(10)
  print *, '   First step = ', rout(1)
  print *, '   Last step = ', rout(2)
  print *, '   Current step = ', rout(3)
  print *, '  '

  ! check final solution against reference values for problem
  call solut(neqn, Tf, ytrue)
  print *, '     y(Tf) = ', y
  print *, '  yref(Tf) = ', ytrue
  print *, 'Error: max = ', maxval(abs(ytrue-y)), &
       '  rms = ', sqrt(sum((ytrue-y)**2)/NEQ)
  print *, ' Oversolve = ', rtol/sqrt(sum((ytrue-y)**2)/NEQ+1.d-20)
  print *, '  '

  ! clean up
  call FARKStopDiagnostics(ier);
  call FARKFree()
  deallocate(y, ytrue, rtols, atols)

end program driver
!-----------------------------------------------------------------



!-----------------------------------------------------------------
! Required subroutines for FARKODE interface
!-----------------------------------------------------------------


subroutine farkifun(t, y, ydot, ipar, rpar, ier)
!-----------------------------------------------------------------
! Implicit portion of the right-hand side of the ODE system
!-----------------------------------------------------------------

  ! Declarations
  implicit none

  ! Arguments
  real*8,  intent(in)  :: t, rpar(1)
  integer, intent(in)  :: ipar(2)
  integer, intent(out) :: ier
  real*8,  intent(in)  :: y(ipar(1))
  real*8,  intent(out) :: ydot(ipar(1))

  ! temporary variables
  real*8 :: tmp(ipar(1))

  ! fully implicit solver only, so just call feval
  call feval(ipar(1), t, y, tmp, ydot, ier, rpar, ipar)

  ier = 0
  
end subroutine farkifun
!-----------------------------------------------------------------


subroutine farkefun(t, y, ydot, ipar, rpar, ier)
!-----------------------------------------------------------------
! Explicit portion of the right-hand side of the ODE system
!-----------------------------------------------------------------

  ! Declarations
  implicit none

  ! Arguments
  real*8,  intent(in)  :: t, rpar(1)
  integer, intent(in)  :: ipar(2)
  integer, intent(out) :: ier
  real*8,  intent(in)  :: y(ipar(1))
  real*8,  intent(out) :: ydot(ipar(1))

  ! temporary variables
  real*8 :: tmp(ipar(1))

  ! fully implicit solver only, return error flag
  ydot = 0.d0
  ier = -1
  
end subroutine farkefun
!-----------------------------------------------------------------


subroutine farkdjac(neq,t,y,fy,DJac,h,ipar,rpar,wk1,wk2,wk3,ier)
!-----------------------------------------------------------------
! Jacobian computation routine
!-----------------------------------------------------------------

  ! Declarations
  implicit none

  ! Arguments
  real*8,  intent(in)  :: t, h, rpar(1)
  integer, intent(in)  :: neq, ipar(2)
  integer, intent(out) :: ier
  real*8,  intent(in), dimension(neq) :: y, fy, wk1, wk2, wk3
  real*8,  intent(out) :: DJac(neq,neq)

  ! fully implicit only, call jeval
  call jeval(ipar(1), ipar(1), t, y, fy, DJac, ier, rpar, ipar)
  ier = 0
  
  
end subroutine farkdjac
!-----------------------------------------------------------------
