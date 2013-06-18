c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Ring Modulator (ODE case)
c        ODE of dimension 15
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/ringmod.f
c
c     This is revision
c     $Id: ringmod.F,v 1.2 2006/10/02 10:29:14 testset Exp $
c
c
c     NOTE: Daniel R. Reynolds, January 2013
c     This file has been modified from the original to create 
c     operator-split versions of the feval and jeval routines, named 
c     feval_e, feval_i and jeval_i corresponding to the explicit and 
c     implicit components of the right-hand side function and the 
c     Jacobian of the implicit RHS function.
c     
c     In addition, the prob() routine has been modified to explicitly
c     specify the length of the input arrays to remove all possibility
c     of segmentation faults.
c
c     All other routines remain untouched.
c
c-----------------------------------------------------------------------
      integer function pidate()
      pidate = 20060828
      return
      end
c-----------------------------------------------------------------------
      subroutine prob(fullnm,problm,type,
     +                neqn,ndisc,t,
     +                numjac,mljac,mujac,
     +                nummas,mlmas,mumas,
     +                ind)
      character fullnm*14, problm*7, type*3
      integer neqn,ndisc,mljac,mujac,mlmas,mumas,ind(1)
      double precision t(0:1)
      logical numjac, nummas

      fullnm = 'Ring Modulator'
      problm = 'ringmod'
      type   = 'ODE'
      neqn   = 15
      ndisc  = 0
      t(0)   = 0d0
      t(1)   = 1d-3
      numjac = .false.
      mljac  = neqn
      mujac  = neqn

      return
      end
c-----------------------------------------------------------------------
      subroutine init(neqn,t,y,yprime,consis)
      integer neqn
      double precision t,y(neqn),yprime(neqn)
      logical consis

      integer i

      do 10 i=1,neqn
         y(i) = 0d0
   10 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine settolerances(neqn,rtol,atol,tolvec)
      integer neqn 
      logical tolvec
      double precision rtol(neqn), atol(neqn)
       
      tolvec  = .false.
      

      return
      end
c-----------------------------------------------------------------------
      subroutine setoutput(neqn,solref,printsolout,
     +                    nindsol,indsol)

      logical solref, printsolout
      integer neqn, nindsol
      integer indsol(neqn)

c the reference solution is available
      solref = .true.  

c output file is required
      printsolout = .true.

c default values if printsolout is .true.
      nindsol = neqn
c only nindsol component of indsol are referenced
      do i=1,nindsol
          indsol(i) = i
      end do  

  

      return
      end


c-----------------------------------------------------------------------
      subroutine feval(neqn,t,y,yprime,f,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)

      double precision c,cs,cp,r,rp,lh,ls1,ls2,ls3,rg1,rg2,rg3,ri,rc,
     +                 gamma,delta,pi,
     +                 uin1,uin2,ud1,ud2,ud3,ud4,qud1,qud2,qud3,qud4
      parameter (c=1.6d-8,cs=2d-12,cp=1d-8,r=25d3,rp=50d0,
     +           lh=4.45d0,ls1=2d-3,ls2=5d-4,ls3=5d-4,
     +           rg1=36.3d0,rg2=17.3d0,rg3=17.3d0,ri=5d1,rc=6d2,
     +           gamma=40.67286402d-9,delta=17.7493332d0,
     +           pi=3.141592653589793238462643383d0)

      uin1   = 0.5d0*sin(2d3*pi*t)
      uin2   = 2d0*sin(2d4*pi*t)
      ud1    = +y(3)-y(5)-y(7)-uin2
      ud2    = -y(4)+y(6)-y(7)-uin2
      ud3    = +y(4)+y(5)+y(7)+uin2
      ud4    = -y(3)-y(6)+y(7)+uin2

c     prevent overflow
c     (native NEC SX double precision .le. 1d75)
c      if (delta*max(ud1,ud2,ud3,ud4).gt.172d0) then
c         ierr = -1
c         return
c      endif
c
c      (double precisione ieee .le. 1d304)
      if (delta*max(ud1,ud2,ud3,ud4).gt.300d0) then
         ierr = -1
         return
      endif

      qud1   = gamma*(exp(delta*ud1)-1d0)
      qud2   = gamma*(exp(delta*ud2)-1d0)
      qud3   = gamma*(exp(delta*ud3)-1d0)
      qud4   = gamma*(exp(delta*ud4)-1d0)

      f(1)  = (y(8)-0.5d0*y(10)+0.5d0*y(11)+y(14)-y(1)/r)/c
      f(2)  = (y(9)-0.5d0*y(12)+0.5d0*y(13)+y(15)-y(2)/r)/c
      f(3)  = (y(10)-qud1+qud4)/cs
      f(4)  = (-y(11)+qud2-qud3)/cs
      f(5)  = (y(12)+qud1-qud3)/cs
      f(6)  = (-y(13)-qud2+qud4)/cs
      f(7)  = (-y(7)/rp+qud1+qud2-qud3-qud4)/cp
      f(8)  = -y(1)/lh
      f(9)  = -y(2)/lh
      f(10) = (0.5d0*y(1)-y(3)-rg2*y(10))/ls2
      f(11) = (-0.5d0*y(1)+y(4)-rg3*y(11))/ls3
      f(12) = (0.5d0*y(2)-y(5)-rg2*y(12))/ls2
      f(13) = (-0.5d0*y(2)+y(6)-rg3*y(13))/ls3
      f(14) = (-y(1)+uin1-(ri+rg1)*y(14))/ls1
      f(15) = (-y(2)-(rc+rg1)*y(15))/ls1

      return
      end
c-----------------------------------------------------------------------
      subroutine feval_e(neqn,t,y,yprime,f,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)

      double precision c,cs,cp,r,rp,lh,ls1,ls2,ls3,rg1,rg2,rg3,ri,rc,
     +                 gamma,delta,pi,
     +                 uin1,uin2,ud1,ud2,ud3,ud4,qud1,qud2,qud3,qud4
      parameter (c=1.6d-8,cs=2d-12,cp=1d-8,r=25d3,rp=50d0,
     +           lh=4.45d0,ls1=2d-3,ls2=5d-4,ls3=5d-4,
     +           rg1=36.3d0,rg2=17.3d0,rg3=17.3d0,ri=5d1,rc=6d2,
     +           gamma=40.67286402d-9,delta=17.7493332d0,
     +           pi=3.141592653589793238462643383d0)

      uin1   = 0.5d0*sin(2d3*pi*t)
      uin2   = 2d0*sin(2d4*pi*t)
      ud1    = +y(3)-y(5)-y(7)-uin2
      ud2    = -y(4)+y(6)-y(7)-uin2
      ud3    = +y(4)+y(5)+y(7)+uin2
      ud4    = -y(3)-y(6)+y(7)+uin2

c     prevent overflow
c      (double precision ieee .le. 1d304)
      if (delta*max(ud1,ud2,ud3,ud4).gt.300d0) then
         ierr = -1
         return
      endif

      qud1   = gamma*(exp(delta*ud1)-1d0)
      qud2   = gamma*(exp(delta*ud2)-1d0)
      qud3   = gamma*(exp(delta*ud3)-1d0)
      qud4   = gamma*(exp(delta*ud4)-1d0)

      f(1)  = 0.d0
      f(2)  = 0.d0
      f(3)  = 0.d0
      f(4)  = 0.d0
      f(5)  = 0.d0
      f(6)  = 0.d0
      f(7)  = 0.d0
      f(8)  = -y(1)/lh
      f(9)  = -y(2)/lh
      f(10) = (0.5d0*y(1)-y(3)-rg2*y(10))/ls2
      f(11) = (-0.5d0*y(1)+y(4)-rg3*y(11))/ls3
      f(12) = (0.5d0*y(2)-y(5)-rg2*y(12))/ls2
      f(13) = (-0.5d0*y(2)+y(6)-rg3*y(13))/ls3
      f(14) = (-y(1)+uin1-(ri+rg1)*y(14))/ls1
      f(15) = (-y(2)-(rc+rg1)*y(15))/ls1

      return
      end
c-----------------------------------------------------------------------
      subroutine feval_i(neqn,t,y,yprime,f,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)

      double precision c,cs,cp,r,rp,lh,ls1,ls2,ls3,rg1,rg2,rg3,ri,rc,
     +                 gamma,delta,pi,
     +                 uin1,uin2,ud1,ud2,ud3,ud4,qud1,qud2,qud3,qud4
      parameter (c=1.6d-8,cs=2d-12,cp=1d-8,r=25d3,rp=50d0,
     +           lh=4.45d0,ls1=2d-3,ls2=5d-4,ls3=5d-4,
     +           rg1=36.3d0,rg2=17.3d0,rg3=17.3d0,ri=5d1,rc=6d2,
     +           gamma=40.67286402d-9,delta=17.7493332d0,
     +           pi=3.141592653589793238462643383d0)

      uin1   = 0.5d0*sin(2d3*pi*t)
      uin2   = 2d0*sin(2d4*pi*t)
      ud1    = +y(3)-y(5)-y(7)-uin2
      ud2    = -y(4)+y(6)-y(7)-uin2
      ud3    = +y(4)+y(5)+y(7)+uin2
      ud4    = -y(3)-y(6)+y(7)+uin2

c     prevent overflow
c      (double precision ieee .le. 1d304)
      if (delta*max(ud1,ud2,ud3,ud4).gt.300d0) then
         ierr = -1
         return
      endif

      qud1  = gamma*(exp(delta*ud1)-1d0)
      qud2  = gamma*(exp(delta*ud2)-1d0)
      qud3  = gamma*(exp(delta*ud3)-1d0)
      qud4  = gamma*(exp(delta*ud4)-1d0)

      f(1)  = (y(8)-0.5d0*y(10)+0.5d0*y(11)+y(14)-y(1)/r)/c
      f(2)  = (y(9)-0.5d0*y(12)+0.5d0*y(13)+y(15)-y(2)/r)/c
      f(3)  = (y(10)-qud1+qud4)/cs
      f(4)  = (-y(11)+qud2-qud3)/cs
      f(5)  = (y(12)+qud1-qud3)/cs
      f(6)  = (-y(13)-qud2+qud4)/cs
      f(7)  = (-y(7)/rp+qud1+qud2-qud3-qud4)/cp
      f(8)  = 0.d0
      f(9)  = 0.d0
      f(10) = 0.d0
      f(11) = 0.d0
      f(12) = 0.d0
      f(13) = 0.d0
      f(14) = 0.d0
      f(15) = 0.d0

      return
      end
c-----------------------------------------------------------------------
      subroutine jeval(ldim,neqn,t,y,yprime,dfdy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfdy(ldim,neqn),rpar(*)

      integer i,j
      double precision c,cs,cp,r,rp,lh,ls1,ls2,ls3,rg1,rg2,rg3,ri,rc,
     +                 gamma,delta,pi,
     +                 uin2,ud1,ud2,ud3,ud4,qpud1,qpud2,qpud3,qpud4
      parameter (c=1.6d-8,cs=2d-12,cp=1d-8,r=25d3,rp=50d0,
     +           lh=4.45d0,ls1=2d-3,ls2=5d-4,ls3=5d-4,
     +           rg1=36.3d0,rg2=17.3d0,rg3=17.3d0,ri=5d1,rc=6d2,
     +           gamma=40.67286402d-9,delta=17.7493332d0,
     +           pi=3.141592653589793238462643383d0)

      do 20 j=1,neqn
         do 10 i=1,neqn
            dfdy(i,j) = 0d0
   10    continue
   20 continue

      uin2  = 2d0*sin(2d4*pi*t)
      ud1   = +y(3)-y(5)-y(7)-uin2
      ud2   = -y(4)+y(6)-y(7)-uin2
      ud3   = +y(4)+y(5)+y(7)+uin2
      ud4   = -y(3)-y(6)+y(7)+uin2
      qpud1 = gamma*delta*exp(delta*ud1)
      qpud2 = gamma*delta*exp(delta*ud2)
      qpud3 = gamma*delta*exp(delta*ud3)
      qpud4 = gamma*delta*exp(delta*ud4)

      dfdy(1,1)   = -1d0/(c*r)
      dfdy(1,8)   = 1d0/c
      dfdy(1,10)  = -0.5d0/c
      dfdy(1,11)  = -dfdy(1,10)
      dfdy(1,14)  = dfdy(1,8)
      dfdy(2,2)   = dfdy(1,1)
      dfdy(2,9)   = dfdy(1,8)
      dfdy(2,12)  = dfdy(1,10)
      dfdy(2,13)  = dfdy(1,11)
      dfdy(2,15)  = dfdy(1,14)
      dfdy(3,3)   = (-qpud1-qpud4)/cs
      dfdy(3,5)   = qpud1/cs
      dfdy(3,6)   = -qpud4/cs
      dfdy(3,7)   = (qpud1+qpud4)/cs
      dfdy(3,10)  = 1d0/cs
      dfdy(4,4)   = (-qpud2-qpud3)/cs
      dfdy(4,5)   = -qpud3/cs
      dfdy(4,6)   = qpud2/cs
      dfdy(4,7)   = (-qpud2-qpud3)/cs
      dfdy(4,11)  = -1d0/cs
      dfdy(5,3)   = qpud1/cs
      dfdy(5,4)   = -qpud3/cs
      dfdy(5,5)   = (-qpud1-qpud3)/cs
      dfdy(5,7)   = (-qpud1-qpud3)/cs
      dfdy(5,12)  = 1d0/cs
      dfdy(6,3)   = -qpud4/cs
      dfdy(6,4)   = qpud2/cs
      dfdy(6,6)   = (-qpud2-qpud4)/cs
      dfdy(6,7)   = (qpud2+qpud4)/cs
      dfdy(6,13)  = -1d0/cs
      dfdy(7,3)   = (qpud1+qpud4)/cp
      dfdy(7,4)   = (-qpud2-qpud3)/cp
      dfdy(7,5)   = (-qpud1-qpud3)/cp
      dfdy(7,6)   = (qpud2+qpud4)/cp
      dfdy(7,7)   = (-qpud1-qpud2-qpud3-qpud4-1d0/rp)/cp
      dfdy(8,1)   = -1d0/lh
      dfdy(9,2)   = dfdy(8,1)
      dfdy(10,1)  = 0.5d0/ls2
      dfdy(10,3)  = -1d0/ls2
      dfdy(10,10) = -rg2/ls2
      dfdy(11,1)  = -0.5d0/ls3
      dfdy(11,4)  = 1d0/ls3
      dfdy(11,11) = -rg3/ls3
      dfdy(12,2)  = dfdy(10,1)
      dfdy(12,5)  = dfdy(10,3)
      dfdy(12,12) = dfdy(10,10)
      dfdy(13,2)  = dfdy(11,1)
      dfdy(13,6)  = dfdy(11,4)
      dfdy(13,13) = dfdy(11,11)
      dfdy(14,1)  = -1d0/ls1
      dfdy(14,14) = -(ri+rg1)/ls1
      dfdy(15,2)  = dfdy(14,1)
      dfdy(15,15) = -(rc+rg1)/ls1

      return
      end
c-----------------------------------------------------------------------
      subroutine jeval_i(ldim,neqn,t,y,yprime,dfdy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfdy(ldim,neqn),rpar(*)

      integer i,j
      double precision c,cs,cp,r,rp,gamma,delta,pi,
     +                 uin2,ud1,ud2,ud3,ud4,qpud1,qpud2,qpud3,qpud4
      parameter (c=1.6d-8,cs=2d-12,cp=1d-8,r=25d3,rp=50d0,
     +           gamma=40.67286402d-9,delta=17.7493332d0,
     +           pi=3.141592653589793238462643383d0)

      do 20 j=1,neqn
         do 10 i=1,neqn
            dfdy(i,j) = 0.d0
   10    continue
   20 continue

      uin2  = 2d0*sin(2d4*pi*t)
      ud1   = +y(3)-y(5)-y(7)-uin2
      ud2   = -y(4)+y(6)-y(7)-uin2
      ud3   = +y(4)+y(5)+y(7)+uin2
      ud4   = -y(3)-y(6)+y(7)+uin2
      qpud1 = gamma*delta*exp(delta*ud1)
      qpud2 = gamma*delta*exp(delta*ud2)
      qpud3 = gamma*delta*exp(delta*ud3)
      qpud4 = gamma*delta*exp(delta*ud4)

      dfdy(1,1)   = -1d0/(c*r)
      dfdy(1,8)   = 1d0/c
      dfdy(1,10)  = -0.5d0/c
      dfdy(1,11)  = -dfdy(1,10)
      dfdy(1,14)  = dfdy(1,8)
      dfdy(2,2)   = dfdy(1,1)
      dfdy(2,9)   = dfdy(1,8)
      dfdy(2,12)  = dfdy(1,10)
      dfdy(2,13)  = dfdy(1,11)
      dfdy(2,15)  = dfdy(1,14)
      dfdy(3,3)   = (-qpud1-qpud4)/cs
      dfdy(3,5)   = qpud1/cs
      dfdy(3,6)   = -qpud4/cs
      dfdy(3,7)   = (qpud1+qpud4)/cs
      dfdy(3,10)  = 1d0/cs
      dfdy(4,4)   = (-qpud2-qpud3)/cs
      dfdy(4,5)   = -qpud3/cs
      dfdy(4,6)   = qpud2/cs
      dfdy(4,7)   = (-qpud2-qpud3)/cs
      dfdy(4,11)  = -1d0/cs
      dfdy(5,3)   = qpud1/cs
      dfdy(5,4)   = -qpud3/cs
      dfdy(5,5)   = (-qpud1-qpud3)/cs
      dfdy(5,7)   = (-qpud1-qpud3)/cs
      dfdy(5,12)  = 1d0/cs
      dfdy(6,3)   = -qpud4/cs
      dfdy(6,4)   = qpud2/cs
      dfdy(6,6)   = (-qpud2-qpud4)/cs
      dfdy(6,7)   = (qpud2+qpud4)/cs
      dfdy(6,13)  = -1d0/cs
      dfdy(7,3)   = (qpud1+qpud4)/cp
      dfdy(7,4)   = (-qpud2-qpud3)/cp
      dfdy(7,5)   = (-qpud1-qpud3)/cp
      dfdy(7,6)   = (qpud2+qpud4)/cp
      dfdy(7,7)   = (-qpud1-qpud2-qpud3-qpud4-1d0/rp)/cp

      return
      end
c-----------------------------------------------------------------------
      subroutine meval(ldim,neqn,t,y,yprime,dfddy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfddy(ldim,neqn),rpar(*)
c
c     dummy subroutine
c
      return
      end
c-----------------------------------------------------------------------
      subroutine solut(neqn,t,y)
      integer neqn
      double precision t,y(neqn)
c
c determined with RADAU
c
c User input:
c
c give relative error tolerance: 1d-12
c give absolute error tolerance: 1d-12
c give initial stepsize: 1d-14
c
c Integration characteristics:
c
c    number of integration steps       27195
c    number of accepted steps          23541
c    number of f evaluations          632631
c    number of Jacobian evaluations     9468
c    number of LU decompositions       22386

      y(  1) = -0.2339057358486745d-001
      y(  2) = -0.7367485485540825d-002
      y(  3) =  0.2582956709291169d+000
      y(  4) = -0.4064465721283450d+000
      y(  5) = -0.4039455665149794d+000
      y(  6) =  0.2607966765422943d+000
      y(  7) =  0.1106761861269975d+000
      y(  8) =  0.2939904342435596d-006
      y(  9) = -0.2840029933642329d-007
      y( 10) =  0.7267198267264553d-003
      y( 11) =  0.7929487196960840d-003
      y( 12) = -0.7255283495698965d-003
      y( 13) = -0.7941401968526521d-003
      y( 14) =  0.7088495416976114d-004
      y( 15) =  0.2390059075236570d-004

      return
      end
