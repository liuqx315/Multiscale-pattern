c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Problem E5
c        ODE of dimension 4
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/e5.f
c
c     This is revision
c     $Id: e5.F,v 1.3 2006/10/03 11:14:46 testset Exp $
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
      character fullnm*23, problm*2, type*3
      integer neqn,ndisc,mljac,mujac,mlmas,mumas,ind(1)
      double precision t(0:1)
      logical numjac, nummas

      fullnm = 'Problem E5 stiff-detest'
      problm = 'e5'
      type   = 'ODE'
      neqn   = 4
      ndisc  = 0
      t(0)   = 0d0
      t(1)   = 1.0d13
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

        y(1)=1.76D-3
        y(2)=0.0D0
        y(3)=0.0D0
        y(4)=0.0D0
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
      double precision prod1, prod2, prod3, prod4

        prod1=7.89D-10*y(1)
        prod2=1.1D7*y(1)*y(3)
        prod3=1.13D9*y(2)*y(3)
        prod4=1.13D3*y(4)
        f(1)=-prod1-prod2
        f(2)=prod1-prod3
        f(4)=prod2-prod4
        f(3)=f(2)-f(4)
    
      return
      end
c-----------------------------------------------------------------------
      subroutine feval_e(neqn,t,y,yprime,f,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)
      double precision prod1, prod4

        prod1=7.89D-10*y(1)
        prod4=1.13D3*y(4)
        f(1)=-prod1
        f(2)=prod1
        f(4)=-prod4
        f(3)=prod1+prod4
    
      return
      end
c-----------------------------------------------------------------------
      subroutine feval_i(neqn,t,y,yprime,f,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)
      double precision prod2, prod3

        prod2=1.1D7*y(1)*y(3)
        prod3=1.13D9*y(2)*y(3)
        f(1)=-prod2
        f(2)=-prod3
        f(4)=prod2
        f(3)=-prod3-prod2
    
      return
      end
c-----------------------------------------------------------------------
      subroutine jeval(ldim,neqn,t,y,yprime,dfdy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfdy(ldim,neqn),rpar(*)
      double precision A,B,CM,C

     

        A=7.89D-10
        B=1.1D7
        CM=1.13D9
        C=1.13D3
        dfdy(1,1)=-A-B*y(3)
        dfdy(1,2)=0.D0
        dfdy(1,3)=-B*y(1)
        dfdy(1,4)=0.D0
        dfdy(2,1)=A
        dfdy(2,2)=-CM*y(3)
        dfdy(2,3)=-CM*y(2)
        dfdy(2,4)=0.D0
        dfdy(3,1)=A-B*y(3)
        dfdy(3,2)=-CM*y(3)
        dfdy(3,3)=-B*y(1)-CM*y(2)
        dfdy(3,4)=C
        dfdy(4,1)=B*y(3)
        dfdy(4,2)=0.D0
        dfdy(4,3)=B*y(1)
        dfdy(4,4)=-C
     
      return
      end
c-----------------------------------------------------------------------
      subroutine jeval_i(ldim,neqn,t,y,yprime,dfdy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfdy(ldim,neqn),rpar(*)
      double precision B,CM

        B=1.1D7
        CM=1.13D9
        dfdy(1,1)=-B*y(3)
        dfdy(1,2)=0.D0
        dfdy(1,3)=-B*y(1)
        dfdy(1,4)=0.D0
        dfdy(2,1)=0.D0
        dfdy(2,2)=-CM*y(3)
        dfdy(2,3)=-CM*y(2)
        dfdy(2,4)=0.D0
        dfdy(3,1)=-B*y(3)
        dfdy(3,2)=-CM*y(3)
        dfdy(3,3)=-B*y(1)-CM*y(2)
        dfdy(3,4)=0.D0
        dfdy(4,1)=B*y(3)
        dfdy(4,2)=0.D0
        dfdy(4,3)=B*y(1)
        dfdy(4,4)=0.D0
     
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
c computed using double precision RADAU on an 
c     Alphaserver DS20E, with a 667 MHz EV67 processor.
c          
c          uround = 1.01d-19
c          rtol  = h0 = 1.1d-18, atol = 1.1d-40
c
c

        y(1) =  0.1152903278711829D-290                     
        y(2) =  0.8867655517642120D-022              
        y(3) =  0.8854814626268838D-022           
        y(4) =  0.0000000000000000D+000
   
      return
      end
