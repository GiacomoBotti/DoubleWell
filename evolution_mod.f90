!**********************************************************************!
! Module containing the Runge-Kutta evolution of the parameters and    !
! coefficients                                                         !
!**********************************************************************!

       module evolution_module

       use normalization_module
       use observable_module
       use BOT_module

       implicit none

       private
       public :: bot_evo

       contains

!......BOT evolution....................................................

       subroutine bot_evo(nd,trj)
       ! nd: bath dimension
       ! trj : trajectory parameters (first step, last step, nstep)
       integer, intent(in) :: nd
       integer*8, dimension(3), intent(in) :: trj

       integer*8 :: i,j,first,last,nstep
       real*8 :: h,time,N,E,a,alpha,xi,zeta,q,p,qold,pold,nover
       real*8,dimension(4) :: hvec = [0.5d0,0.5d0,1.d0,0.d0]
       real*8,dimension(3) :: cond,nparam 
       !real*8,dimension(nd) :: vecxi,veczeta
       
       complex*16,dimension(npar+nd) :: yj,yi
       complex*16,dimension(nd) :: cj,cjsupp 
!       complex*16,dimension(npar) :: lambdaj,lambdai,lambdaold
       real*8,dimension(npar) :: lambdaj,lambdai,lambdaold

       real*8,dimension(npar,4) :: k

       ! trajectory parameters
       first = trj(1)
       last = trj(2)
       nstep = trj(3)

       open(unit=321,file="trajectory_BOT.dat",status="unknown")
       h = dfloat(last-first)/dfloat(nstep)
       hvec = hvec*h

       N = normalization(nd,npar,y0,work)
       ! Numerical:
!       N = normnumtot(nd,npar,y0,work)

       write(321,*) "#Evolution parameters:"
       write(321,*) "#Range: ",first,last
       write(321,*) "#Steps: ",nstep
       write(321,*) "#Timestep: ",h
       write(321,*) "#Normalization constant: ",N
       !write(321,*) "#Square root of N ",dsqrt(N)
  
       lambdaj(1:npar) = y0(1:npar)
       cj(:) = y0(npar+1:npar+nd)/dsqrt(N)

       lambdai(1:npar) = y0(1:npar)

       yj(1:npar) = lambdaj(:)
       yj(npar+1:npar+nd) = cj(:)

       N = normalization(nd,npar,yj,work)
       E = energy(nd,npar,yj,work)
       ! Numerical:
!       N = normnumtot(nd,npar,yj,work)
!       E = energynum(nd,npar,yj,work)

       write(321,*) "#Time ","q ","p ","rec1 ","rec2 ","imc1 ","imc2 ",&
                    &"Norm ","<H> "
       write(321,*) "#0.d0",dreal(y0(1)),dreal(y0(2)),dreal(y0(3)),&
                   &dreal(y0(4)),dimag(y0(3)),dimag(y0(4)),N,E      
       write(321,*) 0.d0,dreal(yj(1)),dreal(yj(2)),dreal(yj(3)),&
                   &dreal(yj(4)),dimag(yj(3)),dimag(yj(4)),N,E     

       k(:,:) = 0.d0
       do j = 1,nstep
          time = j*h
          cj = c_static(nd,npar,yj,work,h)
!          cj = c_stat_num(nd,npar,yj,work,h)
          do i = 1,4
             q = lambdai(1)
             p = lambdai(2)
             k(:,i) = dotlambda(nd,npar,q,p,a,alpha,xi,zeta,cj)
!             k(:,i) = dotlamnum(nd,npar,q,p,a,alpha,xi,zeta,cj)
             lambdai = lambdaj + hvec(i)*k(:,i)
          end do
          lambdaold = lambdaj
          lambdaj=lambdaj+h*(k(:,1)+2.d0*k(:,2)+2.d0*k(:,3)+k(:,4))/6.d0

          cj = c_update(nd,npar,lambdaj,lambdaold,cj,work)
!          cj = c_upd_num(nd,npar,lambdaj,lambdaold,cj,work)

          yj(1:npar) = lambdaj(:)
          yj(npar+1:npar+nd) = cj(:)

          N = normalization(nd,npar,yj,work)
          E = energy(nd,npar,yj,work)
          ! Numerical:
!          N = normnumtot(nd,npar,yj,work)
!          E = energynum(nd,npar,yj,work)

          write(321,*) time,dreal(yj(1)),dreal(yj(2)),dreal(yj(3)),&
                   &dreal(yj(4)),dimag(yj(3)),dimag(yj(4)),N,E      
      
          lambdai = lambdaj
       end do

       close(321)

       end subroutine

       end module
