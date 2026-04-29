!**********************************************************************!
! Module to compare everything with an exact harmonic oscillator       !
!**********************************************************************!

       module coherent_module
 
       use constants

       implicit none

       contains

!......Coherent calculator...............................................

       subroutine coherent_calc(nd,trj,q0,p0,masses,c0,Bcmplx)
       ! nd: bath dimension
       ! trj : trajectory parameters (first step, last step, nstep)
       ! q0 : initial gaussian center (x&y)
       ! p0 : initial gaussian momentum (x&y)
       ! masses : vector of the masses
       ! c0 : initial coefficients 
       ! Bcmplx : initial gaussian width matrix, real & imaginary
       integer, intent(in) :: nd
       integer*8, dimension(3), intent(in) :: trj
       real*8, dimension(nd+1), intent(in) :: q0,p0,masses
       complex*16,dimension(nh), intent(in) :: c0 
       complex*16,dimension(nd+1,nd+1), intent(in) :: Bcmplx 
      
       integer*8 :: i,j,first,last,nstep,k
       real*8 :: h,time,N,E,q,p,E0,Nsq
       real*8,dimension(nh) :: csq 
       real*8,dimension(nd) :: qvec,pvec 
       real*8,dimension(nd+1) :: qtoti,ptoti,qtotj,ptotj,qold,pold 
       real*8,dimension(4) :: hvec = [0.5d0,0.5d0,1.d0,0.d0]
       complex*16,dimension(nd+1,nd+1) :: Bcoh,Bt,num,den
       
       complex*16,dimension(nh) :: cvec,cj,ctemp

       real*8,dimension(nd+1,4) :: kq,kp
       complex*16,dimension(nd+1,nd+1,4) :: kb

       real*8 :: a,Y0
       real*8, dimension(nd) :: avec
       real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
       real*8, dimension(nd+1,nd+1) :: Bmat
       real*8, dimension(nh,nh) :: S00M,invS,X0Mat

       open(unit=777,file='coherent.dat',status='unknown')

       ! trajectory parameters
       first = trj(1)
       last = trj(2)
       nstep = trj(3)

       ! Separate x & y
       q = q0(1)
       qvec = q0(2:nd+1)
       p = p0(1)
       pvec = p0(2:nd+1)

       h = dfloat(last-first)/dfloat(nstep)

       Bcoh(:,:) =0.d0
       do i = 1,nd+1
          Bcoh(i,i) = dsqrt(masses(i))
       end do

       write(777,*) '#Time ','q ', 'p ','B(1,1) ','B(2,2)'

       do j = 1,nstep
          time = j*h

          qtot(:) = q0*dcos(1.d0*time) +p0(:)*dsin(1.d0*time)/masses(:)
          ptot(:) = p0*dcos(1.d0*time) - masses(:)*q0*dsin(1.d0*time)
          num(:,:)=Bcmplx(:,:)*dcos(1.d0*time)+&
                   &iu*Bcoh(:,:)*dsin(1.d0*time)
          den(:,:)=iu*Bcmplx(:,:)*dsin(1.d0*time)+&
                   &Bcoh(:,:)*dcos(1.d0*time)
          Bt(:,:) = Bcoh(:,:)*num(:,:)/den(:,:)

          write(777,*) time,qtot,ptot,Bt(1,1),Bt(2,2)

       end do

       end subroutine

       end module 

 
