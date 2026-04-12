!**********************************************************************!
! Module that contains all things pertaining the evolution of the      !
! Gaussian wavepacket following Coalson & Karplus 1990                 !
!**********************************************************************!

      module eofmotion_module

      use potential_module
      use basisset_module
      use effectivepot_module
      use integrals_module
      use matrix_module

      implicit none

      real*8, dimension(nv+1,nv+1), public :: invMassMat

      private
      public :: MassesMat,KarplusTimeDer 

      contains

!.....Matrix of the masses..............................................

      subroutine MassesMat(masses)
      ! masses: vector of the masses
      ! invMassMat: Matrix of the inverted masses
       real*8, dimension(nv+1), intent(in) :: masses  

       integer :: i

       invMassMat(:,:) = 0.d0

       do i = 1,nv+1
          invMassMat(i,i) = 1.d0/masses(i)
!          write(*,*) invMassMat(i,:)
       end do

       write(*,*) "Matrix of inverse masses generated"

      end subroutine

!.....C&K equations of motion...........................................

      subroutine KarplusTimeDer(nd,cvec,qtot,ptot,Bcmplx,dotq,dotp,dotB)
      ! nd: dimension of bath
      ! cvec: vector of the basis coefficients
      ! qtot: gaussian center position vector (x&y)
      ! ptot: gaussian center momenta vector (x&y)
      ! Bmat: total gaussian width (x&y)
      ! dotq: time-derivative of qtot
      ! dotp: time-derivative of ptot
      ! dotB: time-derivative of Bmat
       integer, intent(in) :: nd
       complex*16, dimension(nh), intent(in) :: cvec
       real*8, dimension(nd+1), intent(in) :: qtot,ptot
       complex*16, dimension(nd+1,nd+1), intent(in) :: Bcmplx
     
       real*8, dimension(nd+1), intent(out) :: dotq, dotp
       complex*16, dimension(nd+1,nd+1), intent(out) :: dotB

       real*8, dimension(nd+1) :: V1
       real*8, dimension(nd+1,nd+1) :: V2,Bmat
       real*8 :: a,Y0
       real*8, dimension(nd) :: avec
       real*8, dimension(nd,nd) :: Amat,Tmat,LambdaMat
       real*8, dimension(nh,nh) :: X3,X2,X1,X0
       complex*16, dimension(nd+1,nd+1) :: MB,prova
 
       Bmat = dreal(Bcmplx)
       call extractA(nd,Bmat,Amat,avec,a)
       call diagonalization(nd,Amat,LambdaMat,Tmat)
 
       dotq = matmul(invMassMat,ptot)

       ! Integrals
       Y0=int_Y0(nd,LambdaMat)
       X3=int_XnMat(nd,3,a,avec,Amat,qtot(1))
       X2=int_XnMat(nd,2,a,avec,Amat,qtot(1))
       X1=int_XnMat(nd,1,a,avec,Amat,qtot(1))
       X0=int_XnMat(nd,0,a,avec,Amat,qtot(1))

       V1 = fun_V1(nd,qtot,cvec,Bmat,Y0,X3,X2,X1,X0)
       dotp = - V1

       V2 = fun_V2(nd,qtot,cvec,Bmat,Y0,X2,X0)
       MB=matmul(invMassMat,Bcmplx)              
!       dotB = -2*matmul(Bcmplx,MB) - V2/2.d0
       prova=matmul(Bcmplx,MB)
       dotB = -(0.d0,1.d0)*matmul(Bcmplx,MB) + (0.d0,1.d0)*V2
!       dotB =0.d0

!       write(*,*) dotB(1,1), V2(1,1), prova(1,1)

      end subroutine

      end module
