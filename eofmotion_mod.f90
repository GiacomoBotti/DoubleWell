!**********************************************************************!
! Module that contains all things pertaining the evolution of the      !
! Gaussian wavepacket following Coalson & Karplus 1990                 !
!**********************************************************************!

      module eofmotion_module

      use potential_module
      use basisset_module
      use effectivepot_module

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
       complex*16, dimension(nd+1,nd+1) :: MB,prova
 
       Bmat = dreal(Bcmplx)
 
       dotq = matmul(invMassMat,ptot)

       V1 = fun_V1(nd,qtot,cvec,Bmat)
       dotp = - V1

       V2 = fun_V2(nd,qtot,cvec,Bmat)
       MB=matmul(invMassMat,Bcmplx)              
!       dotB = -2*matmul(Bcmplx,MB) - V2/2.d0
       prova=matmul(Bcmplx,MB)
       dotB = -(0.d0,1.d0)*matmul(Bcmplx,MB) + (0.d0,1.d0)*V2
!       dotB =0.d0

!       write(*,*) dotB(1,1), V2(1,1), prova(1,1)

      end subroutine

      end module
