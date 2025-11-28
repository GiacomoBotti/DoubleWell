!**********************************************************************!
! Module that computes the effective potential terms aka the           !
! expectation values of V, dV and ddV                                  !
!**********************************************************************!

      module effectivepot_module

      use constants
      use matrix_module
      use basisset_module
      use potential_module
      use integrals_module

      implicit none
 
      private
      public :: fun_V0

      contains

!.....EFFECTIVE POTENTIAL V0 = <V>......................................

      function fun_V0(nd,q,cvec,Bmat) result(V0)
      ! nd: bath dimensions
      ! q: x Gaussian center
      ! cvec: vector of the coefficients
      ! Bmat: total Gaussian width matrix
       integer, intent(in) :: nd
       real*8, intent(in) :: q
       complex*16, dimension(nh), intent(in) :: cvec
       real*8, dimension(nd+1,nd+1), intent(in) :: Bmat

       integer :: i
       real*8 :: Nsq,Y0,a
       complex*16 :: V0,Vx,Vxy,Vy
       real*8, dimension(nd) :: avec
       real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
       real*8, dimension(nh) :: Xcvec
       real*8, dimension(nh,nh) :: X4mat, X2mat, Xtot

       Nsq=fun_Nsq(nd+1,Bmat)

       call extractA(nd,Bmat,Amat,avec,a)
      
       call diagonalization(nd,Amat,LambdaMat,Tmat)

!       do i = 1,nd
!         write(*,*) LambdaMat(i,:)
!       end do
 
       Y0=int_Y0(nd,LambdaMat)

!       write(*,*) Y0

       X4mat=int_XnMat(nd,4,a,avec,Amat,q)
       write(*,*) X4mat(3,3)
       X2mat=int_XnMat(nd,2,a,avec,Amat,q)
       write(*,*) X2mat(3,3)

       Xtot=(X4mat/(16.d0*eta_const)) - X2mat/2.d0      
       
       Xcvec=matmul(Xtot,cvec)

       Vx=Nsq*Y0*dot_product(cvec,Xcvec)
             
       Vxy=0.d0
       Vy=0.d0
       V0 = Vx+Vxy+Vy
        
      end function

      end module
