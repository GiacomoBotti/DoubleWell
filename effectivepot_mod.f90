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
      use inversion_module

      implicit none
 
      private
      public :: fun_V0

      contains

!.....EFFECTIVE POTENTIAL V0 = <V>......................................

      function fun_V0(nd,qtot,cvec,Bmat) result(V0)
      ! nd: bath dimensions
      ! qtot: total Gaussian center vector (x&y)
      ! cvec: vector of the coefficients
      ! Bmat: total Gaussian width matrix (x&y)
       integer, intent(in) :: nd
       real*8, dimension(nd+1), intent(in) :: qtot
       complex*16, dimension(nh), intent(in) :: cvec
       real*8, dimension(nd+1,nd+1), intent(in) :: Bmat

       integer :: i
       real*8 :: Nsq,Y0,a,q
       complex*16 :: V0,Vx,Vxy,Vy
       real*8, dimension(nd) :: avec, qvec, Aa
       real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat,invA
       real*8, dimension(nh) :: Xcvec, XYcvec
       real*8, dimension(nh,nh) :: X4mat, X2mat, X1mat, Xtot, XYtot

       ! Unreavel qtot
       q = qtot(1)
       qvec(:) = qtot(2:nd+1)

       ! Nsq
       Nsq=fun_Nsq(nd+1,Bmat)

       ! Matrix work 
       call extractA(nd,Bmat,Amat,avec,a)
       call diagonalization(nd,Amat,LambdaMat,Tmat)
       invA = invgen_real(nd,Amat)
       Aa = matmul(invA,avec)

       do i = 1,nd
         write(*,*) invA(i,:)
       end do

       write(*,*) Aa(:)

!       do i = 1,nd
!         write(*,*) LambdaMat(i,:)
!       end do
 
       ! Integrals
       Y0=int_Y0(nd,LambdaMat)
!       write(*,*) Y0
       X4mat=int_XnMat(nd,4,a,avec,Amat,q)
!       write(*,*) X4mat(3,3)
       X2mat=int_XnMat(nd,2,a,avec,Amat,q)
!       write(*,*) X2mat(3,3)
       X1mat=int_XnMat(nd,1,a,avec,Amat,q)
!       write(*,*) X1mat(3,3)

       ! Total Hermite Matrices
       Xtot=(X4mat/(16.d0*eta_const)) - X2mat/2.d0      
       Xcvec=matmul(Xtot,cvec)
       XYtot=(qvec(1) +q*Aa(1))*X1mat - Aa(1)*X2mat
       XYcvec=matmul(XYtot,cvec)

       ! V elements
       Vx=Nsq*Y0*dot_product(cvec,Xcvec)
       Vxy=Nsq*Y0*gamma_const*dot_product(cvec,XYcvec)
             
       Vx=0.d0
!       Vxy=0.d0
       Vy=0.d0
       V0 = Vx+Vxy+Vy
        
      end function

      end module
