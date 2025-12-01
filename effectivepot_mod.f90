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
      use quadratic_module

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
       real*8 :: Nsq,Y0,a,q,aAVAa,qVAa,qVq,uWu
       complex*16 :: V0,Vx,Vxy,Vy
       real*8, dimension(nd) :: avec, qvec, Aa, VAa
       real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat,invA
       real*8, dimension(maxorder,nd) :: MomMat
       real*8, dimension(nh) :: Xcvec, XYcvec, Ycvec
       real*8, dimension(nh,nh) :: X4mat,X2mat,X1mat,X0mat
       real*8, dimension(nh,nh) :: Xtot,XYtot,Ytot

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
       VAa = matmul(Vmat,Aa)
       aAVAa = dot_product(Aa,VAa)
       qVAa = dot_product(qvec,VAa)

       MomMat = momenta(nd,LambdaMat)
       qVq = fun_qVq(nd,qvec)
       uWu = fun_uWu(nd,Tmat,MomMat)

!       do i = 1,nd
!         write(*,*) invA(i,:) !matmul(invA(i,:),Amat)
!       end do

       write(*,*) "Matrix of the Momenta:"
       do i = 1,maxorder
         write(*,*) MomMat(i,:) 
       end do


       write(*,*) uWu,qVq!qVAa!aAVAa
!       write(*,*) avec 
!       write(*,*) VAa(:) !Aa(:) !matmul(Amat,Aa(:))

       write(*,*) "Diagonal Matrix"
       do i = 1,nd
         write(*,*) LambdaMat(i,:)
       end do
 
       ! Integrals
       Y0=int_Y0(nd,LambdaMat)
!       write(*,*) "Y0:", Y0
       X4mat=int_XnMat(nd,4,a,avec,Amat,q)
!       write(*,*) "X4mat(3,3)", X4mat(3,3)
       X2mat=int_XnMat(nd,2,a,avec,Amat,q)
!       write(*,*) "X2mat(3,3)", X2mat(3,3)
       X1mat=int_XnMat(nd,1,a,avec,Amat,q)
!       write(*,*) "X1mat(3,3)", X1mat(3,3)
       X0mat=int_XnMat(nd,0,a,avec,Amat,q)
!       write(*,*) "X0mat(3,3)", X0mat(3,3)

       ! Total Hermite Matrices
       Xtot=(X4mat/(16.d0*eta_const)) - X2mat/2.d0      
       Xcvec=matmul(Xtot,cvec)
       XYtot=(qvec(1) +q*Aa(1))*X1mat - Aa(1)*X2mat
       XYcvec=matmul(XYtot,cvec)
       Ytot=X0mat*(uWu+qVq+aAVAa*q**2-q*qVAa)+&
           &X1mat*(2*qVAa+2*q*aAVAa)+&
           &X2mat*aAVAa
       Ycvec=matmul(Ytot,cvec)

       ! V elements
       Vx=Nsq*Y0*dot_product(cvec,Xcvec)
       Vxy=Nsq*Y0*gamma_const*dot_product(cvec,XYcvec)
       Vy=Nsq*Y0*dot_product(cvec,Ycvec)
             
!       Vx=0.d0
!       Vxy=0.d0
       Vy=0.d0
       V0 = Vx+Vxy+Vy
        
      end function

      end module
