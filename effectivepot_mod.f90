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
      public :: fun_V0,fun_V1,fun_V2

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
       real*8 :: V0,Vx,Vxy,Vy
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

!       write(*,*) "Matrix of the Momenta:"
!       do i = 1,maxorder
!         write(*,*) MomMat(i,:) 
!       end do

!       write(*,*) uWu,qVq!qVAa!aAVAa
!       write(*,*) avec 
!       write(*,*) VAa(:) !Aa(:) !matmul(Amat,Aa(:))

!       write(*,*) "Diagonal Matrix"
!       do i = 1,nd
!         write(*,*) LambdaMat(i,:)
!       end do
 
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
!       Vy=0.d0
       V0 = Vx+Vxy+Vy
        
      end function

!.....EFFECTIVE POTENTIAL V1 = <nabla V>................................

      function fun_V1(nd,qtot,cvec,Bmat) result(V1)
      ! nd: bath dimensions
      ! qtot: total Gaussian center vector (x&y)
      ! cvec: vector of the coefficients
      ! Bmat: total Gaussian width matrix (x&y)
       integer, intent(in) :: nd
       real*8, dimension(nd+1), intent(in) :: qtot
       complex*16, dimension(nh), intent(in) :: cvec
       real*8, dimension(nd+1,nd+1), intent(in) :: Bmat

       real*8, dimension(nd+1) :: V1

       integer :: i
       real*8 :: Nsq,Y0,a,q,cX3c,cX1c,cX0c
       real*8 :: dxV, dy1V
       real*8, dimension(nd) :: avec,qvec,Aa,VAa,Vq,VVAa,V1prime
       real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat,invA
       real*8, dimension(nh) :: X3c, X1c, X0c 
       real*8, dimension(nh,nh) :: X3mat,X2mat,X1mat,X0mat

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
       VVAa = matmul(Vmat,VAa)

       Vq = matmul(Vmat,qvec)

!       write(*,*) Aa(:) 
!       write(*,*) Vq(:)
!       write(*,*) VVAa(:)

       ! Integrals
       Y0=int_Y0(nd,LambdaMat)
!       write(*,*) "Y0:", Y0
       X3mat=int_XnMat(nd,3,a,avec,Amat,q)
!       write(*,*) "X3mat(3,3)", X3mat(3,3)
       X1mat=int_XnMat(nd,1,a,avec,Amat,q)
!       write(*,*) "X1mat(3,3)", X1mat(3,3)
       X0mat=int_XnMat(nd,0,a,avec,Amat,q)
!       write(*,*) "X0mat(3,3)", X0mat(3,3)

       ! Total Hermite Matrices
       X3c=matmul(X3mat,cvec)
       X1c=matmul(X1mat,cvec)
       X0c=matmul(X0mat,cvec)
       cX3c=dot_product(cvec,X3c)
       cX1c=dot_product(cvec,X1c)
       cX0c=dot_product(cvec,X0c)
       dxV=(cX3c/(4.d0*eta_const) - (1+gamma_const*Aa(1))*cX1c+&
             &gamma_const*(qvec(1)+Aa(1)*q)*cX0c)*Nsq*Y0

       V1prime=2*Nsq*Y0*(cX0c*(Vq+q*VAa)-cX1c*VVAa)

       dy1V=gamma_const*Nsq*Y0*cX1c + V1prime(1)

       V1(1) = dxV
       V1(2) = dy1V
       V1(3:nd+1) = V1prime(2:nd)

!       write(*,*) "V1", V1(:)
!       write(*,*) "Aa", Aa(:)
!       write(*,*) "VAa", VAa(:)
!       write(*,*) "Nsq", Nsq

      end function

!.....EFFECTIVE POTENTIAL V2= <nabla otimes nabla V>....................

      function fun_V2(nd,qtot,cvec,Bmat) result(V2)
      ! nd: bath dimensions
      ! qtot: total Gaussian center vector (x&y)
      ! cvec: vector of the coefficients
      ! Bmat: total Gaussian width matrix (x&y)
       integer, intent(in) :: nd
       real*8, dimension(nd+1), intent(in) :: qtot
       complex*16, dimension(nh), intent(in) :: cvec
       real*8, dimension(nd+1,nd+1), intent(in) :: Bmat

       real*8, dimension(nd+1,nd+1) :: V2

       integer :: i
       real*8 :: Nsq,Y0,a,q,cX2c,cX0c
       real*8 :: dxV, dy1V
       real*8, dimension(nd) :: avec,qvec,Aa,VAa,Vq,VVAa,V1prime
       real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat,invA
       real*8, dimension(nh) :: X2c, X0c 
       real*8, dimension(nh,nh) :: X2mat,X0mat

       ! Unreavel qtot
       q = qtot(1)
       qvec(:) = qtot(2:nd+1)

       ! Nsq
       Nsq=fun_Nsq(nd+1,Bmat)

       ! Matrix work 
       call extractA(nd,Bmat,Amat,avec,a)
       call diagonalization(nd,Amat,LambdaMat,Tmat)
       !invA = invgen_real(nd,Amat)

       !Aa = matmul(invA,avec)
       !VAa = matmul(Vmat,Aa)
       !VVAa = matmul(Vmat,VAa)

       !Vq = matmul(Vmat,qvec)

!       write(*,*) Aa(:) 
!       write(*,*) Vq(:)
!       write(*,*) VVAa(:)

       ! Integrals
       Y0=int_Y0(nd,LambdaMat)
!       write(*,*) "Y0:", Y0
       X2mat=int_XnMat(nd,2,a,avec,Amat,q)
!       write(*,*) "X2mat(3,3)", X2mat(3,3)
       X0mat=int_XnMat(nd,0,a,avec,Amat,q)
!       write(*,*) "X0mat(3,3)", X0mat(3,3)

       ! Total Hermite Matrices
       X2c=matmul(X2mat,cvec)
       X0c=matmul(X0mat,cvec)
       cX2c=dot_product(cvec,X2c)
       cX0c=dot_product(cvec,X0c)

       V2(:,:) = 0.d0

       V2(1,1) = Nsq*Y0*(3.d0*cX2c/(4.d0*eta_const)-cX0c)
       V2(1,2) = gamma_const*Nsq*Y0*cX0c
       V2(2,1) = V2(1,2) 
       V2(2:nd+1,2:nd+1) = 2.d0*Nsq*Y0*cX0c*Vmat

!       write(*,*) "V2"
!       do i = 1,nd+1
!          write(*,*) V2(i,:)
!       end do

      end function
     

      end module
