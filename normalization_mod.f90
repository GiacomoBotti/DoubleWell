!**********************************************************************!
! Module containing everything pertaining normalization, both          !
! analytical and numerical                                             !
!**********************************************************************!

      module normalization_module

      use integrals_module
      use basisset_module
      use matrix_module

      implicit none
      private
      public :: normalization

      contains

!.....Normalization.....................................................

      function normalization(nd,q,cvec,Bmat) result(Nout)
      ! nd: dimension of the bath
      ! q: center of the well gaussian
      ! cvec: vector of the coefficients
      ! work: work array for I/O
       integer, intent(in) :: nd
       real*8, intent(in) :: q
       complex*16, dimension(nh), intent(in) :: cvec
       real*8, dimension(nd+1,nd+1), intent(in) :: Bmat

       real*8 :: Nout

       integer :: i,j
       real*8 :: Nsq,Y0,a,qq
       complex*16 :: cX0c
       complex*16, dimension(nh) :: X0c
       real*8, dimension(nd) :: avec
       real*8, dimension(nh,nh) :: X0Mat
       real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
       complex*16, dimension(nh) :: test

       test(:) = complex(1.d0,1.d0)

       qq=q

       Nsq=fun_Nsq(nd+1,Bmat)

       call extractA(nd,Bmat,Amat,avec,a)
       call diagonalization(nd,Amat,LambdaMat,Tmat)
     
       Y0=int_Y0(nd,LambdaMat)
!       write(111,*) "Y0: ", Y0
!       write(111,*) "Nsq: ", Nsq

!       write(*,*) cvec

       X0Mat=int_XnMat(nd,0,a,avec,Amat,qq)
       X0c=matmul(X0Mat,cvec)
       cX0c=dot_product(cvec,X0c)
!       X0c=matmul(X0Mat,test)
!       cX0c=dot_product(test,X0c)

       Nout = Nsq*Y0*dreal(cX0c)
!       Nout = dot_product(cvec,cvec) 

!       write(*,*) "S00"
!       do i = 1,nh
!         write(*,*) Y0*X0Mat(i,:)
!       end do

      end function
 
      end module

