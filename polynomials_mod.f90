!**********************************************************************!
! Module that contains the functions that assemble the polynomials of  !
! the momenta, as a step to solve the integrals present in the Karplus !
! effective potential parameters                                       !
!**********************************************************************!

      module polynomials_module

      use quadratic_module
      use matrix_module 
      use potential_module 
      use ypowers_module 

      implicit none

      private
      public :: tildeP00M,fun_Y00

      contains

!.....tildeP00M = uWu + QVQ + (uZQ + QRu)...............................

      function tildeP00M(nd,q,Tmat,MomMat) result(tP00M)
      ! nd: dimensions
      ! q: gaussian variational center vector
      ! Tmat: eigenvector matrix
      ! MomMat: matrix of the momenta
      ! tP00M: matricial part of Y00

       integer,intent(in) :: nd
       real*8, dimension(nd), intent(in) :: q
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 
       real*8 :: tP00M 
 
       real*8 :: uWu,QVQ,uZQ,QRu

       uWu = fun_uWu(nd,Tmat,MomMat)
       QVQ = fun_qVq(nd,q)
       !uZQ = fun_uZQ(nd,q,Tmat,MomMat) 
       !QRu = fun_QRu(nd,q,Tmat,MomMat) 
       uZQ = 0.d0 
       QRu = 0.d0 

       tP00M = uWu + QVQ + uZQ + QRu
 
      end function

!.....Y00 = y4/16eta + tildeP00M...................................

      function fun_Y00(nd,q,Tmat,MomMat) result(Y00)
      ! nd: dimensions
      ! q: gaussian variational center vector
      ! Tmat: eigenvector matrix
      ! MomMat: matrix of the momenta
      ! tP00: solution of Y00

       integer,intent(in) :: nd
       real*8, dimension(nd), intent(in) :: q
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 
       real*8 :: Y00 
 
       real*8 :: y4,tP00M

       !write(*,*) Tmat
       !write(*,*) MomMat

       y4=fun_y4(nd,q,Tmat,MomMat)
       tP00M=tildeP00M(nd,q,Tmat,MomMat)

       Y00 = y4/(16.d0*eta_const) + tP00M 

      end function
 

      end module 
