!**********************************************************************!
! Module that contains the functions that assemble the polynomials of  !
! the momenta, as a step to solve the integrals present in the Karplus !
! effective potential parameters                                       !
!**********************************************************************!

      module polynomials_module

      use quadratic_module
      use diagonal_module 

      implicit none

      private
      public :: tildeP00M

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

      end module 
