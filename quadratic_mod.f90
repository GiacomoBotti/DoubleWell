!**********************************************************************!
! Module that contains the functions that compute the INTEGRALS of the !
! quadratic forms, i.e. the quadratic forms of the polynomials         !
!**********************************************************************!

      module quadratic_module
      
      use diagonal_module
      use potential_module

      implicit none
 
      private
      public :: fun_qVq,fun_uWu

      contains

!.....qVq...............................................................

      function fun_qVq(nd,q) result(qVq)
      ! nd: dimensions
      ! Vmat: potential matrix
      ! q: gaussian variational center vector
      ! qVq: result (x gaussian normalization...)

       integer,intent(in) :: nd
       real*8, dimension(nd), intent(in) :: q

       real*8 :: qVq
       real*8, dimension(nd) :: Vq
      
       Vq = matmul(Vmat,q)
       qVq = dot_product(q,Vq)

      end function 

!.....uWu=uTVTu.........................................................

      function fun_uWu(nd,Tmat,MomMat) result(uWu)
      ! nd: dimensions
      ! Tmat: eigenvector matrix
      ! MomMat: matrix of the momenta
      ! uWu: polynomial of the momenta
       
       integer,intent(in) :: nd
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 
 
       integer :: i
       real*8 :: uWu
       real*8, dimension(nd,nd) :: W1,W2

       W1 = matmul(Vmat,Tmat)
       W2 = matmul(transpose(Tmat),W1)

       uWu=0.d0

       do i = 1,nd
         uWu=uWu + W2(i,i)*MomMat(2,i)
       end do

      end function

      end module
       

