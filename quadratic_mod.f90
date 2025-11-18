!**********************************************************************!
! Module that contains the functions that compute the INTEGRALS of the !
! quadratic forms, i.e. the quadratic forms of the polynomials         !
!**********************************************************************!

      module quadratic_module
      
      use diagonal_module
      use potential_module

      implicit none
 
      private
      public :: fun_qVq,fun_uWu,fun_uZQ,fun_QRu

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

!.....uZQ=uTVQ..........................................................

      function fun_uZQ(nd,q,Tmat,MomMat) result(uZQ)
      ! nd: dimensions
      ! q: gaussian variational center vector
      ! Tmat: eigenvector matrix
      ! MomMat: matrix of the momenta
      ! uZQ: polynomial of the momenta
       
       integer,intent(in) :: nd
       real*8, dimension(nd), intent(in) :: q
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 
 
       integer :: i,j
       real*8 :: uZQ
       real*8, dimension(nd,nd) :: Z

       !Z = matmul(transpose(Tmat),Vmat)

       uZQ=0.d0

       !do i = 1,nd
       !  do j =1,nd
       !    uZQ=uZQ + MomMat(1,i)*Z(i,j)*q(j)
       !  end do
       !end do

      end function

!.....QRu=QVTu..........................................................

      function fun_QRu(nd,q,Tmat,MomMat) result(QRu)
      ! nd: dimensions
      ! q: gaussian variational center vector
      ! Tmat: eigenvector matrix
      ! MomMat: matrix of the momenta
      ! QRu: polynomial of the momenta
       
       integer,intent(in) :: nd
       real*8, dimension(nd), intent(in) :: q
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 
 
       integer :: i,j
       real*8 :: QRu
       real*8, dimension(nd,nd) :: R

       !R = matmul(Vmat,Tmat)

       QRu=0.d0

       !do i = 1,nd
       !  do j =1,nd
       !    QRu=QRu + q(i)*R(i,j)*MomMat(1,j)
       !  end do
       !end do

      end function

      end module
       

