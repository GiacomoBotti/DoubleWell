!**********************************************************************!
! Module containing everything pertaining the Karplus diagonalization  !
! of the multidimensional Gaussian                                     !
!**********************************************************************!

      module diagonal_module

      implicit none
      private
      public :: diagonalization

      contains

!.....Diagonalization with LAPACK.......................................

      subroutine diagonalization(nd,Amat,LambdaMat,Tmat) 
      ! nd: dimension of the matrices
      ! Amat: Gaussian width matrix (precision)
      ! LambdaMat: diagonal matrix
      ! Tmat: eigenvector matrix

       integer, intent(in) :: nd
       complex*16, dimension(nd,nd), intent(in) :: Amat
 
       real*8, dimension(nd,nd), intent(out) :: LambdaMat,Tmat

       integer :: i,n,lwork,info
       real*8, dimension(nd) :: eigenvalues
       real*8, dimension(3*nd-1) :: work

       external DSYEV

       n = nd
       lwork = 3*nd-1
       Tmat = real(Amat)

       call DSYEV('V','L',nd,Tmat,nd,eigenvalues,work,lwork,info) 

       LambdaMat(:,:) = 0.d0

       do i = 1,nd
          LambdaMat(i,i) = eigenvalues(i)
       end do
      
      end subroutine

      end module 
