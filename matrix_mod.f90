!**********************************************************************!
! Module containing everything pertaining matrix manipulation and      !
! creation                                                             !
!**********************************************************************!

      module matrix_module

      implicit none

      integer,public :: maxorder = 8
      private
      public :: diagonalization,momenta,extractA

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

!.....Momenta of a factorized central normal distribution...............

      function momenta(nd,LambdaMat) result(MomMat)
      ! nd: dimension of the matrices
      ! LambdaMat: diagonal matrix
      ! MomMat: matrix of the momenta
      ! Each column of MomMat contains the i-th component momenta

       integer, intent(in) :: nd

       integer :: i,j
       real*8, dimension(nd,nd), intent(in) :: LambdaMat 
       real*8, dimension(maxorder,nd) :: MomMat 
 
       real*8, dimension(8) :: momcoeff
 
       if(size(momcoeff).lt.maxorder) then
         write(*,*) "O-----------------------------------------O"
         write(*,*) "|                  ERROR                  |"
         write(*,*) "+-----------------------------------------+"
         write(*,*) "| I don't have enough coefficients to     |"
         write(*,*) "| generate the momenta, please check      |"
         write(*,*) "| momcoeff in diagonal_mod.f90            |"      
         write(*,*) "O-----------------------------------------O"
         stop
       end if

       momcoeff = [0.d0,1.d0,0.d0,3.d0,0.d0,15.d0,0.d0,105.d0 ]
  
       !write(*,*) "Coefficients for the momenta"
       !write(*,*) momcoeff

       do i = 1,nd
          do j = 1,maxorder,2
             MomMat(j,i) = 0.d0
             MomMat(j+1,i) = momcoeff(j+1)/(2*LambdaMat(i,i))**(j+1) 
          end do
       end do      
 
      end function

!.....Extract Amat and avec.............................................

      subroutine extractA(nd,Bmat,Amat,avec,a)
      ! nd: dimension of the matrices
      ! Bmat: total gaussian width matrix
      ! Amat: bath gaussian width matrix
      ! avec: system-bath gaussian width vector
      
       integer, intent(in) :: nd
       real*8, dimension(nd+1,nd+1), intent(in) :: Bmat

       real*8, intent(out) ::  a
       real*8, dimension(nd), intent(out) ::  avec
       real*8, dimension(nd,nd), intent(out) :: Amat 

       integer :: i,j

       a = Bmat(1,1) 

       do i = 1,nd
         avec(i) = Bmat(1,i+1)
         do j = 1,nd
           Amat(i,j) = Bmat(i+1,j+1)
         end do
       end do

      end subroutine

      end module 
