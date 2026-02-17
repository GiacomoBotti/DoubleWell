!**********************************************************************!
! Module containing everything pertaining matrix manipulation and      !
! creation                                                             !
!**********************************************************************!

      module matrix_module

      implicit none

      integer,public :: maxorder = 8
      private
      public :: diagonalization,momenta,extractA,determinant,trace
      public :: extracttildeA,det_cmplx

      contains

!.....Diagonalization with LAPACK.......................................

      subroutine diagonalization(nd,Amat,LambdaMat,Tmat) 
      ! nd: dimension of the matrices
      ! Amat: Gaussian width matrix (precision)
      ! LambdaMat: diagonal matrix
      ! Tmat: eigenvector matrix

       integer, intent(in) :: nd
       real*8, dimension(nd,nd), intent(in) :: Amat
 
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
             MomMat(j+1,i) = momcoeff(j+1)/(2*LambdaMat(i,i))**((j+1)/2) 
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

!.....Determinant of a matrix...........................................
 
      function determinant(ndim,Bmat) result(Bdet)
      ! nd: dimension of the Bmat matrix CAREFUL!
      ! Bmat: matrix
      ! Bidet: I miss it
       integer, intent(in) :: ndim
       real*8, dimension(ndim,ndim), intent(out) :: Bmat

       integer :: i
       real*8 :: Bdet
       real*8, dimension(ndim,ndim) :: DiagMat,Tmat 

       call diagonalization(ndim,Bmat,DiagMat,Tmat)
        
       Bdet=DiagMat(1,1)
!       write(*,*) DiagMat(1,:)

       do i = 2,ndim
         Bdet = Bdet*DiagMat(i,i)
!       write(*,*) DiagMat(i,:)
       end do
  
      end function
           
!.....Extract complex Amat and avec.....................................

      subroutine extracttildeA(nd,tildeBmat,tildeAmat,tildeavec,tildea)
      ! nd: dimension of the matrices
      ! tildeBmat: total gaussian width matrix (COMPLEX)
      ! tildeAmat: bath gaussian width matrix (COMPLEX)
      ! tildeavec: system-bath gaussian width vector (COMPLEX)
      
       integer, intent(in) :: nd
       complex*16, dimension(nd+1,nd+1), intent(in) :: tildeBmat

       complex*16, intent(out) ::  tildea
       complex*16, dimension(nd), intent(out) ::  tildeavec
       complex*16, dimension(nd,nd), intent(out) :: tildeAmat 

       integer :: i,j

       tildea = tildeBmat(1,1) 

       do i = 1,nd
         tildeavec(i) = tildeBmat(1,i+1)
         do j = 1,nd
           tildeAmat(i,j) = tildeBmat(i+1,j+1)
         end do
       end do

      end subroutine

!.....Compute the trace of a generic matrix.............................

      function trace(nd,Mat) result(traceMat)
      ! nd: matrix dimension
      ! Mat: matrix
      ! traceMat: trace of the matrix
       integer, intent(in) :: nd
       complex*16, dimension(nd,nd), intent(in) :: Mat

       complex*16 :: traceMat

       integer :: i

       traceMat = (0.d0,0.d0)

       do i=1,nd
         traceMat = traceMat + Mat(i,i)
       end do

      end function

!.....Diagonalization with LAPACK (complex).............................

      function det_cmplx(nd,Amat) result(detAmat) 
      ! nd: dimension of the matrices
      ! Amat: complex matrix matrix (precision)
      ! detAmat: determinant

       integer, intent(in) :: nd
       complex*16, dimension(nd,nd), intent(in) :: Amat
 
       complex*16 :: detAmat 


       integer, dimension(nd) :: ipiv   ! pivot indices
       real*8 :: detL,detP
       complex*16 :: detU
       complex*16, dimension(nd) :: work  ! work array for LAPACK
       complex*16, dimension(nd,nd) :: Awork
       integer :: i,n, info

       external ZGETRF

       n = nd
       ! Store A in Ainv to prevent it from being overwritten by LAPACK
       Awork = Amat

       call ZGETRF(n,n,Awork,n,ipiv,info) 

       if (info /= 0) then
         !write(*,*) "DGETRF info : ",info
         stop 'diag_complx Matrix is numerically singular!'
       end if

       ! Determinants of the decomposition 
       detU = complex(1.d0,0.d0)
       detL = 1.d0
       detP = 1.d0

       do i = 1,nd
         detU = detU*Awork(i,i)
         if(ipiv(i).ne.i) then
            detP = - detP
         end if
       end do
 
       ! Total determinant

       detAmat = detP*detU*detL

      end function
      end module 
