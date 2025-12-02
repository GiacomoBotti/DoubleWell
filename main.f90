!***********************************************************************!
! Fortran code to evolve a multiconfigurational wavefunction using      !
! McLachlan variational principle, applied to a double well             !
!***********************************************************************!

      program doublewell

      use constants
      use potential_module
      use basisset_module
      use eofmotion_module
      use check_module

      implicit none

      integer :: i
      real*8,dimension(nv+1) :: masses !Masses vector
      complex*16,dimension(nh) :: c0 !Initial coefficient vector

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "|               MAIN CODE EXECUTION                 |"
      write(*,*) "+---------------------------------------------------+"

! TO BE SURE: GENERATE POTENTIAL MATRIX HERE
      call matrix_pot() 
! TO BE SURE: GENERATE HERMITE COEFFICIENT MATRIX HERE
      call GenHermMat()

!.....Define initial conditions.........................................

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Initial coefficients:"

      c0(:) = 0.d0
      c0(1) = 1.d0!/dsqrt(1.6435395517174083d0)
    
      do i = 1,nh
!        c0(i) = 1.d0/nh
        write(*,*) c0(i) 
      end do

      write(*,*) "+---------------------------------------------------+"

!.....Define masses vector..............................................

      write(*,*) "Masses vector:"

      do i = 1,nv+1
        masses(i) = 1.d0
        write(*,*) masses(i)
      end do

      call MassesMat(masses)

      write(*,*) "+---------------------------------------------------+"

!.....Check Diagonalization.............................................
!      call check_diagonalization(nv)
!.....Check Momenta Matrix..............................................
!      call check_momenta(nv)
!.....Check A matrix extraction.........................................
!      call check_Amat(nv)
!.....Check Matrix Potential............................................
!      call check_vmat()
!.....Check Hermite Matrix..............................................
!      call check_hermmat()
!.....Check Y0..........................................................
!      call check_Y0(nv)
!.....Check XnMat.......................................................
!      call check_XnMat(nv)
!.....Check V0..........................................................
!      call check_V0(nv,c0)
!.....Check Norm........................................................
      call check_norm(nv,c0)

      end program

