!***********************************************************************!
! Fortran code to evolve a multiconfigurational wavefunction using      !
! McLachlan variational principle, applied to a double well             !
!***********************************************************************!

      program doublewell

      use constants
      use potential_module
      use basisset_module
      use check_module

      implicit none

! TO BE SURE: GENERATE POTENTIAL MATRIX HERE
      call matrix_pot() 
! TO BE SURE: GENERATE HERMITE COEFFICIENT MATRIX HERE
      call GenHermMat()
!.....Check Diagonalization.............................................
      call check_diagonalization(nv)
!.....Check Momenta Matrix..............................................
      call check_momenta(nv)
!.....Check A matrix extraction.........................................
      call check_Amat(nv)
!.....Check Matrix Potential............................................
      call check_vmat()
!.....Check Hermite Matrix..............................................
      call check_hermmat()
!.....Check Y0..........................................................
      call check_Y0(nv)
!.....Check XnMat..........................................................
      call check_XnMat(nv)

      end program

