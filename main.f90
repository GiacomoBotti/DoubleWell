!***********************************************************************!
! Fortran code to evolve a multiconfigurational wavefunction using      !
! McLachlan variational principle, applied to a double well             !
!***********************************************************************!

      program doublewell

      use constants
      use check_module

      implicit none

! TO BE SURE: GENERATE POTENTIAL MATRIX HERE
      call matrix_pot() 

!.....Check Diagonalization.............................................
!      call check_diagonalization(5)
!.....Check Momenta Matrix..............................................
!      call check_momenta(5)
!.....Check Matrix Potential............................................
!      call check_vmat()
!.....Check Quadratic Forms.............................................
!      call check_quad(nv)
!.....Check Tu powers...................................................
!      call check_Tupow(nv)
!.....Check y1 powers...................................................
!      call check_ypow(nv)
!.....Check V0 polynomials..............................................
      call check_V0pol(nv)

      end program

