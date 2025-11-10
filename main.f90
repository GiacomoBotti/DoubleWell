!***********************************************************************!
! Fortran code to evolve a multiconfigurational wavefunction using      !
! McLachlan variational principle, applied to a double well             !
!***********************************************************************!

      program doublewell

      use constants
      use check_module

      implicit none

!.....Check Diagonalization.............................................

      call check_diagonalization(5)

      end program

