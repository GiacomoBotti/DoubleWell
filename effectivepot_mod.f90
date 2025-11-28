!**********************************************************************!
! Module that computes the effective potential terms aka the           !
! expectation values of V, dV and ddV                                  !
!**********************************************************************!

      module effectivepot_module

      use constants
      use basisset_module
      use potential_module

      implicit none
 
      private
      public :: fun_V0

!.....EFFECTIVE POTENTIAL V0 = <V>......................................

      function fun_V0(nd,q,cvec,Bmat) result(V0)
      ! nd: bath dimensions
      ! q: x Gaussian center
      ! cvec: vector of the coefficients
      ! Bmat: total Gaussian width matrix
       integer, intent(in) :: nd
       real*8, intent(in) :: q
       complex*16, dimension(nh), intent(in) :: cvec
       real*8, dimension(nd+1,nd+1), intent(in) :: Bmat

       complex*16 :: V0

      end function

      end module
