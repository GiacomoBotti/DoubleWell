!**********************************************************************!
! Module containing all the functions needed to compute the integrals  !
!**********************************************************************!

       module integrals_mod

       use constants

       implicit none

       private
       public :: int_Y0

       contains

!......Y0 integral......................................................
       
       function int_Y0(nd,LambdaMat) result(Y0)
       ! nd: dimensions
       ! LambdaMat: Diagonalized Bath Gaussian width
        integer, intent(in) :: nd
        real*8, dimension(nd,nd), intent(in) :: LambdaMat

        real*8 :: Y0

        integer :: i
        real*8 :: det

        det = LambdaMat(1,1)

        do i = 2,nd
           det = det*LambdaMat(i,i)
        end do

        Y0 = dsqrt(pi/det)
       
       end function

       end module
     
