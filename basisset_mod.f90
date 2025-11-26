!**********************************************************************!
! Module containing all the things needed with to work with the        !
! Hermite polynomials                                                  !
!**********************************************************************!

      module basisset_module

      use constants

      implicit none

      !Number of Hermite polynomials
      integer,parameter,public :: nh = 5  
      !Maximum order of x for the Hermite pol. in database
      integer, parameter, public :: max_x = nh+1  
      !Database of the coefficients
      real*8, dimension(nh,max_x), public :: Mherm

      private
      public :: GenHermMat

      contains

!.....Matrix of Hermite polynomials coefficients........................ 
  
      subroutine GenHermMat()
      ! Mherm: matrix of hermite polynomials coefficients

       integer :: i,j

       Mherm(:,:) = 0.d0

       ! The coefficients are listed for INCREASING ORDER of x power

       Mherm(1,1) = 1.d0
       Mherm(2,2) = 2.d0 

       do i = 3,nh
         Mherm(i,1) =-Mherm(i-1,2)
         do j = 2,max_x-1
           Mherm(i,j) = 2*Mherm(i-1,j-1) -j*Mherm(i-1,j+1)
         end do
       end do

      end subroutine



      end module
