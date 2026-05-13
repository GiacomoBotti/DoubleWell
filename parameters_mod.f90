!**********************************************************************!
! Module containing all the parameters needed for DOUBLEWELL           !
!**********************************************************************!

      module parameters_module

      implicit none
      save

      ! Bath potential matrix dimensions
      integer,parameter,public :: nv = 1   
      !Number of Hermite polynomials
      integer,parameter,public :: nh = 3 
      !Maximum order of x for the Hermite pol. in database
      integer, parameter, public :: max_x = nh+1  
      !Maximum order of y momenta
      integer, parameter, public :: maxorder = 8
      ! Gamma: coupling constant
      real*8, parameter,public  :: gamma_const = 0.5d0
!      real*8, parameter,public  :: gamma_const = 0.0d0
      ! Eta: quartic constant
!      real*8, parameter,public  :: eta_const = 10e+30!1.3544d0 
      real*8, parameter,public  :: eta_const = 1.3544d0 
      ! Sigma: active coordinate quadratic constant
!      real*8, parameter,public  :: sigma_const = 1.d0
      real*8, parameter,public  :: sigma_const = -1.d0
!      real*8, parameter,public  :: sigma_const = 0.d0
      ! Kappa: bath coordinates quadratic constant 
!      real*8, parameter,public  :: kappa_const = 0.d0
      real*8, parameter,public  :: kappa_const = 1.d0

      real*8, public :: coalc
      real*8, dimension(nh), public :: coalvec

      contains 

!.....DYNAMICS SETUP....................................................

      subroutine dynamics_setup(coalson)
      ! coalson: flag for using gaussian average
      integer, intent(in) :: coalson

       coalvec(:) = 1.d0
       coalc = 0.d0

       if (coalson.eq.1) then
         coalvec(:) = 0.d0
         coalc = 1.d0
         write(*,*) "WATCH OUT, YOU OPTED FOR GAUSSIAN AVERAGE"
         write(*,*) "coalc:", coalc
         write(*,*) "coalvec", coalvec
       end if

      end subroutine
      
      end module 
