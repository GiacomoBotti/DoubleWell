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
      real*8, public  :: gamma_const = 0.5d0
!      real*8, parameter,public  :: gamma_const = 0.0d0
      ! Eta: quartic constant
!      real*8, parameter,public  :: eta_const = 10e+30!1.3544d0 
      real*8, public  :: eta_const = 1.3544d0 
      ! Sigma: active coordinate quadratic constant
!      real*8, parameter,public  :: sigma_const = 1.d0
      real*8, public  :: sigma_const = -1.d0
!      real*8, parameter,public  :: sigma_const = 0.d0
      ! Kappa: bath coordinates quadratic constant 
!      real*8, parameter,public  :: kappa_const = 0.d0
      real*8, public  :: kappa_const = 1.d0

      real*8, public :: coalc
      real*8, dimension(nh), public :: coalvec
      real*8, dimension(nv+1), public :: scalvec
      real*8, dimension(nv+1,nv+1), public :: scalmat

      contains 

!.....POTENTIAL SETUP...................................................

      subroutine potential_setup()

      namelist /pot_param/ eta_const,sigma_const,gamma_const,kappa_const


       !read(2222,*) !x^4/16/eta
       !read(2222,*) eta_const
       !read(2222,*) !x^2 coeff
       !read(2222,*) sigma_const
       !read(2222,*) !coupling const
       !read(2222,*) gamma_const
       !read(2222,*) !y^2 coeff
       !read(2222,*) kappa_const
       read(2222,nml=pot_param)
      
      end subroutine
       

!.....DYNAMICS SETUP....................................................

      subroutine dynamics_setup(coalson,scaling,frozen)
      ! coalson: flag for using gaussian average
      ! scaling: flag for scaled parameters dynamics
      ! frozen: flag for frozen gaussian dynamics
       integer, intent(in) :: coalson,scaling,frozen

       integer :: i

       coalvec(:) = 1.d0
       coalc = 0.d0

       scalvec(:) = 1.d0
       scalmat(:,:) = 1.d0

       if (coalson.eq.1) then
         coalvec(:) = 0.d0
         coalc = 1.d0
         write(*,*) "WATCH OUT, YOU OPTED FOR GAUSSIAN AVERAGE"
         write(*,*) "coalc:", coalc
         write(*,*) "coalvec", coalvec
       end if

       if (scaling.eq.1) then
         scalvec(1) = 0.d0
         scalmat(1,:) = 0.d0
         scalmat(:,1) = 0.d0
         write(*,*) "WATCH OUT, YOU OPTED FOR SCALED DYNAMICS"
         write(*,*) "(It only works for SCP propagator)"
         write(*,*) "scalvec:", scalvec
         write(*,*) "scalmat:"
         do i = 1,nv+1
           write(*,*) scalmat(i,:)
         end do
       end if

       if (frozen.eq.1) then
         scalmat(:,:) = 0.d0
         write(*,*) "WATCH OUT, YOU OPTED FOR FROZEN G. DYNAMICS"
         write(*,*) "(It only works for SCP propagator)"
         write(*,*) "            (right now)           "
         write(*,*) "scalmat:"
         do i = 1,nv+1
           write(*,*) scalmat(i,:)
         end do
       end if

      end subroutine
      
      end module 
