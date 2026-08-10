!**********************************************************************!
! Module containing all the parameters needed for DOUBLEWELL           !
!**********************************************************************!

      module parameters_module

      implicit none
      save

      ! Bath potential matrix dimensions
      integer,parameter,public :: nv = 1
      !Number of Hermite polynomials
      integer,parameter,public :: nh = 7 !15
      !Maximum order of x for the Hermite pol. in database
      integer, parameter, public :: max_x = nh+1  
      !Maximum order of y momenta
      integer, parameter, public :: maxorder = 8
      ! Gamma: coupling constant
      real*8, public  :: gamma_const = 0.5d0
      ! Eta: quartic constant
      real*8, public  :: eta_const = 1.3544d0 
      ! Sigma: active coordinate quadratic constant
      real*8, public  :: sigma_const = -1.d0
      ! Kappa: bath coordinates quadratic constant 
      real*8, public  :: kappa_const = 1.d0

      integer, public :: coalmode
      real*8, public :: coalc,ffact
      real*8, dimension(nh), public :: coalvec
      real*8, dimension(nv+1), public :: scalvec
      real*8, dimension(nv+1,nv+1), public :: scalmat
      real*8, public :: lwb=-12.d0
      real*8, public :: hgb=12.d0
      integer*8, public :: gstep=500

      contains 

!.....POTENTIAL SETUP...................................................

      subroutine potential_setup()

      namelist /pot_param/ eta_const,sigma_const,gamma_const,kappa_const

       read(2222,nml=pot_param)
      
      end subroutine

!.....GRID SETUP........................................................

      subroutine grid_setup()

      namelist /grid/ lwb,hgb,gstep

       read(2222,nml=grid)

      end subroutine

!.....DYNAMICS SETUP....................................................

      subroutine dynamics_setup(coalson,scaling,frozen,stationary)
      ! coalson: flag for using gaussian average
      ! scaling: flag for scaled parameters dynamics
      ! frozen: flag for frozen gaussian dynamics
       integer, intent(in) :: coalson,scaling,frozen,stationary

       integer :: i
       real*8, dimension(nv+1) :: scalv
       real*8, dimension(nv+1,nv+1) :: scalm

       namelist /scaling_factors/ scalv,scalm

       scalvec(:) = 1.d0
       scalmat(:,:) = 1.d0

       ! Default
       scalv = scalvec
       scalm = scalmat
       scalv(1) = 0.d0
       scalm(1,1) = 0.d0
       
       read(2222,nml=scaling_factors)

       coalvec(:) = 1.d0
       coalc = 0.d0
       coalmode = 1


       ffact = 0.d0

       if (coalson.ne.0) then
         coalvec(:) = 0.d0
         coalc = 1.d0
         coalmode = coalson
         write(*,*) "WATCH OUT, YOU OPTED FOR GAUSSIAN AVERAGE"
         write(*,*) "coalc:", coalc
         write(*,*) "coalvec", coalvec
         write(*,*) "coalmode", coalmode
       end if

       if (scaling.eq.1) then
         scalvec = scalv
         scalmat = scalm
         write(*,*) "WATCH OUT, YOU OPTED FOR SCALED DYNAMICS"
         write(*,*) "(It works differently for SCP and VTV prop)"
         write(*,*) "scalvec:", scalvec
         write(*,*) "scalmat:"
         do i = 1,nv+1
           write(*,*) scalmat(i,:)
         end do
       end if

       if (frozen.eq.1) then
         scalmat(:,:) = 0.d0
         write(*,*) "WATCH OUT, YOU OPTED FOR FROZEN G. DYNAMICS"
         write(*,*) "scalmat:"
         do i = 1,nv+1
           write(*,*) scalmat(i,:)
         end do
         ffact = 1.d0
       end if

       if (stationary.eq.1) then
         scalvec(:) = 0.d0
         scalmat(:,:) = 0.d0
         write(*,*) "WATCH OUT, YOU OPTED FOR STATIONARY B. DYNAMICS"
         write(*,*) "(It only works for SCP propagator)"
         write(*,*) "            (right now)           "
         write(*,*) "scalmat:"
         do i = 1,nv+1
           write(*,*) scalmat(i,:)
         end do
       end if

      end subroutine
      
      end module 
