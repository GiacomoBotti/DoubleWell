!*********************************************************************!
! Module containing all the things pertaining potential               !
!*********************************************************************!
      module potential_module

      use constants

      implicit none

      !real*8, parameter,public  :: gamma_const = 0.5d0
      real*8, parameter,public  :: gamma_const = 0.0d0
      !real*8, parameter,public  :: eta_const = 10e+30!1.3544d0 
      real*8, parameter,public  :: eta_const = 1.3544d0 
      !real*8, parameter,public  :: sigma_const = 1.d0
      real*8, parameter,public  :: sigma_const = -1.d0
  
      ! Bath potential matrix dimensions
      integer,parameter,public :: nv = 1   
      real*8, dimension(nv,nv),public :: Vmat 

      private
      public :: matrix_pot
      contains

!.....MATRIX POTENTIAL..................................................
     
      subroutine matrix_pot()
      ! Vmat: matrix part of the potential
      
       integer :: i
       real*8, dimension(nv,nv) :: Vharm
       real*8, dimension(nv,nv) :: Vcoupl
      
       Vharm(:,:) = 0.d0
       Vcoupl(:,:) = 0.d0

       Vharm(1,1) = 1.d0
    
       do i = 2,nv
         Vharm(i,i) = 1.d0
         Vcoupl(i,i-1) = 1.d0
         Vcoupl(i-1,i) = 1.d0
       end do

       Vmat = 0.5d0*Vharm + 0.5d0*gamma_const*Vcoupl

       write(*,*) "Matrix potential generated"
 
      end subroutine

      end module
