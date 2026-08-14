!*********************************************************************!
! Module containing all the things pertaining potential               !
!*********************************************************************!
      module potential_module

      use constants
      use parameters_module

      implicit none

      real*8, dimension(nv,nv),public :: Vmat 

      private
      public :: matrix_pot,write_potential2D
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

       Vmat = 0.5d0*kappa_const*Vharm + 0.5d0*bath_const*Vcoupl
!       Vmat = 0.5d0*Vharm + 0.5d0*gamma_const*Vcoupl

       write(*,*) "Matrix potential generated"
 
      end subroutine

!.....PLOT POTENTIAL....................................................

      subroutine write_potential2D()
      
       integer :: i,j
       real*8 :: x,y,grid,V

       x = lwb
       y = lwb

       grid = -2.d0*lwb/100.d0

       open(unit=999,file='pot.dat',status='replace',action='write')

       do i = 1, 101
         y = lwb +(i-1)*grid
         do j = 1, 101
           x = lwb + (j-1)*grid
           V = x**4/16.d0/eta_const+sigma_const*x**2/2.d0 &
             &+ gamma_const*x*y +&
             & kappa_const*y**2/2.d0 
           write(999,'(3ES24.16)') x, y, V
         end do
         write(999,*) 
         flush(999)
       end do
  
      close(999)
 
      end subroutine

      end module
