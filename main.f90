!***********************************************************************!
! Fortran code to evolve a multiconfigurational wavefunction using      !
! McLachlan variational principle, applied to a double well             !
!***********************************************************************!

      program doublewell

      use constants
      use parameters_module
      use potential_module
      use basisset_module
      use eofmotion_module
      use check_module
      use evolution_module
      use BOT_module
      use observable_module

      implicit none

      integer :: i,j,coalson,scaling,frozen,stationary
      real*8 :: t0,t1
      real*8,dimension(nv+1) :: masses !masses vector
      real*8,dimension(nv+1) :: q0 !inital centers vector
      real*8,dimension(nv+1) :: p0 !initial momenta vector
      real*8,dimension(nv+1) :: qeq !equilibrium centers vector
      real*8,dimension(nv+1) :: peq !equilibrium momenta vector
      complex*16,dimension(nh) :: c0 !initial coefficient vector
      complex*16,dimension(nh) :: ceq !equilibrium coefficient vector
      complex*16,dimension(nv+1,nv+1) :: Bcmplx !initial width matrix
      complex*16,dimension(nv+1,nv+1) :: Beq !equilibrium width matrix
      integer*8,dimension(4) :: trj

      namelist /setup/ trj, coalson, scaling, frozen, stationary
      namelist /inp_mass/ masses
      namelist /equilibrium/ qeq, peq, ceq, Beq
      namelist /initial/ q0, p0, Bcmplx 

      call execute_command_line('cat banner.txt')

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "|               MAIN CODE EXECUTION                 |"
      write(*,*) "+---------------------------------------------------+"

      open(unit=111,file='wfx_BOT.dat',status='unknown',action='write')
      open(unit=222,file='wfy_BOT.dat',status='unknown',action='write')
      open(unit=2222,file="input",status="old",action="read")

!.....Input reading.....................................................

      trj = [0, 1, 100, 1]
      coalson = 0
      scaling = 0
      frozen = 0
      stationary = 0
      read(2222,nml=setup)

!.....Potential parameters reading......................................

      call potential_setup() 

!.....SETUP.............................................................

! TO BE SURE: GENERATE POTENTIAL MATRIX HERE
      call matrix_pot() 
! TO BE SURE: GENERATE HERMITE COEFFICIENT MATRIX HERE
      call GenHermMat()
! DYNAMICS SETUP IT'S IMPORTANT
      call dynamics_setup(coalson,scaling,frozen,stationary)

!.....Print potential constants.........................................

      write(*,*) "Potential constants"

      write(*,*) "Eta: ", eta_const
      write(*,*) "Sigma: ", sigma_const
      write(*,*) "Gamma: ", gamma_const
      write(*,*) "Kappa: ", kappa_const

      write(*,*) "+---------------------------------------------------+"

!.....Define masses vector..............................................

      write(*,*) "Masses vector:"

      masses(:) = 1.d0
      masses(1) = 1.5d0
      masses(2) = 4.5d0
      read(2222,nml=inp_mass)
      do i = 1,nv+1
        write(*,*) masses(i)
      end do

      call MassesMat(masses)

!.....Define initial conditions.........................................

      ! Boring default
      ceq(:) =1.d0
      qeq(:) =0.d0
      peq(:) =0.d0
      Beq(:,:) = 0.d0
      do i = 1, nv+1
        Beq(i,i) = dsqrt(masses(i))
      end do

      read(2222,nml=equilibrium)

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Equilibrium Coefficients:"

      do i = 1,nh
        write(*,*) ceq(i) 
      end do

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Initial Gaussian Width Matrix:"

      do i = 1,nv+1
        write(*,*) Beq(i,:)
      end do

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Initial q and p:"
       
      do i = 1,nv+1
        write(*,*) qeq(i), peq(i)
      end do

!.....Basis Projection..................................................

      ! Boring defaults
      q0(:) = qeq(:) 
      p0(:) = peq(:)
      Bcmplx(:,:) = Beq(:,:)
      c0(:) = ceq(:) ! safeguard

      read(2222,nml=initial)

      c0 = c_update(nv,q0,p0,qeq,peq,ceq,Bcmplx,Beq) 
    
      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Projected coefficients"

      do i = 1,nh
        write(*,*) c0(i) 
      end do

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Projected Gaussian Width Matrix:"

      do i = 1,nv+1
        write(*,*) Bcmplx(i,:)
      end do
    
      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Projected q and p:"
       
      do i = 1,nv+1
        write(*,*) q0(i), p0(i)
      end do
      
      call plot_wfn(nv,qeq,peq,ceq,Beq,1.d0)
      call plot_wfn(nv,q0,p0,c0,Bcmplx,1.d0)

!.....Check Diagonalization.............................................
!      call check_diagonalization(nv)
!.....Check Momenta Matrix..............................................
!      call check_momenta(nv)
!.....Check A matrix extraction.........................................
!      call check_Amat(nv)
!.....Check Matrix Potential............................................
!      call check_vmat()
!.....Check Hermite Matrix..............................................
!      call check_hermmat()
!.....Check Y0..........................................................
!      call check_Y0(nv)
!.....Check XnMat.......................................................
!      call check_XnMat(nv)
!.....Check V0..........................................................
!      call check_V0(nv,c0)
!.....Check Norm........................................................
!      call check_norm(nv,c0)
!.....Check NN..........................................................
!      call check_shiftedoverlap(nv)
!.....Check K and S as SG...............................................
!      call check_KSnum(nv)
!.....Check projection..................................................
!      call check_projection(nv,qeq,peq,ceq,Bcmplx)
!.....Evolution.........................................................

      close(2222)
      write(*,*) "WE ARE RUNNING"
      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Start    ", "Stop    ", "Lenght  ", "Print "  
      write(*,*) trj

      call cpu_time(t0)
      call bot_evo(nv,trj,q0,p0,c0,Bcmplx,qeq,peq,ceq,Beq)
      call coherent_calc(nv,trj,q0,p0,masses,c0,Bcmplx)
      call cpu_time(t1)
      write(*,*) "End of a successful run"
      write(*,*) "Have a nice day"
      write(*,*) "I took ", t1-t0, "time"
      write(*,*) "+---------------------------------------------------+"

      close(111)
      close(222)

      end program

