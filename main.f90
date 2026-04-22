!***********************************************************************!
! Fortran code to evolve a multiconfigurational wavefunction using      !
! McLachlan variational principle, applied to a double well             !
!***********************************************************************!

      program doublewell

      use constants
      use potential_module
      use basisset_module
      use eofmotion_module
      use check_module
      use evolution_module
      use BOT_module
      use observable_module

      implicit none

      integer :: i,j
      real*8 :: t0,t1
      real*8,dimension(nv+1) :: masses !masses vector
      real*8,dimension(nv+1) :: q0 !inital centers vector
      real*8,dimension(nv+1) :: p0 !initial momenta vector
      real*8,dimension(nv+1) :: qeq !equilibrium centers vector
      real*8,dimension(nv+1) :: peq !equilibrium momenta vector
      complex*16,dimension(nh) :: c0 !initial coefficient vector
      complex*16,dimension(nh) :: ceq !equilibrium coefficient vector
      complex*16,dimension(nv+1,nv+1) :: Bcmplx !initial width matrix
      integer*8,dimension(3) :: trj

!      call print_double_well_banner()
      call execute_command_line('cat banner.txt')

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "|               MAIN CODE EXECUTION                 |"
      write(*,*) "+---------------------------------------------------+"

! TO BE SURE: GENERATE POTENTIAL MATRIX HERE
      call matrix_pot() 
! TO BE SURE: GENERATE HERMITE COEFFICIENT MATRIX HERE
      call GenHermMat()

!.....Print potential constants.........................................

      write(*,*) "Potential constants"

      write(*,*) "Eta: ", eta_const
      write(*,*) "Sigma: ", sigma_const
      write(*,*) "Gamma: ", gamma_const
      write(*,*) "Kappa: ", kappa_const

      write(*,*) "+---------------------------------------------------+"

!.....Define masses vector..............................................

      write(*,*) "Masses vector:"

      do i = 1,nv+1
        masses(i) = 1.d0
!        masses(i) = 1.1d0*i
        write(*,*) masses(i)
      end do

!      masses(2) = masses(1) !WATCH OUT

      call MassesMat(masses)

!.....Define initial conditions.........................................

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Equilibrium Coefficients:"

      ceq(:) =0.d0
      ceq(1) =1.0d0  !creal
    
      do i = 1,nh
!        ceq(i) = 1.d0/nh
        write(*,*) ceq(i) 
      end do

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Initial Gaussian Width Matrix:"

      do i = 1,nv+1
!         qeq(i) = i*dsqrt(2.d0)/3.d0
!         peq(i) = i*dsqrt(3.d0)/7.d0
!         Bcmplx(i,i) = (i+i)*(1+i/100.d0) + iu*(i+i)*(1+i/40.d0)/10.d0
         Bcmplx(i,i) = dsqrt(masses(i))
!         do j= i+1,nv+1
!           Bcmplx(i,j) = (i+j)/40.d0 + iu*(i+j)/30.d0
!           Bcmplx(j,i) = Bcmplx(i,j)
!         end do   
        write(*,*) Bcmplx(i,:)
      end do

!      qeq(1) = -2.d0*dsqrt(eta_const)
!      qeq(1) = -2.3199137915 
!      qeq(1) = -.23199137915 
!      qeq(1) = 0.d0
      qeq(1) = 1.d0
!      qeq(2) = 0.d0 

      peq(:) = 0.d0
!      peq(1) = 0.005d0

!.....Basis Projection..................................................

!      call check_projection(nv,qeq,peq,ceq,Bcmplx)

!      stop
      !q0(:) = 1.d0
      !q0(1) = 1.d0
      !p0(:) = 0.d0
      !c0 = c_update(nv,qeq,peq,q0,p0,ceq,Bcmplx,Bcmplx) 
      q0(:) = qeq(:) 
      p0(:) = peq(:)
      c0(:) = ceq(:)

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
!.....Evolution.........................................................

      write(*,*) "WE ARE RUNNING"
      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Start       ", "Stop       ", "Lenght     "  
      open(unit=2222,file="input",status="old",action="read")
      read(2222,*)
      read(2222,*) trj
      close(2222)
      !trj = [0,50,100]
      write(*,*) trj

      call cpu_time(t0)
      call bot_evo(nv,trj,q0,p0,c0,Bcmplx)
!      call coherent_calc(nv,trj,q0,p0,masses,c0,Bcmplx)
      call cpu_time(t1)
      write(*,*) "End of a successful run"
      write(*,*) "Have a nice day"
      write(*,*) "I took ", t1-t0, "time"
      write(*,*) "+---------------------------------------------------+"

      end program

