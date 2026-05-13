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
      complex*16,dimension(nv+1,nv+1) :: Beq !equilibrium width matrix
      integer*8,dimension(4) :: trj

      call execute_command_line('cat banner.txt')

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "|               MAIN CODE EXECUTION                 |"
      write(*,*) "+---------------------------------------------------+"

! TO BE SURE: GENERATE POTENTIAL MATRIX HERE
      call matrix_pot() 
! TO BE SURE: GENERATE HERMITE COEFFICIENT MATRIX HERE
      call GenHermMat()

      open(unit=111,file='wfx_BOT.dat',status='unknown',action='write')
      open(unit=222,file='wfy_BOT.dat',status='unknown',action='write')

!.....Print potential constants.........................................

      write(*,*) "Potential constants"

      write(*,*) "Eta: ", eta_const
      write(*,*) "Sigma: ", sigma_const
      write(*,*) "Gamma: ", gamma_const
      write(*,*) "Kappa: ", kappa_const

      write(*,*) "+---------------------------------------------------+"

!.....Define masses vector..............................................

      write(*,*) "Masses vector:"

      masses(1) = 1.5d0
      masses(2) = 4.5d0
      do i = 1,nv+1
!        masses(i) = 1.d0
!        masses(i) = 1.1d0*i
        write(*,*) masses(i)
      end do

!      masses(2) = masses(1) !WATCH OUT

      call MassesMat(masses)

!      write(*,*) InvMassMat

!.....Define initial conditions.........................................

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Equilibrium Coefficients:"

      ceq(:) =0.d0
!      ceq(1) =1.0d0  !creal

      ceq(1) = complex(0.45919157095202662,-8.14834718946430558E-003)
      ceq(2) = complex(-0.61978342037299328,2.46232212017700389E-002)
      ceq(3) = complex(0.48109453241877925,-4.73920543545408093E-002)
!      ceq(1) = complex(0.43502107680599428,-7.71944215178154303E-003)
!      ceq(2) = complex(-0.58715984345742667,2.33271272367176624E-002) 
!      ceq(3) = complex(0.45577113078183779,-4.48974759589270342E-002)

    
      do i = 1,nh
        write(*,*) ceq(i) 
      end do

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Initial Gaussian Width Matrix:"

!      Beq(1,1) = complex(1.d0,0.1d0)
      Beq(1,1) = complex(0.56d0,0.0d0)
      Beq(1,2) = complex(0.2d0,0.3d0)
      Beq(2,1) = complex(0.2d0,0.3d0)
      Beq(2,2) = complex(2.d0,4.d0)

!      Beq(1,1) = 0.5628366657
!      Beq(1,1) = 1.d0*dsqrt(masses(1))
!      Beq(1,2) = 0.d0
!      Beq(2,1) = 0.d0
!      Beq(2,2) = 1.d0*dsqrt(masses(2)) 

      do i = 1,nv+1
!         qeq(i) = i*dsqrt(2.d0)/3.d0
!         peq(i) = i*dsqrt(3.d0)/7.d0
!         Bcmplx(i,i) = (i+i)*(1+i/100.d0) + iu*(i+i)*(1+i/40.d0)/10.d0
!         Beq(i,i) = dsqrt(masses(i))
!         do j= i+1,nv+1
!           Bcmplx(i,j) = (i+j)/40.d0 + iu*(i+j)/30.d0
!           Bcmplx(j,i) = Bcmplx(i,j)
!         end do   
        write(*,*) Beq(i,:)
      end do

      qeq(1) = -2.d0
      qeq(2) = 0.2d0
      peq(1) = 0.1d0
      peq(2) = 0.3d0

!      qeq(1) = 1.d0
!      qeq(2) = 0.d0
!      peq(1) = 0.d0
!      peq(2) = 0.d0

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Initial q and p:"
       
      do i = 1,nv+1
        write(*,*) qeq(i), peq(i)
      end do

!.....Basis Projection..................................................

!      call check_projection(nv,qeq,peq,ceq,Bcmplx)

!      q0(:) = 0.d0
!      p0(:) = 0.d0
      q0(:) = qeq(:) 
      p0(:) = peq(:)
      Bcmplx(:,:) = Beq(:,:)
!      Bcmplx(1,1) = 0.5628366657 !Beq(:,:)
!      Bcmplx(1,1) = 0.56d0 
!      Bcmplx(1,1) = 0.46d0 
!      Bcmplx(1,1) = 1.256d0 
!      Bcmplx(1,1) = 0.856d0 
!      Bcmplx(1,1) = 0.666d0 
!      Beq(1,1) = 0.56d0

      c0 = ceq
      c0 = c_update(nv,q0,p0,qeq,peq,ceq,Bcmplx,Beq) 
    
      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Projected Gaussian Width Matrix:"

      do i = 1,nv+1
        write(*,*) Bcmplx(i,:)
      end do
    
      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Projected coefficients"

      do i = 1,nh
        write(*,*) c0(i) 
      end do
      
      call plot_wfn(nv,qeq,peq,ceq,Beq,1.d0)
      call plot_wfn(nv,q0,p0,c0,Bcmplx,1.d0)

 
!      stop

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
      write(*,*) trj

! FOR DEBUG PURPOSES ONLY, REMOVE IT AFTER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      q0(:) = qeq(:) 
!      p0(:) = peq(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

