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

      implicit none

      integer :: i,j
      real*8,dimension(nv+1) :: masses !Masses vector
      real*8,dimension(nv+1) :: q0 !inital centers vector
      real*8,dimension(nv+1) :: p0 !initial momenta vector
      complex*16,dimension(nh) :: c0 !Initial coefficient vector
      complex*16,dimension(nv+1,nv+1) :: Bcmplx !Initial width matrix
      integer*8,dimension(3) :: trj

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "|               MAIN CODE EXECUTION                 |"
      write(*,*) "+---------------------------------------------------+"

! TO BE SURE: GENERATE POTENTIAL MATRIX HERE
      call matrix_pot() 
! TO BE SURE: GENERATE HERMITE COEFFICIENT MATRIX HERE
      call GenHermMat()

!.....Define initial conditions.........................................

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Initial coefficients:"

      c0(:) = 0.d0
      c0(1) = 1.d0!/dsqrt(1.6435395517174083d0)
    
      do i = 1,nh
!        c0(i) = 1.d0/nh
        write(*,*) c0(i) 
      end do

      q0(:)=0.3d0
      p0(:)=0.3d0

      q0(1) = 1.5d0
      p0(1) = 1.0d0

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Initial Gaussian Width Matrix:"

      do i = 1,nv+1
         Bcmplx(i,i) = i+iu*i
         do j = i+1,nv+1
            Bcmplx(i,j) = (j+iu*j)/20.d0 !Gershgoring circle theorem
            Bcmplx(j,i) = Bcmplx(i,j)
         end do
         write(*,*) Bcmplx(i,:)
      end do

!      Bcmplx(1,1) = complex(1.d0,1.d0)

!.....Define masses vector..............................................

      write(*,*) "Masses vector:"

      do i = 1,nv+1
        masses(i) = 1.d0
        write(*,*) masses(i)
      end do

      call MassesMat(masses)

      write(*,*) "+---------------------------------------------------+"

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
!.....Evolution.........................................................

      write(*,*) "WE ARE RUNNING"
      write(*,*) "+---------------------------------------------------+"
      write(*,*) "Start       ", "Stop       ", "Lenght     "  
      trj = [0,5,1000]
      write(*,*) trj

      call bot_evo(nv,trj,q0,p0,c0,Bcmplx)

      write(*,*) "End of a successful run"
      write(*,*) "Have a nice day"
      write(*,*) "+---------------------------------------------------+"

      end program

