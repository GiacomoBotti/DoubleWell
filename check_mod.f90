!*********************************************************************!
! Module containing checking subroutines                              !
!*********************************************************************!

      module check_module

      use matrix_module 
      use potential_module
      use basisset_module
      use integrals_module
      use effectivepot_module
      use normalization_module

      implicit none

      contains

!.....Check Diagonalization.............................................

      subroutine check_diagonalization(nd)

       integer :: i,j,nd
       real*8 :: harvest
       real*8, dimension(nd,nd) :: LambdaMat,Tmat,RndMat,RecMat,Amat

       do i = 1,nd
          call RANDOM_NUMBER(harvest)
          RndMat(i,i) = harvest
          do j = i+1,nd
             call RANDOM_NUMBER(harvest)
             RndMat(i,j) = harvest
             RndMat(j,i) = RndMat(i,j)
          end do
       end do

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_DIAGONALIZATION"
       write(*,*) "Starting matrix"
     
       do i = 1,nd
          write(*,*) RndMat(i,:)
       end do

       Amat = RndMat 
 
       call diagonalization(nd,Amat,LambdaMat,Tmat)

       write(*,*) "Diagonal matrix"
     
       do i = 1,nd
          write(*,*) LambdaMat(i,:)
       end do
 
       RecMat = matmul(LambdaMat,transpose(Tmat))
       RecMat = matmul(Tmat,RecMat)

       write(*,*) "Recovered matrix"

       do i = 1,nd
          write(*,*) RecMat(i,:)
       end do

       write(*,*) "Error matrix"

       do i = 1,nd
          write(*,*) RecMat(i,:)-RndMat(i,:)
       end do

      end subroutine

!.....Check Momenta.....................................................

      subroutine check_momenta(nd)
      ! nd: dimensions of the matrix

       integer :: i,nd
       real*8, dimension(nd,nd) :: DiagMat
       real*8, dimension(maxorder,nd) :: OutMat

       DiagMat(:,:) = 0.d0

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_MOMENTA"
       write(*,*) "Diagonal Matrix"
       do i = 1,nd
          DiagMat(i,i) = 5.d0*i
          write(*,*) DiagMat(i,:)
       end do

       OutMat = momenta(nd,DiagMat)

       write(*,*) "Matrix of the Momenta"

       do i = 1,maxorder
          write(*,*) Outmat(i,:)
       end do

      end subroutine

!.....Check extractA....................................................

      subroutine check_Amat(nd)
      ! nd : dimension of the matrix

       integer :: i,j,nd
       real*8 :: harvest
       real*8, dimension(nd+1,nd+1) :: RndMat

       real*8 :: a
       real*8, dimension(nd) :: avec 
       real*8, dimension(nd,nd) :: Amat

       do i = 1,nd+1
          call RANDOM_NUMBER(harvest)
          RndMat(i,i) = harvest
          do j = i+1,nd+1
             call RANDOM_NUMBER(harvest)
             RndMat(i,j) = harvest
             RndMat(j,i) = RndMat(i,j)
          end do
       end do

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_AMAT"
       write(*,*) "Starting matrix"
     
       do i = 1,nd+1
          write(*,*) RndMat(i,:)
       end do

       call extractA(nd,RndMat,Amat,avec,a)

       write(*,*) "A matrix"

       do i = 1,nd
          write(*,*) Amat(i,:)
       end do

       write(*,*) "a vector"
       
       write(*,*) avec(:)

       write(*,*) "a: ", a

      end subroutine

!.....Check Vmat........................................................

      subroutine check_vmat()
      
       integer :: i

       write(*,*) "# HELLO I'M CHECK_VMAT"
       write(*,*) "Potential matrix"
       do i = 1,nv
         write(*,*) Vmat(i,:)
       end do

      end subroutine  

!.....Check HermMat.....................................................

      subroutine check_hermmat()
      
       integer :: i,j,k,l
       real*8 :: H
       real*8,dimension(nh,nh) :: Hmat

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_HERMMAT"
       write(*,*) "Hermite Coefficients matrix"
       do i = 1,nh
         write(*,*) Mherm(i,:)
       end do
       
       write(*,*) "------------------------------"
       write(*,*) "The values of some polynomials"
       write(*,*) "------------------------------"
       do i = 1,3
         do j = 1,3
           do l = 1,nh
             H = herm_pol(l,dfloat(i),dfloat(j))
             write(*,*) "npol: ",l,"x: ",i,"q: ",j,"H: ",H
           end do 
         end do
       end do

       write(*,*) "------------------------------"
       write(*,*) "The values of Hmat"
       write(*,*) "------------------------------"

       do i = 1,3
         do j =1,3
           write(*,*) "x: ",i,"q: ",j
           Hmat = fun_Hmat(dfloat(i),dfloat(j))
           do k = 1,nh
             write(*,*) Hmat(k,:)
           end do
         end do
       end do

      end subroutine  

!.....Check Y0..........................................................

      subroutine check_Y0(nd)
      ! nd: dimensions of the matrix

       integer :: i,nd
       real*8 :: Y0
       real*8, dimension(nd,nd) :: DiagMat

       DiagMat(:,:) = 0.d0

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_Y0"
       write(*,*) "------------------------------"
       write(*,*) "Diagonal Matrix"
       do i = 1,nd
          DiagMat(i,i) = 5.d0*i
          write(*,*) DiagMat(i,:)
       end do

       Y0 = int_Y0(nd,DiagMat)

       write(*,*) "------------------------------"
       write(*,*) "Y0: ", Y0

      end subroutine

!.....Check XnMat.......................................................

      subroutine check_XnMat(nd)
      !nd: dimensions of the bath matrix
       integer, intent(in) :: nd

       integer :: i,j
       real*8 :: a,q
       real*8, dimension(nd) :: avec 
       real*8, dimension(nd,nd) :: Amat
       real*8, dimension(nd+1,nd+1) :: Bmat
       real*8, dimension(nh,nh) :: XnMat

       Bmat(:,:) = 0.1d0
       Bmat(1,:) = 0.5d0
       Bmat(:,1) = 0.5d0
       Bmat(1,1) = 1.d0

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_XNMAT"
       write(*,*) "B matrix"
     
       do i = 1,nd+1
          write(*,*) Bmat(i,:)
       end do

       call extractA(nd,Bmat,Amat,avec,a)

       write(*,*) "------------------------------"
       write(*,*) "A matrix"

       do i = 1,nd
          write(*,*) Amat(i,:)
       end do

       write(*,*) "------------------------------"
       write(*,*) "a vector"
       
       write(*,*) avec(:)

       write(*,*) "------------------------------"
       write(*,*) "a: ", a
         
       do i = 1,4
         q = dfloat(i)/2.d0
         write(*,*) "------------------------------"
         write(*,*) "q: ", q
         XnMat = int_XnMat(nd,4,a,avec,Amat,q)
         do j = 1,nh
           write(*,*) XnMat(j,:)
         end do
       end do

      end subroutine

!.....Check V0..........................................................

      subroutine check_V0(nd,cvec)
      ! nd: dimensions of the bath 
      ! cvec: vector of the coefficients
       integer, intent(in) :: nd
       complex*16, dimension(nh),intent(in) :: cvec

       integer :: i,j
       real*8 :: Nsq
       real*8, dimension(nd+1) :: qtot, V1
       real*8, dimension(nd+1,nd+1) :: Bmat,V2
       real*8 :: harvest
       real*8 :: V0 

       do i = 1,nd+1
          Bmat(i,i) = i
          do j = i+1,nd+1
             Bmat(i,j) = j/20.d0 !Gershgoring circle theorem
             Bmat(j,i) = Bmat(i,j)
          end do
       end do

       qtot(1) = 1.d0
       qtot(2) = 2.d0
       qtot(3:nd+1) = 10000.d0 !This way I know if something is wrong


       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_V0"
       write(*,*) "B matrix"
     
       do i = 1,nd+1
       
          write(*,*) Bmat(i,:)
       end do

       write(*,*) "------------------------------"
       Nsq = fun_Nsq(nd+1,Bmat)
       write(*,*) "N squared: ", Nsq

       write(*,*) "------------------------------"
       V0 = fun_V0(nd,qtot,cvec,Bmat)
       write(*,*) "V0 ", V0
       
       write(*,*) "------------------------------"
       V1 = fun_V1(nd,qtot,cvec,Bmat)
       write(*,*) "V1: "
       write(*,*) V1

       write(*,*) "------------------------------"
       V2 = fun_V2(nd,qtot,cvec,Bmat)
       write(*,*) "V2: "
       do i = 1,nd+1
         write(*,*) V2(i,:)
       end do

      end subroutine

!.....Check Normalization...............................................

      subroutine check_norm(nd,cvec)
      ! nd: bath dimensions
      ! cvec: vector of the coefficients
       integer :: nd
       complex*16, dimension(nh),intent(in) :: cvec

       integer :: i,j
       real*8 :: Nout
       real*8, dimension(nd+1,nd+1) :: Bmat

       do i = 1,nd+1
          Bmat(i,i) = i
          do j = i+1,nd+1
             Bmat(i,j) = j/20.d0 !Gershgoring circle theorem
             Bmat(j,i) = Bmat(i,j)
          end do
       end do

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_NORM"
       write(*,*) "B matrix"
     
       do i = 1,nd+1
          write(*,*) Bmat(i,:)
       end do

       write(*,*) "------------------------------"
  
       do i = 1,10
         Nout = normalization(nd,dfloat(i),cvec,Bmat)

         write(*,*) "Norm at q = ", dfloat(i)
         write(*,*) Nout
       end do

      end subroutine


      end module
