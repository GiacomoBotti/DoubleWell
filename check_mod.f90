!*********************************************************************!
! Module containing checking subroutines                              !
!*********************************************************************!

      module check_module

      use matrix_module 
      use potential_module
      use basisset_module

      implicit none

      contains

!.....Check Diagonalization.............................................

      subroutine check_diagonalization(nd)

       integer :: i,j,nd
       real*8 :: harvest
       real*8, dimension(nd,nd) :: LambdaMat,Tmat,RndMat,RecMat

       complex*16, dimension(nd,nd) :: Amat

       do i = 1,nd
          call RANDOM_NUMBER(harvest)
          RndMat(i,i) = harvest
          do j = i+1,nd
             call RANDOM_NUMBER(harvest)
             RndMat(i,j) = harvest
             RndMat(j,i) = RndMat(i,j)
          end do
       end do

       write(*,*) "# HELLO I'M CHECK_DIAGONALIZATION"
       write(*,*) "Starting matrix"
     
       do i = 1,nd
          write(*,*) RndMat(i,:)
       end do

       Amat = cmplx(RndMat) 
 
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

       write(*,*) "# HELLO I'M CHECK_AMAT"
       write(*,*) "Starting matrix"
     
       do i = 1,nd+1
          write(*,*) RndMat(i,:)
       end do

       call extractA(nd,RndMat,Amat,avec)

       write(*,*) "A matrix"

       do i = 1,nd
          write(*,*) Amat(i,:)
       end do

       write(*,*) "a vector"
       
       write(*,*) avec(:)

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

      end module
