!*********************************************************************!
! Module containing checking subroutines                              !
!*********************************************************************!

      module check_module

      use diagonal_module 
      use potential_module
      use quadratic_module
      use polynomials_module

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

!.....Check Vmat........................................................

      subroutine check_vmat()
      
       integer :: i

       call matrix_pot()

       write(*,*) "# HELLO I'M CHECK_VMAT"
       write(*,*) "Potential matrix"
       do i = 1,nv
         write(*,*) Vmat(i,:)
       end do

      end subroutine  

!.....Check quadratic forms.............................................

      subroutine check_quad(nd)
      ! nd : dimensions

       integer :: nd,i,j
       real*8 :: k
       real*8 :: int_qVq,int_uWu,int_uZQ,int_QRu,tP00M
       real*8, dimension(nd) :: qvec
       real*8, dimension(nd,nd) :: DiagMat,Trial
       real*8, dimension(maxorder,nd) :: OutMat
       
       call matrix_pot()

       qvec(:) = 2.d0

       int_qVq = fun_qVq(nd,qvec)

       write(*,*) "# HELLO I'M CHECK_QUAD"
       write(*,*) "Dimensions: ", nd
       write(*,*) "-------------------"
       write(*,*) "qvec:"
       write(*,*) qvec(:)
       write(*,*) "-------------------"
       write(*,*) "qVq:", int_qVq
     
       DiagMat(:,:) = 0.d0
       k=0.d0

       do i = 1,nd
          DiagMat(i,i) = 5.d0
          do j = 1,nd
            k=k+1.25d0
            Trial(i,j) = k
          end do
          !write(*,*) Trial(i,:)
       end do

       OutMat = momenta(nd,DiagMat)

       !write(*,*) "Matrix of the Momenta"

       !do i = 1,maxorder
       !   write(*,*) Outmat(i,:)
       !end do

       int_uWu = fun_uWu(nd,Trial,OutMat)
   
       write(*,*) "-------------------"
       write(*,*) "uWu:", int_uWu

       int_uZQ = fun_uZQ(nd,qvec,Trial,OutMat)
   
       write(*,*) "-------------------"
       write(*,*) "uZQ:", int_uZQ

       int_QRu = fun_QRu(nd,qvec,Trial,OutMat)
   
       write(*,*) "-------------------"
       write(*,*) "QRu:", int_QRu

       tP00M = tildeP00M(nd,qvec,Trial,OutMat)

       write(*,*) "-------------------"
       write(*,*) "tP00M - sum:", tP00M - int_uWu - int_QVQ

      end subroutine


      end module
