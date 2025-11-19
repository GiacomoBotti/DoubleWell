!**********************************************************************!
! Module that contains the functions to compute the powers if Tu,      !
! needed for the powers of y1                                          !
!**********************************************************************!

      module tupowers_module

      use diagonal_module

      implicit none

      private
      public :: fun_Tu1,fun_Tu2,fun_Tu3,fun_Tu4,fun_Tu5,fun_Tu6

      contains

!.....Tu................................................................

      function fun_Tu1(nd,Tmat,MomMat) result(Tu1)
      ! nd: dimension of the system
      ! Tmat: eigenvector matrix
      ! MomMat: momenta matrix

       integer,intent(in) :: nd
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 
      
       integer :: i
       real*8 :: Tu1

       Tu1 = 0.d0

       do i = 1,nd
         Tu1 = Tu1 + Tmat(1,i)*MomMat(1,i)
       end do

      end function

!.....Tu^2..............................................................

      function fun_Tu2(nd,Tmat,MomMat) result(Tu2)
      ! nd: dimension of the system
      ! Tmat: eigenvector matrix
      ! MomMat: momenta matrix

       integer,intent(in) :: nd
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 
      
       integer :: i,j
       real*8 :: Tu2,Tu1

       Tu2=0.d0
       Tu1=0.d0

       do i = 1,nd
         do j = i+1,nd
           Tu1 = Tu1 + Tmat(1,j)*MomMat(1,j)
         end do
         Tu2=Tu2+(Tmat(1,i)**2.d0)*MomMat(2,i)&
             &+2.d0*Tmat(1,i)*MomMat(1,i)*Tu1
         !write(*,*) Tmat(1,i), Tmat(1,i)**2.d0
         Tu1=0.d0
       end do

      end function 

!.....Tu^3..............................................................

      function fun_Tu3(nd,Tmat,MomMat) result(Tu3)
      ! nd: dimension of the system
      ! Tmat: eigenvector matrix
      ! MomMat: momenta matrix

       integer,intent(in) :: nd
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 
      
       integer :: i,j,k
       real*8 :: Tu3,Tu2,Tu1

       Tu3=0.d0
       Tu2=0.d0
       Tu1=0.d0

       do i = 1,nd
         do j = i+1,nd
           do k = j+1,nd
             Tu1 = Tu1 + Tmat(1,i)*Tmat(1,j)*Tmat(1,k)&
                   &*MomMat(1,i)*MomMat(1,j)*MomMat(1,k)
           end do
           Tu2 = Tu2 + Tmat(1,j)*MomMat(1,j)*MomMat(2,i)*(Tmat(1,i)**2)&
                 &+ Tmat(1,i)*MomMat(1,i)*MomMat(2,j)*(Tmat(1,j)**2)
         end do
         Tu3 = Tu3 + MomMat(3,i)*(Tmat(1,i)**3) +3.d0*Tu2+6.d0*Tu3
       end do

      end function

!.....Tu^4..............................................................
 
      function fun_Tu4(nd,Tmat,MomMat) result(Tu4)
      ! nd: dimension of the system
      ! Tmat: eigenvector matrix
      ! MomMat: momenta matrix

       integer,intent(in) :: nd
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 
      
       integer :: i,j,k
       real*8 :: Tu4,Tu22

       Tu22=0.d0
       Tu4=0.d0

       do i = 1,nd
         do j = i+1,nd
           Tu22=Tu22+MomMat(2,i)*MomMat(2,j)*(Tmat(1,i)*Tmat(1,j))**2
         end do
         Tu4=Tu4+MomMat(4,i)*Tmat(1,i)**4 +6.d0*Tu22
         Tu22=0.d0
       end do

      end function

!.....Tu^5..............................................................
 
      function fun_Tu5(nd,Tmat,MomMat) result(Tu5)
      ! nd: dimension of the system
      ! Tmat: eigenvector matrix
      ! MomMat: momenta matrix

       integer,intent(in) :: nd
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 
      
       integer :: i,j,k
       real*8 :: Tu5

       Tu5=0.d0

      end function

!.....Tu^6..............................................................

      function fun_Tu6(nd,Tmat,MomMat) result(Tu6)
      ! nd: dimension of the system
      ! Tmat: eigenvector matrix
      ! MomMat: momenta matrix

       integer,intent(in) :: nd
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 
      
       integer :: i,j,k
       real*8 :: Tu6,Tu42,Tu222

       Tu6=0.d0
       Tu42=0.d0
       Tu222=0.d0

       do i = 1,nd
         do j = i+1,nd
           do k = j+1,nd
             Tu222=Tu222+MomMat(2,i)*MomMat(2,j)*MomMat(2,k)*&
                   &(Tmat(1,i)*Tmat(1,j)*Tmat(1,k))**2
           end do
          Tu42=Tu42+MomMat(4,i)*MomMat(2,j)*(Tmat(1,j)*Tmat(1,i)**2)**2&
               &+MomMat(4,j)*MomMat(2,i)*(Tmat(1,i)*Tmat(1,j)**2)**2
         end do
         Tu6=Tu6+MomMat(6,i)*Tmat(1,i)**6+15.d0*Tu42+90.d0*Tu222
         Tu42=0.d0
         Tu222=0.d0
       end do 

      end function

      end module
       
