!**********************************************************************!
! Module that contains the functions to compute the integrals of the   !
! powers of y1                                                         !
!**********************************************************************!

      module ypowers_module

      use tupowers_module
      use diagonal_module

      implicit none

      private
      public :: fun_y2,fun_y3,fun_y4,fun_y5,fun_y6
    
      contains

!.....y1^2..............................................................

      function fun_y2(nd,qvec,Tmat,MomMat) result(y2)
      ! nd: dimension of the system
      ! qvec: gaussian variational center
      ! Tmat: eigenvector matrix
      ! MomMat: momenta matrix
       
       integer,intent(in) :: nd
       real*8, dimension(nd), intent(in) :: qvec 
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 

       real*8 :: y2,Tu2,Tu1 

       y2 = 0.d0
       Tu2 = fun_Tu2(nd,Tmat,MomMat)
       Tu1 = fun_Tu1(nd,Tmat,MomMat)
       y2 = Tu2 + 2.d0*Tu1*qvec(1) + qvec(1)**2.d0

      end function

!.....y1^3..............................................................

      function fun_y3(nd,qvec,Tmat,MomMat) result(y3)
      ! nd: dimension of the system
      ! qvec: gaussian variational center
      ! Tmat: eigenvector matrix
      ! MomMat: momenta matrix
       
       integer,intent(in) :: nd
       real*8, dimension(nd), intent(in) :: qvec 
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 

       real*8 :: y3,Tu3,Tu2,Tu1 

       y3 = 0.d0
       Tu3 = fun_Tu3(nd,Tmat,MomMat)
       Tu2 = fun_Tu2(nd,Tmat,MomMat)
       Tu1 = fun_Tu1(nd,Tmat,MomMat)
       y3 = Tu3+3.d0*Tu2*qvec(1)+3.d0*Tu1*qvec(1)**2.d0+qvec(1)**3.d0

      end function

!.....y1^4..............................................................

      function fun_y4(nd,qvec,Tmat,MomMat) result(y4)
      ! nd: dimension of the system
      ! qvec: gaussian variational center
      ! Tmat: eigenvector matrix
      ! MomMat: momenta matrix
       
       integer,intent(in) :: nd
       real*8, dimension(nd), intent(in) :: qvec 
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 

       real*8 :: y4,Tu4,Tu3,Tu2,Tu1 

       y4 = 0.d0
       Tu4 = fun_Tu4(nd,Tmat,MomMat)
       Tu3 = fun_Tu3(nd,Tmat,MomMat)
       Tu2 = fun_Tu2(nd,Tmat,MomMat)
       Tu1 = fun_Tu1(nd,Tmat,MomMat)
       y4 = Tu4+4.d0*Tu3*qvec(1)+6.d0*Tu2*qvec(1)**2.d0+&
           &4.d0*Tu1*qvec(1)**3.d0+qvec(1)**4.d0

      end function

!.....y1^5..............................................................

      function fun_y5(nd,qvec,Tmat,MomMat) result(y5)
      ! nd: dimension of the system
      ! qvec: gaussian variational center
      ! Tmat: eigenvector matrix
      ! MomMat: momenta matrix
       
       integer,intent(in) :: nd
       real*8, dimension(nd), intent(in) :: qvec 
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 

       real*8 :: y5,Tu5,Tu4,Tu3,Tu2,Tu1 

       y5 = 0.d0
       Tu5 = fun_Tu5(nd,Tmat,MomMat)
       Tu4 = fun_Tu4(nd,Tmat,MomMat)
       Tu3 = fun_Tu3(nd,Tmat,MomMat)
       Tu2 = fun_Tu2(nd,Tmat,MomMat)
       Tu1 = fun_Tu1(nd,Tmat,MomMat)
       y5 = Tu5+5.d0*Tu4*qvec(1)+10.d0*Tu3*qvec(1)**2.d0+&
           &10.d0*Tu2*qvec(1)**3.d0+5.d0*Tu1*qvec(1)**4.d0+&
           &qvec(1)**5.d0

      end function
  
!.....y1^6..............................................................

      function fun_y6(nd,qvec,Tmat,MomMat) result(y6)
      ! nd: dimension of the system
      ! qvec: gaussian variational center
      ! Tmat: eigenvector matrix
      ! MomMat: momenta matrix
       
       integer,intent(in) :: nd
       real*8, dimension(nd), intent(in) :: qvec 
       real*8, dimension(nd,nd), intent(in) :: Tmat
       real*8, dimension(maxorder,nd), intent(in) :: MomMat 

       real*8 :: y6,Tu6,Tu5,Tu4,Tu3,Tu2,Tu1 

       y6 = 0.d0
       Tu6 = fun_Tu6(nd,Tmat,MomMat)
       Tu5 = fun_Tu5(nd,Tmat,MomMat)
       Tu4 = fun_Tu4(nd,Tmat,MomMat)
       Tu3 = fun_Tu3(nd,Tmat,MomMat)
       Tu2 = fun_Tu2(nd,Tmat,MomMat)
       Tu1 = fun_Tu1(nd,Tmat,MomMat)
       y6 = Tu6+6.d0*Tu5*qvec(1)+15.d0*Tu4*qvec(1)**2.d0+&
           &20.d0*Tu3*qvec(1)**3.d0+15.d0*Tu2*qvec(1)**4.d0+&
           &6.d0*Tu1*qvec(1)**5.d0+qvec(1)**6.d0

      end function
  
      end module
