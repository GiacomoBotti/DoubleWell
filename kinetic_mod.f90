!**********************************************************************!
! Module containing everything pertaining kinetic energy               !
!**********************************************************************!

       module kinetic_module

       use constants
       use integrals_module
       use basisset_module

       implicit none

       private  
!       public  

       contains

!......dx Hi dx Hj integral.............................................

       function dxHdxH(nd,LambdaMat,a,avec,Amat,q) result(intdHdH)
       ! nd: dimensions
       ! LambdaMat: Diagonalized Bath Gaussian width
       ! a: x gaussian width
       ! avec: xy gaussian width vector
       ! Amat: y gaussian width matrix
       ! q: x variational parameter
       ! intdHdH: first integral of T
        integer, intent(in) :: nd
        real*8, intent(in) :: a,q
        real*8, dimension(nd), intent(in) :: avec
        real*8, dimension(nd,nd), intent(in) :: Amat
        real*8, dimension(nd,nd), intent(in) :: LambdaMat
        
        real*8, dimension(nh,nh) :: intdHdH

        integer :: i,j
        real*8 :: Y0
        real*8, dimension(nh,nh) :: X0mat

        Y0 = int_Y0(nd,LambdaMat)
        X0mat = int_XnMat(nd,0,a,avec,Amat,q)
        
        intdHdH(:,:) = 0.d0

        do i = 2,nh
           do j = 2,nh
              intdHdH(i,j) = 4.d0*i*j*Y0*X0mat(i-1,j-1)
           end do
        end do
       
       end function

!......dx ln G dx Hi integral...........................................

!......|dx ln G |^2 integral............................................

!......dy ln G dy ln G integral.........................................

!......Kinetic energy...................................................

       end module
