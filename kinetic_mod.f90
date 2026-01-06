!**********************************************************************!
! Module containing everything pertaining kinetic energy               !
!**********************************************************************!

       module kinetic_module

       use constants
       use integrals_module
       use basisset_module
       use matrix_module
       use inversion_module

       implicit none

       private  
!       public  

       contains

!......dx Hi dx Hj integral.............................................

       function dxHdxH(nd,Bmat,q) result(intdHdH)
       ! nd: dimensions
       ! Bmat: total gaussian width (real)
       ! q: x variational parameter
       ! intdHdH: first integral of T
        integer, intent(in) :: nd
        real*8, dimension(nd+1,nd+1), intent(in) :: Bmat
        real*8, intent(in) :: q

        real*8 :: a
        real*8, dimension(nd) :: avec
        real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
        
        real*8, dimension(nh,nh) :: intdHdH

        integer :: i,j
        real*8 :: Y0
        real*8, dimension(nh,nh) :: X0mat

        call extractA(nd,Bmat,Amat,avec,a)
        call diagonalization(nd,Amat,LambdaMat,Tmat)

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

       function dxlnGdxHi(nd,q,p,qvec,tildeBmat) result(intdlnG)
       ! nd : dimensions of the bath
       ! q : active mode position
       ! p : active mode momentum
       ! qvec : bath position vector
       ! tildeBmat : complex total gaussian width
        integer, intent(in) :: nd
        real*8, intent(in) :: q,p
        real*8, dimension(nd), intent(in) :: qvec
        complex*16, dimension(nd+1,nd+1), intent(in) :: tildeBmat

        complex*16, dimension(nh,nh) :: part,intdlnG

        integer :: i,j
        real*8 :: Y0
        real*8, dimension(nh,nh) :: X0mat,X1mat

        real*8 :: a
        real*8, dimension(nd) :: avec,invAa
        real*8, dimension(nd,nd) :: Amat,invA,LambdaMat,Tmat
        real*8, dimension(nd+1,nd+1) :: Bmat

        complex*16 :: tildea,tildeaq,tildeainvAa
        complex*16, dimension(nd) :: tildeavec
        complex*16, dimension(nd,nd) :: tildeAmat

        Bmat = real(tildeBmat)

        call extractA(nd,Bmat,Amat,avec,a)
        call extracttildeA(nd,tildeBmat,tildeAmat,tildeavec,tildea)
        call diagonalization(nd,Amat,LambdaMat,Tmat)

        invA = invgen_real(nd,Amat)

        Y0 = int_Y0(nd,LambdaMat)
        X0mat = int_XnMat(nd,0,a,avec,Amat,q)
        X1mat = int_XnMat(nd,1,a,avec,Amat,q)

        tildeaq=dot_product(dconjg(tildeavec),qvec)
        invAa = matmul(invA,avec)
        tildeAinvAa= dot_product(dconjg(tildeavec),invAa)
        

        intdlnG(:,:) = (0.d0,0.d0) 
        part(:,:) = (0.d0,0.d0) 

        do i = 2,nh
           do j = 1,nh
              part(i,j) = -tildea*X1mat(i-1,j)&
                       &+tildea*q*X0mat(i-1,j)&
                       &-tildeaq*X0mat(i-1,j)&
                       &+tildeainvAa*X1mat(i-1,j)&
                       &-tildeainvAa*q*X0mat(i-1,j)&
                       &+tildeaq*X0mat(i-1,j)&
                       &+(0.d0,1.d0)*p*X0mat(i-1,j)
              intdlnG(i,j) = 2.d0*i*Y0*part(i,j)
           end do
        end do
        
       end function 

!......|dx ln G |^2 integral............................................

!......dy ln G dy ln G integral.........................................

!......Kinetic energy...................................................

       end module
