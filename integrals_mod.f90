!**********************************************************************!
! Module containing all the functions needed to compute the integrals  !
!**********************************************************************!

       module integrals_module

       use constants
       use basisset_module
       use matrix_module
       use inversion_module

       implicit none

       real*8, parameter :: lwb=-10.d0
       real*8, parameter :: hgb=10.d0
       integer*8, parameter :: nstep=500

       private
       public :: int_Y0,int_XnMat,fun_Nsq,fun_NiNj

       contains

!......Y0 integral......................................................
       
       function int_Y0(nd,LambdaMat) result(Y0)
       ! nd: dimensions
       ! LambdaMat: Diagonalized Bath Gaussian width
        integer, intent(in) :: nd
        real*8, dimension(nd,nd), intent(in) :: LambdaMat

        real*8 :: Y0

        integer :: i
        real*8 :: det

        det = LambdaMat(1,1)

        do i = 2,nd
           det = det*LambdaMat(i,i)
        end do

        Y0 = dsqrt((pi**nd)/det)
       
       end function

!......Gx function......................................................

       function fun_Gx(nd,a,avec,Amat,x,q) result(Gx)
       ! nd: dimension of y
       ! a: x gaussian width
       ! avec: xy gaussian width vector
       ! Amat: y gaussian width matrix
       ! x: x coordinate value
       ! q: x variational parameter
       ! Gx: Gx(x-q)
        integer, intent(in) :: nd
        real*8, intent(in) :: a,x,q
        real*8, dimension(nd), intent(in) :: avec
        real*8, dimension(nd,nd), intent(in) :: Amat

        real*8 :: Gx,aAa
        real*8, dimension(nd) :: vec1 
        real*8, dimension(nd,nd) :: invA

        invA=invgen_real(nd,Amat)

        vec1=matmul(invA,avec)

        aAa=dot_product(avec,vec1)

        Gx=dexp((-a+aAa)*(x-q)**2.d0)

       end function

!......Xnij integral....................................................

       function int_XnMat(nd,pow,a,avec,Amat,q) result(XnMat)
       ! nd: dimension of y
       ! pow: power of x
       ! a: x gaussian width
       ! avec: xy gaussian width vector
       ! Amat: y gaussian width matrix
       ! q: x variational parameter
       ! int_X0ij: matrix of integrals 
        integer, intent(in) :: nd,pow
        real*8, intent(in) :: a,q
        real*8, dimension(nd), intent(in) :: avec
        real*8, dimension(nd,nd), intent(in) :: Amat

        integer :: i
        real*8 :: x,Gx,h
        real*8, dimension(nh,nh) :: integral,integrand,s,Hmat,XnMat

        
        h = (hgb-lwb)/dfloat(nstep)
 
        ! Compute integral in boundaries
        integral(:,:) = 0.d0
        ! Lower bound
        x = lwb
        Gx=fun_Gx(nd,a,avec,Amat,x,q)
        Hmat=fun_Hmat(x,q)
        integral=Gx*Hmat*x**pow
        ! Higher bound
        x = hgb
        Gx=fun_Gx(nd,a,avec,Amat,x,q)
        Hmat=fun_Hmat(x,q)
        integral=integral+Gx*Hmat*x**pow
        ! First step
        x = lwb+h
        Gx=fun_Gx(nd,a,avec,Amat,x,q)
        Hmat=fun_Hmat(x,q)
        integral=integral+Gx*Hmat*x**pow

        s(:,:) = 0.d0
        do i = 2, nstep-2, 2 !only even
           x = lwb + i*h
           Gx=fun_Gx(nd,a,avec,Amat,x,q)
           Hmat=fun_Hmat(x,q)
           integrand = Gx*Hmat*x**pow 
           s = s + 2.d0*integrand ! even
           x = x + h
           Gx=fun_Gx(nd,a,avec,Amat,x,q)
           Hmat=fun_Hmat(x,q)
           integrand = Gx*Hmat*x**pow
           s = s + 4.d0*integrand ! odd
        end do

        integral = (integral + s)*h/3.d0
 
        XnMat = integral

!        write(*,*) "power:", pow
!        write(*,*) "Xn11:", XnMat(1,1)
     
       end function

!......N^2 factor.......................................................

       function fun_Nsq(ndim,Bmat) result(Nsq)
       ! ndim: dimension of the Bmat
       ! Bmat: fullD Gaussian width matrix (x & y)
       ! Nsq: square of the normalization factor
        integer, intent(in) :: ndim
        real*8, dimension(ndim,ndim), intent(in) :: Bmat

        real*8 :: Nsq,Bdet
        real*8, dimension(ndim,ndim) :: Support

        Support = Bmat
        Bdet = determinant(ndim,Support)
        !write(*,*) "Bdet", Bdet
        Nsq = dsqrt(Bdet/(pi**ndim))
        !write(*,*) "Nsq", Nsq

       end function

!......Gx function......................................................

       function fun_Gx2(nd,a,avec,Amat,x,q) result(Gx)
       ! nd: dimension of y
       ! a: x gaussian width
       ! avec: xy gaussian width vector
       ! Amat: y gaussian width matrix
       ! x: x coordinate value
       ! q: x variational parameter
       ! Gx: Gx(x-q)
        integer, intent(in) :: nd
        real*8, intent(in) :: a,x,q
        real*8, dimension(nd), intent(in) :: avec
        real*8, dimension(nd,nd), intent(in) :: Amat

        integer :: i
        real*8 :: Gx,aLa
        real*8, dimension(nd) :: vec1 
        real*8, dimension(nd,nd) :: LambdaMat,Tmat,invLambda

        call diagonalization(nd,Amat,LambdaMat,Tmat)
        
        invLambda(:,:) = 0.d0

        do i = 1,nd
          invLambda(i,i) = 1.d0/LambdaMat(i,i)
        end do

        vec1=matmul(invLambda,avec)

        aLa=dot_product(avec,vec1)

        Gx=dexp(-(a-aLa)*(x-q)**2)

       end function

!......NiNj factor......................................................

       function fun_NiNj(ndim,Bimat,Bjmat) result(NiNj)
       ! ndim: dimension of the Bmat
       ! Bimat: fullD Gaussian width matrix (x & y) BRA
       ! Bjmat: fullD Gaussian width matrix (x & y) KET
       ! NiNj: square of the normalization factor
        integer, intent(in) :: ndim
        real*8, dimension(ndim,ndim), intent(in) :: Bimat
        real*8, dimension(ndim,ndim), intent(in) :: Bjmat

        real*8 :: NiNj,Bidet,Bjdet
        real*8, dimension(ndim,ndim) :: Supporti,Supportj

        Supporti = Bimat
        Bidet = determinant(ndim,Supporti)
        Supportj = Bjmat
        Bjdet = determinant(ndim,Supportj)
        !write(*,*) "Bdet", Bdet
        NiNj =((Bjdet/pi**ndim)*(Bidet/pi**ndim))**(1.d0/4.d0) 
        !write(*,*) "NiNj", NiNj
       end function
       end module
     
