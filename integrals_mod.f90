!**********************************************************************!
! Module containing all the functions needed to compute the integrals  !
!**********************************************************************!

       module integrals_module

       use constants
       use basisset_module
       use matrix_module
       use inversion_module

       implicit none

       !real*8, parameter :: lwb=-5.d0
       !real*8, parameter :: hgb=5.d0
       !integer*8, parameter :: nstep=200
       ! HUGE GRID
       real*8, parameter :: lwb=-10.d0
       real*8, parameter :: hgb=10.d0
       integer*8, parameter :: nstep=1000

       private
       public :: int_Y0,int_XnMat,fun_Nsq,fun_NiNj,fun_Sb,int_TauMat
       public :: int_PMat

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
        Hmat=fun_Hmat(x,q,a)
        integral=Gx*Hmat*x**pow
        ! Higher bound
        x = hgb
        Gx=fun_Gx(nd,a,avec,Amat,x,q)
        Hmat=fun_Hmat(x,q,a)
        integral=integral+Gx*Hmat*x**pow
        ! First step
        x = lwb+h
        Gx=fun_Gx(nd,a,avec,Amat,x,q)
        Hmat=fun_Hmat(x,q,a)
        integral=integral+Gx*Hmat*x**pow

        s(:,:) = 0.d0
        do i = 2, nstep-2, 2 !only even
           x = lwb + i*h
           Gx=fun_Gx(nd,a,avec,Amat,x,q)
           Hmat=fun_Hmat(x,q,a)
           integrand = Gx*Hmat*x**pow 
           s = s + 2.d0*integrand ! even
           x = x + h
           Gx=fun_Gx(nd,a,avec,Amat,x,q)
           Hmat=fun_Hmat(x,q,a)
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
!        write(*,*) "Nsq", Nsq

!        Nsq = 1.d0

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
!        write(*,*) "NiNj", NiNj

!        NiNj = 1.d0

       end function

!......Sb function......................................................

       function fun_Sb(nd,x,qi,qj,pi,pj,Bimat,Bjmat) result(Sb)
       ! nd: bath dimensions
       ! x: active mode coordinate
       ! qi: total gaussian center vector (bra)
       ! qj: total gaussian center vector (ket)
       ! pi: total gaussian momenta vector (bra)
       ! pj: total gaussian momenta vector (ket)
       ! Bimat: complex gaussian width matrix (bra)
       ! Bjmat: complex gaussian width matrix (ket)
        integer, intent(in) :: nd
        real*8, intent(in) :: x
        real*8, dimension(nd+1), intent(in) :: qi,qj,pi,pj
        complex*16, dimension(nd+1,nd+1), intent(in) :: Bimat,Bjmat

        integer :: i,j
        integer*8 :: ndouble
        real*8 :: q,p
        real*8, dimension(nd) :: pbi,pbj,qbi,qbj
        complex*16 :: alphai,alphaj,Sb,gAg,qiAiqi,qjAjqj,prodi,prodj
        complex*16 :: h0i, h0j
        complex*16, dimension(nd) :: aveci,avecj,gveci,gvecj,gtot,qiAi
        complex*16, dimension(nd) :: qjAj,Ag,veci,vecj
        complex*16, dimension(nd,nd) :: Ai,Aj,Atot,invAtot

        pbi=pi(2:nd+1)
        pbj=pj(2:nd+1)
        qbi=qi(2:nd+1)
        qbj=qj(2:nd+1)

        call extracttildeA(nd,Bimat,Ai,aveci,alphai) 
        call extracttildeA(nd,Bjmat,Aj,avecj,alphaj) 

        Atot = Aj + transpose(dconjg(Ai))
        invAtot = invgen(nd,Atot)
         
        qiAi=matmul(transpose(dconjg(Ai)),qbi)
        qiAiqi = dot_product(qbi,qiAi)
        qjAj=matmul(Aj,qbj)
        qjAjqj = dot_product(qbj,qjAj)

        gveci=-(x-qi(1))*dconjg(aveci)-iu*pbi+qiAi
        gvecj=-(x-qj(1))*avecj+iu*pbj+qjAj
        gtot=gveci+gvecj

        Ag = matmul(invAtot,gtot)
        gAg = dot_product(dconjg(gtot),Ag)

        veci = (x-qi(1))*dconjg(aveci) + iu*pbi
        prodi = dot_product(qbi,veci)

        h0i=-0.5d0*dconjg(alphai)*(x-qi(1))**2-iu*pi(1)*(x-qi(1))&
            -0.5d0*qiAiqi+prodi

        vecj = (x-qj(1))*avecj - iu*pbj
        prodj = dot_product(qbj,vecj)

        h0j=-0.5d0*alphaj*(x-qj(1))**2+iu*pj(1)*(x-qj(1))&
            -0.5d0*qjAjqj+prodj

        Sb=exp(+0.5d0*gAg+h0i+h0j)

       end function

!......Time-shifted overlap.............................................

       function int_TauMat(nd,qi,qj,ppi,pj,Bimat,Bjmat) result(TauMat)
       ! nd: bath dimensions
       ! x: active mode coordinate
       ! qi: total gaussian center vector (bra)
       ! qj: total gaussian center vector (ket)
       ! ppi: total gaussian momenta vector (bra)
       ! pj: total gaussian momenta vector (ket)
       ! Bimat: complex gaussian width matrix (bra)
       ! Bjmat: complex gaussian width matrix (ket)
        integer, intent(in) :: nd
        real*8, dimension(nd+1), intent(in) :: qi,qj,ppi,pj
        complex*16, dimension(nd+1,nd+1), intent(in) :: Bimat,Bjmat

        integer :: i
        real*8 :: x,h,NiNj,ralphai,ralphaj
        complex*16 :: Sb,norm,alphai,alphaj,detAtot
        real*8, dimension(nh,nh) :: HiHjmat
        complex*16, dimension(nh,nh) :: integral,integrand,s,TauMat
        complex*16, dimension(nd) :: aveci,avecj 
        complex*16, dimension(nd,nd) :: Atot,Ai,Aj 

        NiNj = fun_NiNj(nd+1,real(Bimat),real(Bjmat))

        call extracttildeA(nd,Bimat,Ai,aveci,alphai) 
        call extracttildeA(nd,Bjmat,Aj,avecj,alphaj) 
        ralphai = real(alphai)
        ralphaj = real(alphaj)
        Atot = Aj + transpose(dconjg(Ai))
        detAtot = det_cmplx(nd,Atot)
!        write(*,*) "detAtot: ", detAtot
        norm = zsqrt(((2*pi)**nd)/detAtot)

        h = (hgb-lwb)/dfloat(nstep)
 
        ! Compute integral in boundaries
        integral(:,:) = 0.d0
        ! Lower bound
        x = lwb
        Sb=fun_Sb(nd,x,qi,qj,ppi,pj,Bimat,Bjmat)
        HiHjmat=fun_HmatShift(x,qi(1),qj(1),ralphai,ralphaj)
        integral=Sb*HiHjmat
        ! Higher bound
        x = hgb
        Sb=fun_Sb(nd,x,qi,qj,ppi,pj,Bimat,Bjmat)
        HiHjmat=fun_HmatShift(x,qi(1),qj(1),ralphai,ralphaj)
        integral=integral+Sb*HiHjmat
        ! First step
        x = lwb+h
        Sb=fun_Sb(nd,x,qi,qj,ppi,pj,Bimat,Bjmat)
        HiHjmat=fun_HmatShift(x,qi(1),qj(1),ralphai,ralphaj)
        integral=integral+Sb*HiHjmat

        s(:,:) = 0.d0
        do i = 2, nstep-2, 2 !only even
           x = lwb + i*h
           Sb=fun_Sb(nd,x,qi,qj,ppi,pj,Bimat,Bjmat)
           HiHjmat=fun_HmatShift(x,qi(1),qj(1),ralphai,ralphaj)
           integrand=Sb*HiHjmat
           s = s + 2.d0*integrand ! even
           x = x + h
           Sb=fun_Sb(nd,x,qi,qj,ppi,pj,Bimat,Bjmat)
           HiHjmat=fun_HmatShift(x,qi(1),qj(1),ralphai,ralphaj)
           integrand=Sb*HiHjmat
           s = s + 4.d0*integrand ! odd
        end do

        integral = (integral + s)*h/3.d0
 
        TauMat = integral*NiNj*norm

!        write(*,*) "power:", pow
!        write(*,*) "Xn11:", XnMat(1,1)
     
       end function

!......Batista reaction probability.....................................

       function int_PMat(nd,q,Bmat) result(PMat)
       ! nd: dimension of y
       ! q: x variational parameter
       ! Bmat: complex gaussian width
       ! Pmat: reaction probability 
        integer, intent(in) :: nd
        real*8, intent(in) :: q
        complex*16, dimension(nd+1,nd+1), intent(in) :: Bmat

        integer :: i,pow
        real*8 :: x,Gx,h,a,zero,Nsq,Y0
        real*8, dimension(nd) :: avec
        real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
        real*8, dimension(nh,nh) :: integral,integrand,s,Hmat,PMat

        Nsq=fun_Nsq(nd+1,real(Bmat))
        call extractA(nd,real(Bmat),Amat,avec,a) 
        call diagonalization(nd,Amat,LambdaMat,Tmat)
        Y0=int_Y0(nd,LambdaMat)

        !write(*,*) "MATRICES"
        !write(*,*) Bmat(1,1), Bmat(2,2)
        !write(*,*) a, Amat(1,1), avec(1)

        zero = 0.d0
        
        h = (hgb-zero)/dfloat(nstep)
 
        ! Compute integral in boundaries
        integral(:,:) = 0.d0
        ! Lower bound
        x = zero 
        Gx=fun_Gx(nd,a,avec,Amat,x,q)
        Hmat=fun_Hmat(x,q,a)
        integral=Gx*Hmat
        ! Higher bound
        x = hgb
        Gx=fun_Gx(nd,a,avec,Amat,x,q)
        Hmat=fun_Hmat(x,q,a)
        integral=integral+Gx*Hmat
        ! First step
        x = zero+h
        Gx=fun_Gx(nd,a,avec,Amat,x,q)
        Hmat=fun_Hmat(x,q,a)
        integral=integral+Gx*Hmat

        s(:,:) = 0.d0
        do i = 2, nstep-2, 2 !only even
           x = zero + i*h
           Gx=fun_Gx(nd,a,avec,Amat,x,q)
           Hmat=fun_Hmat(x,q,a)
           integrand = Gx*Hmat 
           s = s + 2.d0*integrand ! even
           x = x + h
           Gx=fun_Gx(nd,a,avec,Amat,x,q)
           Hmat=fun_Hmat(x,q,a)
           integrand = Gx*Hmat
           s = s + 4.d0*integrand ! odd
        end do

        integral = (integral + s)*h/3.d0
 
        Pmat = Nsq*Y0*integral

!        write(*,*) "power:", pow
!        write(*,*) "Xn11:", XnMat(1,1)
     
       end function
       end module
     
