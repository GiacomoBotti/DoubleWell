!**********************************************************************!
! Module containing all the functions needed to compute the integrals  !
!**********************************************************************!

       module integrals_module

       use constants
       use parameters_module
       use basisset_module
       use matrix_module
       use inversion_module

       implicit none


       private
       public :: int_Y0,int_XnMat,fun_Nsq,fun_NiNj,fun_Sb,int_TauMat
       public :: int_PMat,gauss_hermite,int_XnMat_CS,fun_Sb_coeff
       public :: int_TauMat_CS,int_PMat_og

       contains

!......Y0 integral......................................................
       
       function int_Y0(nd,LambdaMat) result(Y0)
       ! Computes Y0 integral

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
       ! Computes tempered x Gaussian distribution

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

       function int_XnMat_CS(nd,pow,a,avec,Amat,q) result(XnMat)
       ! Computes the expectation value of x^n

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

        
        h = (hgb-lwb)/dfloat(gstep)
 
        ! Compute integral in boundaries
        integral(:,:) = 0.d0
        ! Lower bound
        x =  lwb
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
        do i = 2, gstep-2, 2 !only even
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

        where (abs(XnMat).le.1d-15)
           XnMat = complex(0.d0,0.d0)
        end where

!        write(*,*) "power:", pow
!        write(*,*) "Xn11:", XnMat(1,1)
     
       end function

!......N^2 factor.......................................................

       function fun_Nsq(ndim,Bmat) result(Nsq)
       ! Computes normalization factor
 
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
! SG NORMALIZATION
!        Nsq = dsqrt(Bmat(1,1)/pi) 

       end function

!......Gx function......................................................

       function fun_Gx2(nd,a,avec,Amat,x,q) result(Gx)
       ! Another way to compute the tempered x Gaussian distribution

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
       ! Computes the normalization factor for the time-shifted overlap

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
! SG normalization
!        NiNj = (Bimat(1,1)*Bjmat(1,1)/(pi*pi))**(1.d0/4.d0)

       end function

!......Sb function......................................................

       function fun_Sb(nd,x,qi,qj,ppi,pj,Bimat,Bjmat) result(Sb)
       ! Computes the exponential function in the time-shifted overlap

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
        real*8, dimension(nd+1), intent(in) :: qi,qj,ppi,pj
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

        pbi=ppi(2:nd+1)
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

        h0i=-0.5d0*dconjg(alphai)*(x-qi(1))**2-iu*ppi(1)*(x-qi(1))&
            -0.5d0*qiAiqi+prodi

        vecj = (x-qj(1))*avecj - iu*pbj
        prodj = dot_product(qbj,vecj)

        h0j=-0.5d0*alphaj*(x-qj(1))**2+iu*pj(1)*(x-qj(1))&
            -0.5d0*qjAjqj+prodj

        Sb=zexp(+0.5d0*gAg+h0i+h0j)
        !Sb=complex(dexp(-(313*x**2)/700)*dsqrt(15.d0)/(5.d0*dsqrt(pi)),&
        !      &0.d0) 

       end function

!......Time-shifted overlap (optimized CS)............................

      function int_TauMat_CS(nd,qi,qj,ppi,pj,Bimat,Bjmat) result(TauMat)
       ! Computes the time-shifted overlap using Cavalieri-Simpson
       ! with precomputed Gaussian coefficients

       ! nd: bath dimensions
       ! qi,qj: Gaussian center vectors
       ! ppi,pj: Gaussian momentum vectors
       ! Bimat,Bjmat: complex Gaussian width matrices
       ! TauMat: time-shifted overlap matrix

       integer, intent(in) :: nd
       real*8, dimension(nd+1), intent(in) :: qi,qj,ppi,pj
       complex*16, dimension(nd+1,nd+1), intent(in) :: Bimat,Bjmat

       integer :: i
       real*8 :: x,h,NiNj,ralphai,ralphaj
       complex*16 :: Sb,norm,alphai,alphaj,detAtot
       complex*16 :: Cc,Dc,Ec
       real*8, dimension(nh,nh) :: HiHjmat
       complex*16, dimension(nh,nh) :: integral,integrand,s,TauMat
       complex*16, dimension(nd) :: aveci,avecj
       complex*16, dimension(nd,nd) :: Atot,Ai,Aj

       ! Gaussian normalization
       NiNj = fun_NiNj(nd+1,real(Bimat,kind=8), &
                      real(Bjmat,kind=8))

       call extracttildeA(nd,Bimat,Ai,aveci,alphai)
       call extracttildeA(nd,Bjmat,Aj,avecj,alphaj)

       ralphai = real(alphai,kind=8)
       ralphaj = real(alphaj,kind=8)

       Atot = Aj + transpose(dconjg(Ai))
       detAtot = det_cmplx(nd,Atot)

       norm = zsqrt(((2*pi)**nd)/detAtot)

       ! Gaussian coefficients
       call fun_Sb_coeff(nd,qi,qj,ppi,pj,Bimat,Bjmat,Cc,Dc,Ec)

       ! Cavalieri-Simpson integration
       h = (hgb-lwb)/dfloat(gstep)

       integral(:,:) = (0.d0,0.d0)

       ! Lower bound
       x = lwb
       Sb = exp(-Cc*x**2+Dc*x+Ec)
       HiHjmat = fun_HmatShift(x,qi(1),qj(1),ralphai,ralphaj)
       integral = Sb*HiHjmat

       ! Higher bound
       x = hgb
       Sb = exp(-Cc*x**2+Dc*x+Ec)
       HiHjmat = fun_HmatShift(x,qi(1),qj(1),ralphai,ralphaj)
       integral = integral + Sb*HiHjmat

       ! First step
       x = lwb+h
       Sb = exp(-Cc*x**2+Dc*x+Ec)
       HiHjmat = fun_HmatShift(x,qi(1),qj(1),ralphai,ralphaj)
       integral = integral + 4.d0*Sb*HiHjmat

       s(:,:) = (0.d0,0.d0)

       do i = 2,gstep-2,2

          ! Even point
          x = lwb+i*h
          Sb = exp(-Cc*x**2+Dc*x+Ec)
          HiHjmat = fun_HmatShift(x,qi(1),qj(1),ralphai,ralphaj)
          integrand = Sb*HiHjmat
          s = s + 2.d0*integrand

          ! Odd point
          x = x+h
          Sb = exp(-Cc*x**2+Dc*x+Ec)
          HiHjmat = fun_HmatShift(x,qi(1),qj(1),ralphai,ralphaj)
          integrand = Sb*HiHjmat
          s = s + 4.d0*integrand

       end do

       integral = (integral+s)*h/3.d0

       TauMat = integral*NiNj*norm

       where (abs(TauMat).le.1.d-15)
          TauMat = (0.d0,0.d0)
       end where

       end function

!......Batista reaction probability.....................................

       function int_PMat_og(nd,q,Bmat) result(PMat)
       ! Computes reaction probability (x>0)

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
        
        h = (hgb-zero)/dfloat(gstep)
 
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
        integral=integral+4.d0*Gx*Hmat

        s(:,:) = 0.d0
        do i = 2, gstep-2, 2 !only even
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

!......Xnij integral (Gauss-Hermite)..................................

       function int_XnMat(nd,pow,a,avec,Amat,q) result(XnMat)
       ! Computes the expectation value of x^n
       ! using Gauss-Hermite quadrature

       ! nd: dimension of y
       ! pow: power of x
       ! a: x gaussian width
       ! avec: xy gaussian width vector
       ! Amat: y gaussian width matrix
       ! q: x variational parameter
       ! XnMat: matrix of integrals

       integer, intent(in) :: nd,pow
       real*8, intent(in) :: a,q
       real*8, dimension(nd), intent(in) :: avec
       real*8, dimension(nd,nd), intent(in) :: Amat

       integer :: i
       real*8 :: x,beta
       real*8, dimension(nd) :: vec1
       real*8, dimension(nd,nd) :: invA
       real*8, dimension(nh,nh) :: XnMat,Hmat

       ! Effective Gaussian width
       invA = invgen_real(nd,Amat)
       vec1 = matmul(invA,avec)
       beta = a-dot_product(avec,vec1)

       !if (beta.le.0.d0) then
       !   write(*,*) "Error: nonpositive Gaussian width",beta
       !   stop
       !end if

       ! Gauss-Hermite integration
       XnMat(:,:) = 0.d0

       do i = 1,ngh
          x = q + gh_nodes(i)/sqrt(beta)

          Hmat = fun_Hmat(x,q,a)

          XnMat = XnMat + gh_weights(i)*Hmat*x**pow
       end do

       XnMat = XnMat/sqrt(beta)

       where (abs(XnMat).le.1.d-15)
          XnMat = 0.d0
       end where

       end function

!......Sb quadratic coefficients......................................

       subroutine fun_Sb_coeff(nd,qi,qj,ppi,pj,Bimat,Bjmat,Cc,Dc,Ec)
       ! Computes coefficients of Sb(x)=exp(-C*x**2+D*x+E)

       ! nd: bath dimensions
       ! qi,qj: Gaussian center vectors
       ! ppi,pj: Gaussian momentum vectors
       ! Bimat,Bjmat: complex Gaussian width matrices
       ! Ccoef,Dcoef,Ecoef: quadratic exponent coefficients

       integer, intent(in) :: nd
       real*8, dimension(nd+1), intent(in) :: qi,qj,ppi,pj
       complex*16, dimension(nd+1,nd+1), intent(in) :: Bimat,Bjmat
       complex*16, intent(out) :: Cc,Dc,Ec

       complex*16 :: alphai,alphaj,ei,ej
       complex*16, dimension(nd) :: aveci,avecj
       complex*16, dimension(nd) :: uvec,vvec,qiAi,qjAj
       complex*16, dimension(nd) :: Mu,Mv
       complex*16, dimension(nd,nd) :: Ai,Aj,Atot,invAtot

       ! Extract Gaussian width blocks
       call extracttildeA(nd,Bimat,Ai,aveci,alphai)
       call extracttildeA(nd,Bjmat,Aj,avecj,alphaj)

       Atot = Aj + transpose(dconjg(Ai))
       invAtot = invgen(nd,Atot)

       qiAi = matmul(transpose(dconjg(Ai)),qi(2:nd+1))
       qjAj = matmul(Aj,qj(2:nd+1))

       ! g(x) = vvec - x*uvec
       uvec = dconjg(aveci) + avecj

       vvec = qi(1)*dconjg(aveci) + qj(1)*avecj &
            - iu*ppi(2:nd+1) + iu*pj(2:nd+1) &
            + qiAi + qjAj

       Mu = matmul(invAtot,uvec)
       Mv = matmul(invAtot,vvec)

       ! Constant contributions from bra and ket
       ei = -0.5d0*dconjg(alphai)*qi(1)**2 &
            + iu*ppi(1)*qi(1) &
            - 0.5d0*dot_product(qi(2:nd+1),qiAi) &
            - qi(1)*dot_product(qi(2:nd+1),dconjg(aveci)) &
            + iu*dot_product(qi(2:nd+1),ppi(2:nd+1))

       ej = -0.5d0*alphaj*qj(1)**2 &
            - iu*pj(1)*qj(1) &
            - 0.5d0*dot_product(qj(2:nd+1),qjAj) &
            - qj(1)*dot_product(qj(2:nd+1),avecj) &
            - iu*dot_product(qj(2:nd+1),pj(2:nd+1))

       ! Quadratic coefficient
       Cc = 0.5d0*(dconjg(alphai)+alphaj &
                     - dot_product(dconjg(uvec),Mu))

       ! Linear coefficient
       Dc = -dot_product(dconjg(uvec),Mv) &
               + dconjg(alphai)*qi(1) + alphaj*qj(1) &
               + iu*(pj(1)-ppi(1)) &
               + dot_product(qi(2:nd+1),dconjg(aveci)) &
               + dot_product(qj(2:nd+1),avecj)

       ! Constant coefficient
       Ec = 0.5d0*dot_product(dconjg(vvec),Mv) + ei + ej

       end subroutine

!......Time-shifted overlap (Gauss-Hermite)..........................

      function int_TauMat(nd,qi,qj,ppi,pj,Bimat,Bjmat) result(TauMat)
       ! Computes the time-shifted overlap using real Gauss-Hermite
       ! quadrature

       ! nd: bath dimensions
       ! qi,qj: Gaussian center vectors
       ! ppi,pj: Gaussian momentum vectors
       ! Bimat,Bjmat: complex Gaussian width matrices
       ! TauMat: time-shifted overlap matrix

       integer, intent(in) :: nd
       real*8, dimension(nd+1), intent(in) :: qi,qj,ppi,pj
       complex*16, dimension(nd+1,nd+1), intent(in) :: Bimat,Bjmat

       integer :: k
       real*8 :: x,c,d,mu,phase,fac,NiNj,ralphai,ralphaj
       complex*16 :: Cc,Dc,Ec,norm,alphai,alphaj,detAtot,Sb
       real*8, dimension(nh,nh) :: HiHjmat
       complex*16, dimension(nh,nh) :: TauMat
       complex*16, dimension(nd) :: aveci,avecj
       complex*16, dimension(nd,nd) :: Atot,Ai,Aj

       ! Gaussian normalization
       NiNj = fun_NiNj(nd+1,real(Bimat),real(Bjmat))

       call extracttildeA(nd,Bimat,Ai,aveci,alphai)
       call extracttildeA(nd,Bjmat,Aj,avecj,alphaj)

       ralphai = real(alphai,kind=8)
       ralphaj = real(alphaj,kind=8)

       Atot = Aj + transpose(dconjg(Ai))
       detAtot = det_cmplx(nd,Atot)

       norm = zsqrt(((2*pi)**nd)/detAtot)

       ! Quadratic Gaussian coefficients
       call fun_Sb_coeff(nd,qi,qj,ppi,pj,Bimat,Bjmat,Cc,Dc,Ec)

       ! Real Gaussian envelope
       c = real(Cc,kind=8)
       d = real(Dc,kind=8)

       if (c.le.0.d0) then
          write(*,*) "Error: nonpositive Gaussian width in TauMat"
          write(*,*) "c = ",c
          stop
       end if

       mu = d/(2.d0*c)

       ! Gaussian prefactor and Jacobian
       fac = exp(real(Ec,kind=8)+d**2/(4.d0*c))/sqrt(c)

       ! Gauss-Hermite integration
       TauMat(:,:) = (0.d0,0.d0)

       do k = 1,ngh

          x = mu + gh_nodes(k)/sqrt(c)

          ! Residual complex phase
          phase = -aimag(Cc)*x**2 + aimag(Dc)*x + aimag(Ec)

          Sb = exp(iu*phase)

          ! Hermite polynomial products
          HiHjmat = fun_HmatShift(x,qi(1),qj(1),ralphai,ralphaj)

          TauMat = TauMat + gh_weights(k)*Sb*HiHjmat

       end do

       TauMat = TauMat*fac*NiNj*norm

       where (abs(TauMat).le.1.d-15)
          TauMat = (0.d0,0.d0)
       end where

       end function

!......PMat integral (optimized CS).................................

       function int_PMat(nd,q,Bmat) result(PMat)
       ! Computes reaction probability (x>0)

       ! nd: dimension of y
       ! q: x variational parameter
       ! Bmat: complex gaussian width
       ! PMat: reaction probability

        integer, intent(in) :: nd
        real*8, intent(in) :: q
        complex*16, dimension(nd+1,nd+1), intent(in) :: Bmat

        integer :: i
        real*8 :: x,Gx,h,a,zero,Nsq,Y0,beta
        real*8, dimension(nd) :: avec,vec1
        real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat,invA
        real*8, dimension(nh,nh) :: integral,integrand,s,Hmat,PMat

        Nsq=fun_Nsq(nd+1,real(Bmat,kind=8))
        call extractA(nd,real(Bmat,kind=8),Amat,avec,a)
        call diagonalization(nd,Amat,LambdaMat,Tmat)
        Y0=int_Y0(nd,LambdaMat)

        ! Effective Gaussian width
        invA=invgen_real(nd,Amat)
        vec1=matmul(invA,avec)
        beta=a-dot_product(avec,vec1)

        if (beta.le.0.d0) then
           write(*,*) "Error: nonpositive beta in int_PMat"
           write(*,*) "beta = ",beta
           stop
        end if

        zero = 0.d0

        h = (hgb-zero)/dfloat(gstep)

        ! Compute integral in boundaries
        integral(:,:) = 0.d0

        ! Lower bound
        x = zero
        Gx=exp(-beta*(x-q)**2)
        Hmat=fun_Hmat(x,q,a)
        integral=Gx*Hmat

        ! Higher bound
        x = hgb
        Gx=exp(-beta*(x-q)**2)
        Hmat=fun_Hmat(x,q,a)
        integral=integral+Gx*Hmat

        ! First step
        x = zero+h
        Gx=exp(-beta*(x-q)**2)
        Hmat=fun_Hmat(x,q,a)
        integral=integral+4.d0*Gx*Hmat

        s(:,:) = 0.d0
        do i = 2, gstep-2, 2 !only even
           x = zero + i*h
           Gx=exp(-beta*(x-q)**2)
           Hmat=fun_Hmat(x,q,a)
           integrand = Gx*Hmat
           s = s + 2.d0*integrand ! even

           x = x + h
           Gx=exp(-beta*(x-q)**2)
           Hmat=fun_Hmat(x,q,a)
           integrand = Gx*Hmat
           s = s + 4.d0*integrand ! odd
        end do

        integral = (integral + s)*h/3.d0

        PMat = Nsq*Y0*integral

       end function

       end module
      
