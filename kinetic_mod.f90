!**********************************************************************!
! Module containing everything pertaining kinetic energy               !
!**********************************************************************!

       module kinetic_module

       use constants
       use eofmotion_module
       use integrals_module
       use basisset_module
       use matrix_module
       use inversion_module

       implicit none

       private  
       public :: kin_energy 

       contains

!......K bath...........................................................

       function Kbath(nd,q,p,qvec,pvec,tildeBmat) result(intKb)
       ! nd : dimensions of the bath
       ! q : active mode position
       ! p : active mode momentum
       ! qvec : bath position vector
       ! pvec : bath momentum vector
       ! tildeBmat : complex total gaussian width
        integer, intent(in) :: nd
        real*8, intent(in) :: q,p
        real*8, dimension(nd), intent(in) :: qvec,pvec
        complex*16, dimension(nd+1,nd+1), intent(in) :: tildeBmat

        complex*16, dimension(nh,nh) :: intKb

        integer :: i,j
        real*8 :: Y0
        real*8, dimension(nh,nh) :: X0mat,X1mat,X2mat
        real*8, dimension(nh,nh) :: lin, sqr 

        real*8 :: a,pMp,Nsq
        real*8, dimension(nd) :: avec,invAa,kvec,Mp
        real*8, dimension(nd,nd) :: Amat,invA,LambdaMat,Tmat,invMy
        real*8, dimension(nd,nd) :: qqMat,invAainvAaMat,qinvAa
        real*8, dimension(nd+1,nd+1) :: Bmat

        complex*16 :: tildea,aMa,qtAMtAq,iAatAMtAq,qtAMtAiAa,pMa,aMtAq
        complex*16 :: aMtAiAa,Tr0,Tr1,Tr2,Tr3,Tr4,pMtAa
        complex*16 :: coeff0,coeff1,coeff2 
        complex*16, dimension(nd) :: Ma,tAMtAq,tAMtAiAa,aMtA,MtAa 
        complex*16, dimension(nd,nd) :: MtA,tAMtA,tAMtAqq,tAMtAiAaiAa 
        complex*16, dimension(nd,nd) :: tAMtAqiAa,tAMtAiA

        complex*16, dimension(nd) :: tildeavec 
        complex*16, dimension(nd,nd) :: tildeAmat 

        write(111,*) "KINETIC BATH"
        Bmat = real(tildeBmat)
        Nsq=fun_Nsq(nd+1,Bmat)

        call extractA(nd,Bmat,Amat,avec,a)
        call extracttildeA(nd,tildeBmat,tildeAmat,tildeavec,tildea)
        call diagonalization(nd,Amat,LambdaMat,Tmat)

        invA = invgen_real(nd,Amat)
        invAa = matmul(invA,avec)

        invMy = invMassMat(2:nd+1,2:nd+1)

        Y0 = int_Y0(nd,LambdaMat)
        X0mat = int_XnMat(nd,0,a,avec,Amat,q)
        X1mat = int_XnMat(nd,1,a,avec,Amat,q)
        X2mat = int_XnMat(nd,2,a,avec,Amat,q)
        lin = X1mat -q*X0mat
        sqr = X2mat -2*q*X1mat +q*q*X0mat

        !write(111,*) "X0mat(1,1): ", X0mat(1,1)
        !write(111,*) "X1mat(1,1): ", X1mat(1,1)
        !write(111,*) "X2mat(1,1): ", X2mat(1,1)

        !\mathbf{p}^{T}\mathbb{M}_{y}^{-1}\mathbf{p}
        Mp = matmul(invMy,pvec)
        pMp = dot_product(pvec,Mp)
        !\tilde{\mathbf{a}}^{\dagger}\mathbb{M}_{y}^{-1}\tilde{\mathbf{a}}
        Ma = matmul(invMy,tildeavec)
        aMa = dot_product(conjg(tildeavec),Ma)
        !\tilde{\mathbb{A}}^{\dagger}\mathbb{M}_{y}^{-1}\tilde{\mathbb{A}}
        MtA = matmul(invMy,tildeAmat)
        tAMtA = matmul(tildeAmat,MtA)
        !\mathbf{q}^{T} (above) \mathbf{q}
        tAMtAq = matmul(tAMtA,qvec)
        qtAMtAq = dot_product(qvec,tAMtAq)        
        !(\mathbb{A}^{-1}\mathbf{a})^{T} (AMA) \mathbf{q}
        iAatAMtAq = dot_product(invAa,tAMtAq)
        !\mathbf{q}^{T} (above) \mathbb{A}^{-1}\mathbf{a}
        tAMtAiAa = matmul(tAMtA,invAa)
        qtAMtAiAa = dot_product(qvec,tAMtAiAa)
        !\mathbf{p}^{T}\mathbb{M}_{y}^{-1}\Im(\tilde{\mathbf{a}})
        !pMa = dot_product(pvec,aimag(Ma))
        pMa = dot_product(pvec,Ma)
        !\Re(\tilde{\mathbf{a}}^{\dagger}\mathbb{M}_{y}^{-1}\tilde{\mathbb{A}})
        aMtA = matmul(tildeavec,MtA)       
        ! (above)\mathbf{q}
        !aMtAq = dot_product(aMtA,qvec)
        aMtAq = dot_product(qvec,aMtA)
        ! (above)\mathbb{A}^{-1}\mathbf{a}
        !aMtAiAa = dot_product(aMtA,invAa)
        aMtAiAa = dot_product(invAa,aMtA)
        ! \mathbf{p} \mathbb{M}^{-1}\tilde{\mathbb{A}}\mathbb{A}^{-1}\mathbf{a}
        MtAa = matmul(MtA,avec)
        pMtAa = dot_product(pvec,MtAa)

        do i =1,nd
          do j = 1,nd
            qqMat(i,j) = qvec(i)*qvec(j)
            invAainvAaMat(i,j) = invAa(i)*invAa(j)
            qinvAa(i,j) = qvec(i)*invAa(j)
          end do
        end do

        ! Tr[ MA ]
        Tr0 = trace(nd,MtA)
!        write(111,*) "Tr0: ", Tr0
        ! Tr[AMA qq]
        tAMtAqq = matmul(transpose(tAMtA),qqMat)
        Tr1 = trace(nd,tAMtAqq)
!        write(111,*) "Tr1: ", Tr1
        ! Tr[AMA AaAa]
        tAMtAiAaiAa = matmul(transpose(tAMtA),invAainvAaMat)
        Tr2 = trace(nd,tAMtAiAaiAa)
!        write(111,*) "Tr2: ", Tr2
        ! Tr[AMA qAa]
        tAMtAqiAa = matmul(transpose(tAMtA),qinvAa)
        Tr3 = trace(nd,tAMtAqiAa)
!        write(111,*) "Tr3: ", Tr3
        ! Tr[AMA A-1]
        tAMtAiA = matmul(transpose(tAMtA),invA)
        Tr4 = trace(nd,tAMtAiA)
!        write(111,*) "Tr4: ", Tr4

        coeff0=-Tr0+0.5*Tr4+Tr1-qtAMtAq-pMp
!        write(111,*) "coeff0: ", coeff0
        coeff1=-2*Tr3+2*qtAMtAiAa-2*pMtAa+2*pMa
!        write(111,*) "coeff1: ", coeff1
        coeff2=Tr2-2*aMtAiAa+aMa
!        write(111,*) "coeff2: ", coeff2
 
        intKb = coeff0*X0mat+coeff1*lin+coeff2*sqr

        intKb = intKb*Y0*Nsq

        write(111,*) "intKb: ", intKb(1,1)

       end function

!......K active.........................................................

       function Kact(nd,q,p,qvec,tildeBmat) result(intKa)
       ! nd : dimensions of the bath
       ! q : active mode position
       ! p : active mode momentum
       ! qvec : bath position vector
       ! tildeBmat : complex total gaussian width
        integer, intent(in) :: nd
        real*8, intent(in) :: q,p
        real*8, dimension(nd), intent(in) :: qvec
        complex*16, dimension(nd+1,nd+1), intent(in) :: tildeBmat

        complex*16, dimension(nh,nh) :: intKa,intKa1,intKa2,intKa3

        integer :: i,j
        real*8 :: Y0
        real*8, dimension(nh,nh) :: X0mat,X1mat,X2mat,lin,sqr

        real*8 :: a,Nsq
        real*8, dimension(nd) :: avec,invAa,kvec
        real*8, dimension(nd,nd) :: Amat,invA,LambdaMat,Tmat,qqMat
        real*8, dimension(nd,nd) :: invAainvAaMat,qinvAa 
        real*8, dimension(nd+1,nd+1) :: Bmat

        complex*16 :: Ka2c0,Ka2c1,Ka3c0,Ka3c1,Ka3c2
        complex*16 :: tildea,tildeaq,tildeainvAa,tildeacq
        complex*16 :: tildeasq
        complex*16 :: Tr0,Tr1,Tr2,Tr3
        complex*16, dimension(nd) :: tildeavec
        complex*16, dimension(nd,nd) :: tildeAmat,alpha
        complex*16, dimension(nd,nd) :: alphaqq,alphaInvAaInvAa
        complex*16, dimension(nd,nd) :: alphaqinvAa,alphainvA

        write(111,*) "KINETIC ACTIVE"
        Bmat = real(tildeBmat)
        Nsq=fun_Nsq(nd+1,Bmat)

        call extractA(nd,Bmat,Amat,avec,a)
        call extracttildeA(nd,tildeBmat,tildeAmat,tildeavec,tildea)
        call diagonalization(nd,Amat,LambdaMat,Tmat)

        invA = invgen_real(nd,Amat)

        Y0 = int_Y0(nd,LambdaMat)
        !write(111,*) "Y0: ", Y0
        X0mat = int_XnMat(nd,0,a,avec,Amat,q)
        X1mat = int_XnMat(nd,1,a,avec,Amat,q)
        X2mat = int_XnMat(nd,2,a,avec,Amat,q)

        lin = X1mat -q*X0mat
        sqr = X2mat -2*q*X1mat +q*q*X0mat

        invAa = matmul(invA,avec)
!        write(111,*) invAa
        tildeAinvAa= dot_product(invAa,tildeavec)
        !tildeAinvAa= dot_product(tildeavec,invAa)
!        write(111,*) tildeavec 

        ! |G|^2 Hi dxdxHj

        intKa1(:,:) = 0.d0

        do i = 1,nh
           do j = 3,nh
              intKa1(i,j) = 4.d0*(j-1)*(j-2)*Y0*X0mat(i,j-2)*Nsq
           end do
        end do

        ! |G|^2 Hi dxHj dxlnG

        intKa2(:,:) = (0.d0,0.d0) 

        Ka2c0=iu*p
!        write(111,*) Ka2c0
        Ka2c1=-tildea+tildeAinvAa
!        write(111,*) tildea, tildeAinvAa 

        do i = 1,nh
           do j = 2,nh
             intKa2(i,j)=2*(j-1)*(Ka2c0*X0mat(i,j-1)+Ka2c1*lin(i,j-1))
           end do
        end do

        intKa2=intKa2*Y0*Nsq

        ! |G|^2 HiHj dxdxlnG

        tildeasq=tildea*tildea
        tildeaq=dot_product(qvec,tildeavec)
        tildeacq=dot_product(tildeavec,qvec)

        do i =1,nd
          do j = 1,nd
            alpha(i,j) = tildeavec(i)*tildeavec(j)
            qqMat(i,j) = qvec(i)*qvec(j)
            invAainvAaMat(i,j) = invAa(i)*invAa(j)
            qinvAa(i,j) = qvec(i)*invAa(j)
          end do
        end do
        !Tr[alpha A]
        alphainvA = matmul(transpose(alpha),invA)
        Tr0 = trace(nd,alphainvA)
        !write(111,*) "Trace 0: ", Tr0
        !Tr[alpha qq]
        alphaqq = matmul(transpose(alpha),qqMat)
        Tr1 = trace(nd,alphaqq)
        !write(111,*) "Trace 1: ", Tr1
        !Tr[alpha InvAa InvAa]
        alphaInvAaInvAa = matmul(transpose(alpha),invAainvAaMat)
        Tr2 = trace(nd,alphaInvAaInvAa)
        !write(111,*) "Trace 2: ", Tr2
        !Tr[alpha q invAa]
        alphaqinvAa = matmul(transpose(alpha),qinvAa)
        Tr3 = trace(nd,alphaqinvAa)
        !write(111,*) "Trace 3: ", Tr3

        Ka3c0=-tildea-p*p-tildeaq*tildeaq+0.5*Tr0+Tr1
!        write(111,*) "Ka3c0: ", Ka3c0
        Ka3c1=-2*iu*p*tildea+2*iu*p*tildeainvAa-2*Tr3&
             &+2*tildeaq*tildeainvAa
!        write(111,*) "Ka3c1: ", Ka3c1
        Ka3c2=+tildeasq-2*tildea*tildeainvAa+Tr2
!        write(111,*) "Ka3c2: ", Ka3c2
 
        intKa3=(Ka3c0*X0mat+Ka3c1*lin+Ka3c2*sqr)*Y0*Nsq

        write(111,*) "intKa1(5,3): ", intKa1(5,3)
        write(111,*) "intKa2(5,3): ", intKa2(5,3)
        write(111,*) "intKa3: ", intKa3(1,1)

        intKa=intKa1 + 2*intKa2 + intKa3
       
       end function
 
!......Kinetic energy...................................................

       function kin_energy(nd,q,p,qvec,pvec,tildeBmat) result(K00)
       ! nd : dimensions of the bath
       ! q : active mode position
       ! p : active mode momentum
       ! qvec : bath position vector
       ! pvec : bath momentum vector
       ! tildeBmat : complex total gaussian width
        integer, intent(in) :: nd
        real*8, intent(in) :: q,p
        real*8, dimension(nd), intent(in) :: qvec,pvec
        complex*16, dimension(nd+1,nd+1), intent(in) :: tildeBmat

        complex*16, dimension(nh,nh) :: K00 !Complex for debugging

        real*8 :: mx
        complex*16, dimension(nh,nh) :: intKb,intKa
        complex*16, dimension(nh,nh) :: intdHdH,intdlnG
        complex*16, dimension(nh,nh) :: intdlnGsq,intdyln

        !intdHdH = dxHdxH(nd,real(tildeBmat),q) 
        !intdlnG = dxlnGdxHi(nd,q,p,qvec,tildeBmat)
        !intdlnGsq = dxlnGsq(nd,q,p,qvec,tildeBmat)  
        !intdyln = dylnGdylnG(nd,q,p,qvec,pvec,tildeBmat)  
        intKb = Kbath(nd,q,p,qvec,pvec,tildeBmat)  
        intKa = Kact(nd,q,p,qvec,tildeBmat)  

        mx = invMassMat(1,1)

        !K00 = -(mx*(intdlnGsq+intdHdH+transpose(dconjg(intdlnG))+&
        !      &intdlnG) + intdyln)/2.d0
        K00 = -0.5d0*(mx*intKa+intKb)

        write(111,*) "Ka: ", real(intKa(1,1)), aimag(intKa(1,1))
        write(111,*) "Kb: ", real(intKb(1,1)), aimag(intKb(1,1))
        write(111,*) "K00: ", real(K00(1,1)), aimag(K00(1,1))
        !write(111,*) "K00: ", K00(3,2)

       end function
       end module
