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
       
       !write(111,*) "intdHdH: ", intdHdH(1,1)        

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

       !write(111,*) "intdlnG: ", intdlnG(1,1)        

       end function 

!......|dx ln G |^2 integral............................................

       function dxlnGsq(nd,q,p,qvec,tildeBmat) result(intdlnGsq)
       ! nd : dimensions of the bath
       ! q : active mode position
       ! p : active mode momentum
       ! qvec : bath position vector
       ! tildeBmat : complex total gaussian width
        integer, intent(in) :: nd
        real*8, intent(in) :: q,p
        real*8, dimension(nd), intent(in) :: qvec
        complex*16, dimension(nd+1,nd+1), intent(in) :: tildeBmat

        complex*16, dimension(nh,nh) :: intdlnGsq

        integer :: i,j
        real*8 :: Y0
        real*8, dimension(nh,nh) :: X0mat,X1mat,X2mat

        real*8 :: a,qk,tildeasq,kinvAa,Nsq
        real*8, dimension(nd) :: avec,invAa,kvec
        real*8, dimension(nd,nd) :: Amat,invA,LambdaMat,Tmat,qqMat
        real*8, dimension(nd,nd) :: invAainvAaMat,qinvAa 
        real*8, dimension(nd+1,nd+1) :: Bmat

        complex*16 :: tildea,tildeaq,tildeainvAa,tildeacq,talphainvA
        complex*16 :: Talphaqq,TalphaInvAaInvAa,TalphaqinvAa
        complex*16, dimension(nd) :: tildeavec
        complex*16, dimension(nd,nd) :: tildeAmat,alpha,alphainvA
        complex*16, dimension(nd,nd) :: alphaqq,alphaInvAaInvAa
        complex*16, dimension(nd,nd) :: alphaqinvAa

        Bmat = real(tildeBmat)
        Nsq=fun_Nsq(nd+1,Bmat)

        call extractA(nd,Bmat,Amat,avec,a)
        call extracttildeA(nd,tildeBmat,tildeAmat,tildeavec,tildea)
        call diagonalization(nd,Amat,LambdaMat,Tmat)

        invA = invgen_real(nd,Amat)

        Y0 = int_Y0(nd,LambdaMat)
        X0mat = int_XnMat(nd,0,a,avec,Amat,q)
        X1mat = int_XnMat(nd,1,a,avec,Amat,q)
        X2mat = int_XnMat(nd,2,a,avec,Amat,q)

        tildeaq=dot_product(dconjg(tildeavec),qvec)
        tildeacq=dot_product(tildeavec,qvec)
        invAa = matmul(invA,avec)
        tildeAinvAa= dot_product(dconjg(tildeavec),invAa)

        kvec=real(tildea)*real(tildeavec)+aimag(tildea)*aimag(tildeavec)

        qk=dot_product(qvec,kvec)
        kinvAa=dot_product(kvec,invAa)

        tildeasq=tildea*dconjg(tildea)

        do i =1,nd
          do j = 1,nd
            alpha(i,j) = dconjg(tildeavec(i))*tildeavec(j)
            qqMat(i,j) = qvec(i)*qvec(j)
            invAainvAaMat(i,j) = invAa(i)*invAa(j)
            qinvAa(i,j) = qvec(i)*invAa(j)
          end do
        end do

        alphainvA = matmul(transpose(alpha),invA)
        alphaqq = matmul(transpose(alpha),qqMat)
        alphaInvAaInvAa = matmul(transpose(alpha),invAainvAaMat)
        alphaqinvAa = matmul(transpose(alpha),qinvAa)

        TalphainvA = trace(nd,alphainvA)
        Talphaqq = trace(nd,alphaqq)
        TalphaInvAaInvAa = trace(nd,alphaInvAaInvAa)
        TalphaqinvAa = trace(nd,alphaqinvAa)

        intdlnGsq=tildeasq*X0mat&
                  &-(2*tildeasq*q+2*p*aimag(tildea)+2*qk)*X1mat&
                  &+(tildeasq*q*q+p*p+2*p*aimag(tildea)*q&
                  &+2*p*aimag(tildeaq)+tildeacq*tildeaq+2*q*qk)*X0mat&
                  &-(2*p*aimag(tildeaq)+tildeacq*tildeaq+2*q*qk)*X0mat&
                  &-(2*p*aimag(tildeAinvAa)+tildeacq*tildeAinvAa&
                  &+2*q*kinvAa)*(-X1mat+q*X0mat)-2*kinvAa*X2mat&
                  &+(2*qk+kinvAa*q)*X1mat + 0.5d0*TalphainvA*X0mat&
                  &+Talphaqq*X0mat+TalphaInvAaInvAa*&
                  &(X2mat-2*q*X1mat + X0mat*q**2)&
                  &-2*TalphaqinvAa*(X1mat-q*X0mat)

        intdlnGsq=Y0*intdlnGsq*Nsq

        write(111,*) "intdlnGsq: ", intdlnGsq(1,1)

       end function

!......dy ln G dy ln G integral.........................................

       function dylnGdylnG(nd,q,p,qvec,pvec,tildeBmat) result(intdyln) 
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

        complex*16, dimension(nh,nh) :: intdyln

        integer :: i,j
        real*8 :: Y0
        real*8, dimension(nh,nh) :: X0mat,X1mat,X2mat

        real*8 :: a,pMp,Nsq
        real*8, dimension(nd) :: avec,invAa,kvec,Mp
        real*8, dimension(nd,nd) :: Amat,invA,LambdaMat,Tmat,invMy
        real*8, dimension(nd,nd) :: qqMat,invAainvAaMat,qinvAa
        real*8, dimension(nd+1,nd+1) :: Bmat

        complex*16 :: tildea,aMa,qtAMtAq,iAatAMtAq,qtAMtAiAa,pMa,aMtAq
        complex*16 :: aMtAiAa,Tr1,Tr2,Tr3,Tr4
        complex*16, dimension(nd) :: Ma,tAMtAq,tAMtAiAa,aMtA 
        complex*16, dimension(nd,nd) :: MtA,tAMtA,tAMtAqq,tAMtAiAaiAa 
        complex*16, dimension(nd,nd) :: tAMtAqiAa,tAMtAiA

        complex*16, dimension(nd) :: tildeavec 
        complex*16, dimension(nd,nd) :: tildeAmat 

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

        !\mathbf{p}^{T}\mathbb{M}_{y}^{-1}\mathbf{p}
        Mp = matmul(invMy,pvec)
        pMp = dot_product(pvec,Mp)
        !\tilde{\mathbf{a}}^{\dagger}\mathbb{M}_{y}^{-1}\tilde{\mathbf{a}}
        Ma = matmul(invMy,tildeavec)
        aMa = dot_product(tildeavec,Ma)
        !\tilde{\mathbb{A}}^{\dagger}\mathbb{M}_{y}^{-1}\tilde{\mathbb{A}}
        MtA = matmul(invMy,tildeAmat)
        tAMtA = matmul(dconjg(tildeAmat),MtA)
        !\mathbf{q}^{T} (above) \mathbf{q}
        tAMtAq = matmul(tAMtA,qvec)
        qtAMtAq = dot_product(qvec,tAMtAq)        
        !(\mathbb{A}^{-1}\mathbf{a})^{T} (AMA) \mathbf{q}
        iAatAMtAq = dot_product(invAa,tAMtAq)
        !\mathbf{q}^{T} (above) \mathbb{A}^{-1}\mathbf{a}
        tAMtAiAa = matmul(tAMtA,invAa)
        qtAMtAiAa = dot_product(qvec,tAMtAiAa)
        !\mathbf{p}^{T}\mathbb{M}_{y}^{-1}\Im(\tilde{\mathbf{a}})
        pMa = dot_product(pvec,dimag(Ma))
        !\Re(\tilde{\mathbf{a}}^{\dagger}\mathbb{M}_{y}^{-1}\tilde{\mathbb{A}})
        aMtA = real(matmul(tildeavec,MtA))       
        ! (above)\mathbf{q}
        aMtAq = dot_product(aMtA,qvec)
        ! (above)\mathbb{A}^{-1}\mathbf{a}
        aMtAiAa = dot_product(aMtA,invAa)
 
        do i =1,nd
          do j = 1,nd
            qqMat(i,j) = qvec(i)*qvec(j)
            invAainvAaMat(i,j) = invAa(i)*invAa(j)
            qinvAa(i,j) = qvec(i)*invAa(j)
          end do
        end do

        ! Tr[AMA qq]
        tAMtAqq = matmul(transpose(tAMtA),qqMat)
        Tr1 = trace(nd,tAMtAqq)
        ! Tr[AMA AaAa]
        tAMtAiAaiAa = matmul(transpose(tAMtA),invAainvAaMat)
        Tr2 = trace(nd,tAMtAiAaiAa)
        ! Tr[AMA qAa]
        tAMtAqiAa = matmul(transpose(tAMtA),qinvAa)
        Tr3 = trace(nd,tAMtAqiAa)
        ! Tr[AMA A-1]
        tAMtAiA = matmul(transpose(tAMtA),invA)
        Tr4 = trace(nd,tAMtAiA)

       ! No y terms
        intdyln=(X2mat-2*q*X1mat+q*q*X0mat)*aMa+X0mat*pMp+X0mat*qtAMtAq&
               &-(2*pMa +2*aMtAq)*(X1mat-q*X0mat)

       ! y1 terms
        intdyln=intdyln-qtAMtAq*X0mat+iAatAMtAq*(X1mat-q*X0mat)&
               &-X0mat*qtAMtAq+(X1mat-q*X0mat)*(2*aMtAq+qtAMtAiAa)&
               &-2*(X2mat-2*X1mat+q*q*X0mat)*aMtAiAa
 
       ! y2 terms
        intdyln=intdyln+0.5d0*Tr4*X0mat+Tr1*X0mat-2*Tr3*(X1mat-q*X0mat)&
               &+Tr2*(X2mat-2*q*X1mat+q*q*X0mat)

        intdyln = intdyln*Y0*Nsq

        write(111,*) "intdyln: ", intdyln(1,1)

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
        complex*16, dimension(nh,nh) :: intdHdH,intdlnG
        complex*16, dimension(nh,nh) :: intdlnGsq,intdyln

        intdHdH = dxHdxH(nd,real(tildeBmat),q) 
        intdlnG = dxlnGdxHi(nd,q,p,qvec,tildeBmat)
        intdlnGsq = dxlnGsq(nd,q,p,qvec,tildeBmat)  
        intdyln = dylnGdylnG(nd,q,p,qvec,pvec,tildeBmat)  

        write(111,*) real(intdHdH(1,1)),real(intdlnG(1,1)),&
                   &real(intdlnGsq(1,1)),real(intdyln(1,1))

        mx = invMassMat(1,1)

        K00 = (mx*(intdlnGsq+intdHdH+transpose(dconjg(intdlnG))+&
              &intdlnG) + intdyln)/2.d0

       end function
       end module
