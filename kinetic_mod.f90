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

        real*8 :: a,qk,tildeasq,kinvAa
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

        intdlnGsq=Y0*intdlnGsq

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

        real*8 :: a,p
        real*8, dimension(nd) :: avec,invAa,kvec,Mp
        real*8, dimension(nd,nd) :: Amat,invA,LambdaMat,Tmat,invMy
        real*8, dimension(nd,nd) :: qqMat,invAainvAaMat,qinvAa
        real*8, dimension(nd+1,nd+1) :: Bmat

        complex*16 :: tildea,TtAMtAinvA,TtAMtAqq,TtAMtAiAaiAa
        complex*16 :: TtAMtAqiAa,aMtAq,qtAMtAiAa,qtAMtAq,aMtAiAa
        complex*16 :: qtAMa,iAaAMa,aMa,pMa 
        complex*16, dimension(nd) :: tildeavec,MtAq,tAMtAiAa,tAMtAq
        complex*16, dimension(nd) :: MtAiAa,Ma,tAMa 
        complex*16, dimension(nd,nd) :: tildeAmat,MtA,tAMtA,tAMtAinvA
        complex*16, dimension(nd,nd) :: tAMtAqq,tAMtAiAaiAa,tAMtAqiAa

        Bmat = real(tildeBmat)

        call extractA(nd,Bmat,Amat,avec,a)
        call extracttildeA(nd,tildeBmat,tildeAmat,tildeavec,tildea)
        call diagonalization(nd,Amat,LambdaMat,Tmat)

        invA = invgen_real(nd,Amat)

        invMy = invMassMat(2:nd,2:nd)

        Y0 = int_Y0(nd,LambdaMat)
        X0mat = int_XnMat(nd,0,a,avec,Amat,q)
        X1mat = int_XnMat(nd,1,a,avec,Amat,q)
        X2mat = int_XnMat(nd,2,a,avec,Amat,q)

        Mp = matmul(invMy,pvec)
        pMp = dot_product(pvec,Mp)
        MtA = matmul(invMy,tildeA)
        MtAq = matmul(MtA,qvec)
        aMtAq = dot_product(tildeavec,MtAq)
        tAMtA = matmul(tildeA,MtA)
        tAMtAq = matmul(tAMtA,qvec)
        qtAMtAq = dot_product(qvec,tAMtAq)        
        tAMtAinvA = matmul(transpose(tAMtA),invA)

        TtAMtAinvA = trace(nd,tAMtAinvA)

        invAa = matmul(invA,avec)
        MtAiAa = matmul(MtA,invAa)
        aMtAiAa = dot_product(tildeavec,MtAiAa)
        tAMtAiAa = matmul(tAMtA,invAa)
        qtAMtAiAa = dot_product(qvec,tAMtAiAa)
 
        do i =1,nd
          do j = 1,nd
            qqMat(i,j) = qvec(i)*qvec(j)
            invAainvAaMat(i,j) = invAa(i)*invAa(j)
            qinvAa(i,j) = qvec(i)*invAa(j)
          end do
        end do

        tAMtAqq = matmul(transpose(tAMtA),qqMat)
        tAMtAiAaiAa = matmul(transpose(tAMtA),invAainvAaMat)
        tAMtAqiAa = matmul(transpose(tAMtA),qinvAa)

        TtAMtAqq = trace(nd,tAMtAqq)
        TtAMtAiAaiAa = trace(nd,tAMtAiAaiAa)
        TtAMtAqiAa = trace(nd,tAMtAqiAa)

        Ma = matmul(invMy,tildeavec)
        aMa = dot_product(tildeavec,Ma)
        pMa = dot_product(pvec,Ma)
        tAMa = matmul(tildeAmat,Ma)
        qtAMa = dot_product(qvec,tAMa)
        iAaAMa = dot_product(invAa,tAMa)

        intdyln = (0.5d0*TtAMtAinvA +TtAMtAqq)*X0mat&
             &+TtAMtAiAaiAa*(X2mat-2*q*X1mat+q*q*X0mat)&
             &-2*TtAMtAqiAa*(X1mat-q*X0mat)&
             &+(aMtAq + 2*qtAMtAiAa)*(X1mat - q*X0mat)&
             &-2*qtAMtAq*X0mat-aMtAiAa*(X2mat-2*q*X1mat+q*q*X0mat)&
             &+qtAMa*(X1mat-q*X0mat)-iAaAMa*(X2mat-2*q*X1mat+q*q*X0mat)&
             &+aMa*X2mat+X1mat*(-qtAMa-aMtAq-2*q*aMa-dimag(pMa))&
             &-X0mat*(qtAMtAq+qtAMa*q+q*aMtAq+q*q*aMa+2*dimag(pMa)+pMp)

        intdyln = intdyln*Y0


!......Kinetic energy...................................................

       end module
