!**********************************************************************!
! Module containing the functions required by the BOT_evo subroutine   !
! in evolution_mod.f90                                                 !
!**********************************************************************!

       module BOT_module

       use parameters_module
       use basisset_module
       use kinetic_module
       use effectivepot_module
       use integrals_module
       use inversion_module
       use matrix_module
       use observable_module
       use constants

       implicit none

       integer*8 :: info

       private
       public ::c_static,c_update,c_update_fb,c_update_fbs,c_update_full

       contains

!......Coefficient update with static basis (analytical)................

       function c_static(nd,h,cvec,S00M,H00M) result(csout)
       ! nd: bath dimension 
       ! h : time-step size
       ! S00M : overlap matrix
       ! H00M : hamiltonian matrix
       ! q0: full position vector
       ! p0: full momenta vector
       ! tildeBmat: full complex gaussian width matrix
        implicit none
        integer, intent(in) :: nd
        real*8, intent(in) :: h
        complex*16, dimension(nh), intent(in) :: cvec
        real*8, dimension(nh,nh) :: S00M
        complex*16, dimension(nh,nh) :: H00M

        integer*8 :: i,lwork,nh8
        real*8 :: q,p,a,Nsq,Y0
        real*8, dimension(nd) :: qvec,pvec,avec
        real*8, dimension(nh) :: eigenv
        real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
        real*8, dimension(nd+1,nd+1) :: Bmat
        real*8, dimension(nh,nh) :: X4,X3,X2,X1,X0
        complex*16, dimension(nh) :: c,csout,expvec
        complex*16, dimension(nh,nh) :: T00M,V00M
        complex*16, dimension(nh,nh) :: Z,adjZ,B

        complex*16, dimension(2*nh-1) :: lapwork
        real*8, dimension(3*nh-1) :: rwork

        ! External procedures defined in LAPACK
        external ZHEGV
        lwork = 2*nh-1

!        write(*,*) "I AM  STATIC"

!        write(*,*) qtot
!        write(*,*) ptot
!        write(*,*) cvec
!        do i = 1,nd+1
!          write(*,*) tildeBmat(i,:)
!        end do

        ! Copy H00M so LAPACK can overwrite
        Z = H00M
        B = S00M
!        B(:,:) = 0.d0
!        do i = 1,nh
!          B(i,i) = 1.d0
!        end do

!        write(*,*) "T00M:"
!        do i = 1,nh
!          write(*,*) T00M(i,:)
!        end do

!        write(*,*) "V00M:"
!        do i = 1,nh
!          write(*,*) V00M(i,:)
!        end do

!        write(*,*) "H00M:"
!        do i = 1,nh
!          write(*,*) H00M(i,:)
!        end do
!        write(*,*) "????????????"
!        write(*,*) "H00M:"
!        do i = 1,nh
!          write(*,*) Z(i,:)
!        end do


        call ZHEGV(1,'V','U',nh,Z,nh,B,nh,eigenv,&
                  &lapwork,lwork,rwork,info)

        adjZ = dconjg(transpose(Z))
         
!        call test_static(nh,S00M,H00M,Z,eigenv)

        c = cvec
        c = matmul(S00M,c)
        c = matmul(adjZ,c) 
        do i = 1,nh
          expvec(i) = zexp(-iu*eigenv(i)*h)*c(i)
        end do

        csout = matmul(Z,expvec)

       end function


!......TEST STATIC EVOLUTION............................................

       subroutine test_static(nd,S,H,Z,eigenv)
       ! nd : matrices dimensions
       ! S : overlap matrix
       ! H : hamiltonian matrix
       ! Z : eigenvectors of HZ = LSZ
        implicit none
        integer, intent(in) :: nd

        integer*8 :: i
        real*8, dimension(nd) :: eigenv
        complex*16, dimension(nd,nd) :: S,H
        complex*16, dimension(nd,nd) :: Z,adjZ,TEST1,TEST2,TEST6
        complex*16, dimension(nd,nd) :: TEST3,TEST4,TEST5,invS
        character(len=100) :: formato
      
        formato="(E10.1,E10.1,E10.1,E10.1,E10.1,E10.1)"
  
        write(*,*) "S:"
        do i = 1,nd
           write(*,*) S(i,:)
        end do
  
        write(*,*) "H:"
        do i = 1,nd
           write(*,*) H(i,:)
        end do
  
        write(*,*) "Z:"
        do i = 1,nd
           write(*,*) Z(i,:)
        end do

        adjZ = dconjg(transpose(Z))

        TEST4 = matmul(adjZ,matmul(H,Z))

        TEST5(:,:) = complex(0.d0,0.d0)
        write(*,*) "Eigenvalues:"
        do i = 1,nd
           TEST5(i,i) = eigenv(i)
           write(*,*) eigenv(i)
        end do
        
        TEST6 = matmul(S,matmul(Z,TEST5))

        ! TEST TEST TEST TEST TEST
        TEST1 = matmul(S,Z)
        TEST2 = matmul(adjZ,TEST1)
        write(*,*) "Z^H*S*Z = 1"
        write(*,formato) TEST2(1,1)-1.d0, TEST2(1,2), TEST2(1,3)
        write(*,formato) TEST2(2,1), TEST2(2,2)-1.d0, TEST2(2,3)
        write(*,formato) TEST2(3,1), TEST2(3,2), TEST2(3,3)-1.d0
        write(*,*) " "
        invS = invgen(nd,S)
        TEST3 = matmul(Z,adjZ)
        write(*,*) "ZZ^H=S-1"
        write(*,formato) TEST3(1,1)-invS(1,1), TEST3(1,2)-invS(1,2)
        write(*,formato) TEST3(2,1)-invS(2,1), TEST3(2,2)-invS(2,2)
        write(*,formato) TEST3(3,1)-invS(3,1), TEST3(3,2)-invS(3,2)
        write(*,*) " "
        write(*,*) "Z^H*H*Z=L"
        write(*,formato) TEST4(1,1)-TEST5(1,1), TEST4(1,2)-TEST5(1,2)
        write(*,formato) TEST4(2,1)-TEST5(2,1), TEST4(2,2)-TEST5(2,2)
        write(*,formato) TEST4(3,1)-TEST5(3,1), TEST4(3,2)-TEST5(3,2)
        write(*,*) " "

       end subroutine

!......Analytical update of electronic coefficients.....................

       function c_update(nd,qb,pb,qk,pk,c,Bb,Bk) result(cout)
       ! nd : bath dimension
       ! qb : full position vector (Bra)
       ! pb : full momentum vector (Bra)
       ! qk : full position vector (Ket)
       ! pk : full momentum vector (Ket)
       ! c : exponentially-evolved coefficients
       ! Bb : full width matrix (Bra)
       ! Bk : full width matrix (Ket)
        integer, intent(in) :: nd
        real*8, dimension(nd+1), intent(in) :: qb,pb,qk,pk
        complex*16, dimension(nh), intent(in) :: c
        complex*16, dimension(nd+1,nd+1), intent(in) :: Bb,Bk

        integer :: i,j,info
        real*8 :: a,q,Nsq,Y0,reS,imS
        real*8, dimension(nd) :: avec 
        real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
        real*8, dimension(nd+1,nd+1) :: Bmat
        complex*16, dimension(nh) :: csupp,cout
        complex*16, dimension(nh,nh) :: Aux,S00M,Tt0M,invS,X0Mat

        external ZPOSV

        q = qb(1)

!        write(*,*) "I AM C UPDATE"

!        write(*,*) qb
!        write(*,*) pb
!        write(*,*) qk
!        write(*,*) pk

        Bmat = real(Bb)

        Nsq=fun_Nsq(nd+1,Bmat)
        call extractA(nd,Bmat,Amat,avec,a)
        call diagonalization(nd,Amat,LambdaMat,Tmat)
        Y0=int_Y0(nd,LambdaMat)
        X0Mat=int_XnMat(nd,0,a,avec,Amat,q)

        Aux = X0Mat*Y0*Nsq 
        S00M = Aux

! SG qtag.f line 1046 : cleaning new overlap (?)

        do i = 1, nh
          S00M(i,i) = Aux(i,i)*complex(1.d0,0.d0)
          do j = i+1,nh
            reS = dreal(Aux(i,j)+Aux(j,i))/2.d0
            imS = dimag(Aux(i,j)-Aux(j,i))/2.d0
            S00M(i,j) = reS+iu*imS
            S00M(j,i) = conjg(S00M(i,j))
          end do
        end do
            
        Tt0M = int_TauMat(nd,qb,qk,pb,pk,Bb,Bk) 
        csupp = matmul(Tt0M,c)
 
!        write(*,*) "TtOMc"
!        write(*,*) csupp

!        write(*,*) "Tt0M:"
!        do i = 1,nh
!          write(*,*) Tt0M(i,:)
!        end do

!        write(*,*) "S00M:"
!        do i = 1,nh
!          write(*,*) S00M(i,:)
!        end do

!        cout = linsys(nh,S00M,csupp) 
         
         call zposv('U',nh,1,S00M,nh,csupp,nh,info)

         cout = csupp

       end function

!......Analytical update of electronic coefficients (FB).................

       function c_update_fb(nd,qb,pb,qk,pk,c,Bb,Bk) result(cout)
       ! nd : bath dimension
       ! qb : full position vector (Bra)
       ! pb : full momentum vector (Bra)
       ! qk : full position vector (Ket)
       ! pk : full momentum vector (Ket)
       ! c : exponentially-evolved coefficients
       ! Bb : full width matrix (Bra)
       ! Bk : full width matrix (Ket)
        integer, intent(in) :: nd
        real*8, dimension(nd+1), intent(in) :: qb,pb,qk,pk
        complex*16, dimension(nh), intent(in) :: c
        complex*16, dimension(nd+1,nd+1), intent(in) :: Bb,Bk

        integer :: i
        real*8 :: a,q,Nsq,Y0
        real*8, dimension(nd) :: avec 
        real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
        real*8, dimension(nd+1,nd+1) :: Bmat
        complex*16, dimension(nh) :: csupp,cout,csupptest,csuppdiff
        complex*16, dimension(nh,nh) :: SttM,St0M,invS,X0Mat,S0tM
        complex*16, dimension(nh,nh) :: prod1,prod2,summa,S00M

!        write(*,*) "I AM C UPDATE FB"
        q = qb(1)
        Bmat = real(Bb)

        Nsq=fun_Nsq(nd+1,Bmat)
        call extractA(nd,Bmat,Amat,avec,a)
        call diagonalization(nd,Amat,LambdaMat,Tmat)
        Y0=int_Y0(nd,LambdaMat)
        X0Mat=int_XnMat(nd,0,a,avec,Amat,q)

        SttM = X0Mat*Y0*Nsq 

!        write(*,*) "SttM"
 
!        do i = 1,nh
!          write(*,*) SttM(i,:)
!        end do


        invS = invgen(nh,SttM) 
!        write(*,*) "SttM^{-1}"
 
!        do i = 1,nh
!          write(*,*) invS(i,:)
!        end do

        q = qk(1)
        Bmat = real(Bk)

!        write(*,*) q, Bmat(1,1)

        Nsq=fun_Nsq(nd+1,Bmat)
        call extractA(nd,Bmat,Amat,avec,a)
        call diagonalization(nd,Amat,LambdaMat,Tmat)
        Y0=int_Y0(nd,LambdaMat)
        X0Mat=int_XnMat(nd,0,a,avec,Amat,q)

        S00M = X0Mat*Y0*Nsq 
!       write(*,*) "S00M"
 
!        do i = 1,nh
!          write(*,*) S00M(i,:)
!        end do

        St0M = int_TauMat(nd,qb,qk,pb,pk,Bb,Bk)
!        write(*,*) "St0M"
 
!        do i = 1,nh
!          write(*,*) St0M(i,:)
!        end do
        S0tM = conjg(transpose(St0M))
!        write(*,*) "S0tM"
 
!        do i = 1,nh
!          write(*,*) S0tM(i,:)
!        end do

        prod1=matmul(invS,St0M)
        prod2=matmul(S0tM,prod1)

!        write(*,*) "prod2"
 
!        do i = 1,nh
!          write(*,*) prod2(i,:)
!        end do

        summa = 0.5*(prod2+S00M)
 
        !write(*,*) "summa"
 
!       do i = 1,nh
!          write(*,*) summa(i,:)
!        end do

        csupp = matmul(summa,c)
 
        cout = linsys(nh,S0tM,csupp) 

       end function

!......Analytical update of electronic coefficients (FB - Schr).........

       function c_update_fbs(nd,qb,pb,qk,pk,c,Bb,Bk) result(cout)
       ! nd : bath dimension
       ! qb : full position vector (Bra)
       ! pb : full momentum vector (Bra)
       ! qk : full position vector (Ket)
       ! pk : full momentum vector (Ket)
       ! c : exponentially-evolved coefficients
       ! Bb : full width matrix (Bra)
       ! Bk : full width matrix (Ket)
        integer, intent(in) :: nd
        real*8, dimension(nd+1), intent(in) :: qb,pb,qk,pk
        complex*16, dimension(nh), intent(in) :: c
        complex*16, dimension(nd+1,nd+1), intent(in) :: Bb,Bk

        integer :: i
        real*8 :: a,q,Nsq,Y0
        real*8, dimension(nd) :: avec 
        real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
        real*8, dimension(nd+1,nd+1) :: Bmat
        complex*16, dimension(nh) :: csupp,cout,csupptest,csuppdiff
        complex*16, dimension(nh,nh) :: SttM,St0M,invS,X0Mat,S0tM
        complex*16, dimension(nh,nh) :: prod1,prod2,summa,S00M

!        write(*,*) "I AM C UPDATE FB"
        q = qb(1)
        Bmat = real(Bb)

!        write(*,*) q, Bmat(1,1)

        Nsq=fun_Nsq(nd+1,Bmat)
        call extractA(nd,Bmat,Amat,avec,a)
        call diagonalization(nd,Amat,LambdaMat,Tmat)
        Y0=int_Y0(nd,LambdaMat)
        X0Mat=int_XnMat(nd,0,a,avec,Amat,q)

        SttM = X0Mat*Y0*Nsq 

!        write(*,*) "SttM"
 
!        do i = 1,nh
!          write(*,*) SttM(i,:)
!        end do


!        write(*,*) "SttM^{-1}"
 
!        do i = 1,nh
!          write(*,*) invS(i,:)
!        end do

        q = qk(1)
        Bmat = real(Bk)

!        write(*,*) q, Bmat(1,1)

        Nsq=fun_Nsq(nd+1,Bmat)
        call extractA(nd,Bmat,Amat,avec,a)
        call diagonalization(nd,Amat,LambdaMat,Tmat)
        Y0=int_Y0(nd,LambdaMat)
        X0Mat=int_XnMat(nd,0,a,avec,Amat,q)

        S00M = X0Mat*Y0*Nsq 
        invS = invgen(nh,S00M) 
!       write(*,*) "S00M"
 
!        do i = 1,nh
!          write(*,*) S00M(i,:)
!        end do

        St0M = int_TauMat(nd,qb,qk,pb,pk,Bb,Bk)
!        write(*,*) "St0M"
 
!        do i = 1,nh
!          write(*,*) St0M(i,:)
!        end do
        S0tM = conjg(transpose(St0M))
!        write(*,*) "S0tM"
 
!        do i = 1,nh
!          write(*,*) S0tM(i,:)
!        end do

        prod1=matmul(invS,S0tM)
        prod2=matmul(St0M,prod1)

!        write(*,*) "prod2"
 
!        do i = 1,nh
!          write(*,*) prod2(i,:)
!        end do

        summa = 0.5*(prod2+SttM)
 
        !write(*,*) "summa"
 
!       do i = 1,nh
!          write(*,*) summa(i,:)
!        end do

        csupp = matmul(summa,c)
 
        cout = linsys(nh,S0tM,csupp) 

       end function

!......Numerical update of c............................................

       function c_update_full(nd,h,q0,qt,tq,p0,pt,tp,B0,Bt,tB,c0,tc) &
        & result(cout)
        ! nd: dimensions of the bath
        ! h: timestep
        ! q0,qt,tq: gaussian center at t, t+dt, t-dt
        ! p0,pt,tp: gaussian momentum at t, t+dt, t-dt
        ! B0,Bt,tB: gaussian width at t, t+dt, t-dt
        ! c0,tc: coefficients at t, t-dt
        integer, intent(in) :: nd
        real*8 :: h
        real*8, dimension(nd+1), intent(in) :: q0,qt,tq,p0,pt,tp 
        complex*16, dimension(nh), intent(in) :: c0,tc
        complex*16, dimension(nd+1,nd+1), intent(in) :: B0,Bt,tB

        integer :: i,j
        real*8 :: a,q,Nsq,Y0,reS,imS
        real*8, dimension(nd) :: avec 
        real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
        real*8, dimension(nd+1,nd+1) :: Bmat
        complex*16, dimension(nh) :: csupp,cout,cexpo
        complex*16, dimension(nh,nh) :: S00M,St0M,sumS,X0Mat,S0tM
        complex*16, dimension(nh,nh) :: prod1,prod2,summa,H00M,Aux
        real*8 :: E0,Mx
        real*8, dimension(nd) :: My

        ! S00M
        q = q0(1)
        Bmat = real(B0)

        Nsq=fun_Nsq(nd+1,Bmat)
        call extractA(nd,Bmat,Amat,avec,a)
        call diagonalization(nd,Amat,LambdaMat,Tmat)
        Y0=int_Y0(nd,LambdaMat)
        X0Mat=int_XnMat(nd,0,a,avec,Amat,q)

        S00M = X0Mat*Y0*Nsq 

        ! S0t
        St0M = int_TauMat(nd,q0,qt,p0,pt,B0,Bt)
        !St0M = int_TauMat(nd,qt,qt,pt,p0,Bt,B0)
        ! S0-t
        S0tM = int_TauMat(nd,q0,tq,p0,tp,B0,tB)
        ! H00M
        call energy(nd,q0,p0,c0,B0,H00M,E0,Mx,My)

! SG qtag.f line 1046 : cleaning new overlap (?)
        Aux(:,:) =- iu*(St0M-S0tM)

        do i = 1, nh
          sumS(i,i) = Aux(i,i)*complex(1.d0,0.d0) 
          do j = i+1,nh
            reS = dreal(Aux(i,j)+Aux(j,i))/2.d0
            imS = dimag(Aux(i,j)-Aux(j,i))/2.d0
            sumS(i,j) = reS+iu*imS
            sumS(j,i) = conjg(sumS(i,j))
          end do
        end do

! Hermitizing H00M too, why not 
        Aux(:,:) = H00M(:,:)

        do i = 1, nh
          H00M(i,i) = Aux(i,i)*complex(1.d0,0.d0) 
          do j = i+1,nh
            reS = dreal(Aux(i,j)+Aux(j,i))/2.d0
            imS = dimag(Aux(i,j)-Aux(j,i))/2.d0
            H00M(i,j) = reS+iu*imS
            H00M(j,i) = conjg(H00M(i,j))
          end do
        end do

       ! cexpo = c_static(nd,h,c0,dreal(S00M),H00M)
        summa = S0tM - St0M -2*h*iu*H00M
       ! summa = -iu*sumS -2*h*iu*H00M
       ! summa = S0tM - St0M 
        csupp = matmul(summa,c0) + matmul(S00M,tc)! +2.d0*(cexpo - c0)
 
        cout = linsys(nh,S00M,csupp) 

       end function
       end module
