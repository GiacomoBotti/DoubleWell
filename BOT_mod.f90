!**********************************************************************!
! Module containing the functions required by the BOT_evo subroutine   !
! in evolution_mod.f90                                                 !
!**********************************************************************!

       module BOT_module

       use basisset_module
       use kinetic_module
       use effectivepot_module
       use integrals_module
       use inversion_module
       use matrix_module
       use constants

       implicit none

       integer*8 :: info

       private
       public :: c_static,c_update,c_update_fb,c_update_fbs

       contains

!......Coefficient update with static basis (analytical)................

       function c_static(nd,h,qtot,ptot,cvec,tildeBmat) result(csout)
       ! nd: bath dimension 
       ! h : time-step size
       ! q0: full position vector
       ! p0: full momenta vector
       ! tildeBmat: full complex gaussian width matrix
        implicit none
        integer, intent(in) :: nd
        real*8, intent(in) :: h
        real*8, dimension(nd+1), intent(in) :: qtot,ptot
        complex*16, dimension(nh), intent(in) :: cvec
        complex*16, dimension(nd+1,nd+1), intent(in) :: tildeBmat

        integer*8 :: i,lwork,nh8
        real*8 :: q,p,a,Nsq,Y0
        real*8, dimension(nd) :: qvec,pvec,avec
        real*8, dimension(nh) :: eigenv
        real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
        real*8, dimension(nd+1,nd+1) :: Bmat
        real*8, dimension(nh,nh) :: X4,X3,X2,X1,X0
        complex*16, dimension(nh) :: c,csout,expvec
        complex*16, dimension(nh,nh) :: S00M,T00M,V00M,H00M
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

        q = qtot(1)
        p = ptot(1)
        qvec = qtot(2:nd+1)
        pvec = ptot(2:nd+1)

        Bmat = real(tildeBmat)

        Nsq=fun_Nsq(nd+1,Bmat)
        call extractA(nd,Bmat,Amat,avec,a)
        call diagonalization(nd,Amat,LambdaMat,Tmat)
        Y0=int_Y0(nd,LambdaMat)
        X4=int_XnMat(nd,4,a,avec,Amat,qtot(1))
        X3=int_XnMat(nd,3,a,avec,Amat,qtot(1))
        X2=int_XnMat(nd,2,a,avec,Amat,qtot(1))
        X1=int_XnMat(nd,1,a,avec,Amat,qtot(1))
        X0=int_XnMat(nd,0,a,avec,Amat,q)

        S00M = X0*Y0*Nsq 
        T00M = kin_energy(nd,q,p,qvec,pvec,tildeBmat,Y0,X2,X1,X0)  
        V00M = fun_V0(nd,qtot,cvec,Bmat,Y0,X4,X2,X1,X0) 
        H00M = T00M + V00M
        ! Copy H00M so LAPACK can overwrite
        Z = H00M
        B = S00M
!        B = invgen(nh,S00M)

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

        call ZHEGV(1,'V','U',nh,Z,nh,B,nh,eigenv,&
                  &lapwork,lwork,rwork,info)

!         write(*,*) info
        adjZ = dconjg(transpose(Z))
         
!        call test_static(nh,S00M,H00M,Z,eigenv)

        c = cvec
!        write(*,*) "IN" 
!        write(*,*) c
        c = matmul(S00M,c)
!        write(*,*) "Sc" 
!        write(*,*) c
        c = matmul(adjZ,c) 
       
!        write(*,*) "Z^HSc" 
!        write(*,*) c

!        B(:,:) = 0.d0
        do i = 1,nh
          expvec(i) = zexp(-iu*eigenv(i)*h)*c(i)
!          B(i,i) = zexp(-iu*eigenv(i)*h)
        end do

!        B = matmul(B,adjZ)
!        B = matmul(Z,B)

!        write(*,*) "ZexpZ"
!        do i = 1,nh
!          write(*,*) B(i,:)
!        end do

!        write(*,*) "exp Z^HSc" 
!        write(*,*) expvec 
        csout = matmul(Z,expvec)
        
!        csout = expvec
!        write(*,*) "Zexp Z^HSc" 
!        write(*,*) "c final"
!        write(*,*) csout(:)

!.......1H DEBUGGING....................................................

!        write(*,*) "I am doing a 1H evolution!"
       ! csout(1) = cvec(1)*zexp(-iu*H00M(1,1)*h)/S00M(1,1)
      

!        csout(:) = 0.d0
!        do i= 1,nh
!           csout(i) = cvec(i)*zexp(-iu*H00M(i,i)*h)!/S00M(i,i)
!           csout(i) = cvec(i)*zexp(-iu*eigenv(i)*h)!/S00M(i,i)
!        end do

!        write(*,*) zexp(-iu*eigenv(1)*h), zexp(-iu*eigenv(2)*h)
!        write(*,*) zexp(-iu*H00M(1,1)*h), zexp(-iu*H00M(2,2)*h)

!        write(*,*) csout(:)

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

        integer :: i
        real*8 :: a,q,Nsq,Y0
        real*8, dimension(nd) :: avec 
        real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
        real*8, dimension(nd+1,nd+1) :: Bmat
        complex*16, dimension(nh) :: csupp,cout
        complex*16, dimension(nh,nh) :: S00M,Tt0M,invS,X0Mat

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

        S00M = X0Mat*Y0*Nsq 

!        invS = invgen(nh,S00M) 

!        write(*,*) "IN UPDATE"
!        write(*,*) c

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

        cout = linsys(nh,S00M,csupp) 

!        write(*,*) "cout"
!        write(*,*) cout
 
!        cout = matmul(invS,csupp)

!.......1H DEBUGGING....................................................

       ! write(*,*) "I am doing a 1H evolution!"
       ! cout(:) = 0.d0
       ! cout(1) = Tt0M(1,1)*c(1)/S00M(1,1)
       ! write(*,*) cout(nh)

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
        !csupptest = matmul(St0M,c)

       ! csuppdiff = csupp-csupptest 

       ! write(*,*) "csuppdiff"
       ! write(*,*) csuppdiff(:)
 
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
        !csupptest = matmul(St0M,c)

       ! csuppdiff = csupp-csupptest 

       ! write(*,*) "csuppdiff"
       ! write(*,*) csuppdiff(:)
 
        cout = linsys(nh,S0tM,csupp) 

       end function
       end module
