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
       public :: c_static,c_update

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

        integer*8 :: i,lwork
        real*8 :: q,p,a,Nsq,Y0
        real*8, dimension(nd) :: qvec,pvec,avec
        real*8, dimension(nh) :: eigenv
        real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
        real*8, dimension(nd+1,nd+1) :: Bmat
        real*8, dimension(nh,nh) :: X0Mat
        complex*16, dimension(nh) :: c,csout,expvec
        complex*16, dimension(nh,nh) :: S00M,T00M,V00M,H00M
        complex*16, dimension(nh,nh) :: Z,adjZ,B

        complex*16, dimension(2*nh-1) :: lapwork
        real*8, dimension(3*nh-1) :: rwork

        ! External procedures defined in LAPACK
        external ZHEGV
        lwork = 2*nh-1

        write(*,*) "I AM  STATIC"

        write(*,*) qtot
        write(*,*) ptot
        write(*,*) cvec
        do i = 1,nd+1
          write(*,*) tildeBmat(i,:)
        end do

        q = qtot(1)
        p = ptot(1)
        qvec = qtot(2:nd+1)
        pvec = ptot(2:nd+1)

        Bmat = real(tildeBmat)

        Nsq=fun_Nsq(nd+1,Bmat)
        call extractA(nd,Bmat,Amat,avec,a)
        call diagonalization(nd,Amat,LambdaMat,Tmat)
        Y0=int_Y0(nd,LambdaMat)
        X0Mat=int_XnMat(nd,0,a,avec,Amat,q)

        S00M = X0Mat*Y0*Nsq 
        T00M = kin_energy(nd,q,p,qvec,pvec,tildeBmat)  
        V00M = fun_V0(nd,qtot,cvec,Bmat) 
        H00M = T00M + V00M
        ! Copy H00M so LAPACK can overwrite
        Z = H00M
        B = S00M

!        write(*,*) "I AM STATIC"

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
         
        call test_static(nh,S00M,H00M,Z,eigenv)

        c = cvec
        c = matmul(S00M,c)
        c = matmul(adjZ,c) 
       
        write(*,*) nd, nh
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

        integer*8 :: i,nd_double
        real*8, dimension(nd) :: eigenv
        complex*16, dimension(nd,nd) :: S,H
        complex*16, dimension(nd,nd) :: Z,adjZ,TEST1,TEST2,TEST6
        complex*16, dimension(nd,nd) :: TEST3,TEST4,TEST5,invS
        character(len=100) :: formato
      
        nd_double = nd
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
        invS = invgen(nd_double,S)
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

        write(*,*) "I AM C UPDATE"

        write(*,*) qb
        write(*,*) pb
        write(*,*) qk
        write(*,*) pk

        Bmat = real(Bb)

        Nsq=fun_Nsq(nd+1,Bmat)
        call extractA(nd,Bmat,Amat,avec,a)
        call diagonalization(nd,Amat,LambdaMat,Tmat)
        Y0=int_Y0(nd,LambdaMat)
        X0Mat=int_XnMat(nd,0,a,avec,Amat,q)

        S00M = X0Mat*Y0*Nsq 

        Tt0M = int_TauMat(nd,qb,qk,pb,pk,Bb,Bk) 
        csupp = matmul(Tt0M,c)

        write(*,*) "Tt0M:"
        do i = 1,nh
          write(*,*) Tt0M(i,:)
        end do

        cout = linsys(nh,S00M,csupp) 

       end function
       end module
