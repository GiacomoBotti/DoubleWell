!**********************************************************************!
! Module containing the functions required by the BOT_evo subroutine   !
! in evolution_mod.f90                                                 !
!**********************************************************************!

       module BOT_module

       use integral_module
       use smatrix_module
       use hmatrix_module
       use overlap_module
       use hamiltonian_module
       use inversion_module

       implicit none

       integer*8 :: info

       private
       public :: c_static,c_stat_num,c_update,c_upd_num

       contains

!......Coefficient update with static basis (analytical)................

       function c_static(nd,npar,y,work,h) result(csout)
       ! nd: basis set dimension
       ! npar: parameter space dimension
       ! y: runge-kutta function - [parameters - coefficients]
       ! work: work array for I/O
       ! h : time-step size
        implicit none
        integer*8, intent(in) :: nd,npar
        real*8, intent(in) :: h
        complex*16, dimension(npar+nd), intent(in) :: y
        real*8, dimension(4), intent(in) :: work

        integer*8 :: i,lwork
        real*8 :: q,p,a,alpha,xi,zeta
        real*8, dimension(nd) :: eigenv
        complex*16, dimension(nd) :: c,csout,expvec
        complex*16, dimension(nd,nd) :: S00M,T00M,V00M,H00M,B
        complex*16, dimension(nd,nd) :: Z,adjZ

        complex*16, dimension(2*nd-1) :: lapwork
        complex*16, dimension(3*nd-1) :: rwork

        ! External procedures defined in LAPACK
        external ZHEGV
        lwork = 2*nd-1

        ! extract work
        a = work(1)
        alpha = work(2)
        xi = work(3)
        zeta = work(4)

        ! extract y
        q = dreal(y(1))
        p = dreal(y(2))
        c(:) = y(npar+1:npar+nd)

        S00M = S00maple(alpha,q,p,xi,zeta)
        T00M = T00maple(a,alpha,q,p,xi,zeta) 
        V00M = V00maple(a,alpha,q,p,xi,zeta) 
        H00M = T00M + V00M
        ! Copy H00M so LAPACK can overwrite
        Z = H00M
        B = S00M

        call ZHEGV(1,'V','U',nd,Z,nd,B,nd,eigenv,&
                  &lapwork,lwork,rwork,info)

!         write(*,*) info
        adjZ = dconjg(transpose(Z))
         
!        call test_static(nd,S00M,H00M,Z,eigenv)

        c = matmul(S00M,c)
        c = matmul(adjZ,c) 
       
        do i = 1,nd
           expvec(i) = zexp(-(0,1.d0)*eigenv(i)*h)*c(i)
        end do

        csout = matmul(Z,expvec)

       end function

!......Coefficient update with static basis (numerical).................

       function c_stat_num(nd,npar,y,work,h) result(csout)
       ! nd: basis set dimension
       ! npar: parameter space dimension
       ! y: runge-kutta function - [parameters - coefficients]
       ! work: work array for I/O
       ! h : time-step size
        implicit none
        integer*8, intent(in) :: nd,npar
        real*8, intent(in) :: h
        complex*16, dimension(npar+nd), intent(in) :: y
        real*8, dimension(4), intent(in) :: work

        integer*8 :: i,lwork
        real*8 :: q,p,a,alpha,xi,zeta,Tn
        ! PUT IT IN A MODULEEEEEEEEEE
        real*8 :: hm = 10.d0 !heavy mass 
        real*8, dimension(nd) :: eigenv
        complex*16, dimension(nd) :: c,csout,expvec
        complex*16, dimension(nd,nd) :: S00mat,T00mat,V00mat,H00mat
        complex*16, dimension(nd,nd) :: Z,adjZ,B
        complex*16 :: nover,ndyyover
        !real*8, dimension(nd) :: vecxi,veczeta
        real*8, dimension(3) :: cond,nparam 

        complex*16, dimension(2*nd-1) :: lapwork
        complex*16, dimension(3*nd-1) :: rwork

        ! External procedures defined in LAPACK
        external ZHEGV
        lwork = 2*nd-1

        ! extract work
        a = work(1)
        alpha = work(2)
        xi = work(3)
        zeta = work(4)

        ! extract y
        q = dreal(y(1))
        p = dreal(y(2))
        c(:) = y(npar+1:npar+nd)
        
        ! build vecxi and veczeta
        !vecxi(1) = xi
        !vecxi(2) = -1.d0*xi
        !veczeta(1) = zeta
        !veczeta(2) = -1.d0*zeta
        
        ! NUCLEAR OVERLAPS
        cond = [-10.d0,10.d0,1000.d0]
        nparam = [a,q,p]
        nover = nbraket1D(nwfn,nwfn,cond,nparam)

        ndyyover =nbraket1D(nwfn,dyynwfn,cond,nparam)
        Tn = -dreal(ndyyover/(2.d0*hm)) 

        ! ELECTRONIC OVERLAP
        S00mat = S00(nd,q,p,alpha,vecxi,veczeta)

        ! ELECTRONIC KINETIC MATRIX
        T00mat = T00(nd,q,p,alpha,vecxi,veczeta)

        ! TOTAL OVERLAP MATRICES
        ! WATCH OUT TO USE THE ELECTRONIC ONES
        ! SO FIRST COMPUTE THE ONES THAT YOU
        ! ARE NOT GOING TO USE
        T00mat = Tn*S00mat + nover*T00mat
        T00mat = dreal(T00mat) 
        S00mat = S00mat*nover 

        V00mat = V00(nd,q,p,a,alpha,vecxi,veczeta)
        V00mat = dreal(V00mat)

        H00mat = T00mat + V00mat
        ! Copy H00M so LAPACK can overwrite
        Z = H00mat
        B = S00mat

        call ZHEGV(1,'V','U',nd,Z,nd,B,nd,eigenv,&
                  &lapwork,lwork,rwork,info)

        adjZ = dconjg(transpose(Z))

!        call test_static(nd,S00mat,H00mat,Z,eigenv)

        c = matmul(S00mat,c)
        c = matmul(adjZ,c) 
       
        do i = 1,nd
           expvec(i) = zexp(-(0,1.d0)*eigenv(i)*h)*c(i)
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
        integer*8, intent(in) :: nd

        integer*8 :: i
        real*8, dimension(nd) :: eigenv
        complex*16, dimension(nd,nd) :: S,H
        complex*16, dimension(nd,nd) :: Z,adjZ,TEST1,TEST2,TEST6
        complex*16, dimension(nd,nd) :: TEST3,TEST4,TEST5,invS
        character(len=100) :: formato
      
        formato="(E10.1,E10.1,E10.1,E10.1)"
  
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
        write(*,formato) TEST2(1,1)-1.d0, TEST2(1,2)
        write(*,formato) TEST2(2,1), TEST2(2,2)-1.d0
        write(*,*) " "
        invS = invgen(nd,S)
        TEST3 = matmul(Z,adjZ)
        write(*,*) "ZZ*=S-1"
        write(*,formato) TEST3(1,1)-invS(1,1), TEST3(1,2)-invS(1,2)
        write(*,formato) TEST3(2,1)-invS(2,1), TEST3(2,2)-invS(2,2)
        write(*,*) " "
        write(*,*) "Z^H*H*Z=L"
        write(*,formato) TEST4(1,1)-TEST5(1,1), TEST4(1,2)-TEST5(1,2)
        write(*,formato) TEST4(2,1)-TEST5(2,1), TEST4(2,2)-TEST5(2,2)
        write(*,*) " "

       end subroutine

!......Analytical update of electronic coefficients.....................

       function c_update(nd,npar,lambdat,lambda0,c,work) result(cout)
       ! nd : basis set dimension
       ! npar : parameters space dimension
       ! lambdat : parameters at new step
       ! lambda0 : parameters at old step
       ! c : exponentially-evolved coefficienst
       ! work : a, alpha, xi, zeta
        integer*8, intent(in) :: nd,npar
        real*8, dimension(npar), intent(in) :: lambdat, lambda0
        real*8, dimension(4), intent(in) :: work
        complex*16, dimension(nd), intent(in) :: c

        real*8 :: a,alpha,xi,zeta,q,p,qold,pold
        complex*16, dimension(nd) :: csupp,cout
        complex*16, dimension(nd,nd) :: S00M,Tt0M,invS

        ! extract work
        a = work(1)
        alpha = work(2)
        xi = work(3)
        zeta = work(4)

        ! extract lambdas
        q = lambdat(1)
        p = lambdat(2)
        qold = lambda0(1)
        pold = lambda0(2)

        S00M = S00maple(alpha,q,p,xi,zeta)
        Tt0M = Tt0maple(a,alpha,q,p,qold,pold,xi,zeta)
        csupp = matmul(Tt0M,c)

        cout = linsys(npar,S00M,csupp) 

       end function

!......Numerical update of electronic coefficients......................

       function c_upd_num(nd,npar,lambdat,lambda0,c,work) result(cout)
       ! nd : basis set dimension
       ! npar : parameters space dimension
       ! lambdat : parameters at new step
       ! lambda0 : parameters at old step
       ! c : exponentially-evolved coefficienst
       ! work : a, alpha, xi, zeta
        integer*8, intent(in) :: nd,npar
        real*8, dimension(npar), intent(in) :: lambdat, lambda0
        real*8, dimension(4), intent(in) :: work
        complex*16, dimension(nd), intent(in) :: c

        real*8 :: a,alpha,xi,zeta,q,p,qold,pold,nover
        complex*16 :: tover
        real*8, dimension(3) :: cond, nparam
        real*8, dimension(5) :: tparam
        !real*8, dimension(nd) :: vecxi,veczeta
        complex*16, dimension(nd) :: csupp,cout
        complex*16, dimension(nd,nd) :: S00M,Tt0M,invS

        ! extract work
        a = work(1)
        alpha = work(2)
        xi = work(3)
        zeta = work(4)

        ! extract lambdas
        q = lambdat(1)
        p = lambdat(2)
        qold = lambda0(1)
        pold = lambda0(2)
        
        ! build vecxi and veczeta
        !vecxi(1) = xi
        !vecxi(2) = -1.d0*xi
        !veczeta(1) = zeta
        !veczeta(2) = -1.d0*zeta
        
        ! NUCLEAR OVERLAPS
        cond = [-10.d0,10.d0,1000.d0]
        nparam = [a,q,p]
        tparam = [a,q,p,qold,pold]
        nover = nbraket1D(nwfn,nwfn,cond,nparam)
        tover = ntbraket1D(nwfn,nwfn,cond,tparam)

        S00M = S00(nd,q,p,alpha,vecxi,veczeta)
        S00M = S00M*nover
        !S00M = S00maple(alpha,q,p,xi,zeta)
        Tt0M = Tt0(nd,q,p,qold,pold,alpha,vecxi,veczeta)
        Tt0M = Tt0M*tover
        !Tt0M = Tt0maple(a,alpha,q,p,qold,pold,xi,zeta)
        csupp = matmul(Tt0M,c)

        cout = linsys(npar,S00M,csupp) 

       end function




       end module


