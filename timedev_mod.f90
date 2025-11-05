!**********************************************************************!
! Module containing the time derivatives equations (analytical)        !
!**********************************************************************!

      module timedev_module

      use constants
      use hamiltonian_module
      use overlap_module
      use potential_module 
      use inversion_module 

      implicit none

      private
      public :: dotcoeff,inv2D,dvector,Bmatrix,dotlambda

      contains

!.....Coefficients time-derivative......................................

      function dotcoeff(nd,npar,q,p,a,alpha,xi,zeta,c) result(dotcout)
      ! nd: electronic bs dimensions
      ! npar: parameter space dimensions
      ! q,p: nucl gau phase space center
      ! alpha: elec gau width
      ! xi: pos scaling factors
      ! zeta: mom scaling factors
        implicit none
        integer*8, intent(in) :: nd,npar
        real*8, intent(in) :: q,p,a,alpha,xi,zeta
        complex*16, dimension(nd), intent(in) :: c

        integer*8 :: i,j
        ! A module for fuck's sake
        real*8 :: hm = 10.d0
        real*8, dimension(npar) :: dotlam
        real*8, dimension(2) :: grad
        complex*16, dimension(nd) :: dotcout,c_conj,sysvec
        complex*16, dimension(nd,nd) :: T00M,V00M,H00M
        complex*16, dimension(nd,nd) :: S00M,S0pM,S0qM,MatSum

        T00M = T00maple(a,alpha,q,p,xi,zeta)
        V00M = V00maple(a,alpha,q,p,xi,zeta)
        H00M = T00M + V00M
        
        S0pM = S0pmaple(alpha,q,p,xi,zeta)
        S0qM = S0qmaple(alpha,q,p,xi,zeta)
        S00M = S00maple(alpha,q,p,xi,zeta)

        ! dotlambda from classical evolution
        !dotlam(1) = p/hm
        !grad = -1.d0*vgrad(0.d0,q)
        !dotlam(2) = grad(2) 
        
        dotlam = dotlambda(nd,npar,q,p,a,alpha,xi,zeta,c)

        MatSum = H00M -(0,1)*(dotlam(1)*S0qM+dotlam(2)*S0pM)
         
        sysvec = (0,-1)*matmul(MatSum,c) 

        dotcout = linsys(nd,S00M,sysvec)

      end function

!.....d vector..........................................................

      function dvector(nd,npar,q,p,a,alpha,xi,zeta,c) result(dout)
      ! nd: electronic bs dimensions
      ! npar: parameter space dimensions
      ! q,p: nucl gau phase space center
      ! a,alpha: nucl,elec gau width
      ! xi: pos scaling factors
      ! zeta: mom scaling factors
      ! c: coefficient vector
        implicit none
        integer*8, intent(in) :: nd,npar
        real*8, intent(in) :: q,p,a,alpha,xi,zeta
        complex*16, dimension(nd), intent(in) :: c

        integer*8 :: i,j
        real*8 :: hm = 10.d0
        complex*16, dimension(npar) :: dout
        complex*16 :: dq,dp
        complex*16, dimension(nd) :: dqvec,dpvec 
        complex*16, dimension(nd,nd) :: T00M,V00M,H00M
        complex*16, dimension(nd,nd) :: Tq0M,Vq0M,Hq0M
        complex*16, dimension(nd,nd) :: Tp0M,Vp0M,Hp0M
        complex*16, dimension(nd,nd) :: Sp0M,Sq0M,S00M,invS00M
        complex*16, dimension(nd,nd) :: Dqmat,Dpmat 

        T00M = T00maple(a,alpha,q,p,xi,zeta)
        V00M = V00maple(a,alpha,q,p,xi,zeta)
        H00M = T00M + V00M

        Tq0M = Tq0maple(a,alpha,q,p,xi,zeta)
        Vq0M = Vq0maple(a,alpha,q,p,xi,zeta)
        Hq0M = Tq0M + Vq0M

        Tp0M = Tp0maple(a,alpha,q,p,xi,zeta)
        Vp0M = Vp0maple(a,alpha,q,p,xi,zeta)
        Hp0M = Tp0M + Vp0M

        S00M = S00maple(alpha,q,p,xi,zeta)
        Sq0M = Sq0maple(alpha,q,p,xi,zeta)
        Sp0M = Sp0maple(alpha,q,p,xi,zeta)
        !invS00M = inv2D(S00M)
        invS00M = invgen(npar,S00M)

        Dqmat = Hq0M - matmul(Sq0M,matmul(invS00M,H00M))
        Dpmat = Hp0M - matmul(Sp0M,matmul(invS00M,H00M))

        dqvec = matmul(Dqmat,c)
        dpvec = matmul(Dpmat,c)

        ! if vectors are complex, the first one
        ! is automatically conjugated
        dq = dot_product(c,dqvec) 
        dp = dot_product(c,dpvec) 

        dout(1) = dq
        dout(2) = dp

      end function

!.....B matrix..........................................................

      function Bmatrix(nd,npar,q,p,a,alpha,xi,zeta,c) result(Bout)
      ! nd: electronic bs dimensions
      ! npar: parameter space dimensions
      ! q,p: nucl gau phase space center
      ! a,alpha: nucl,elec gau width
      ! xi: pos scaling factors
      ! zeta: mom scaling factors
      ! c: coefficient vector
        implicit none
        integer*8, intent(in) :: nd,npar
        real*8, intent(in) :: q,p,a,alpha,xi,zeta
        complex*16, dimension(nd), intent(in) :: c

        integer*8 :: i,j
        real*8 :: hm = 10.d0
        complex*16, dimension(npar,npar) :: Bout
        complex*16 :: Bqq,Bqp,Bpq,Bpp
        complex*16, dimension(nd) :: Bqqvec,Bqpvec,Bpqvec,Bppvec 
        complex*16, dimension(nd,nd) :: Sp0M,Sq0M,S00M,invS00M
        complex*16, dimension(nd,nd) :: SpqM,SqpM,SppM,SqqM
        complex*16, dimension(nd,nd) :: S0qM,S0pM
        complex*16, dimension(nd,nd) :: Bqqmat,Bqpmat,Bpqmat,Bppmat


        S00M = S00maple(alpha,q,p,xi,zeta)
        Sq0M = Sq0maple(alpha,q,p,xi,zeta)
        S0qM = S0qmaple(alpha,q,p,xi,zeta)
        SqqM = Sqqmaple(a,alpha,q,p,xi,zeta)
        Sp0M = Sp0maple(alpha,q,p,xi,zeta)
        S0pM = S0pmaple(alpha,q,p,xi,zeta)
        SppM = Sppmaple(a,alpha,q,p,xi,zeta)
        SqpM = Sqpmaple(alpha,q,p,xi,zeta)
        SpqM = Spqmaple(alpha,q,p,xi,zeta)
        !invS00M = inv2D(S00M)
        invS00M = invgen(npar,S00M)

        Bqqmat = SqqM - matmul(Sq0M,matmul(invS00M,S0qM))
        Bqpmat = SqpM - matmul(Sq0M,matmul(invS00M,S0pM))
        Bpqmat = SpqM - matmul(Sp0M,matmul(invS00M,S0qM))
        Bppmat = SppM - matmul(Sp0M,matmul(invS00M,S0pM))
 
        Bqqvec = matmul(Bqqmat,c)
        Bqpvec = matmul(Bqpmat,c)
        Bpqvec = matmul(Bpqmat,c)
        Bppvec = matmul(Bppmat,c)

        Bqq = dot_product(c,Bqqvec)
        Bqp = dot_product(c,Bqpvec)
        Bpq = dot_product(c,Bpqvec)
        Bpp = dot_product(c,Bppvec)

        Bout(1,1) = Bqq
        Bout(1,2) = Bqp
        Bout(2,1) = Bpq
        Bout(2,2) = Bpp

      end function 

!.....Parameters time derivative........................................

      function dotlambda(nd,npar,q,p,a,alpha,xi,zeta,c) result(lamout)
      ! nd: electronic bs dimensions
      ! npar: parameter space dimensions
      ! q,p: nucl gau phase space center
      ! a,alpha: nucl,elec gau width
      ! xi: pos scaling factors
      ! zeta: mom scaling factors
      ! c: coefficient vector
        implicit none
        integer*8, intent(in) :: nd,npar
        real*8, intent(in) :: q,p,a,alpha,xi,zeta
        complex*16, dimension(nd), intent(in) :: c

        real*8, dimension(npar) :: lamout
        complex*16, dimension(npar) :: dvec,imdvec,Redvec
        complex*16, dimension(npar,npar) :: Bmat,ReBmat,invBmat,imBmat

        dvec = dvector(nd,npar,q,p,a,alpha,xi,zeta,c)
        Bmat = Bmatrix(nd,npar,q,p,a,alpha,xi,zeta,c)

        imdvec = dimag(dvec)
        ReBmat = dreal(Bmat)
        Redvec = -dreal(dvec)
        imBmat = dimag(Bmat)
        !invBmat = inv2D(ReBmat)
        !invBmat = invgen(npar,ReBmat)

        !lamout = matmul(invBmat,imdvec)
        lamout = linsys(npar,ReBmat,imdvec) 
        !lamout = linsys(npar,imBmat,Redvec) 

      end function        

      end module
