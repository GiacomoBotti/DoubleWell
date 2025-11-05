!***********************************************************************!
! Fortran code to evolve a multiconfigurational wavefunction using      !
! McLachlan variational principle, applied to a double well             !
!***********************************************************************!

      program doublewell

      use potential_module 
      use basisset_module 
      use constants
      use check_module
      use integral_module
      use smatrix_module
      use hmatrix_module
      use evolution_module

      implicit none

      integer*8 :: i,j,nd,npar
      integer*8,dimension(3) :: trj
      real*8 :: V,N_el,N_nucl,rho,NN_el
      complex*16 :: N,elec_g,elec_wfn,elec_g1,elec_g2,psi
      real*8 :: ifloat,jfloat

      real*8 :: q,p
      real*8 :: a,alpha
      real*8 :: xi,zeta

      real*8, dimension(5) :: eparam
      real*8, dimension(4) :: work 
      complex*16, dimension(4) :: y 
      real*8, dimension(3) :: init,nparam
      real*8, dimension(2,2) :: test_H

      complex*16, dimension(2) :: c_par,c_conj

!.....Check Potential Hessian...........................................

      call check_hessian(dsqrt(dlog(40.d0)),a)

!.....Useful parameters.................................................
      write(*,*) "vecxi = ",vecxi

      alpha = 1.d0
      xi = vecxi(1) !(x+xi*q)^2
      q = dsqrt(dlog(40.d0))
      zeta = -1.d0*vecxi(1)*alpha/a !+i*zeta*p

      call zetakahler(a,alpha)

      write(*,*) "zeta = ",zeta
      write(*,*) "veczeta = ",veczeta
      p = 0.5d0

      !c_par(1) =complex(1.d0/dsqrt(2.d0),1.d0/dsqrt(2.d0))
      !c_par(2) =complex(1.d0/dsqrt(2.d0),1.d0/dsqrt(2.d0))
      c_par(1) =complex(1.d0,0.d0)
      c_par(2) =complex(0.d0,1.d0)

!.....Check Potential function..........................................
!      call plot_potential()
!.....Check Nuclear wfn function........................................
!      call plot_nwfn(a,q,p)
!.....Check Potential Gradient..........................................
!      call plot_grad()
!.....Check Electronic gaussian function................................
!      call plot_egau(alpha,vecxi(1),q,zeta,1.d0)
!.....Check Electronic wavefunction.....................................
!      call plot_2d_ewfn(alpha,vecxi,q,zeta,1.d0)
!.....Density plot q=4..................................................
!      call density(4.d0,c_par,a,p,alpha,xi,zeta) 
!.....Density plot q=1..................................................
!      call density(1.d0,c_par,a,p,alpha,xi,zeta)
!.....Nuclear derivatives plot..........................................
!      call plot_nucl_wfn_der(a,q,1.d0) 
!.....Electronic derivatives plot.......................................
!      call plot_egau_der(alpha,xi,q,zeta,1.d0)
!.....Normalization at start check......................................
!      call normalization(c_par,q,p,a,alpha,xi,zeta) 
!.....Check S...........................................................
!      call S00_check(a,alpha,xi,zeta)
!.....Check H...........................................................
!      call H00_check(a,alpha,xi,zeta)
!.....Check 2D integral.................................................
!      call check_2D(a,alpha,xi,zeta)
!.....Check coeff evolution.............................................
!      call check_coeff(q,1.d0,c_par,a,alpha,xi,zeta)
!.....Check 2x2 inversion...............................................
!      call check_inv2D(q,p,alpha,xi,zeta)
!.....Check d vector....................................................
!      call check_d_vec(a,alpha,xi,zeta,c_par)
!.....Check hamiltonian.................................................     
!      call hamexptval(q,1.d0,a,alpha,xi,zeta,c_par)
!.....Check homogeneus term.............................................
!      call homo_check(a,alpha,xi,zeta,c_par)
!.....Check BOT.........................................................
!      call check_bot(a,alpha,xi,zeta,c_par)

!.....Dynamics????......................................................
      nd = 2
      npar = 2

      work(1) = a
      work(2) = alpha
      work(3) = xi
      work(4) = zeta

      y(1) = q 
      y(2) = p 
      y(3:4) = c_par(:)

      trj = [0,4,200]
       
      call rungekutta(nd,npar,trj,y,work)
!      call bot_evo(nd,npar,trj,y,work)
!      call scprop(nd,npar,trj,y,work)
!      call bot_scp(nd,npar,trj,y,work)

      end program

