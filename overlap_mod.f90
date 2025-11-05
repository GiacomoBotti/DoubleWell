!**********************************************************************!
! Module containing all the analytical integrals extracted from MAPLE  !
! worksheet h2plus_param_gencode.mw                                    !
!                                                                      !
! It assumes 2g ewfn w/                                                !
! vecxi = [xi, -xi], xi <0                                             !
! veczeta = [zeta, -zeta], zeta >0                                     !
!                                                                      !
! Right now: used to debug numerical integrals                         !
! In the future: used instead of the numerical integrals               !
!**********************************************************************!

      module overlap_module
      
      use constants

      implicit none

      integer*8 :: nd =2
      private
      public :: S00maple,S0qmaple,Sq0maple,S0pmaple,Sp0maple
      public :: Sppmaple,Sqpmaple,Spqmaple,Sqqmaple,Tt0maple

      contains

!.....S00 matrix........................................................
      
      function S00maple(alpha,q,p,xi,zeta) result(S00out)
      ! alpha: elec gau width
      ! q,p: nucl gau phase space center
      ! xi,zeta: pos,mom scaling factors
        implicit none
        real*8, intent(in) :: alpha,q,p,xi,zeta

        complex*16, dimension(nd,nd) :: S00out

        S00out(1,1) = 1.d0
        S00out(1,2) = dexp((-(alpha*q*xi)**2.d0-(p*zeta)**2.d0)/alpha)
        S00out(2,1) = dexp((-(alpha*q*xi)**2.d0-(p*zeta)**2.d0)/alpha)
        S00out(2,2) = 1.d0
 
      end function

!.....S0q matrix........................................................
       
      function S0qmaple(alpha,q,p,xi,zeta) result(S0qout)
      ! alpha: elec gau width
      ! q,p: nucl gau phase space center
      ! xi,zeta: pos,mom scaling factors
        implicit none
        real*8, intent(in) :: alpha,q,p,xi,zeta

        complex*16, dimension(nd,nd) :: S0qout

        S0qout(1,1)=(0.d0,-1.d0)*p 
        S0qout(1,2)=-dexp((-(alpha*q*xi)**2.d0-(p*zeta)**2.d0)/alpha)
        S0qout(1,2)=S0qout(1,2)*((0,1)*p*zeta*xi+alpha*q*xi**2+(0,1)*p)
        S0qout(2,1)=-dexp((-(alpha*q*xi)**2.d0-(p*zeta)**2.d0)/alpha)
        S0qout(2,1)=S0qout(2,1)*((0,1)*p*zeta*xi+alpha*q*xi**2+(0,1)*p)
        S0qout(2,2)=(0.d0,-1.d0)*p

      end function

!.....Sq0 matrix........................................................

      function Sq0maple(alpha,q,p,xi,zeta) result(Sq0out)
      ! alpha: elec gau width
      ! q,p: nucl gau phase space center
      ! xi,zeta: pos,mom scaling factors
        implicit none
        real*8, intent(in) :: alpha,q,p,xi,zeta

        complex*16, dimension(nd,nd) :: Sq0out

        Sq0out(1,1)=(0.d0,1.d0)*p 
        Sq0out(1,2)=dexp((-(alpha*q*xi)**2.d0-(p*zeta)**2.d0)/alpha)
        Sq0out(1,2)=Sq0out(1,2)*((0,1)*p*zeta*xi-alpha*q*xi**2+(0,1)*p)
        Sq0out(2,1)=dexp((-(alpha*q*xi)**2.d0-(p*zeta)**2.d0)/alpha)
        Sq0out(2,1)=Sq0out(2,1)*((0,1)*p*zeta*xi-alpha*q*xi**2+(0,1)*p)
        Sq0out(2,2)=(0.d0,1.d0)*p

      end function

!.....S0p matrix........................................................

      function S0pmaple(alpha,q,p,xi,zeta) result(S0pout)
      ! alpha: elec gau width
      ! q,p: nucl gau phase space center
      ! xi,zeta: pos,mom scaling factors
        implicit none
        real*8, intent(in) :: alpha,q,p,xi,zeta

        complex*16, dimension(nd,nd) :: S0pout

        S0pout(1,1)=(0.d0,-1.d0)*q*xi*zeta 
        S0pout(1,2)=-dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)
        S0pout(1,2)=S0pout(1,2)*p*zeta**2/alpha
        S0pout(2,1)=-dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)
        S0pout(2,1)=S0pout(2,1)*p*zeta**2/alpha
        S0pout(2,2)=(0.d0,-1.d0)*q*xi*zeta
      
      end function

!.....Sp0 matrix........................................................

      function Sp0maple(alpha,q,p,xi,zeta) result(Sp0out)
      ! alpha: elec gau width
      ! q,p: nucl gau phase space center
      ! xi,zeta: pos,mom scaling factors
        implicit none
        real*8, intent(in) :: alpha,q,p,xi,zeta

        complex*16, dimension(nd,nd) :: Sp0out

        Sp0out(1,1)=(0.d0,1.d0)*q*xi*zeta 
        Sp0out(1,2)=-dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)
        Sp0out(1,2)=Sp0out(1,2)*p*zeta**2/alpha
        Sp0out(2,1)=-dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)
        Sp0out(2,1)=Sp0out(2,1)*p*zeta**2/alpha
        Sp0out(2,2)=(0.d0,1.d0)*q*xi*zeta
      
      end function

!.....Spp matrix........................................................

      function Sppmaple(a,alpha,q,p,xi,zeta) result(Sppout)
      ! a,alpha: nucl,elec gau width
      ! q,p: nucl gau phase space center
      ! xi,zeta: pos,mom scaling factors
        implicit none
        real*8, intent(in) :: a,alpha,q,p,xi,zeta

        complex*16, dimension(nd,nd) :: Sppout

        Sppout(1,1)=zeta**2*(2*alpha*q**2*xi**2+1)/alpha/2.d0+1.d0/a/2.d0
        Sppout(1,2)=dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)&
                    &*(2.d0*a*p**2*zeta**4-alpha*a*zeta**2+alpha**2)&
                    &/alpha**2/a/2.d0
        Sppout(2,1)=dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)&
                    &*(2.d0*a*p**2*zeta**4-alpha*a*zeta**2+alpha**2)&
                    &/alpha**2/a/2.d0
        Sppout(2,2)=zeta**2*(2*alpha*q**2*xi**2+1)/alpha/2.d0+1.d0/a/2.d0
      
      end function

!.....Sqp matrix........................................................

      function Sqpmaple(alpha,q,p,xi,zeta) result(Sqpout)
      ! alpha: elec gau width
      ! q,p: nucl gau phase space center
      ! xi,zeta: pos,mom scaling factors
        implicit none
        real*8, intent(in) :: alpha,q,p,xi,zeta

        complex*16, dimension(nd,nd) :: Sqpout

        Sqpout(1,1)=(0,1.d0)/2.d0-(-(2*p*q)+(0,1))*xi*zeta/2.d0
        Sqpout(1,2)=-dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)&
                    &*((-p*q*(xi*zeta)**2+(0, -1.d0)/2.d0*xi*zeta&
                    &+(0,-1.d0)/2.d0)*alpha+(0,1)*(xi*zeta+1)*&
                    &(zeta*p)**2)/alpha
        Sqpout(2,1)=-dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)&
                    &*((-p*q*(xi*zeta)**2+(0, -1.d0)/2.d0*xi*zeta&
                    &+(0,-1.d0)/2.d0)*alpha+(0,1)*(xi*zeta+1)*&
                    &(zeta*p)**2)/alpha
        Sqpout(2,2)=(0,1.d0)/2.d0-(-(2*p*q)+(0,1))*xi*zeta/2.d0
      
      end function

!.....Spq matrix........................................................

      function Spqmaple(alpha,q,p,xi,zeta) result(Spqout)
      ! alpha: elec gau width
      ! q,p: nucl gau phase space center
      ! xi,zeta: pos,mom scaling factors
        implicit none
        real*8, intent(in) :: alpha,q,p,xi,zeta

        complex*16, dimension(nd,nd) :: Spqout

        Spqout(1,1)=(0,-1.d0)/2.d0+((2*p*q)+(0,1))*xi*zeta/2.d0
        Spqout(1,2)=dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)&
                    &*((p*q*(xi*zeta)**2+(0, -1.d0)/2.d0*xi*zeta&
                    &+(0,-1.d0)/2.d0)*alpha+(0,1)*(xi*zeta+1)*&
                    &(zeta*p)**2)/alpha
        Spqout(2,1)=dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)&
                    &*((p*q*(xi*zeta)**2+(0, -1.d0)/2.d0*xi*zeta&
                    &+(0,-1.d0)/2.d0)*alpha+(0,1)*(xi*zeta+1)*&
                    &(zeta*p)**2)/alpha
        Spqout(2,2)=(0,-1.d0)/2.d0+((2*p*q)+(0,1))*xi*zeta/2.d0
      
      end function

!.....Sqq matrix........................................................

      function Sqqmaple(a,alpha,q,p,xi,zeta) result(Sqqout)
      ! a,alpha: nucl,elec gau width
      ! q,p: nucl gau phase space center
      ! xi,zeta: pos,mom scaling factors
        implicit none
        real*8, intent(in) :: a,alpha,q,p,xi,zeta

        complex*16, dimension(nd,nd) :: Sqqout

        Sqqout(1,1)=alpha*xi**2/2.d0+p**2+a/2.d0
        Sqqout(1,2)=dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)&
                    &*((alpha*q)**2*xi**4+((p*zeta)**2-alpha/2.d0)&
                    &*xi**2+2.d0*p**2*xi*zeta+p**2+a/2.d0)
        Sqqout(2,1)=dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)&
                    &*((alpha*q)**2*xi**4+((p*zeta)**2-alpha/2.d0)&
                    &*xi**2+2.d0*p**2*xi*zeta+p**2+a/2.d0)
        Sqqout(2,2)=alpha*xi**2/2.d0+p**2+a/2.d0
      
      end function

!.....Tt0 matrix........................................................

      function Tt0maple(a,alpha,q,p,qold,pold,xi,zeta) result(Tt0out)
      ! a,alpha: nucl,elec gau width
      ! q,p: new nucl gau phase space center
      ! qold,pold: old nucl gau phase space center
      ! xi,zeta: pos,mom scaling factors
        implicit none
        real*8, intent(in) :: a,alpha,q,p,qold,pold,xi,zeta

        complex*16, dimension(nd,nd) :: Tt0out

        Tt0out(1,1) = zexp((-a*zeta**2*(p-pold)**2+2.d0*(0,1)*&
                      &(q+qold)*alpha*xi*a*(p-pold)*zeta+2.d0*alpha*&
                      &(-(q-qold)**2*a**2/2+(-alpha*(q-qold)*xi**2/2+&
                      &(0,1)*p+(0,1)*pold)*(q-qold)*a-(p-pold)**2/2))&
                      &/a/alpha/4.d0)
        Tt0out(1,2) = zexp((-a*xi**2*(q+qold)**2*alpha**2+(-(q-qold)**2&
                      &*a**2+2.d0*(0,1)*(p+pold)*(q-qold)*(zeta*xi+1)&
                      &*a-(p-pold)**2)*alpha-a*zeta**2*(p+pold)**2)&
                      &/a/alpha/4.d0)
        Tt0out(2,1) = zexp((-a*xi**2*(q+qold)**2*alpha**2+(-(q-qold)**2&
                      &*a**2+2.d0*(0,1)*(p+pold)*(q-qold)*(zeta*xi+1)&
                      &*a-(p-pold)**2)*alpha-a*zeta**2*(p+pold)**2)&
                      &/a/alpha/4.d0)
        Tt0out(2,2) = zexp((-a*zeta**2*(p-pold)**2+2.d0*(0,1)*&
                      &(q+qold)*alpha*xi*a*(p-pold)*zeta+2.d0*alpha*&
                      &(-(q-qold)**2*a**2/2+(-alpha*(q-qold)*xi**2/2+&
                      &(0,1)*p+(0,1)*pold)*(q-qold)*a-(p-pold)**2/2))&
                      &/a/alpha/4.d0)
      end function 

      end module

