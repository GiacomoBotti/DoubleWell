!**********************************************************************!
! Module containing all the analytical integrals needed for the        !
! hamiltonian matrix, extracted from the MAPLE worksheet               !
! h2plus_param_gencode.mw                                              !
!                                                                      !
! It assumes 2g ewfn w/                                                !
! vecxi = [xi, -xi], xi <0                                             !
! veczeta = [zeta, -zeta], zeta >0                                     !
!                                                                      !
! Right now: used to debug numerical integrals                         !
! In the future: used instead of the numerical integrals               !
!**********************************************************************!

      module hamiltonian_module
      use constants
      use potential_module 

      implicit none

      integer*8 :: nd = 2
      real*8 :: lm = 1.d0 !light particle mass
      real*8 :: hm = 10.d0 !heavy particle mass
      private
      public :: T00maple,Tq0maple,Tp0maple
      public :: V00maple,Vq0maple,Vp0maple

      contains

!.....T00 matrix........................................................

      function T00maple(a,alpha,q,p,xi,zeta) result(Tout)
      !a,alpha: nucl,elec gau width
      !q,p: nucl gaus phase space center
      !xi,zeta: pos, mom scaling factor
        implicit none
        real *8, intent(in) :: a,alpha,q,p,xi,zeta 

        real*8 :: exponential
        complex*16, dimension(nd,nd) :: Tout

        exponential=-dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)

        Tout(1,1)=(2*hm*(p*zeta)**2+2*lm*p**2+hm*alpha+a*lm)&
                  &/(lm*hm*4.d0)
        Tout(1,2)=exponential*(2*hm*(alpha*q*xi)**2-2*lm*p**2-hm*alpha&
                  &-a*lm)/(lm*hm*4.d0)
        Tout(2,1)=Tout(1,2)
        Tout(2,2)=Tout(1,1)

      end function

!.....Tq0 matrix........................................................

      function Tq0maple(a,alpha,q,p,xi,zeta) result(Tqout)
      !a,alpha: nucl,elec gau width
      !q,p: nucl gaus phase space center
      !xi,zeta: pos, mom scaling factor
        implicit none
        real *8, intent(in) :: a,alpha,q,p,xi,zeta 

        real*8 :: exponential
        complex*16, dimension(nd,nd) :: Tqout

        exponential=-dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)

        Tqout(1,1)=(0.d0,1.d0)/2.d0*p*(((p*zeta)**2-alpha*xi*zeta+alpha&
                   &/2.d0)*hm+lm*(p**2+3.d0/2.d0*a))/(lm*hm)
        Tqout(1,2)=exponential*(-((q*alpha)**3*xi**4*hm)+(0,1)*zeta*p*& 
                  &(q*alpha)**2*xi**3*hm+((0,1)*q*alpha*p*hm+(lm*p**2)&
                  &+(3.d0/2.d0*hm*alpha)+(a*lm/2.d0))*q*alpha*(xi**2)-&
                  &(0,1)/2.d0*(2.d0*lm*p**2+hm*alpha+a*lm)*p*zeta*&
                  &xi-(0,1)/2.d0*p*(2.d0*lm*p**2+hm*alpha+3.d0*a*&
                  &lm))/(lm*hm*2.d0)
        Tqout(2,1)=Tqout(1,2)
        Tqout(2,2)=Tqout(1,1)

      end function

!.....Tp0 matrix........................................................

      function Tp0maple(a,alpha,q,p,xi,zeta) result(Tpout)
      !a,alpha: nucl,elec gau width
      !q,p: nucl gaus phase space center
      !xi,zeta: pos, mom scaling factor
        implicit none
        real*8, intent(in) :: a,alpha,q,p,xi,zeta 

        real*8 :: exponential
        complex*16, dimension(nd,nd) :: Tpout

        exponential=-dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)

        Tpout(1,1)=(2.d0*(0.d0,1.d0)*hm*p**2*q*xi*zeta**3+&
                   &2.d0*hm*p*zeta**2+(0.d0,1.d0)*((2.d0*p**2+a)*lm+&
                   &hm*alpha)*xi*q*zeta+2.d0*lm*p)/(hm*lm*4.d0)
        Tpout(1,2)=exponential*(hm*zeta*xi*q*(-p*q*xi*zeta+(0.d0,1.d0))&
                   &*alpha**2+p*(hm*zeta**2-2.d0*lm)*alpha/2.d0+&
                   &lm*p*zeta**2*(p**2+a/2.d0))/(alpha*hm*lm*2.d0)
        Tpout(2,1)=Tpout(1,2)
        Tpout(2,2)=Tpout(1,1)

      end function

!.....V00 matrix........................................................

      function V00maple(a,alpha,q,p,xi,zeta) result(V00out)
      !a,alpha: nucl,elec gau width
      !q,p: nucl gaus phase space center
      !xi,zeta: pos, mom scaling factor
        implicit none
        real*8, intent(in) :: a,alpha,q,p,xi,zeta 

        real*8 :: exponential
        complex*16, dimension(nd,nd) :: V00out

        exponential=-dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)
 
        V00out(1,1)=((a+(q*kappa)**2-kappa/2.d0)*V0*alpha*&
                    &dexp(-kappa*q**2)+((1.d0/8.d0+q**2*&
                    &(xi**2+1.d0/4.d0)*a)*alpha+a/2.d0)*k)/(alpha*a)
        V00out(1,2)=exponential*(-(a+(q*kappa)**2-kappa/2.d0)&
                    &*V0*alpha**2*dexp(-kappa*q**2)+&
                    &((-a*q**2/4.d0-1.d0/8.d0)*alpha**2-alpha*a/2.d0&
                    &+p**2*a*zeta**2)*k)/(a*alpha**2)
        V00out(2,1)=V00out(1,2)
        V00out(2,2)=V00out(1,1)
      end function

!.....Vq0 matrix........................................................

      function Vq0maple(a,alpha,q,p,xi,zeta) result(Vq0out)
      !a,alpha: nucl,elec gau width
      !q,p: nucl gaus phase space center
      !xi,zeta: pos, mom scaling factor
        implicit none
        real*8, intent(in) :: a,alpha,q,p,xi,zeta 

        real*8 :: exponential,exponential1
        complex*16, dimension(nd,nd) :: Vq0out

        exponential=-dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)
        exponential1=dexp((-(alpha*q*xi)**2-alpha*kappa*q**2-&
                     &(p*zeta)**2)/alpha)
 
        Vq0out(1,1)=(alpha*V0*(((0, 1)*p-q*kappa)*a+(0, 1)*&
                    &(kappa*q**2-(1.d0/2.d0))*p*kappa)*&
                    &dexp(-kappa*q**2)+((q*(xi**2+1.d0/4.d0)&
                    &*((0,1)*p*q+1.d0)*a+(0,1)*p/8.d0)*alpha+(0,1)&
                    &*p*a/2.d0)*k)/(a*alpha)
        Vq0out(1,2)=(alpha**2*V0*(-q*(a+(q*kappa)**2-kappa/2.d0)&
                    &*(xi**2)*alpha+((0,1)*(xi*zeta+1.d0)*p-q*kappa)&
                    &*a+(0,1)*(kappa*q**2-1.d0/2.d0)*p*kappa&
                    &*(xi*zeta+1.d0))*exponential1+exponential*&
                    &((a*q**2+(1.d0/2.d0))*q*(xi**2)*alpha**3/4.d0&
                    &+(-q*((0,1)*q*(xi*zeta+1.d0)*p-(2.d0*xi**2)+1.d0)&
                    &*a/4.d0-(0,1)*(xi*zeta+1.d0)*p/8.d0)*alpha**2&
                    &-(3.d0/2.d0)*a*p*((2.d0/3.d0)*p*q*(xi*zeta)**2&
                    &+(0,1)*xi*zeta+(0,1)/3.d0)*alpha+(0,1)*&
                    &(xi*zeta+1.d0)*a*p**3*zeta**2)*k)/(a*alpha** 2)
        Vq0out(2,1)=Vq0out(1,2)
        Vq0out(2,2)=Vq0out(1,1)

      end function

!.....Vp0 matrix........................................................

      function Vp0maple(a,alpha,q,p,xi,zeta) result(Vp0out)
      !a,alpha: nucl,elec gau width
      !q,p: nucl gaus phase space center
      !xi,zeta: pos, mom scaling factor
        implicit none
        real*8, intent(in) :: a,alpha,q,p,xi,zeta 

        real*8 :: exponential,exponential1
        complex*16, dimension(nd,nd) :: Vp0out

        exponential=-dexp((-(alpha*q*xi)**2-(p*zeta)**2)/alpha)
        exponential1=dexp((-(alpha*q*xi)**2-alpha*kappa*q**2-&
                     &(p*zeta)**2)/alpha)
 
        Vp0out(1,1)=(0,1)*((xi*((kappa*q)**2+a-kappa/2.d0)*zeta+kappa)&
                    &*V0*alpha*dexp(-kappa*q**2)+k*((-1.d0/4.d0+xi*&
                    &((q*xi)**2*a+a*q**2/4.d0+1.d0/8.d0)*zeta)*alpha&
                    &+3.d0/2.d0*zeta*xi*a))*q/(alpha*a)        
        Vp0out(1,2)=(8.d0*((0,1)*alpha*q*kappa-(((kappa*q)**2&
                    &+a-kappa/2.d0)*zeta**2*p))*(alpha**2)*V0*&
                    &exponential1+exponential*2*k*((0, 1)*q*(alpha**3)&
                    &+(p*zeta**2*(a*q**2+1.d0/2.d0)*alpha**2)+&
                    &(6.d0*zeta**2*alpha*p*a)-(4.d0*zeta**4*p**3*a)))&
                    &/(a*(alpha**3)*8)
        Vp0out(2,1)=Vp0out(1,2)
        Vp0out(2,2)=Vp0out(1,1)

      end function

      end module
