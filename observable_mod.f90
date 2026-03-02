!*********************************************************************!
! Module containing the observables functions                         !
!*********************************************************************!

       module observable_module

       use kinetic_module
       use effectivepot_module
       use basisset_module

       implicit none

       private
       public :: energy

       contains

!......Analytical Hamiltonian...........................................

       function energy(nd,qtot,ptot,cvec,tildeBmat) result(Hout)
       ! nd: bath dimension
       ! qtot: total position vector
       ! ptot: total momentum vector
       ! cvec: basis set coefficients vector
       ! tildeBmat: total complex gaussian width
         integer, intent(in) :: nd
         real*8, dimension(nd+1), intent(in) :: qtot,ptot
         complex*16, dimension(nh), intent(in) :: cvec
         complex*16, dimension(nd+1,nd+1), intent(in) :: tildeBmat

         real*8 :: V0,T0,Hout,q,p
         real*8, dimension(nd) :: qvec,pvec
         complex*16, dimension(nh) :: Tc,Vc
         complex*16, dimension(nh,nh) :: H00M,T00M,V00M

         q=qtot(1)
         p=ptot(1)
         qvec=qtot(2:nd+1)
         pvec=ptot(2:nd+1)

         T00M = kin_energy(nd,q,p,qvec,pvec,tildeBmat)
         V00M = fun_V0(nd,qtot,cvec,real(tildeBmat)) 
   
         Tc = matmul(T00M,cvec)
         T0 = dreal(dot_product(cvec,Tc))

         Vc = matmul(V00M,cvec)
         V0 = dreal(dot_product(cvec,Vc))

         Hout = T0 + V0

         write(324,*) Hout, T0, V0
         
       end function
       end module
