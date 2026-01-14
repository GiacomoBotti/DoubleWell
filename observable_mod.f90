!*********************************************************************!
! Module containing the observables functions                         !
!*********************************************************************!

       module observable_module

       use kinetic_module
       use effectivepot_module

       implicit none

       private
       public :: energynum

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
         real*8, dimension(nh), intent(in) :: cvec
         complex*16, dimension(nd+1,nd+1), intent(in) :: tildeBmat

         real*8 :: V0,Hout
         complex*16, dimension(nh) :: c, Hc
         complex*16, dimension(nh,nh) :: H00M,T00M

         T00M = kin_energy(nd,q,p,qvec,pvec,tildeBmat)
         V0 = fun_V0(nd,qtot,cvec,real(tildeBmat)) 
         H00M = T00M 
   
         Hc = matmul(H00M,c)
         Hout = dreal(dot_product(c,Hc)) + V0
         
       end function
       end module
