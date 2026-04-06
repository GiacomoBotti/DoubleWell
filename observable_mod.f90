!*********************************************************************!
! Module containing the observables functions                         !
!*********************************************************************!

       module observable_module

       use constants
       use integrals_module
       use kinetic_module
       use effectivepot_module
       use basisset_module

       implicit none

       private
       public :: energy,plot_wfn

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

!......Plot wavefunction................................................

       subroutine plot_wfn(nd,qtot,ptot,cvec,tildeBmat,uuunit)
       ! nd: bath dimension
       ! qtot: total position vector
       ! ptot: total momentum vector
       ! cvec: basis set coefficients vector
       ! tildeBmat: total complex gaussian width
       ! uuunit: unit in which to print => fort.unit 
         integer, intent(in) :: nd,uuunit
         real*8, dimension(nd+1), intent(in) :: qtot,ptot
         complex*16, dimension(nh), intent(in) :: cvec
         complex*16, dimension(nd+1,nd+1), intent(in) :: tildeBmat

         integer :: i,j
         real*8 :: Nsq,N,phase,pol,herm
         real*8, dimension(nd+1) :: rvec,xvec
         complex*16 :: sqr,img,psi
         complex*16, dimension(nd+1) :: sqrvec

         write(*,*) "Plotting gaussian wfn on the x=y cut"
   
         !Nsq = fun_Nsq(nd,real(tildeBmat))
         !N = dsqrt(Nsq)

         phase = datan((aimag(cvec(1))/real(cvec(1))))

         xvec(:) = -5.d0
         do i = 1,100
            xvec(:) = xvec(:)+0.1d0 
            rvec(:) = xvec(:) - qtot(:)
            sqrvec = matmul(tildeBmat,rvec)
            sqr = dot_product(rvec,sqrvec)
            img = iu*dot_product(ptot,rvec)
            pol = 0.d0
            do j = 1,nh
              herm = herm_pol(j,xvec(1),qtot(1),real(tildeBmat(1,1)))
              pol = pol + cvec(j)*herm
            end do
            psi = zexp(-0.5d0*sqr + img)*pol
            write(uuunit,*) xvec, real(psi), aimag(psi), &
                & real(psi*dconjg(psi))
         end do

       end subroutine

       end module
