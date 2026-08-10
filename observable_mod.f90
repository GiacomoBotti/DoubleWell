!*********************************************************************!
! Module containing the observables functions                         !
!*********************************************************************!

       module observable_module

       use constants
       use parameters_module
       use integrals_module
       use kinetic_module
       use effectivepot_module
       use basisset_module
       use matrix_module
       use inversion_module

       implicit none

       private
       public :: energy,plot_wfn

       contains

!......Analytical Hamiltonian...........................................

       subroutine energy(nd,qtot,ptot,cvec,tildeBmat,H00M,Hout,&
                  &Mx,Mx2,My)
       ! nd: bath dimension
       ! qtot: total position vector
       ! ptot: total momentum vector
       ! cvec: basis set coefficients vector
       ! tildeBmat: total complex gaussian width
         integer, intent(in) :: nd
         real*8, dimension(nd+1), intent(in) :: qtot,ptot
         complex*16, dimension(nh), intent(in) :: cvec
         complex*16, dimension(nd+1,nd+1), intent(in) :: tildeBmat

         real*8, intent(out) :: Hout,Mx,Mx2
         real*8, dimension(nd) :: My
         complex*16, dimension(nh,nh), intent(out) :: H00M

         integer :: i
         real*8 :: V0,T0,q,p,Nsq,cX1c,cX0c,cX2c
         real*8, dimension(nd) :: qvec,pvec,Aa
         complex*16, dimension(nh) :: Tc,Vc,X1c,X0c,X2c
         complex*16, dimension(nh,nh) :: T00M,V00M
         real*8 :: a,Y0
         real*8, dimension(nd) :: avec
         real*8, dimension(nd,nd) :: Amat,Tmat,LambdaMat,invA
         real*8, dimension(nh,nh) :: X4,X3,X2,X1,X0

         q=qtot(1)
         p=ptot(1)
         qvec=qtot(2:nd+1)
         pvec=ptot(2:nd+1)

         call extractA(nd,real(tildeBmat),Amat,avec,a)
         call diagonalization(nd,Amat,LambdaMat,Tmat)
         invA = invgen_real(nd,Amat)
         Nsq=fun_Nsq(nd+1,real(tildeBmat))
      
         ! Integrals
         Y0=int_Y0(nd,LambdaMat)

         !write(324,*) Y0 

         X4=int_XnMat(nd,4,a,avec,Amat,qtot(1))
         X3=int_XnMat(nd,3,a,avec,Amat,qtot(1))
         X2=int_XnMat(nd,2,a,avec,Amat,qtot(1))
         X1=int_XnMat(nd,1,a,avec,Amat,qtot(1))
         X0=int_XnMat(nd,0,a,avec,Amat,qtot(1))
         !write(324,*) "X4mat out"
         !write(324,*) X4(1,1)

         T00M = kin_energy(nd,q,p,qvec,pvec,tildeBmat,Y0,X2,X1,X0)
         V00M = fun_V0(nd,qtot,cvec,real(tildeBmat),Y0,X4,X2,X1,X0) 
         H00M = T00M + V00M
   
         Tc = matmul(T00M,cvec)
         T0 = dreal(dot_product(cvec,Tc))

         Vc = matmul(V00M,cvec)
         V0 = dreal(dot_product(cvec,Vc))

         Hout = T0 + V0

!         write(324,*) "+++++++++++++++++++++++++++++"
         write(324,*) Hout, T0, V0, dreal(H00M(1,1))

!         write(324,*) "+++++++++++++++++++++++++++++"
!         do i = 1, nh
!           write(324,*) H00M(i,:) 
!         end do 
!         write(324,*) "+++++++++++++++++++++++++++++"
!         do i = 1, nh
!           write(324,*) T00M(i,:)
!         end do 
!         write(324,*) "+++++++++++++++++++++++++++++"
!         do i = 1, nh
!           write(324,*) V00M(i,:)
!         end do 
         
!........Momenta

         X0c = matmul(X0,cvec)
         cX0c = dot_product(cvec,X0c)

         X1c = matmul(X1,cvec)
         cX1c = dot_product(cvec,X1c)

         X2c = matmul(X2,cvec)
         cX2c = dot_product(cvec,X2c)

         Aa = matmul(invA,avec)

         Mx = Nsq*Y0*cX1c
         Mx2 = Nsq*Y0*cX2c

         My(:) = Nsq*Y0*(cX0c*qvec(:)+q*cX0c*Aa(:)-cX1c*Aa(:))


       end subroutine 

!......Plot wavefunction................................................

       subroutine plot_wfn(nd,qtot,ptot,cvec,tildeBmat,norm)
       ! nd: bath dimension
       ! qtot: total position vector
       ! ptot: total momentum vector
       ! cvec: basis set coefficients vector
       ! tildeBmat: total complex gaussian width
         integer, intent(in) :: nd
         real*8, intent(in) :: norm
         real*8, dimension(nd+1), intent(in) :: qtot,ptot
         complex*16, dimension(nh), intent(in) :: cvec
         complex*16, dimension(nd+1,nd+1), intent(in) :: tildeBmat

         integer :: i,j
         real*8 :: Nsq,N,phase,herm,step
         real*8, dimension(nd+1) :: rvec,xvec
         complex*16 :: sqr,img,psi,pol
         complex*16, dimension(nd+1) :: sqrvec

         !write(*,*) "Plotting gaussian wfn on the x=y cut"
   

         !phase = datan((aimag(cvec(1))/real(cvec(1))))
         phase = 0.d0
         Nsq = fun_Nsq(nd+1,real(tildeBmat))
         N = dsqrt(Nsq)

        ! write(*,*) "Plotting gaussian wfn on the y=0 cut"
         xvec(:) = qtot(:) 
         write(111,*) "#", qtot(:)
         write(111,*) "#", tildeBmat(1,1),tildeBmat(1,2),tildeBmat(2,2) 
         write(111,*) "#", cvec(:) 
         rvec(:) = 0.d0
         xvec(1) = -10.d0  
         step = -2*xvec(1)/200
        ! xvec(:) = -5.d0 
         do i = 1,200
            xvec(1) = xvec(1)+step 
            rvec(:) = xvec(:) - qtot(:)
            sqrvec = matmul(tildeBmat,rvec)
            sqr = dot_product(rvec,sqrvec)
            img = iu*dot_product(ptot,rvec)
            pol = 0.d0
            do j = 1,nh
              herm = herm_pol(j,xvec(1),qtot(1),real(tildeBmat(1,1)))
              pol = pol + cvec(j)*herm
            end do
            psi = zexp(-0.5d0*sqr + img)*pol*N!*dsqrt(norm)
            write(111,*) xvec(1), real(psi), aimag(psi), &
                & real(psi*dconjg(psi))
         end do

         xvec(:) = qtot(:) 
         write(222,*) "#", qtot(:)
         write(222,*) "#", tildeBmat(1,1),tildeBmat(1,2),tildeBmat(2,2) 
         write(222,*) "#", cvec(:) 
         rvec(:) = 0.d0
         xvec(2) = -10.d0 
         step = -2*xvec(2)/200
         do i = 1,200
            xvec(2) = xvec(2)+step 
            !xvec(:) = xvec(:)+0.1d0 
            rvec(:) = xvec(:) - qtot(:)
            sqrvec = matmul(tildeBmat,rvec)
            sqr = dot_product(rvec,sqrvec)
            img = iu*dot_product(ptot,rvec)
            pol = 0.d0
            do j = 1,nh
              herm = herm_pol(j,xvec(1),qtot(1),real(tildeBmat(1,1)))
              pol = pol + cvec(j)*herm
            end do
            psi = zexp(-0.5d0*sqr + img)*pol*N!*dsqrt(norm)
            write(222,*) xvec(2), real(psi), aimag(psi), &
                & real(psi*dconjg(psi))
         end do

         write(111,*) " " 
         write(111,*) " " 
         write(222,*) " " 
         write(222,*) " " 

       end subroutine

       end module
