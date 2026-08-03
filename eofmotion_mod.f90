!**********************************************************************!
! Module that contains all things pertaining the evolution of the      !
! Gaussian wavepacket following Coalson & Karplus 1990                 !
!**********************************************************************!

      module eofmotion_module

      use constants
      use parameters_module
      use potential_module
      use basisset_module
      use effectivepot_module
      use integrals_module
      use matrix_module
      use inversion_module
      !use observable_module 

      implicit none

      integer :: i
      real*8, dimension(nv+1,nv+1), public :: invMassMat

      private
      public :: MassesMat,KarplusTimeDer,rungekutta,scprop,vtvprop 
      public :: scprop_der,pece_param,scpece_param

      contains

!.....Matrix of the masses..............................................

      subroutine MassesMat(masses)
      ! masses: vector of the masses
      ! invMassMat: Matrix of the inverted masses
       real*8, dimension(nv+1), intent(in) :: masses  

       invMassMat(:,:) = 0.d0

       do i = 1,nv+1
          invMassMat(i,i) = 1.d0/masses(i)
!          write(*,*) invMassMat(i,:)
       end do

       write(*,*) "Matrix of inverse masses generated"

      end subroutine

!.....C&K equations of motion...........................................

      subroutine KarplusTimeDer(nd,cvec,qtot,ptot,Bcmplx,dotq,dotp,dotB)
      ! nd: dimension of bath
      ! cvec: vector of the basis coefficients
      ! qtot: gaussian center position vector (x&y)
      ! ptot: gaussian center momenta vector (x&y)
      ! Bmat: total gaussian width (x&y)
      ! dotq: time-derivative of qtot
      ! dotp: time-derivative of ptot
      ! dotB: time-derivative of Bmat
       integer, intent(in) :: nd
       complex*16, dimension(nh), intent(in) :: cvec
       real*8, dimension(nd+1), intent(in) :: qtot,ptot
       complex*16, dimension(nd+1,nd+1), intent(in) :: Bcmplx
     
       real*8, dimension(nd+1), intent(out) :: dotq, dotp
       complex*16, dimension(nd+1,nd+1), intent(out) :: dotB

       real*8, dimension(nd+1) :: V1
       real*8, dimension(nd+1,nd+1) :: V2,Bmat
       real*8 :: a,Y0
       real*8, dimension(nd) :: avec
       real*8, dimension(nd,nd) :: Amat,Tmat,LambdaMat
       complex*16, dimension(nh) :: cwork
       real*8, dimension(nh,nh) :: X3,X2,X1,X0
       complex*16, dimension(nd+1,nd+1) :: MB,prova
 
       ! Allows for Gaussian Average only
       ! or average only on coalmode basis function
       cwork(:) = coalvec(:)*cvec(:) 
       cwork(coalmode) = cwork(coalmode) + coalc

       Bmat = dreal(Bcmplx)
       call extractA(nd,Bmat,Amat,avec,a)
       call diagonalization(nd,Amat,LambdaMat,Tmat)
 
       dotq = matmul(invMassMat,ptot)

       ! Integrals
       Y0=int_Y0(nd,LambdaMat)
       X3=int_XnMat(nd,3,a,avec,Amat,qtot(1))
       X2=int_XnMat(nd,2,a,avec,Amat,qtot(1))
       X1=int_XnMat(nd,1,a,avec,Amat,qtot(1))
       X0=int_XnMat(nd,0,a,avec,Amat,qtot(1))

       V1 = fun_V1(nd,qtot,cwork,Bmat,Y0,X3,X2,X1,X0)
       dotp = - V1

       V2 = fun_V2(nd,qtot,cwork,Bmat,Y0,X2,X0)
       MB=matmul(invMassMat,Bcmplx)              
       ! confirmed with debug_karplus_width.mw
       dotB = -(0.d0,1.d0)*matmul(Bcmplx,MB) + (0.d0,1.d0)*V2

       dotq(:) = scalvec(:)*dotq(:)
       dotp(:) = scalvec(:)*dotp(:)
       dotB(:,:) = scalmat(:,:)*dotB(:,:)

! TEST TEST TEST TEST TEST TEST TEST
       
       !dotB(1,2) =complex(0.d0,aimag(dotB(1,2)))
       !dotB(2,1) =complex(0.d0,aimag(dotB(2,1)))

!       write(*,*) dotq
!       write(*,*) dotp
!       write(*,*) dotB(1,1), dotB(1,2), dotB(2,2)

      end subroutine

!.....Runge-Kutta 4 propagator..........................................

      subroutine rungekutta(nd,h,cj,qj,pj,Bj)
      ! Does one Runge-Kutta 4 step
       integer, intent(in) :: nd
       real*8, intent(in) :: h
       complex*16, dimension(nh), intent(in) :: cj
       real*8, dimension(nd+1) :: qj,pj,qi,ppi
       complex*16, dimension(nd+1,nd+1) :: Bj,Bi

       real*8,dimension(4) :: hvec = [0.5d0,0.5d0,1.d0,0.d0]

       real*8,dimension(nd+1,4) :: kq,kp
       complex*16,dimension(nd+1,nd+1,4) :: kb

       hvec = hvec*h
       qi = qj
       ppi = pj
       Bi = Bj

       do i = 1,4
         call KarplusTimeDer(nd,cj,qi,ppi,Bi,&
              &kq(:,i),kp(:,i),kb(:,:,i))             
         qi = qj + hvec(i)*kq(:,i)
         ppi = pj + hvec(i)*kp(:,i)
         Bi = Bj + hvec(i)*kb(:,:,i)
       end do
       qj = qj&
          &+h*(kq(:,1)+2.d0*kq(:,2)+2.d0*kq(:,3)+kq(:,4))/6.d0
       pj = pj&
          &+h*(kp(:,1)+2.d0*kp(:,2)+2.d0*kp(:,3)+kp(:,4))/6.d0
       Bj = Bj&
        &+h*(kb(:,:,1)+2.d0*kb(:,:,2)+2.d0*kb(:,:,3)+kb(:,:,4))/6.d0

      end subroutine

!.....Self-consistent propagator........................................

      subroutine scprop(nd,h,cj,qj,pj,Bj)
      ! Does one self-consistent progator step
       integer, intent(in) :: nd
       real*8, intent(in) :: h
       complex*16, dimension(nh), intent(in) :: cj
       real*8, dimension(nd+1) :: qj,pj,qi,ppi,qiold,piold,qav,pav
       real*8, dimension(nd+1) :: sqq,sqp 
       complex*16, dimension(nd+1,nd+1) :: Bj,Bi,Biold,Bav
       real*8, dimension(nd+1,nd+1) :: sqB

       integer*8 :: maxcycle=20  
       real*8 :: thr = 1.d-5
       real*8 :: error
       real*8,dimension(nd+1,4) :: kq,kp
       complex*16,dimension(nd+1,nd+1,4) :: kb

       qi = qj
       ppi = pj
       Bi = Bj

       do i = 1,maxcycle
         qiold = qi
         piold = ppi
         Biold = Bi
          
         qav = (qi+qj)/2 
         pav = (ppi+pj)/2 
         Bav = (Bi+Bj)/2 

         call KarplusTimeDer(nd,cj,qav,pav,Bav,&
              &kq(:,1),kp(:,1),kb(:,:,1))             

         qi = qj + h*kq(:,1)
         ppi = pj + h*kp(:,1)
         Bi = Bj + h*kb(:,:,1)
  
         sqq = (qiold-qi)**2
         sqp = (piold-ppi)**2
         sqB = (Biold-Bi)*conjg(Biold-Bi)

         error = dsqrt(sum(sqq) + sum(sqp) + sum(sqB))

         if(error.le.thr) then
           exit
         end if
       end do
      
       qj=qi
       pj=ppi
       Bj=Bi

      end subroutine

!.....VTV propagator...................................................

      subroutine vtvprop(nd,h,cj,qj,pj,Bj)
      ! computes one step of VTV integrator, see J.L.Vanicek2023
       integer, intent(in) :: nd
       real*8, intent(in) :: h
       complex*16, dimension(nh), intent(in) :: cj
       real*8, dimension(nd+1) :: qj,pj,qi,ppi
       complex*16, dimension(nd+1,nd+1) :: Bj,Bi,invB,invBi

       real*8, dimension(nd+1) :: V1,dotq
       real*8, dimension(nd+1,nd+1) :: V2,Bmat
       real*8 :: a,Y0
       real*8, dimension(nd) :: avec
       real*8, dimension(nd,nd) :: Amat,Tmat,LambdaMat
       complex*16, dimension(nh) :: cwork
       real*8, dimension(nh,nh) :: X3,X2,X1,X0

       ! Allows for Gaussian Average only
       !cwork(:) = coalvec(:)*cj(:) 
       !cwork(1) = cwork(1) + coalc
       cwork(:) = coalvec(:)*cj(:) 
       cwork(coalmode) = cwork(coalmode) + coalc
       !write(*,*) cwork
       !Half V step 
       call extractA(nd,dreal(Bj),Amat,avec,a)
       call diagonalization(nd,Amat,LambdaMat,Tmat)
       Y0=int_Y0(nd,LambdaMat)
       X3=int_XnMat(nd,3,a,avec,Amat,qj(1))
       X2=int_XnMat(nd,2,a,avec,Amat,qj(1))
       X1=int_XnMat(nd,1,a,avec,Amat,qj(1))
       X0=int_XnMat(nd,0,a,avec,Amat,qj(1))
       ! KARPLUS
       V1 = fun_V1(nd,qj,cwork,dreal(Bj),Y0,X3,X2,X1,X0)
       V2 = fun_V2(nd,qj,cwork,dreal(Bj),Y0,X2,X0)
       ! HELLER
       !V1 =qj(1)**3/(4.d0*eta_const)+sigma_const*qj(1)+gamma_const*qj(2)
       !V2 = 3.d0*qj(1)**2/(4.d0*eta_const)+sigma_const*qj(1)
       ppi(:) = pj(:) - 0.5d0*h*V1(:)!*scalvec(:)
       Bi = Bj +0.5d0*iu*h*V2*(1.d0-ffact)
       !Bi(:,:) = Bj(:,:) +0.5d0*iu*h*V2(:,:)*scalmat(:,:)
       !Full T step
       dotq = matmul(invMassMat,ppi)
       qi(:) = qj(:) + h*dotq(:)!*scalvec(:) 
       invBi=invgen(nd+1,Bi)
       !invBi=Bi
       invB(:,:) = invBi(:,:) + iu*h*invMassMat(:,:)*(1.d0 -ffact)
       !invB(:,:) = invBi(:,:) + iu*h*invMassMat(:,:)*scalmat(:,:)
       Bi = invgen(nd+1,invB) 
       !Half V step
       call extractA(nd,dreal(Bi),Amat,avec,a)
       call diagonalization(nd,Amat,LambdaMat,Tmat)
       Y0=int_Y0(nd,LambdaMat)
       X3=int_XnMat(nd,3,a,avec,Amat,qi(1))
       X2=int_XnMat(nd,2,a,avec,Amat,qi(1))
       X1=int_XnMat(nd,1,a,avec,Amat,qi(1))
       X0=int_XnMat(nd,0,a,avec,Amat,qi(1))
       ! KARPLUS
       V1 = fun_V1(nd,qi,cwork,dreal(Bi),Y0,X3,X2,X1,X0)
       V2 = fun_V2(nd,qi,cwork,dreal(Bi),Y0,X2,X0)
       ! HELLER
       !V1 =qj(1)**3/(4.d0*eta_const)+sigma_const*qj(1)+gamma_const*qj(2)
       !V2 = 3.d0*qj(1)**2/(4.d0*eta_const)+sigma_const*qj(1)
       ppi(:) = ppi(:) - 0.5d0*h*V1(:)!*scalvec(:)
       Bi = Bi +0.5d0*iu*h*V2*(1.d0 -ffact)
       !Bi(:,:) = Bi(:,:) +0.5d0*iu*h*V2(:,:)*scalmat(:,:)
       ! Finish
       qj = qi
       pj = ppi
       Bj = Bi
       
      end subroutine

!.....Self-consistent propagator with derivatives convergence...........

      subroutine scprop_der(nd,h,cj,qj,pj,Bj)
      ! Does one self-consistent progator step with derivatives conv
       integer, intent(in) :: nd
       real*8, intent(in) :: h
       complex*16, dimension(nh), intent(in) :: cj
       real*8, dimension(nd+1) :: qj,pj,qi,ppi,qiold,piold,qav,pav
       real*8, dimension(nd+1) :: sqq,sqp 
       complex*16, dimension(nd+1,nd+1) :: Bj,Bi,Biold,Bav
       real*8, dimension(nd+1,nd+1) :: sqB

       integer*8 :: maxcycle=20  
       real*8 :: thr = 1.d-5
       real*8 :: error
       real*8,dimension(nd+1,4) :: kq,kp
       complex*16,dimension(nd+1,nd+1,4) :: kb

       qi = qj
       ppi = pj
       Bi = Bj

       call KarplusTimeDer(nd,cj,qj,pj,Bj,&
            &kq(:,1),kp(:,1),kb(:,:,1))             

       kq(:,2) = kq(:,1)
       kp(:,2) = kp(:,1)
       kb(:,:,2) = kb(:,:,1)

       do i = 1,maxcycle

        ! qi = qj + h*(kq(:,2)+kq(:,1))*0.5d0
        ! ppi = pj + h*(kp(:,2)+kp(:,1))*0.5d0
        ! Bi = Bj + h*(kb(:,:,2)+kb(:,:,1))*0.5d0
         ! derivatives at midpoint 
         qi = qj + h*(kq(:,2))*0.5d0
         ppi = pj + h*(kp(:,2))*0.5d0
         Bi = Bj + h*(kb(:,:,2))*0.5d0

         call KarplusTimeDer(nd,cj,qi,ppi,Bi,&
              &kq(:,3),kp(:,3),kb(:,:,3))             
 
         !Convergence on parameters 
         !sqq = (qiold-qi)**2
         !sqp = (piold-ppi)**2
         !sqB = (Biold-Bi)*conjg(Biold-Bi)
         !Convergence on derivatives
         sqq = (kq(:,2)-kq(:,3))**2
         sqp = (kp(:,2)-kp(:,3))**2
         sqB = (kb(:,:,2)-kb(:,:,3))*conjg(kb(:,:,2)-kb(:,:,3))

         error = dsqrt(sum(sqq) + sum(sqp) + sum(sqB))
 
         !write(*,*) i
         kq(:,2) = kq(:,3)
         kp(:,2) = kp(:,3)
         kb(:,:,2) = kb(:,:,3)

         if(error.le.thr) then
           qi = qj + h*kq(:,2)
           ppi = pj + h*kp(:,2)
           Bi = Bj + h*kb(:,:,2)
           !qi = qj + h*(kq(:,2)+kq(:,1))*0.5d0
           !ppi = pj + h*(kp(:,2)+kp(:,1))*0.5d0
           !Bi = Bj + h*(kb(:,:,2)+kb(:,:,1))*0.5d0
          ! write(*,*) "I get out at", i
           exit
         end if
       end do
      
       qj=qi
       pj=ppi
       Bj=Bi

      end subroutine

!.....PECE propagator for the parameters only..........................

      subroutine pece_param(nd,h,cj,qj,pj,Bj)
      ! Does one PECE progator step for the parameters
       integer, intent(in) :: nd
       real*8, intent(in) :: h
       complex*16, dimension(nh), intent(in) :: cj
       real*8, dimension(nd+1) :: qj,pj,qi,ppi,qiold,piold,qav,pav
       real*8, dimension(nd+1) :: sqq,sqp 
       complex*16, dimension(nd+1,nd+1) :: Bj,Bi,Biold,Bav
       real*8, dimension(nd+1,nd+1) :: sqB

       integer*8 :: maxcycle=20  
       real*8 :: thr = 1.d-5
       real*8 :: error
       real*8,dimension(nd+1,4) :: kq,kp
       complex*16,dimension(nd+1,nd+1,4) :: kb

       qi = qj
       ppi = pj
       Bi = Bj

       call KarplusTimeDer(nd,cj,qj,pj,Bj,&
            &kq(:,1),kp(:,1),kb(:,:,1))             

       ! Predictor (Lambda tilde) 
       qi = qj + h*(kq(:,1))
       ppi = pj + h*(kp(:,1))
       Bi = Bj + h*(kb(:,:,1))

       call KarplusTimeDer(nd,cj,qi,ppi,Bi,&
            &kq(:,2),kp(:,2),kb(:,:,2))             

       ! Corrector (Lambda hat)
       qi = qj + h*(kq(:,1)+kq(:,2))*0.5d0
       ppi = pj + h*(kp(:,1)+kp(:,2))*0.5d0
       Bi = Bj + h*(kb(:,:,1)+kb(:,:,2))*0.5d0

       call KarplusTimeDer(nd,cj,qi,ppi,Bi,&
            &kq(:,3),kp(:,3),kb(:,:,3))             

       ! Evaluator (Lambda t)
       qi = qj + h*(kq(:,1)+kq(:,3))*0.5d0
       ppi = pj + h*(kp(:,1)+kp(:,3))*0.5d0
       Bi = Bj + h*(kb(:,:,1)+kb(:,:,3))*0.5d0

       qj=qi
       pj=ppi
       Bj=Bi

      end subroutine

!.....SC PECE propagator for the parameters only........................

      subroutine scpece_param(nd,h,cj,qj,pj,Bj)
      ! Does a self-consisten PECE progator step for the parameters
       integer, intent(in) :: nd
       real*8, intent(in) :: h
       complex*16, dimension(nh), intent(in) :: cj
       real*8, dimension(nd+1) :: qj,pj,qi,ppi,qiold,piold,qav,pav
       real*8, dimension(nd+1) :: sqq,sqp 
       complex*16, dimension(nd+1,nd+1) :: Bj,Bi,Biold,Bav
       real*8, dimension(nd+1,nd+1) :: sqB

       integer :: k
       integer*8 :: maxcycle=20  
       real*8 :: thr = 1.d-10
       real*8 :: error
       real*8,dimension(nd+1,4) :: kq,kp
       complex*16,dimension(nd+1,nd+1,4) :: kb

       qi = qj
       ppi = pj
       Bi = Bj

       call KarplusTimeDer(nd,cj,qj,pj,Bj,&
            &kq(:,1),kp(:,1),kb(:,:,1))             

       ! Predictor (Lambda tilde) 
       qi = qj + h*(kq(:,1))
       ppi = pj + h*(kp(:,1))
       Bi = Bj + h*(kb(:,:,1))

       kq(:,4) = kq(:,1)
       kp(:,4) = kp(:,1)
       kb(:,:,4) = kb(:,:,1)

       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
       ! Corrector - Evaluator SC cycle !
       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

       do k = 1,maxcycle
       
          call KarplusTimeDer(nd,cj,qi,ppi,Bi,&
               &kq(:,2),kp(:,2),kb(:,:,2))             

          ! Corrector (Lambda hat)
          qi = qj + h*(kq(:,1)+kq(:,2))*0.5d0
          ppi = pj + h*(kp(:,1)+kp(:,2))*0.5d0
          Bi = Bj + h*(kb(:,:,1)+kb(:,:,2))*0.5d0

          call KarplusTimeDer(nd,cj,qi,ppi,Bi,&
               &kq(:,3),kp(:,3),kb(:,:,3))             

          ! Evaluator (Lambda t)
          qi = qj + h*(kq(:,1)+kq(:,3))*0.5d0
          ppi = pj + h*(kp(:,1)+kp(:,3))*0.5d0
          Bi = Bj + h*(kb(:,:,1)+kb(:,:,3))*0.5d0

          error = abs(sum(kq(:,4)-kq(:,3))) &
                & + abs(sum(kp(:,4)-kp(:,3))) &
                & + abs(sum(kb(:,:,4)-kb(:,:,3))) 

          if(error.le.thr) then
!             write(2345,*) "I exit at ", k
             exit
          end if

          kq(:,4) = kq(:,3)
          kp(:,4) = kp(:,3)
          kb(:,:,4) = kb(:,:,3)

       end do

       qj=qi
       pj=ppi
       Bj=Bi

      end subroutine

      end module
