!**********************************************************************!
! Module containing the Runge-Kutta evolution of the parameters and    !
! coefficients                                                         !
!**********************************************************************!

       module evolution_module

       use constants
       use parameters_module
       use basisset_module
       use normalization_module
       use observable_module
       use BOT_module
       use eofmotion_module
       use integrals_module
       use inversion_module
       use matrix_module
       use effectivepot_module

       implicit none

       private
       public :: bot_evo

       contains

!......BOT evolution....................................................

       subroutine bot_evo(nd,trj,q0,p0,c0,Bcmplx,qeq,peq,ceq,Beq)
       ! nd: bath dimension
       ! trj : trajectory parameters (first step, last step, nstep)
       ! q0 : initial gaussian center (x&y)
       ! p0 : initial gaussian momentum (x&y)
       ! c0 : initial coefficients 
       ! Bcmplx : initial gaussian width matrix, real & imaginary
       integer, intent(in) :: nd
       integer*8, dimension(4), intent(in) :: trj
       real*8, dimension(nd+1), intent(in) :: q0,p0,qeq,peq
       complex*16,dimension(nh), intent(in) :: c0,ceq 
       complex*16,dimension(nd+1,nd+1), intent(in) :: Bcmplx,Beq 

       integer*8 :: i,j,first,last,nstep,k,nprint
       real*8 :: h,time,N,E,q,p,E0,Nsq,Neq,Mx
       real*8,dimension(nh) :: csq,phase 
       real*8,dimension(nd) :: qvec,pvec,My 
       real*8,dimension(nd+1) :: qtoti,ptoti,qtotj,ptotj,qold,pold,qeqX 
       complex*16,dimension(nd+1,nd+1) :: Bcmplxi,Bcmplxj,Bold
       
       complex*16,dimension(nh) :: cvec,cj,ctemp,ceqN

       real*8,dimension(nd+1,4) :: kq,kp
       real*8,dimension(4) :: hvec = [0.5d0,0.5d0,1.d0,0.d0]
       complex*16,dimension(nd+1,nd+1,4) :: kb

       real*8 :: a,Y0,cPc
       real*8, dimension(nd) :: avec
       real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
       real*8, dimension(nd+1,nd+1) :: Bmat
       real*8, dimension(nh,nh) :: S00M,invS,X0Mat
       complex*16, dimension(nh,nh) :: H00M

       complex*16 :: cTau0c,cTauXc,ckg
       complex*16, dimension(nh) :: Tau0c,Pc,TauXc
       complex*16, dimension(nh,nh) :: Tau0,Prob,TauX

       integer :: bar_width,pos
       real :: frac
       character(len=50) :: bar

       bar_width = 50
       hvec = hvec*h

       ! trajectory parameters
       first = trj(1)
       last = trj(2)
       nstep = trj(3)

       ! Separate x & y
       q = q0(1)
       qvec = q0(2:nd+1)
       p = p0(1)
       pvec = p0(2:nd+1)

      write(*,*) "+---------------------------------------------------+"
       write(*,*) "Initial conditions"
       write(*,*) "q          ","p"
       write(*,*) q,p
       do i = 1,nd
         write(*,*) qvec(i), pvec(i)
       end do

       open(unit=321,file="trajectory_BOT.dat",status="unknown")
       open(unit=322,file="coefficients_BOT.dat",status="unknown")
       open(unit=323,file="qbath_BOT.dat",status="unknown")
       open(unit=324,file="energy_BOT.dat",status="unknown")
       open(unit=325,file="pbath_BOT.dat",status="unknown")
       open(unit=326,file="phase_BOT.dat",status="unknown")
       open(unit=327,file="correlation_BOT.dat",status="unknown")
       open(unit=328,file="reaction_BOT.dat",status="unknown")
       open(unit=329,file="crosscorr_BOT.dat",status="unknown")
       open(unit=330,file="momenta_BOT.dat",status="unknown")

       h = dfloat(last-first)/dfloat(nstep)

       call normalization(nd,qeq(1),ceq,dreal(Beq),S00M,Neq)
       ceqN = ceq/dsqrt(Neq)
       call normalization(nd,q,c0,dreal(Bcmplx),S00M,N)
       cvec= c0/dsqrt(N)
       csq(:) = conjg(cvec(:))*cvec(:)
       phase(:) = datan((aimag(cvec(:))/real(cvec(:))))

       call print_banners(first,last,nstep,h,N)

       call normalization(nd,q,cvec,dreal(Bcmplx),S00M,N)
       call energy(nd,q0,p0,cvec,Bcmplx,H00M,E0,Mx,My)
       E=E0

       write(*,*) N, E, q0(1), p0(1), real(Bcmplx(1,1)),&
                  &aimag(Bcmplx(1,1))!,real(Bcmplxj(1,3))
       write(*,*) csq
       write(*,*) q0(2:nd+1) 

       write(321,*) 0.d0,N,E,q0(1), p0(1), real(Bcmplx(1,1)),&
                    &aimag(Bcmplx(1,1)),real(Bcmplx(2,2)),&
                    &aimag(Bcmplx(2,2))

       write(322,*) 0.d0, csq, dreal(cvec), dimag(cvec)

       write(323,*) 0.d0, q0(2:nd+1) 

       write(325,*) 0.d0, p0(2:nd+1) 

       write(326,*) 0.d0, phase

       qtoti = q0
       ptoti = p0
       Bcmplxi = Bcmplx
       qtotj = q0
       ptotj = p0
       Bcmplxj = Bcmplx

       cj = cvec

       Tau0 = int_TauMat(nd,qtotj,qeq,ptotj,peq,Bcmplxj,Beq)
       Tau0c = matmul(Tau0,ceqN)
       cTau0c = dot_product(cj,Tau0c)
       write(327,*) 0.d0, real(cTau0c),aimag(cTau0c),&
                       real(cTau0c*conjg(cTau0c))/N,&
                       dsqrt(real(cTau0c*conjg(cTau0c)))

       Prob = int_PMat(nd,qtotj(1),Bcmplxj)
       Pc = matmul(Prob,cj)
       cPc = dot_product(cj,Pc)
       write(328,*) 0.d0, cPc/dsqrt(N) 

       qeqX(:) = qeq(:)
       !qeqX(1) = -qeq(1)
       qeqX(:) = -qeq(:)

       TauX = int_TauMat(nd,qtotj,qeqX,ptotj,peq,Bcmplxj,Beq)
       TauXc = matmul(TauX,ceqN)
       cTauXc = dot_product(cj,TauXc)
       write(329,*) 0.d0, real(cTauXc),aimag(cTauXc),&
                       real(cTauXc*conjg(cTauXc))/N,&
                       dsqrt(real(cTauXc*conjg(cTauXc)))
       write(330,*) 0.d0, Mx, My

      write(*,*) "+---------------------------------------------------+"
       write(*,*) "Initial projection overlap (|C(0)|^2)" 
       write(*,*) real(cTau0c*conjg(cTau0c))/N


       kq(:,:) = 0.d0
       kp(:,:) = 0.d0
       kb(:,:,:) = 0.d0
       
       write(111,*) "# starting wfn"
       write(222,*) "# starting wfn"
       call plot_wfn(nd,qtotj,ptotj,cj,Bcmplxj,N)

       nprint=nstep/trj(4)
       time = dfloat(first)
       ckg = complex(0.d0,0.d0)
       qold = qtotj
       pold = ptotj
       Bold = Bcmplxj

! BEGIN TRAJECTORY CYCLE
       do k = 1,nprint
       frac = real(k) / real(nprint)
       pos = int(bar_width * frac)
       if (mod(k, nprint/100) == 0 .or. j == nprint) then
       bar = repeat('#', pos) // repeat('-', bar_width - pos)
       end if
       do j = 1,trj(4)
          time = time + h
          ! Static evolution of coefficents
          cj = c_static(nd,h,cj,S00M,H00M)
          ! Variational evolution of parameters
          qold = qtotj
          pold = ptotj
          Bold = Bcmplxj
!          call rungekutta(nd,h,cj,qtotj,ptotj,Bcmplxj)
          call scprop(nd,h,cj,qtotj,ptotj,Bcmplxj)
!          call vtvprop(nd,h,cj,qtotj,ptotj,Bcmplxj)
          ! Selective freezing
          qtotj(:) = scalvec(:)*qtotj(:)+(1.d0-scalvec(:))*qold
          ptotj(:) = scalvec(:)*ptotj(:)+(1.d0-scalvec(:))*pold
         Bcmplxj(:,:)=scalmat(:,:)*Bcmplxj+(1.d0-scalmat(:,:))*Bold(:,:)
          ! Projection of the coefficients
          cj = c_update(nd,qtotj,ptotj,qold,pold,cj,Bcmplxj,Bold)
!          cj = c_update_fb(nd,qtotj,ptotj,qold,pold,cj,Bcmplxj,Bold)
!          cj = c_update_fbs(nd,qtotj,ptotj,qold,pold,cj,Bcmplxj,Bold)
          csq(:) = conjg(cj(:))*cj(:)
          ! Normalization and energy
          call normalization(nd,qtotj(1),cj,dreal(Bcmplxj),S00M,N)
          call energy(nd,qtotj,ptotj,cj,Bcmplxj,H00M,E,Mx,My)
       end do !j

       phase(:) = datan((aimag(cj(:))/real(cj(:))))
       Tau0 = int_TauMat(nd,qtotj,qeq,ptotj,peq,Bcmplxj,Beq)
       Tau0c = matmul(Tau0,ceqN)
       cTau0c = dot_product(cj,Tau0c)

       Prob = int_PMat(nd,qtotj(1),Bcmplxj)
       Pc = matmul(Prob,cj)
       cPc = dot_product(cj,Pc)

       TauX = int_TauMat(nd,qtotj,qeqX,ptotj,peq,Bcmplxj,Beq)
       TauXc = matmul(TauX,ceqN)
       cTauXc = dot_product(cj,TauXc)

      write(321,*) time,N,E,qtotj(1),ptotj(1),real(Bcmplxj(1,1))&
                  &,aimag(Bcmplxj(1,1)),real(Bcmplxj(2,2))&
                  &,aimag(Bcmplxj(2,2))
      write(322,*) time, csq, dreal(cj), dimag(cj)
      write(323,*) time, qtotj(2:nd+1) 
      write(325,*) time, ptotj(2:nd+1) 
      write(326,*) time, phase 
      write(327,*) time,real(cTau0c),aimag(cTau0c),&
                   real(cTau0c*conjg(cTau0c))/N,&
                   dsqrt(real(cTau0c*conjg(cTau0c))/N)
      write(328,*) time, cPc/dsqrt(N) 
      write(329,*) time, real(cTauXc),aimag(cTauXc),&
                   real(cTauXc*conjg(cTauXc))/N,&
                   dsqrt(real(cTauXc*conjg(cTauXc)))
      write(330,*) time, Mx, My
      
      write(111,*) "#", time 
      write(222,*) "#", time 
      call plot_wfn(nd,qtotj,ptotj,cj,Bcmplxj,N)
      write(*,'(A,F6.2,A)', advance='no') char(13)//'['//bar//'] ', &
                 frac*100.0, '%'
          
       end do !k

       write(*,*) "Last step:"
       write(*,*) N,E, qtotj(1), ptotj(1),real(Bcmplxj(1,1)),&  
                  &aimag(Bcmplxj(1,1))
       write(*,*) time, csq

       close(321)
       close(322)
       close(323)
       close(324)
       close(325)
       close(326)
       close(327)
       close(328)
       close(329)
       close(330)

       end subroutine

       end module
