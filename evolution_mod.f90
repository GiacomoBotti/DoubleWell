!**********************************************************************!
! Module containing the Runge-Kutta evolution of the parameters and    !
! coefficients                                                         !
!**********************************************************************!

       module evolution_module

       use basisset_module
       use normalization_module
       use observable_module
!       use BOT_module
       use eofmotion_module
       use integrals_module

       implicit none

       private
       public :: bot_evo

       contains

!......BOT evolution....................................................

       subroutine bot_evo(nd,trj,q0,p0,c0,Bcmplx)
       ! nd: bath dimension
       ! trj : trajectory parameters (first step, last step, nstep)
       ! q0 : initial gaussian center (x&y)
       ! p0 : initial gaussian momentum (x&y)
       ! c0 : initial coefficients 
       ! Bcmplx : initial gaussian width matrix, real & imaginary
       integer, intent(in) :: nd
       integer*8, dimension(3), intent(in) :: trj
       real*8, dimension(nd+1), intent(in) :: q0,p0
       complex*16,dimension(nh), intent(in) :: c0 
       complex*16,dimension(nd+1,nd+1), intent(in) :: Bcmplx 

       integer*8 :: i,j,first,last,nstep,k
       real*8 :: h,time,N,E,q,p,E0,Nsq
       real*8,dimension(nd) :: qvec,pvec 
       real*8,dimension(nd+1) :: qtoti,ptoti,qtotj,ptotj,qold,pold 
       real*8,dimension(4) :: hvec = [0.5d0,0.5d0,1.d0,0.d0]
       complex*16,dimension(nd+1,nd+1) :: Bcmplxi,Bcmplxj,Bold
       
       complex*16,dimension(nh) :: cvec,cj 

       real*8,dimension(nd+1,4) :: kq,kp
       complex*16,dimension(nd+1,nd+1,4) :: kb

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
      write(*,*) "+---------------------------------------------------+"
       write(*,*) "Writing trajectory output on trajectory_BOT.dat"
       write(*,*) "Writing coefficients output on coefficients_BOT.dat"
       write(*,*) "Writing bath positions on bath_BOT.dat"
       write(*,*) "Writing energy components on enegy_BOT.dat"
      write(*,*) "+---------------------------------------------------+"
       h = dfloat(last-first)/dfloat(nstep)
       hvec = hvec*h

       N = normalization(nd,q,c0,dreal(Bcmplx))

       write(321,*) "#Evolution parameters:"
       write(321,*) "#Range: ",first,last
       write(321,*) "#Steps: ",nstep
       write(321,*) "#Timestep: ",h
       write(321,*) "#Normalization constant: ",N

       write(322,*) "#Evolution parameters:"
       write(322,*) "#Range: ",first,last
       write(322,*) "#Steps: ",nstep
       write(322,*) "#Timestep: ",h
       write(322,*) "#Normalization constant: ",N

       write(323,*) "#Evolution parameters:"
       write(323,*) "#Range: ",first,last
       write(323,*) "#Steps: ",nstep
       write(323,*) "#Timestep: ",h
       write(323,*) "#Normalization constant: ",N
  
       write(324,*) "#Evolution parameters:"
       write(324,*) "#Range: ",first,last
       write(324,*) "#Steps: ",nstep
       write(324,*) "#Timestep: ",h
       write(324,*) "#Normalization constant: ",N
       write(324,*) "#H ", "T ", "V "
  
       write(325,*) "#Evolution parameters:"
       write(325,*) "#Range: ",first,last
       write(325,*) "#Steps: ",nstep
       write(325,*) "#Timestep: ",h
       write(325,*) "#Normalization constant: ",N
       write(325,*) "#H ", "T ", "V "

       cvec= c0/dsqrt(N)

       N = normalization(nd,q,cvec,dreal(Bcmplx))
       E0 = energy(nd,q0,p0,c0,Bcmplx)
       E=E0

       write(*,*) "First step:"
       write(*,*) N, E/E0, q0(1), p0(1), real(Bcmplx(1,1)),&
                  &real(Bcmplxj(3,3)),real(Bcmplxj(1,3))
       write(*,*) time, dreal(cvec), dimag(cvec)
       write(*,*) time, q0(2:nd+1) 

       write(321,*) "#Time ","N ","E ","q ","p ","B(1,1) ",&
                    &"B(3,3) ", "B(1,3)"
       write(321,*) 0.d0,N,E/E0,q0(1), p0(1), real(Bcmplx(1,1)),&
                    &real(Bcmplx(3,3)),real(Bcmplx(1,3))

       write(322,*) "#Time ","Real c", "Immaginary c"
       write(322,*) 0.d0, dreal(c0), dimag(c0)

       write(323,*) "#Time ","qbath"
       write(323,*) 0.d0, q0(2:nd+1) 

       write(323,*) "#Time ","pbath"
       write(323,*) 0.d0, p0(2:nd+1) 

       qtoti = q0
       ptoti = p0
       Bcmplxi = Bcmplx
       qtotj = q0
       ptotj = p0
       Bcmplxj = Bcmplx
 
       cj = cvec

       kq(:,:) = 0.d0
       kp(:,:) = 0.d0
       kb(:,:,:) = 0.d0
       do j = 1,nstep
          time = j*h
!          cj = c_static(nd,npar,yj,work,h)
          do i = 1,4
             call KarplusTimeDer(nd,cj,qtoti,ptoti,Bcmplxi,&
                  &kq(:,i),kp(:,i),kb(:,:,i))             
             qtoti = qtotj + hvec(i)*kq(:,i)
             ptoti = ptotj + hvec(i)*kp(:,i)
             Bcmplxi = Bcmplxj + hvec(i)*kb(:,:,i)
          end do
          qold = qtotj
          pold = ptotj
          Bold = Bcmplxj
          qtotj = qtotj&
             &+h*(kq(:,1)+2.d0*kq(:,2)+2.d0*kq(:,3)+kq(:,4))/6.d0
          ptotj = ptotj&
             &+h*(kp(:,1)+2.d0*kp(:,2)+2.d0*kp(:,3)+kp(:,4))/6.d0
          Bcmplxj = Bcmplxj&
          &+h*(kb(:,:,1)+2.d0*kb(:,:,2)+2.d0*kb(:,:,3)+kb(:,:,4))/6.d0

!          cj = c_update(nd,npar,lambdaj,lambdaold,cj,work)

          N = normalization(nd,qtotj(1),cj,dreal(Bcmplxj))
          E = energy(nd,qtotj,ptotj,cj,Bcmplxj)

          write(321,*) time,N,E/E0,qtotj(1),ptotj(1),real(Bcmplxj(1,1))&
                       &,real(Bcmplxj(3,3)),real(Bcmplxj(1,3))
          write(322,*) time, dreal(cj), dimag(cj)
          write(323,*) time, qtotj(2:nd+1) 
          write(325,*) time, ptotj(2:nd+1) 
        
          Nsq=fun_Nsq(nd+1,real(Bcmplxj))

          write(444,*) "Time: ", time, "Nsq: ", Nsq
          do k = 1,nd+1
             write(444,*) Bcmplxj(k,:)
          end do

          qtoti = qtotj  
          ptoti = ptotj
          Bcmplxi = Bcmplxj
       end do

       write(*,*) "Last step:"
       write(*,*) N,E/E0, qtotj(1), ptotj(1),real(Bcmplxj(1,1)),&  
                  &real(Bcmplxj(3,3)),real(Bcmplxj(1,3))
       write(*,*) time, dreal(cj), dimag(cj)
       write(*,*) time, qtotj(2:nd+1) 

       write(222,*) qtotj
       write(222,*) ptotj
       do i = 1,nd+1
         write(222,*) Bcmplxj(i,:)
       end do

       close(321)
       close(322)
       close(323)
       close(324)

       end subroutine

       end module
