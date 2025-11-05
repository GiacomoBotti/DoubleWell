!**********************************************************************!
! Module containing the Runge-Kutta evolution of the parameters and    !
! coefficients                                                         !
!**********************************************************************!

       module evolution_module

       use timedev_module
       use dtmatrix_module
       use normalization_module
       use observable_module
       use BOT_module

       implicit none

       private
       public :: homogen,rungekutta,bot_evo,scprop,bot_scp

       contains

!......Homogeneus term function.........................................

       function homogen(nd,npar,y,work) result(f)
       ! nd: basis set dimension
       ! npar: parameter space dimension
       ! y: runge-kutta function - [parameters - coefficients]
       ! work: work array for I/O
         implicit none
         integer*8, intent(in) :: nd,npar
         complex*16, dimension(npar+nd), intent(in) :: y
         real*8, dimension(4), intent(in) :: work
         
         real*8 :: q,p,a,alpha,xi,zeta
         real*8, dimension(npar) :: dotlam
         complex*16, dimension(nd) :: c,dotc
         complex*16, dimension(npar+nd) :: f 
         
         ! extract work
         a = work(1)
         alpha = work(2)
         xi = work(3)
         zeta = work(4)

         ! extract y
         q = dreal(y(1))
         p = dreal(y(2))
         c(:) = y(npar+1:npar+nd)
  
         ! coefficients time derivative
!         dotc = dotcoeff(nd,npar,q,p,a,alpha,xi,zeta,c)
         ! For numerical matrices
          dotc = dotcoeffnum(nd,npar,q,p,a,alpha,xi,zeta,c)

         ! parameters time derivative
!         dotlam = dotlambda(nd,npar,q,p,a,alpha,xi,zeta,c)   
         ! For numerical matrices
          dotlam = dotlamnum(nd,npar,q,p,a,alpha,xi,zeta,c)   

         f(:) = complex(0.d0,0.d0)
         f(1:npar) = dotlam(:)
         f(npar+1:npar+nd) = dotc(:)

       end function

!......Four Steps Runge-Kutta...........................................

       subroutine rungekutta(nd,npar,trj,y0,work)
       ! trj : trajectory parameters (first step, last step, nstep)
       ! y0 : initial condition vector (q,p,c)
       ! work : array containing all non-variational parameters
       integer*8, intent(in) :: nd,npar
       integer*8, dimension(3), intent(in) :: trj
       real*8,dimension(4),intent(in) :: work
       complex*16, dimension(npar+nd), intent(in) :: y0

       integer*8 :: i,j,first,last,nstep
       real*8 :: h,time,N,E
       real*8,dimension(4) :: hvec = [0.5d0,0.5d0,1.d0,0.d0]
       
       complex*16,dimension(npar+nd) :: yj,yi
       complex*16,dimension(npar+nd,4) :: k

       first = trj(1)
       last = trj(2)
       nstep = trj(3)

       open(unit=123,file="trajectory_RK.dat",status="unknown")
       h = dfloat(last-first)/dfloat(nstep)
       hvec = hvec*h

!       N = normalization(nd,npar,y0,work)
       ! Numerical:
       N = normnumtot(nd,npar,y0,work)

       write(123,*) "#Evolution parameters:"
       write(123,*) "#Range: ",first,last
       write(123,*) "#Steps: ",nstep
       write(123,*) "#Timestep: ",h
       write(123,*) "#Normalization constant: ",N
       !write(123,*) "#Square root of N ",dsqrt(N)
  
       yj(1:npar) = y0(1:npar)
       yj(npar+1:npar+nd) = y0(npar+1:npar+nd)/dsqrt(N)

       yi(1:npar) = y0(1:npar)
       yi(npar+1:npar+nd) = y0(npar+1:npar+nd)/dsqrt(N)

!       N = normalization(nd,npar,yj,work)
!       E = energy(nd,npar,yj,work)
       ! Numerical:
       N = normnumtot(nd,npar,yj,work)
       E = energynum(nd,npar,yj,work)

       write(123,*) "#Time ","q ","p ","rec1 ","rec2 ","imc1 ","imc2 ",&
                    &"Norm ","<H> "
       write(123,*) "#0.d0",dreal(y0(1)),dreal(y0(2)),dreal(y0(3)),&
                   &dreal(y0(4)),dimag(y0(3)),dimag(y0(4)),N,E      
       write(123,*) 0.d0,dreal(yj(1)),dreal(yj(2)),dreal(yj(3)),&
                   &dreal(yj(4)),dimag(yj(3)),dimag(yj(4)),N,E     
       write(123,*) 0.d0,dreal(yi(1)),dreal(yi(2)),dreal(yi(3)),&
                   &dreal(yi(4)),dimag(yi(3)),dimag(yi(4)),N,E      

       k(:,:) = complex(0.d0,0.d0)
       do j = 1,nstep
          time = j*h
          do i = 1,4
             k(:,i) = homogen(nd,npar,yi,work)
             yi = yj + hvec(i)*k(:,i)
          end do
          yj = yj + h*(k(:,1)+2.d0*k(:,2)+2.d0*k(:,3)+k(:,4))/6.d0
!          N = normalization(nd,npar,yj,work) 
!          E = energy(nd,npar,yj,work)
          ! Numerical:
          N = normnumtot(nd,npar,yj,work)
          E = energynum(nd,npar,yj,work)
          write(123,*) time,dreal(yj(1)),dreal(yj(2)),dreal(yj(3)),&
                   &dreal(yj(4)),dimag(yj(3)),dimag(yj(4)),N,E      
          yi = yj 
       end do
       
       close(123)  
     
       end subroutine

!......BOT evolution....................................................

       subroutine bot_evo(nd,npar,trj,y0,work)
       ! trj : trajectory parameters (first step, last step, nstep)
       ! y0 : initial condition vector (q,p,c)
       ! work : array containing all non-variational parameters
       integer*8, intent(in) :: nd,npar
       integer*8, dimension(3), intent(in) :: trj
       real*8,dimension(4),intent(in) :: work
       complex*16, dimension(npar+nd), intent(in) :: y0

       integer*8 :: i,j,first,last,nstep
       real*8 :: h,time,N,E,a,alpha,xi,zeta,q,p,qold,pold,nover
       real*8,dimension(4) :: hvec = [0.5d0,0.5d0,1.d0,0.d0]
       real*8,dimension(3) :: cond,nparam 
       !real*8,dimension(nd) :: vecxi,veczeta
       
       complex*16,dimension(npar+nd) :: yj,yi
       complex*16,dimension(nd) :: cj,cjsupp 
!       complex*16,dimension(npar) :: lambdaj,lambdai,lambdaold
       real*8,dimension(npar) :: lambdaj,lambdai,lambdaold

       real*8,dimension(npar,4) :: k
       complex*16,dimension(nd,nd) :: S00M,S00mat,Tt0M,Tt0mat,invS 

       ! trajectory parameters
       first = trj(1)
       last = trj(2)
       nstep = trj(3)

       ! extract work
       a = work(1)
       alpha = work(2)
       xi = work(3)
       zeta = work(4)
       
       ! A module for Christ's sake
       !vecxi(1) = xi
       !vecxi(2) = -1.d0*xi
       !veczeta(1) = zeta
       !veczeta(2) = -1.d0*zeta
       cond = [-10.d0,10.d0,1000.d0]

       open(unit=321,file="trajectory_BOT.dat",status="unknown")
       h = dfloat(last-first)/dfloat(nstep)
       hvec = hvec*h

       N = normalization(nd,npar,y0,work)
       ! Numerical:
!       N = normnumtot(nd,npar,y0,work)

       write(321,*) "#Evolution parameters:"
       write(321,*) "#Range: ",first,last
       write(321,*) "#Steps: ",nstep
       write(321,*) "#Timestep: ",h
       write(321,*) "#Normalization constant: ",N
       !write(321,*) "#Square root of N ",dsqrt(N)
  
       lambdaj(1:npar) = y0(1:npar)
       cj(:) = y0(npar+1:npar+nd)/dsqrt(N)

       lambdai(1:npar) = y0(1:npar)

       yj(1:npar) = lambdaj(:)
       yj(npar+1:npar+nd) = cj(:)

       N = normalization(nd,npar,yj,work)
       E = energy(nd,npar,yj,work)
       ! Numerical:
!       N = normnumtot(nd,npar,yj,work)
!       E = energynum(nd,npar,yj,work)

       write(321,*) "#Time ","q ","p ","rec1 ","rec2 ","imc1 ","imc2 ",&
                    &"Norm ","<H> "
       write(321,*) "#0.d0",dreal(y0(1)),dreal(y0(2)),dreal(y0(3)),&
                   &dreal(y0(4)),dimag(y0(3)),dimag(y0(4)),N,E      
       write(321,*) 0.d0,dreal(yj(1)),dreal(yj(2)),dreal(yj(3)),&
                   &dreal(yj(4)),dimag(yj(3)),dimag(yj(4)),N,E     

       k(:,:) = 0.d0
       do j = 1,nstep
          time = j*h
          cj = c_static(nd,npar,yj,work,h)
!          cj = c_stat_num(nd,npar,yj,work,h)
          do i = 1,4
             q = lambdai(1)
             p = lambdai(2)
             k(:,i) = dotlambda(nd,npar,q,p,a,alpha,xi,zeta,cj)
!             k(:,i) = dotlamnum(nd,npar,q,p,a,alpha,xi,zeta,cj)
             lambdai = lambdaj + hvec(i)*k(:,i)
          end do
          lambdaold = lambdaj
          lambdaj=lambdaj+h*(k(:,1)+2.d0*k(:,2)+2.d0*k(:,3)+k(:,4))/6.d0

          cj = c_update(nd,npar,lambdaj,lambdaold,cj,work)
!          cj = c_upd_num(nd,npar,lambdaj,lambdaold,cj,work)

          yj(1:npar) = lambdaj(:)
          yj(npar+1:npar+nd) = cj(:)

          N = normalization(nd,npar,yj,work)
          E = energy(nd,npar,yj,work)
          ! Numerical:
!          N = normnumtot(nd,npar,yj,work)
!          E = energynum(nd,npar,yj,work)

          write(321,*) time,dreal(yj(1)),dreal(yj(2)),dreal(yj(3)),&
                   &dreal(yj(4)),dimag(yj(3)),dimag(yj(4)),N,E      
      
          lambdai = lambdaj
       end do

       close(321)

       end subroutine

!......Self Consistent Propagation......................................

       subroutine scprop(nd,npar,trj,y0,work)
       ! trj : trajectory parameters (first step, last step, nstep)
       ! y0 : initial condition vector (q,p,c)
       ! work : array containing all non-variational parameters
       integer*8, intent(in) :: nd,npar
       integer*8, dimension(3), intent(in) :: trj
       real*8,dimension(4),intent(in) :: work
       complex*16, dimension(npar+nd), intent(in) :: y0

       integer*8 :: i,j,first,last,nstep
       integer*8 :: maxcycle=10  
       real*8 :: h,time,N,E,error
       real*8 :: thr = 1.d-5
       real*8,dimension(4) :: hvec = [0.5d0,0.5d0,1.d0,0.d0]
       
       complex*16,dimension(npar+nd) :: yj,yi,yiold,yav
       complex*16,dimension(npar+nd) :: k
       real*8,dimension(npar+nd) :: square

       first = trj(1)
       last = trj(2)
       nstep = trj(3)

       open(unit=132,file="trajectory_scp.dat",status="unknown")
       open(unit=99,file="scp_evo_log.dat",status="unknown")

       h = dfloat(last-first)/dfloat(nstep)
       hvec = hvec*h

       N = normalization(nd,npar,y0,work)
       ! Numerical:
!       N = normnumtot(nd,npar,y0,work)

       write(132,*) "#Evolution parameters:"
       write(132,*) "#Range: ",first,last
       write(132,*) "#Steps: ",nstep
       write(132,*) "#Timestep: ",h
       write(132,*) "#Normalization constant: ",N
       !write(132,*) "#Square root of N ",dsqrt(N)
  
       write(99,*) "#Self consistent propagation data"
       write(99,*) "# maxcycle: ", maxcycle
       write(99,*) "# threshold: ", thr
 
       yj(1:npar) = y0(1:npar)
       yj(npar+1:npar+nd) = y0(npar+1:npar+nd)/dsqrt(N)

       yi(1:npar) = y0(1:npar)
       yi(npar+1:npar+nd) = y0(npar+1:npar+nd)/dsqrt(N)


       N = normalization(nd,npar,yj,work)
       E = energy(nd,npar,yj,work)
       ! Numerical:
!       N = normnumtot(nd,npar,yj,work)
!       E = energynum(nd,npar,yj,work)

       write(132,*) "#Time ","q ","p ","rec1 ","rec2 ","imc1 ","imc2 ",&
                    &"Norm ","<H> "
       write(132,*) "#0.d0",dreal(y0(1)),dreal(y0(2)),dreal(y0(3)),&
                   &dreal(y0(4)),dimag(y0(3)),dimag(y0(4)),N,E      
       write(132,*) 0.d0,dreal(yj(1)),dreal(yj(2)),dreal(yj(3)),&
                   &dreal(yj(4)),dimag(yj(3)),dimag(yj(4)),N,E     
       write(132,*) 0.d0,dreal(yi(1)),dreal(yi(2)),dreal(yi(3)),&
                   &dreal(yi(4)),dimag(yi(3)),dimag(yi(4)),N,E      

       k(:) = complex(0.d0,0.d0)
       do j = 1,nstep
          time = j*h
          do i = 1,maxcycle
             yiold = yi
             yav = (yi + yj)/2.d0
             k(:) = homogen(nd,npar,yav,work)
             yi = yj + h*k
             square = (yiold - yi)*dconjg(yiold - yi)
             error = sum(square)
             if (dsqrt(error).le.thr.or.i.eq.maxcycle) then
                write(99,*) time, i
                exit
             end if
          end do
          yj = yi
          N = normalization(nd,npar,yj,work) 
          E = energy(nd,npar,yj,work)
          ! Numerical:
!          N = normnumtot(nd,npar,yj,work)
!          E = energynum(nd,npar,yj,work)
          write(132,*) time,dreal(yj(1)),dreal(yj(2)),dreal(yj(3)),&
                   &dreal(yj(4)),dimag(yj(3)),dimag(yj(4)),N,E      
       end do
       
       close(132)  
       close(99)  
     
       end subroutine

!......BOT +SCP evolution...............................................

       subroutine bot_scp(nd,npar,trj,y0,work)
       ! trj : trajectory parameters (first step, last step, nstep)
       ! y0 : initial condition vector (q,p,c)
       ! work : array containing all non-variational parameters
       integer*8, intent(in) :: nd,npar
       integer*8, dimension(3), intent(in) :: trj
       real*8,dimension(4),intent(in) :: work
       complex*16, dimension(npar+nd), intent(in) :: y0

       integer*8 :: i,j,first,last,nstep
       integer*8 :: maxcycle=10  
       real*8 :: h,time,N,E,a,alpha,xi,zeta,q,p,qold,pold,nover,error
       real*8 :: thr = 1.d-5
       real*8,dimension(4) :: hvec = [0.5d0,0.5d0,1.d0,0.d0]
       real*8,dimension(3) :: cond,nparam 
       !real*8,dimension(nd) :: vecxi,veczeta
       
       complex*16,dimension(npar+nd) :: yj,yi
       complex*16,dimension(nd) :: cj,cjsupp 
!       complex*16,dimension(npar) :: lambdaj,lambdai,lambdaold
       real*8,dimension(npar) :: lambdaj,lambdai,lambdaold,lambdaav

       real*8,dimension(npar) :: k,lambdaiold,square
       complex*16,dimension(nd,nd) :: S00M,S00mat,Tt0M,Tt0mat,invS 

       ! trajectory parameters
       first = trj(1)
       last = trj(2)
       nstep = trj(3)

       ! extract work
       a = work(1)
       alpha = work(2)
       xi = work(3)
       zeta = work(4)
       
       ! A module for Christ's sake
       !vecxi(1) = xi
       !vecxi(2) = -1.d0*xi
       !veczeta(1) = zeta
       !veczeta(2) = -1.d0*zeta
       cond = [-10.d0,10.d0,1000.d0]

       open(unit=231,file="trajectory_BOTSCP.dat",status="unknown")
       open(unit=88,file="scp_bot_log.dat",status="unknown")
       h = dfloat(last-first)/dfloat(nstep)
       hvec = hvec*h

       N = normalization(nd,npar,y0,work)
       ! Numerical:
!       N = normnumtot(nd,npar,y0,work)

       write(231,*) "#Evolution parameters:"
       write(231,*) "#Range: ",first,last
       write(231,*) "#Steps: ",nstep
       write(231,*) "#Timestep: ",h
       write(231,*) "#Normalization constant: ",N
       !write(231,*) "#Square root of N ",dsqrt(N)
  
       write(88,*) "#Self consistent propagation data"
       write(88,*) "# maxcycle: ", maxcycle
       write(88,*) "# threshold: ", thr

       lambdaj(1:npar) = y0(1:npar)
       cj(:) = y0(npar+1:npar+nd)/dsqrt(N)

       lambdai(1:npar) = y0(1:npar)

       yj(1:npar) = lambdaj(:)
       yj(npar+1:npar+nd) = cj(:)

       N = normalization(nd,npar,yj,work)
       E = energy(nd,npar,yj,work)
       ! Numerical:
!       N = normnumtot(nd,npar,yj,work)
!       E = energynum(nd,npar,yj,work)

       write(231,*) "#Time ","q ","p ","rec1 ","rec2 ","imc1 ","imc2 ",&
                    &"Norm ","<H> "
       write(231,*) "#0.d0",dreal(y0(1)),dreal(y0(2)),dreal(y0(3)),&
                   &dreal(y0(4)),dimag(y0(3)),dimag(y0(4)),N,E      
       write(231,*) 0.d0,dreal(yj(1)),dreal(yj(2)),dreal(yj(3)),&
                   &dreal(yj(4)),dimag(yj(3)),dimag(yj(4)),N,E     

       k(:) =0.d0
       do j = 1,nstep
          time = j*h
          cj = c_static(nd,npar,yj,work,h)
!          cj = c_stat_num(nd,npar,yj,work,h)
          do i = 1,maxcycle
             lambdaiold = lambdai
             lambdaav = (lambdai+lambdaj)/2.d0
             q = lambdaav(1)
             p = lambdaav(2)
             k(:) = dotlambda(nd,npar,q,p,a,alpha,xi,zeta,cj)
!             k(:) = dotlamnum(nd,npar,q,p,a,alpha,xi,zeta,cj)
             lambdai = lambdaj + h*k(:)
             square = (lambdaiold-lambdai)**2
             error = sum(square)
             if (dsqrt(error).le.thr.or.i.eq.maxcycle) then
                write(88,*) time,i
                exit
             end if
          end do
          lambdaold = lambdaj
          lambdaj=lambdai

          cj = c_update(nd,npar,lambdaj,lambdaold,cj,work)
!          cj = c_upd_num(nd,npar,lambdaj,lambdaold,cj,work)

          yj(1:npar) = lambdaj(:)
          yj(npar+1:npar+nd) = cj(:)

          N = normalization(nd,npar,yj,work)
          E = energy(nd,npar,yj,work)
          ! Numerical:
!          N = normnumtot(nd,npar,yj,work)
!          E = energynum(nd,npar,yj,work)

          write(231,*) time,dreal(yj(1)),dreal(yj(2)),dreal(yj(3)),&
                   &dreal(yj(4)),dimag(yj(3)),dimag(yj(4)),N,E      
      
       end do

       close(231)
       close(88)  

       end subroutine

       end module
