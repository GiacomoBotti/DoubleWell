!*********************************************************************!
! Module containing checking subroutines                              !
!*********************************************************************!

      module check_module

      use constants
      use matrix_module 
      use potential_module
      use basisset_module
      use integrals_module
      use effectivepot_module
      use normalization_module
      use kinetic_module
      use BOT_module
      use observable_module

      implicit none

      contains

!.....Check Diagonalization.............................................

      subroutine check_diagonalization(nd)

       integer :: i,j,nd
       real*8 :: harvest
       real*8, dimension(nd,nd) :: LambdaMat,Tmat,RndMat,RecMat,Amat

       do i = 1,nd
          call RANDOM_NUMBER(harvest)
          RndMat(i,i) = harvest
          do j = i+1,nd
             call RANDOM_NUMBER(harvest)
             RndMat(i,j) = harvest
             RndMat(j,i) = RndMat(i,j)
          end do
       end do

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_DIAGONALIZATION"
       write(*,*) "Starting matrix"
     
       do i = 1,nd
          write(*,*) RndMat(i,:)
       end do

       Amat = RndMat 
 
       call diagonalization(nd,Amat,LambdaMat,Tmat)

       write(*,*) "Diagonal matrix"
     
       do i = 1,nd
          write(*,*) LambdaMat(i,:)
       end do
 
       RecMat = matmul(LambdaMat,transpose(Tmat))
       RecMat = matmul(Tmat,RecMat)

       write(*,*) "Recovered matrix"

       do i = 1,nd
          write(*,*) RecMat(i,:)
       end do

       write(*,*) "Error matrix"

       do i = 1,nd
          write(*,*) RecMat(i,:)-RndMat(i,:)
       end do

      end subroutine

!.....Check Momenta.....................................................

      subroutine check_momenta(nd)
      ! nd: dimensions of the matrix

       integer :: i,nd
       real*8, dimension(nd,nd) :: DiagMat
       real*8, dimension(maxorder,nd) :: OutMat

       DiagMat(:,:) = 0.d0

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_MOMENTA"
       write(*,*) "Diagonal Matrix"
       do i = 1,nd
          DiagMat(i,i) = 5.d0*i
          write(*,*) DiagMat(i,:)
       end do

       OutMat = momenta(nd,DiagMat)

       write(*,*) "Matrix of the Momenta"

       do i = 1,maxorder
          write(*,*) Outmat(i,:)
       end do

      end subroutine

!.....Check extractA....................................................

      subroutine check_Amat(nd)
      ! nd : dimension of the matrix

       integer :: i,j,nd
       real*8 :: harvest
       real*8, dimension(nd+1,nd+1) :: RndMat

       real*8 :: a
       real*8, dimension(nd) :: avec 
       real*8, dimension(nd,nd) :: Amat

       do i = 1,nd+1
          call RANDOM_NUMBER(harvest)
          RndMat(i,i) = harvest
          do j = i+1,nd+1
             call RANDOM_NUMBER(harvest)
             RndMat(i,j) = harvest
             RndMat(j,i) = RndMat(i,j)
          end do
       end do

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_AMAT"
       write(*,*) "Starting matrix"
     
       do i = 1,nd+1
          write(*,*) RndMat(i,:)
       end do

       call extractA(nd,RndMat,Amat,avec,a)

       write(*,*) "A matrix"

       do i = 1,nd
          write(*,*) Amat(i,:)
       end do

       write(*,*) "a vector"
       
       write(*,*) avec(:)

       write(*,*) "a: ", a

      end subroutine

!.....Check Vmat........................................................

      subroutine check_vmat()
      
       integer :: i

       write(*,*) "# HELLO I'M CHECK_VMAT"
       write(*,*) "Potential matrix"
       do i = 1,nv
         write(*,*) Vmat(i,:)
       end do

      end subroutine  

!.....Check HermMat.....................................................

      subroutine check_hermmat()
      
       integer :: i,j,k,l
       real*8 :: H
       real*8,dimension(nh,nh) :: Hmat

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_HERMMAT"
       write(*,*) "Hermite Coefficients matrix"
       do i = 1,nh
         write(*,*) Mherm(i,:)
       end do
       
       write(*,*) "------------------------------"
       write(*,*) "The values of some polynomials"
       write(*,*) "------------------------------"
       do i = 1,3
         do j = 1,3
           do l = 1,nh
             H = herm_pol(l,dfloat(i),dfloat(j),3.d0)
             write(*,*) "npol: ",l,"x: ",i,"q: ",j,"H: ",H
           end do 
         end do
       end do

       write(*,*) "------------------------------"
       write(*,*) "The values of Hmat"
       write(*,*) "------------------------------"

       do i = 1,3
         do j =1,3
           write(*,*) "x: ",i,"q: ",j
           Hmat = fun_Hmat(dfloat(i),dfloat(j),3.d0)
           do k = 1,nh
             write(*,*) Hmat(k,:)
           end do
         end do
       end do

       write(*,*) "------------------------------"
       write(*,*) "The values of Hmat Shifted"
       write(*,*) "------------------------------"

       do i = 1,3
         do j =1,3
           write(*,*) "x: ",1.d0,"qi: ",i,"qj: ",j
           Hmat = fun_HmatShift(1.d0,dfloat(i),dfloat(j),3.d0,4.d0)
           do k = 1,nh
             write(*,*) Hmat(k,:)
           end do
         end do
       end do

      end subroutine  

!.....Check Y0..........................................................

      subroutine check_Y0(nd)
      ! nd: dimensions of the matrix

       integer :: i,nd
       real*8 :: Y0
       real*8, dimension(nd,nd) :: DiagMat

       DiagMat(:,:) = 0.d0

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_Y0"
       write(*,*) "------------------------------"
       write(*,*) "Diagonal Matrix"
       do i = 1,nd
          DiagMat(i,i) = 5.d0*i
          write(*,*) DiagMat(i,:)
       end do

       Y0 = int_Y0(nd,DiagMat)

       write(*,*) "------------------------------"
       write(*,*) "Y0: ", Y0

      end subroutine

!.....Check XnMat.......................................................

      subroutine check_XnMat(nd)
      !nd: dimensions of the bath matrix
       integer, intent(in) :: nd

       integer :: i,j
       real*8 :: a,q
       real*8, dimension(nd) :: avec 
       real*8, dimension(nd,nd) :: Amat
       real*8, dimension(nd+1,nd+1) :: Bmat
       real*8, dimension(nh,nh) :: XnMat

       Bmat(:,:) = 0.1d0
       Bmat(1,:) = 0.5d0
       Bmat(:,1) = 0.5d0
       Bmat(1,1) = 1.d0

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_XNMAT"
       write(*,*) "B matrix"
     
       do i = 1,nd+1
          write(*,*) Bmat(i,:)
       end do

       call extractA(nd,Bmat,Amat,avec,a)

       write(*,*) "------------------------------"
       write(*,*) "A matrix"

       do i = 1,nd
          write(*,*) Amat(i,:)
       end do

       write(*,*) "------------------------------"
       write(*,*) "a vector"
       
       write(*,*) avec(:)

       write(*,*) "------------------------------"
       write(*,*) "a: ", a
         
       do i = 1,4
         q = dfloat(i)/2.d0
         write(*,*) "------------------------------"
         write(*,*) "q: ", q
         XnMat = int_XnMat(nd,4,a,avec,Amat,q)
         do j = 1,nh
           write(*,*) XnMat(j,:)
         end do
       end do

      end subroutine

!.....Check V0..........................................................

      subroutine check_V0(nd,cvec)
      ! nd: dimensions of the bath 
      ! cvec: vector of the coefficients
       integer, intent(in) :: nd
       complex*16, dimension(nh),intent(in) :: cvec

       integer :: i,j
       real*8 :: Nsq
       real*8, dimension(nd+1) :: qtot, V1
       real*8, dimension(nd+1,nd+1) :: Bmat,V2
       real*8 :: harvest
       real*8, dimension(nh,nh) :: V0
       real*8 :: a,Y0
       real*8, dimension(nd) :: avec
       real*8, dimension(nd,nd) :: Amat,Tmat,LambdaMat
       real*8, dimension(nh,nh) :: X4,X3,X2,X1,X0

       do i = 1,nd+1
          Bmat(i,i) = i
          do j = i+1,nd+1
             Bmat(i,j) = j/20.d0 !Gershgoring circle theorem
             Bmat(j,i) = Bmat(i,j)
          end do
       end do

       qtot(1) = 1.d0
       qtot(2) = 2.d0
       qtot(3:nd+1) = 10000.d0 !This way I know if something is wrong


       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_V0"
       write(*,*) "B matrix"
     
       do i = 1,nd+1
       
          write(*,*) Bmat(i,:)
       end do
       call extractA(nd,Bmat,Amat,avec,a)
       call diagonalization(nd,Amat,LambdaMat,Tmat)

       ! Integrals
       Y0=int_Y0(nd,LambdaMat)
       X4=int_XnMat(nd,4,a,avec,Amat,qtot(1))
       X3=int_XnMat(nd,3,a,avec,Amat,qtot(1))
       X2=int_XnMat(nd,2,a,avec,Amat,qtot(1))
       X1=int_XnMat(nd,1,a,avec,Amat,qtot(1))
       X0=int_XnMat(nd,0,a,avec,Amat,qtot(1))

       write(*,*) "------------------------------"
       Nsq = fun_Nsq(nd+1,Bmat)
       write(*,*) "N squared: ", Nsq

       write(*,*) "------------------------------"
       V0 = fun_V0(nd,qtot,cvec,Bmat,Y0,X4,X2,X1,X0)
       write(*,*) "V0(1,1) ", V0(1,1)
       
       write(*,*) "------------------------------"
       V1 = fun_V1(nd,qtot,cvec,Bmat,Y0,X3,X2,X1,X0)
       write(*,*) "V1: "
       write(*,*) V1

       write(*,*) "------------------------------"
       V2 = fun_V2(nd,qtot,cvec,Bmat,Y0,X2,X0)
       write(*,*) "V2: "
       do i = 1,nd+1
         write(*,*) V2(i,:)
       end do

      end subroutine

!.....Check Normalization...............................................

      subroutine check_norm(nd,cvec)
      ! nd: bath dimensions
      ! cvec: vector of the coefficients
       integer :: nd
       complex*16, dimension(nh),intent(in) :: cvec

       integer :: i,j
       real*8 :: Nout
       real*8, dimension(nd+1,nd+1) :: Bmat
       real*8, dimension(nh,nh) :: S00M

       do i = 1,nd+1
          Bmat(i,i) = i
          do j = i+1,nd+1
             Bmat(i,j) = j/20.d0 !Gershgoring circle theorem
             Bmat(j,i) = Bmat(i,j)
          end do
       end do

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_NORM"
       write(*,*) "B matrix"
     
       do i = 1,nd+1
          write(*,*) Bmat(i,:)
       end do

       write(*,*) "------------------------------"
  
       do i = 1,10
         call normalization(nd,dfloat(i),cvec,Bmat,S00M,Nout)

         write(*,*) "Norm at q = ", dfloat(i)
         write(*,*) Nout
       end do

      end subroutine

!.....Check shifted overlap.............................................

      subroutine check_shiftedoverlap(nd)
      ! nd: dimension of y
       integer, intent(in) :: nd

       integer :: i,j
       real*8 :: Nsq,NiNj
       real*8,dimension(nd+1) :: qi,qj,ppi,pj 
       complex*16,dimension(nd+1,nd+1) :: Bimat,Bjmat
       complex*16,dimension(nh,nh) :: TauMat

       complex*16 :: Sb

       do i = 1,nd+1
          Bimat(i,i) = i+iu*i
          do j = i+1,nd+1
             Bimat(i,j) = (j+iu*j)/20.d0 !Gershgoring circle theorem
             Bimat(j,i) = Bimat(i,j)
          end do
       end do

       do i = 1,nd+1
          Bjmat(i,i) = (i+2)*(1+iu)
          do j = i+1,nd+1
             Bjmat(i,j) = (j+iu*j)/5.d0 !Gershgoring circle theorem
             Bjmat(j,i) = Bjmat(i,j)
          end do
       end do

       write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
       write(*,*) "# HELLO I'M CHECK_NORM"
       write(*,*) "Bi matrix"
     
       do i = 1,nd+1
          write(*,*) Bimat(i,:)
       end do

       write(*,*) "Bj matrix"
     
       do i = 1,nd+1
          write(*,*) Bjmat(i,:)
       end do

       write(*,*) "------------------------------"

       Nsq = fun_Nsq(nd+1,real(Bimat))
       NiNj = fun_NiNj(nd+1,real(Bimat),real(Bimat))
       write(*,*) "Nsq: ",Nsq,"NiNj: ",NiNj
       NiNj = fun_NiNj(nd+1,real(Bimat),real(Bjmat))
       write(*,*) "Nsq: ",Nsq,"NiNj: ",NiNj

       qj(:) = 1.5d0
       pj(:) = 2.d0
       qi(:) = 1.7d0
       ppi(:) = 1.5d0

       qj(1) = 1.d0
       qi(1) = 0.7d0 

       write(*,*) "------------------------------"
       
       Sb = fun_Sb(nd,1.5d0,qi,qj,ppi,pj,Bimat,Bjmat)
       write(*,*) "Sb: ", Sb

       write(*,*) "------------------------------"

       TauMat = int_TauMat(nd,qi,qi,ppi,ppi,Bimat,Bimat)
 
       write(*,*) "TauMat(1,1): ", TauMat(1,1)

       TauMat = int_TauMat(nd,qi,qj,ppi,pj,Bimat,Bjmat)
 
       write(*,*) "TauMat(1,1): ", TauMat(1,1)
    

      end subroutine

!.....Check SG maple....................................................
!.....maple file: TGWP-K-and-S-numeric_GB.mw............................
      subroutine check_KSnum(nd)
      ! nd: dimension of y
       integer, intent(in) :: nd

       integer :: i,j
       real*8 :: q,p,H,fact1,fact2,norm,normdev,pre,Nout
       real*8, dimension(nd) :: qvec,pvec
       real*8, dimension(nd+1) :: qtot,ptot,dq,dp,qi,ppi
       complex*16, dimension(nd+1) :: cvec
       complex*16, dimension(nd+1,nd+1) :: Bmat,Bimat
       real*8 :: a,Y0
       real*8, dimension(nd) :: avec
       real*8, dimension(nd,nd) :: Amat,Tmat,LambdaMat
       real*8, dimension(nh,nh) :: X4,X3,X2,X1,X0,S00M

       complex*16, dimension(nh,nh) :: Hmat,Kmat,Tau
 
       do i = 1,nd+1
          qtot(i) = i*dsqrt(2.d0)/3.d0
          ptot(i) = i*dsqrt(3.d0)/7.d0
          Bmat(i,i) = (i+i)*(1+i/100.d0) + iu*(i+i)*(1+i/40.d0)/10.d0
          do j= i+1,nd+1
            Bmat(i,j) = (i+j)/40.d0 + iu*(i+j)/30.d0
            Bmat(j,i) = Bmat(i,j)
          end do   
       end do

       q = qtot(1)
       p = ptot(1)
       qvec = qtot(2:nd+1)
       pvec = ptot(2:nd+1)

       write(*,*) "qtot: ", qtot
       write(*,*) "ptot: ", ptot
       write(*,*) "Bmat: "
       do i = 1,nd+1
         write(*,*) Bmat(i,:)
       end do
       call extractA(nd,real(Bmat),Amat,avec,a)
       call diagonalization(nd,Amat,LambdaMat,Tmat)

       ! Integrals
       Y0=int_Y0(nd,LambdaMat)
       X4=int_XnMat(nd,4,a,avec,Amat,qtot(1))
       X3=int_XnMat(nd,3,a,avec,Amat,qtot(1))
       X2=int_XnMat(nd,2,a,avec,Amat,qtot(1))
       X1=int_XnMat(nd,1,a,avec,Amat,qtot(1))
       X0=int_XnMat(nd,0,a,avec,Amat,qtot(1))


       write(*,*) "Hermite pol. in x=1"
       do i = 1,nh
         H = herm_pol(i,1.d0,q,real(Bmat(1,1)))
         write(*,*) H
       end do

       write(*,*) "------------------------------"
       write(*,*) "The values of Hmat in x=1"
       write(*,*) "------------------------------"

       Hmat = fun_Hmat(1.d0,q,real(Bmat(1,1)))
       do i = 1,nh
         write(*,*) Hmat(i,:)
       end do

       write(*,*) "------------------------------"
       write(*,*) "The values of H derv in x=1"
       write(*,*) "------------------------------"

       write(*,*) "SECOND DERIVATIVE of H4"

       fact1 = factorial(4) ! H4
       fact2 = factorial(2) ! H2
       norm = 1/dsqrt(fact1*2**4) ! H4
       normdev = 1/dsqrt(fact2*2**2) ! H2
       pre = real(Bmat(1,1))*norm/normdev
       H = herm_pol(5-2,1.d0,q,real(Bmat(1,1)))
       write(*,*) "No function (H and pre)"
       write(*,*) 4.d0*(5-1)*(5-2)*H*pre, pre
       write(*,*) "Function (H and pre)"
       pre = der_pre(5,2,real(Bmat(1,1)))
       write(*,*) 4.d0*(5-1)*(5-2)*H*pre, pre

       write(*,*) "FIRST DERIVATIVE of H2"

       fact1 = factorial(2) ! H2
       fact2 = factorial(1) ! H1
       norm = 1/dsqrt(fact1*2**2) ! H2
       normdev = 1/dsqrt(fact2*2**1) ! H1
       pre = dsqrt(real(Bmat(1,1)))*norm/normdev
       H = herm_pol(3-1,1.d0,q,real(Bmat(1,1)))
       write(*,*) "No function (H and pre)"
       write(*,*) 2.d0*(3-1)*H*pre,pre
       write(*,*) "Function (H and pre):"
       pre = der_pre(3,1,real(Bmat(1,1)))
       write(*,*) 2.d0*(3-1)*H*pre,pre
      
       Kmat = kin_energy(nd,q,p,qvec,pvec,Bmat,Y0,X2,X1,X0)

       write(*,*) "qtot: ", qtot
       write(*,*) "ptot: ", ptot
       write(*,*) "Bmat: "
       do i = 1,nd+1
         write(*,*) Bmat(i,:)
       end do
       write(*,*) "Kmat: "
       do i = 1,nd+1
         write(*,*) Kmat(i,:)
       end do

       dq = [0.01,0.02,0.03]
       dp = [0.04,0.05,0.06]

       qi = qtot + dq
       ppi = ptot + dp

       write(*,*) dq
       write(*,*) dp
       write(*,*) sin(0.1)
       write(*,*) qi 
       write(*,*) ppi

       Bimat(:,:) = Bmat(:,:)*sin(0.1)

       write(*,*) "Bimat: "
       do i = 1,nd+1
         write(*,*) Bimat(i,:)
       end do

       Tau=int_TauMat(nd,qi,qtot,ppi,ptot,Bimat,Bmat) 

       write(*,*) "Tau: "
       do i = 1,nd+1
         write(*,*) Tau(i,:)
       end do
       
       cvec(:) = complex(1.d0,0.d0)
       
       call normalization(nd,q,cvec,real(Bmat),S00M,Nout)

      end subroutine  

!.....Coherent dynamics calculator......................................

      subroutine coherent_calc(nd,trj,q0,p0,masses,c0,Bcmplx)
      ! nd: bath dimension
      ! trj : trajectory parameters (first step, last step, nstep)
      ! q0 : initial gaussian center (x&y)
      ! p0 : initial gaussian momentum (x&y)
      ! masses : vector of the masses
      ! c0 : initial coefficients 
      ! Bcmplx : initial gaussian width matrix, real & imaginary
      integer, intent(in) :: nd
      integer*8, dimension(3), intent(in) :: trj
      real*8, dimension(nd+1), intent(in) :: q0,p0,masses
      complex*16,dimension(nh), intent(in) :: c0 
      complex*16,dimension(nd+1,nd+1), intent(in) :: Bcmplx 
      
      integer*8 :: i,j,first,last,nstep,k
      real*8 :: h,time,N,E,q,p,E0,Nsq
      real*8,dimension(nh) :: csq 
      real*8,dimension(nd) :: qvec,pvec 
      real*8,dimension(nd+1) :: qtot,ptot,w,phase,phase_c,phase_0 
      real*8,dimension(4) :: hvec = [0.5d0,0.5d0,1.d0,0.d0]
      complex*16,dimension(nd+1,nd+1) :: Bcoh,Bt,num,den
      
      complex*16,dimension(nh) :: cvec,cj,ctemp

      real*8,dimension(nd+1,4) :: kq,kp
      complex*16,dimension(nd+1,nd+1,4) :: kb

      real*8 :: energy,gbot,x,y,rx,ry
      real*8, dimension(nd) :: avec
      real*8, dimension(nd,nd) :: Amat,LambdaMat,Tmat
      real*8, dimension(nd+1,nd+1) :: Bmat
      real*8, dimension(nh,nh) :: S00M,invS,X0Mat

      complex*16 :: img,psi,totcorr
      complex*16, dimension(nd+1) :: A,B,C,correlation

      open(unit=777,file='coherent.dat',status='unknown')
      open(unit=778,file='final_coh_wfn.dat',status='unknown')
      open(unit=888,file='energy_BOT.dat',status='old')

      read(888,*)
      read(888,*)
      read(888,*)
      read(888,*)
      read(888,*)
      read(888,*)

      ! trajectory parameters
      first = trj(1)
      last = trj(2)
      nstep = trj(3)

      ! Separate x & y
      q = q0(1)
      qvec = q0(2:nd+1)
      p = p0(1)
      pvec = p0(2:nd+1)

      h = dfloat(last-first)/dfloat(nstep)

      Bcoh(:,:) =0.d0
      do i = 1,nd+1
         Bcoh(i,i) = dsqrt(masses(i))
      end do

      w(:) = dsqrt(1.d0/masses(:))

      write(777,*) '#Time ','q(1) ','B(1,1) ','|c|^2 ','phase'
      write(778,*) '#x ','y ', 'Re(Psi) ','Im(Psi) '

      read(888,*) energy

      !write(*,*) "I AM COHERENT"
 
     ! write(*,*) energy, 0.5d0*(w(1) + w(2))

      qtot(:) = q0(:)
      ptot(:) = p0(:)
      do i=1,nd+1
      phase_0(i)=0.5d0*(ptot(i)*qtot(i)-p0(i)*q0(i))+&
               &0.5d0*iu*log((iu*Bcmplx(i,i)*dsin(w(i)*time)+&
               &Bcoh(i,i)*dcos(w(i)*time))/Bcoh(i,i))
      end do

      !write(*,*) phase_0

      num(:,:) = 0.d0
      den(:,:) = 0.d0
      do j = 1,nstep
         time = j*h
         read(888,*) energy
         do i = 1,nd+1
           qtot(i)=q0(i)*dcos(w(i)*time)&
                   &+p0(i)*dsin(w(i)*time)/(masses(i)*w(i))
           ptot(i)=p0(i)*dcos(w(i)*time)&
                   &-masses(i)*w(i)*q0(i)*dsin(w(i)*time)
           num(i,i)=Bcmplx(i,i)*dcos(w(i)*time)+&
                    &iu*Bcoh(i,i)*dsin(w(i)*time)
           den(i,i)=iu*Bcmplx(i,i)*dsin(w(i)*time)+&
                    &Bcoh(i,i)*dcos(w(i)*time)
           phase(i)=0.5d0*(ptot(i)*qtot(i)-p0(i)*q0(i))+&
                    &0.5d0*iu*log((iu*Bcmplx(i,i)*dsin(w(i)*time)+&
                    &Bcoh(i,i)*dcos(w(i)*time))/Bcoh(i,i))
           phase_c(i)=0.5d0*(ptot(i)*qtot(i)-p0(i)*q0(i))-&
                    &0.5d0*w(i)*time
           Bt(i,i) = Bcoh(i,i)*num(i,i)/den(i,i)
           A(i) = (conjg(Bt(i,i)) + Bcmplx(i,i))/2.d0
           B(i)=Bcmplx(i,i)*q0(i)-iu*ptot(i)+iu*p0(i)+&
                &conjg(Bt(i,i))*qtot(i)
           C(i)=-iu*phase(i)+iu*phase_0(i)-iu*p0(i)*q0(i)&
           !C(i)=-iu*p0(i)*q0(i)&
                &+iu*ptot(i)*qtot(i)-(Bcmplx(i,i)*q0(i)*q0(i))/2.d0&
                &-(conjg(Bt(i,i))*qtot(i)*qtot(i))/2.d0
         correlation(i)=zsqrt(pi/A(i))*zexp(B(i)*B(i)/(4.d0*A(i))+C(i))&
                        &*(real(Bt(i,i))*real(Bcmplx(i,i))/pi**2)**0.25d0
         ! Check if A(i) is const with q0=p0=0
         !correlation(i)=zsqrt(pi/A(i))!*zexp(B(i)*B(i)/(4.d0*A(i))+C(i))&
                        !&*(real(Bt(i,i))*real(Bcmplx(i,i))/pi**2)**0.25d0
         ! Check if B(i) is const with q0=p0=0
         !correlation(i)=zexp(B(i)*B(i)/(4.d0*A(i)))!+C(i))&
         ! Check if C(i) is const with q0=p0=0
         !correlation(i)=zexp(C(i))
         end do
         gbot = -energy*time!&
      !       &+0.5d0*dot_product(ptot,qtot)-0.5d0*dot_product(p0,q0)
         totcorr =correlation(1)*correlation(2)
      !write(777,*) time,qtot,ptot,real(Bt(1,1)),real(Bt(2,2)),&
      !             phase(1) + phase(2),&!gbot,phase_c(1)+phase_c(2),&
                   !-0.5d0*(w(1)+w(2))*time,real(totcorr),aimag(totcorr)
      !             real(totcorr*conjg(totcorr))
      write(777,*) time,qtot(1),real(Bt(1,1)),&
                   real(totcorr*conjg(totcorr)),&
                   phase(1) + phase(2)!gbot,phase_c(1)+phase_c(2),&
                   !-0.5d0*(w(1)+w(2))*time,real(totcorr),aimag(totcorr)

      end do

      
      x = -5.d0
      y = -5.d0
      do i = 1,100
         x = x+0.1d0
         y = y+0.1d0
         rx = x - qtot(1)
         ry = y - qtot(2)
         img = iu*(ptot(1)*rx+ptot(2)*ry+phase(1)+phase(2))
         psi = zexp(-0.5d0*Bt(1,1)*rx*rx-0.5d0*Bt(2,2)*ry*ry+img)
         write(778,*) x,y,real(psi),aimag(psi)
      end do

      close(777)
      close(888)
      close(778)

      end subroutine

!.....Check projection scheme...........................................

      subroutine check_projection(nv,qeq,peq,ceq,Bcmplx)
       integer, intent(in) :: nv 
       real*8,dimension(nv+1) :: qeq !inital centers vector
       real*8,dimension(nv+1) :: peq !initial momenta vector
       complex*16,dimension(nh) :: ceq !initial coefficient vector
       complex*16,dimension(nv+1,nv+1),intent(in) :: Bcmplx !initial width matrix

       integer :: i,j,uuunit
       real*8 :: increment
       real*8,dimension(nv+1) :: q0 !equilibrium centers vector
       real*8,dimension(nv+1) :: p0 !equilibrium momenta vector
       complex*16,dimension(nh) :: c0 !equilibrium coefficient vector
       complex*16,dimension(nh) :: cref !reference coefficient vector

       cref(:) = ceq(:)
  
       increment = 0.05d0
       uuunit = 1000

       call plot_wfn(nv,qeq,peq,ceq,Bcmplx,1.d0,uuunit)
       uuunit=uuunit+1
 
       open(unit=2000,file="mae_BOT.dat",status='unknown')
       open(unit=2001,file="mae_BOT_fb.dat",status='unknown')
       open(unit=2002,file="mae_BOT_fbs.dat",status='unknown')
 
       do j = 1,10
         q0(:) = qeq(:) + j*increment
         p0(:) = peq(:) + j*increment
         write(*,*) "+------------------------------------------------+"
         write(*,*) "BOT Coefficients forward:"
         c0 = c_update(nv,qeq,peq,q0,p0,ceq,Bcmplx,Bcmplx) 
         do i = 1,nh
           write(*,*) c0(i) 
         end do
         write(uuunit,*) "#BOT forward", j*increment
         call plot_wfn(nv,q0,p0,c0,Bcmplx,1.d0,uuunit)
         uuunit=uuunit+1
         write(*,*) "+------------------------------------------------+"
         write(*,*) "BOT fb Coefficients forward:"
         c0 = c_update_fb(nv,qeq,peq,q0,p0,ceq,Bcmplx,Bcmplx) 
         do i = 1,nh
           write(*,*) c0(i) 
         end do
         write(uuunit,*) "#BOT fb forward", j*increment
         call plot_wfn(nv,q0,p0,c0,Bcmplx,1.d0,uuunit)
         uuunit=uuunit+1
         write(*,*) "+------------------------------------------------+"
         write(*,*) "BOT fbs Coefficients forward:"
         c0 = c_update_fbs(nv,qeq,peq,q0,p0,ceq,Bcmplx,Bcmplx) 
         do i = 1,nh
           write(*,*) c0(i) 
         end do
         write(uuunit,*) "#BOT fbs forward", j*increment
         call plot_wfn(nv,q0,p0,c0,Bcmplx,1.d0,uuunit)
         uuunit=uuunit+1
         write(*,*) "+------------------------------------------------+"
         write(*,*) "BOT Coefficients backward:"
         ceq = c_update(nv,q0,p0,qeq,peq,c0,Bcmplx,Bcmplx) 
         do i = 1,nh
           write(*,*) ceq(i) 
         end do
         !write(*,*) "Increment MAE"
         write(uuunit,*) "#BOT backward", j*increment
         write(2000,*) j*increment, abs(sum(cref-ceq))/dfloat(nh)
         call plot_wfn(nv,qeq,peq,ceq,Bcmplx,1.d0,uuunit)
         uuunit=uuunit+1
         write(*,*) "+------------------------------------------------+"
         write(*,*) "BOT fb Coefficients backward:"
         ceq = c_update_fb(nv,q0,p0,qeq,peq,c0,Bcmplx,Bcmplx) 
         do i = 1,nh
           write(*,*) ceq(i) 
         end do
         !write(*,*) "Increment MAE"
         write(uuunit,*) "#BOT fb backward", j*increment
         write(2001,*) j*increment, abs(sum(cref-ceq))/dfloat(nh)
         call plot_wfn(nv,qeq,peq,ceq,Bcmplx,1.d0,uuunit)
         uuunit=uuunit+1
         write(*,*) "+------------------------------------------------+"
         write(*,*) "BOT fbs Coefficients backward:"
         ceq = c_update_fbs(nv,q0,p0,qeq,peq,c0,Bcmplx,Bcmplx) 
         do i = 1,nh
           write(*,*) ceq(i) 
         end do
         !write(*,*) "Increment MAE"
         write(uuunit,*) "#BOT fbs backward", j*increment
         write(2002,*) j*increment, abs(sum(cref-ceq))/dfloat(nh)
         call plot_wfn(nv,qeq,peq,ceq,Bcmplx,1.d0,uuunit)
         uuunit=uuunit+1
       end do
       close(2000)
       close(2001)
       close(2002)
      end subroutine
              
     

      end module
