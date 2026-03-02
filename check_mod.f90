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

       write(*,*) "------------------------------"
       Nsq = fun_Nsq(nd+1,Bmat)
       write(*,*) "N squared: ", Nsq

       write(*,*) "------------------------------"
       V0 = fun_V0(nd,qtot,cvec,Bmat)
       write(*,*) "V0(1,1) ", V0(1,1)
       
       write(*,*) "------------------------------"
       V1 = fun_V1(nd,qtot,cvec,Bmat)
       write(*,*) "V1: "
       write(*,*) V1

       write(*,*) "------------------------------"
       V2 = fun_V2(nd,qtot,cvec,Bmat)
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
         Nout = normalization(nd,dfloat(i),cvec,Bmat)

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
      
       Kmat = kin_energy(nd,q,p,qvec,pvec,Bmat)

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
       
       Nout = normalization(nd,q,cvec,real(Bmat))

        

      end subroutine  
      end module
