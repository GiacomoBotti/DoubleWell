!*********************************************************************!
! Module containing the observables functions                         !
!*********************************************************************!

       module observable_module

       use hamiltonian_module
       use hmatrix_module
       use smatrix_module
       use integral_module

       implicit none

       private
       public :: energy,energynum

       contains

!......Analytical Hamiltonian...........................................

       function energy(nd,npar,y,work) result(Hout)
       ! nd: basis set dimension
       ! npar: parameter space dimension
       ! y: runge-kutta function - [parameters - coefficients]
       ! work: work array for I/O
         implicit none
         integer*8, intent(in) :: nd,npar
         complex*16, dimension(npar+nd), intent(in) :: y
         real*8, dimension(4), intent(in) :: work

         integer*8 :: i,j
         real*8 :: q,p,a,alpha,xi,zeta,Hout
         complex*16, dimension(nd) :: c, Hc
         complex*16, dimension(nd,nd) :: H00M,T00M,V00M

         ! extract work
         a = work(1)
         alpha = work(2)
         xi = work(3)
         zeta = work(4)

         ! extract y
         q = dreal(y(1))
         p = dreal(y(2))
         c(:) = y(npar+1:npar+nd)

         T00M = T00maple(a,alpha,q,p,xi,zeta)
         V00M = V00maple(a,alpha,q,p,xi,zeta)
         H00M = T00M + V00M
   
         Hc = matmul(H00M,c)
         Hout = dreal(dot_product(c,Hc))
         
       end function

!......Numerical Hamiltonian............................................

       function energynum(nd,npar,y,work) result(Hout)
       ! nd: basis set dimension
       ! npar: parameter space dimension
       ! y: runge-kutta function - [parameters - coefficients]
       ! work: work array for I/O
         implicit none
         integer*8, intent(in) :: nd,npar
         complex*16, dimension(npar+nd), intent(in) :: y
         real*8, dimension(4), intent(in) :: work

         integer*8 :: i,j
         real*8 :: q,p,a,alpha,xi,zeta,Hout,Tn
         real*8 :: hm = 10.d0 !heavy mass 
         complex*16, dimension(nd) :: c, Hc
         complex*16, dimension(nd,nd) :: H00mat,T00mat,V00mat,S00mat
         complex*16 :: nover,ndyyover
         !real*8, dimension(nd) :: vecxi,veczeta
         !real*8, dimension(nd) :: veczeta
         real*8, dimension(3) :: cond,nparam 

         ! extract work
         a = work(1)
         alpha = work(2)
         xi = work(3)
         zeta = work(4)

         ! extract y
         q = dreal(y(1))
         p = dreal(y(2))
         c(:) = y(npar+1:npar+nd)
         
         ! build vecxi and veczeta
         !vecxi(1) = xi
         !vecxi(2) = -1.d0*xi
         !veczeta(1) = -1.d0*vecxi(1)*alpha/a 
         !veczeta(2) = -1.d0*vecxi(2)*alpha/a
         
         ! NUCLEAR OVERLAPS
         cond = [-10.d0,10.d0,1000.d0]
         nparam = [a,q,p]
         nover = nbraket1D(nwfn,nwfn,cond,nparam)

         ndyyover =nbraket1D(nwfn,dyynwfn,cond,nparam)
         Tn = -dreal(ndyyover/(2.d0*hm)) 

         ! ELECTRONIC OVERLAP
         S00mat = S00(nd,q,p,alpha,vecxi,veczeta)

         ! ELECTRONIC KINETIC MATRIX
         T00mat = T00(nd,q,p,alpha,vecxi,veczeta)

         ! TOTAL OVERLAP MATRICES
         ! WATCH OUT TO USE THE ELECTRONIC ONES
         ! SO FIRST COMPUTE THE ONES THAT YOU
         ! ARE NOT GOING TO USE
         T00mat = Tn*S00mat + nover*T00mat
         T00mat = dreal(T00mat) 

         V00mat = V00(nd,q,p,a,alpha,vecxi,veczeta)
         V00mat = dreal(V00mat)

         H00mat = T00mat + V00mat
   
         Hc = matmul(H00mat,c)
         Hout = dreal(dot_product(c,Hc))
         
       end function
       end module

