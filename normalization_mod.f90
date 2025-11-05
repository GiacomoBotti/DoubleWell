!**********************************************************************!
! Module containing everything pertaining normalization, both          !
! analytical and numerical                                             !
!**********************************************************************!

       module normalization_module

       use overlap_module
       use smatrix_module
       use integral_module

       implicit none
       private
       public :: normalization,normnum,normnumtot

       contains

!......Analitical Total Normalization...................................

       function normalization(nd,npar,y,work) result(Nout)
       ! nd: basis set dimension
       ! npar: parameter space dimension
       ! y: runge-kutta function - [parameters - coefficients]
       ! work: work array for I/O
         implicit none
         integer*8, intent(in) :: nd,npar
         complex*16, dimension(npar+nd), intent(in) :: y
         real*8, dimension(4), intent(in) :: work

         integer*8 :: i,j
         real*8 :: q,p,a,alpha,xi,zeta,Nout
         complex*16, dimension(nd) :: c, Sc
         complex*16, dimension(nd,nd) :: S00M

         ! extract work
         a = work(1)
         alpha = work(2)
         xi = work(3)
         zeta = work(4)

         ! extract y
         q = dreal(y(1))
         p = dreal(y(2))
         c(:) = y(npar+1:npar+nd)

         S00M = S00maple(alpha,q,p,xi,zeta)
   
         Sc = matmul(S00M,c)
         Nout = dreal(dot_product(c,Sc))
         
       end function

!......Numerical Electronic Normalization...............................

       function normnum(nd,npar,y,work) result(Nout)
       ! nd: basis set dimension
       ! npar: parameter space dimension
       ! y: runge-kutta function - [parameters - coefficients]
       ! work: work array for I/O
         implicit none
         integer*8, intent(in) :: nd,npar
         complex*16, dimension(npar+nd), intent(in) :: y
         real*8, dimension(4), intent(in) :: work

         integer*8 :: i,j
         real*8 :: q,p,a,alpha,xi,zeta,Nout
         !real*8, dimension(nd) :: vecxi,veczeta
         !real*8, dimension(nd) :: veczeta
         complex*16, dimension(nd) :: c,Sc 
         complex*16, dimension(nd,nd) :: S00mat

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
         !veczeta(1) = zeta
         !veczeta(2) = -1.d0*zeta
         
         S00mat = S00(nd,q,p,alpha,vecxi,veczeta)
 
         Sc = matmul(S00mat,c)
         Nout = dreal(dot_product(c,Sc))

       end function

!......Numerical Total Normalization....................................

       function normnumtot(nd,npar,y,work) result(Nout)
       ! nd: basis set dimension
       ! npar: parameter space dimension
       ! y: runge-kutta function - [parameters - coefficients]
       ! work: work array for I/O
         implicit none
         integer*8, intent(in) :: nd,npar
         complex*16, dimension(npar+nd), intent(in) :: y
         real*8, dimension(4), intent(in) :: work

         integer*8 :: i,j
         real*8 :: q,p,a,alpha,xi,zeta,Nout,Nel,Nnucl
         real*8,dimension(3) :: cond,nparam
         complex*16, dimension(nd) :: c
         complex*16, dimension(nd,nd) :: S00mat

         ! extract work
         a = work(1)
         alpha = work(2)
         xi = work(3)
         zeta = work(4)

         ! extract y
         q = dreal(y(1))
         p = dreal(y(2))
         c(:) = y(npar+1:npar+nd)

         Nel = normnum(nd,npar,y,work)

         cond = [-10.d0,10.d0,1000.d0]
         nparam = [a,q,p]

         Nnucl = nbraket1D(nwfn,nwfn,cond,nparam)
 
         Nout = Nnucl*Nel

       end function
       
       end module
