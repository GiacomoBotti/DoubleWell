!**********************************************************************!
! Module containing all the things needed with to work with the        !
! Hermite polynomials                                                  !
!**********************************************************************!

      module basisset_module

      use constants

      implicit none

      !Number of Hermite polynomials
      integer,parameter,public :: nh = 12 
      !Maximum order of x for the Hermite pol. in database
      integer, parameter, public :: max_x = nh+1  
      !Database of the coefficients
      real*8, dimension(nh,max_x), public :: Mherm

      private
      public :: GenHermMat,herm_pol,fun_Hmat,fun_HmatShift,factorial

      contains

!.....Factorial of a number.............................................

      function factorial(n) result(fact)
      ! n: input number
      ! fact: factorial of n
       implicit none
       integer, intent(in) :: n
       real*8 :: fact
       integer :: i

       fact = 1
       do i = 2, n
          fact = fact * i
       end do

      end function 

!.....Matrix of Hermite polynomials coefficients........................ 
  
      subroutine GenHermMat()
      ! Mherm: matrix of hermite polynomials coefficients
       integer :: i,j

       Mherm(:,:) = 0.d0

       ! The coefficients are listed for INCREASING ORDER of x power

       Mherm(1,1) = 1.d0
       if(nh.ne.1) then
       write(*,*) "Matrix of the Hermite coefficients generated"
          Mherm(2,2) = 2.d0 
       else
       write(*,*) "Only one hermite polynomial in use"
       end if

       do i = 3,nh
         Mherm(i,1) =-Mherm(i-1,2)
         do j = 2,max_x-1
           Mherm(i,j) = 2*Mherm(i-1,j-1) -j*Mherm(i-1,j+1)
         end do
       end do


      end subroutine

!.....Hermite Polynomial................................................
     
      function herm_pol(npol,x,q,alpha) result(Hofx)
      ! npol: order of the polynomial, from 1 to nh
      ! x: variable of the polynomial
      ! q: variable of the polynomial
      ! alpha: x coord gaussian width (REAL)
      ! Hofx: value of the polynomial in x-q
       integer, intent(in) :: npol
       real*8, intent(in) :: x,q,alpha
       
       integer :: i
       real*8 :: Hofx,fact,Hnorm

       Hofx = 0.d0
       do i = 1,max_x
         Hofx = Hofx + Mherm(npol,i)*(dsqrt(alpha)*(x-q))**(i-1)
         !Hofx = Hofx + Mherm(npol,i)*((x-q))**(i-1)
       end do

       fact = factorial(npol-1)

       Hnorm = 1/dsqrt(fact*2**(npol-1))
       !Hnorm = 1

       !write(*,*) "npol: ", npol, "factorial: ", fact

       Hofx = Hofx*Hnorm

      end function

!.....Matrix of Hermite Polynomial Products.............................
     
      function fun_Hmat(x,q,alpha) result(Hmat)
      ! x: variable of the polynomial
      ! q: variable of the polynomial
      ! alpha: x coord gaussian width (REAL)
      ! Hmat: matrix of the polynomial products in x-q
       real*8, intent(in) :: x,q,alpha
       
       integer :: i,j
       real*8 :: H1,H2
       real*8, dimension(nh,nh) :: Hmat

       Hmat(:,:) = 0.d0

       do i =1,nh
         do j =1,nh
           H1 = herm_pol(i,x,q,alpha)
           H2 = herm_pol(j,x,q,alpha)
           Hmat(i,j) = H1*H2
         end do
       end do
 
      end function 

!.....Matrix of Hermite Polynomial Products (Shifted)...................
     
      function fun_HmatShift(x,qi,qj,alphai,alphaj) result(HmatS)
      ! x: variable of the polynomial
      ! qi: variable of the polynomial
      ! qj: variable of the polynomial
      ! alphai: x coord gaussian width (REAL)
      ! alphaj: x coord gaussian width (REAL)
      ! Hmat: matrix of the polynomial products in x-q1 and x-q2
       real*8, intent(in) :: x,qi,qj,alphai,alphaj
       
       integer :: i,j
       real*8 :: H1,H2
       real*8, dimension(nh,nh) :: HmatS

       HmatS(:,:) = 0.d0

       do i =1,nh
         do j =1,nh
           H1 = herm_pol(i,x,qi,alphai)
           H2 = herm_pol(j,x,qj,alphaj)
           HmatS(i,j) = H1*H2
         end do
       end do
 
      end function 

      end module
