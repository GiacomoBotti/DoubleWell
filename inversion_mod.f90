!**********************************************************************!
! Module containing the matrix inversion functions                     !
!**********************************************************************!

      module inversion_module

      implicit none

      private
      public :: inv2D,invgen,linsys,invgen_real

      contains

!.....Inversion for a <2,2> matrix......................................

      function inv2D(A) result(invA)
      ! A : <2,2> matrix
        implicit none
        complex*16, intent(in), dimension(2,2) :: A

        complex*16, dimension(2,2) :: invA
        complex*16 :: determ

        determ = A(1,1)*A(2,2) - A(1,2)*A(2,1)

        if (determ.eq.0.d0) then
           stop "2x2 Matrix is singular"
        end if

        invA(1,1) = A(2,2)
        invA(2,2) = A(1,1)
        invA(1,2) = - A(1,2)
        invA(2,1) = - A(2,1)

        invA = invA/determ

      end function

!.....Inversion for a general COMPLEX matrix............................
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
      function invgen(npar,A) result(Ainv)
        implicit none
        integer*8, intent(in) :: npar
        complex*16, dimension(npar,npar), intent(in) :: A
        complex*16, dimension(npar,npar) :: Ainv

        complex*16, dimension(npar) :: work  ! work array for LAPACK
        integer, dimension(npar) :: ipiv   ! pivot indices
        integer :: n, info

        ! External procedures defined in LAPACK
        external ZGETRF
        external ZGETRI

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = npar

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call ZGETRF(n, n, Ainv, n, ipiv, info)

        if (info /= 0) then
          !write(*,*) "DGETRF info : ",info
          stop 'Invgen Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call ZGETRI(n, Ainv, n, ipiv, work, n, info)

        if (info /= 0) then
          !write(*,*) "DGETRI info : ",info
          stop 'Invgen Matrix inversion failed!'
        end if
      end function 

!.....Linear system for a general matrix................................
! Solve the linear system Ax=B by using the LU
! decomposition.  Depends on LAPACK.
      function linsys(npar,A,B) result(xout)
        implicit none
        integer, intent(in) :: npar
        complex*16, dimension(npar,npar), intent(in) :: A
        complex*16, dimension(npar),intent(in) :: B

        complex*16, dimension(npar) :: work,xout   !work array for LAPACK
        complex*16, dimension(npar,npar) :: Ainv
        integer, dimension(npar) :: ipiv   ! pivot indices
        integer :: n, info

        ! External procedures defined in LAPACK
        external ZGETRF
        external ZGETRS

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        ! Store B in xout to prevent it from being overwritten by LAPACK
        xout = B
        n = npar

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call ZGETRF(n, n, Ainv, n, ipiv, info)

        if (info /= 0) then
          !write(*,*) "DGETRF info : ",info
          stop 'Linsys Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call ZGETRS('N',n,1,Ainv,n,ipiv,xout,n, info)

        if (info /= 0) then
          !write(*,*) "DGETRS info : ",info
          stop 'Linsys Matrix inversion failed!'
        end if
      end function 

!.....Inversion for a general REAL matrix..............................
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
      function invgen_real(npar,A) result(Ainv)
        implicit none
        integer, intent(in) :: npar
        real*8, dimension(npar,npar), intent(in) :: A
        real*8, dimension(npar,npar) :: Ainv

        real*8, dimension(npar) :: work  ! work array for LAPACK
        integer, dimension(npar) :: ipiv   ! pivot indices
        integer :: n, info

        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = npar

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)

        if (info /= 0) then
          !write(*,*) "DGETRF info : ",info
          stop 'Invgen_real Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)

        if (info /= 0) then
          !write(*,*) "DGETRI info : ",info
          stop 'Invgen_real Matrix inversion failed!'
        end if
      end function 

      end module
