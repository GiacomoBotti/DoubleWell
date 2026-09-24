! You need some constants in your life

MODULE constants

     implicit none
     save

     real*8, parameter :: pi = 4.d0*datan(1.d0)

     complex*16,parameter :: iu = complex(0.d0,1.d0)
     
END MODULE constants
