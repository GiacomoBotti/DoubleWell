      subroutine print_banners(first,last,nstep,h,N)
      ! Prints banners
       implicit none
 
       integer, intent(in) :: first,last,nstep
       real*8, intent(in) :: h,N

      write(*,*) "+---------------------------------------------------+"
       write(*,*) "Writing trajectory output on trajectory_BOT.dat"
       write(*,*) "Writing coefficients output on coefficients_BOT.dat"
       write(*,*) "Writing bath positions on qbath_BOT.dat"
       write(*,*) "Writing bath momenta on pbath_BOT.dat"
       write(*,*) "Writing energy components on energy_BOT.dat"
       write(*,*) "Writing total phase on phase_BOT.dat"
       write(*,*) "Writing correlation function on correlation_BOT.dat"
       write(*,*) "Writing reaction prob on reaction_BOT.dat"
       write(*,*) "Writing cross correlation fun on crosscorr_BOT.dat"
       write(*,*) "Writing first momenta on momenta_BOT.dat"
      write(*,*) "+---------------------------------------------------+"

       write(321,*) "#Evolution parameters:"
       write(321,*) "#Range: ",first,last
       write(321,*) "#Steps: ",nstep
       write(321,*) "#Timestep: ",h
       write(321,*) "#Normalization constant: ",N
       write(321,*) "#Time ","N ","E ","q ","p ","Re B(1,1) ",&
                    &"Im B(1,1) ","Re B(2,2) ","Im B(2,2) ",&
                    & "|B(1,2)|^2 ","S00M(nh,nh) ","H00M(nh,nh) "

       write(322,*) "#Evolution parameters:"
       write(322,*) "#Range: ",first,last
       write(322,*) "#Steps: ",nstep
       write(322,*) "#Timestep: ",h
       write(322,*) "#Normalization constant: ",N
       write(322,*) "#Time ", "|c|^2"

       write(323,*) "#Evolution parameters:"
       write(323,*) "#Range: ",first,last
       write(323,*) "#Steps: ",nstep
       write(323,*) "#Timestep: ",h
       write(323,*) "#Normalization constant: ",N
       write(323,*) "#Time ", "qbath"
  
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
       write(325,*) "#Time ", "pbath"
  
       write(326,*) "#Evolution parameters:"
       write(326,*) "#Range: ",first,last
       write(326,*) "#Steps: ",nstep
       write(326,*) "#Timestep: ",h
       write(326,*) "#Normalization constant: ",N
       write(326,*) "#Time ", "phase"
  
       write(327,*) "#Evolution parameters:"
       write(327,*) "#Range: ",first,last
       write(327,*) "#Steps: ",nstep
       write(327,*) "#Timestep: ",h
       write(327,*) "#Normalization constant: ",N
       write(327,*) "#Time ", "correlation: re(C) & im(C) & |C|^2 & |C|"
  
       write(328,*) "#Evolution parameters:"
       write(328,*) "#Range: ",first,last
       write(328,*) "#Steps: ",nstep
       write(328,*) "#Timestep: ",h
       write(328,*) "#Normalization constant: ",N
       write(328,*) "#Time ", "reaction probability"
  
       write(329,*) "#Evolution parameters:"
       write(329,*) "#Range: ",first,last
       write(329,*) "#Steps: ",nstep
       write(329,*) "#Timestep: ",h
       write(329,*) "#Normalization constant: ",N
       write(329,*) "#Time ", "reaction probability"

       write(330,*) "#Evolution parameters:"
       write(330,*) "#Range: ",first,last
       write(330,*) "#Steps: ",nstep
       write(330,*) "#Timestep: ",h
       write(330,*) "#Normalization constant: ",N
       write(330,*) "#Time ", "<x> ", "<x>**2 - <x**2> ", "<y> ..."

       write(*,*) "First step:"
       write(*,*) "N ","E ","q ","p ","B(1,1) ",&
                    &"B(2,2) ", "B(1,3)"


      end subroutine print_banners
