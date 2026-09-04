fc=gfortran
myflgs = -fno-automatic -O3 -ftree-vectorize -fcheck=bounds -ffast-math
lib = -llapack -lblas

obj= constants.o parameters_mod.o inversion_mod.o basisset_mod.o potential_mod.o matrix_mod.o\
     integrals_mod.o quadratic_mod.o effectivepot_mod.o normalization_mod.o eofmotion_mod.o kinetic_mod.o\
     observable_mod.o BOT_mod.o evolution_mod.o check_mod.o print.o main.o

compile: $(obj)
	$(fc) $(myflgs) $(obj) -o doublewell.x $(lib)

%.o: %.f
	$(fc) $(myflgs) -c $<

%.o: %.f90
	$(fc) $(myflgs) -c $<

.clean: 
	rm -f *.o *.x *.mod 
