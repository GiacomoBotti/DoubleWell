fc=gfortran
myflgs = -fno-automatic -O3 -fcheck=bounds
lib = -llapack -lblas

#obj= constants.o inversion_mod.o potential_mod.o basisset_mod.o integral_mod.o integral_lha_mod.o\
     smatrix_mod.o overlap_mod.o hamiltonian_mod.o hmatrix_mod.o normalization_mod.o\
     timedev_mod.o dtmatrix_mod.o observable_mod.o BOT_mod.o evolution_mod.o check_mod.o main.o

obj= constants.o potential_mod.o diagonal_mod.o quadratic_mod.o polynomials_mod.o check_mod.o main.o

compile: $(obj)
	$(fc) $(myflgs) $(obj) -o doublewell.x $(lib)

%.o: %.f
	$(fc) $(myflgs) -c $<

%.o: %.f90
	$(fc) $(myflgs) -c $<

.clean: 
	rm -f *.o *.x *.mod 
