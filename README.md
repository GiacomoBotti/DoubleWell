# Double Well 

## Objective
The objective of this repository is to write a code for a MacLachlan variational principle dynamics of a double well potential ($x$ coordinate) coupled with an harmonic bath ($\mathbf{y}$ coordinate):

$$V = \frac{x^{4}}{16\eta} - \frac{x^{2}}{2} + V_{c} + \sum_{i} \frac{y_{i}^{2}}{2}$$

The coupling potential can be either a **Well and Chain** coupling,

$$V_{c} = \gamma xy_{1} + \gamma \sum_{i=2} y_{i}y_{i-1}$$

or an **Hub and Spoke** coupling,

$$V_{c} = \frac{\gamma}{2}\sum_{i=1}^{F_{B}} xy_{i}^{2}$$

The idea is to represent the double well mode with two Hermite-Gauss basis functions and the bath with a single multidimensional Gaussian. The evolution will be performed with the Heller-Karplus equations of motions. More informations of the Report USC, Chapt. V.

## NOTICE!

Since the matrix potential is initialized with the public parameter `nv`, if this parameter is changed a
```
make .clean
```
is required before recompiling with
```
make compile
```

## TO DO
- [ ] `check_mod.f90` at Report
- [ ] `quadratic_mod.f90` at Report
- [ ] `polynomials_mod.f90` at Report
- [ ] `ypowers_mod.f90` at Report
- [ ] $V_{0}$ polynomials subroutines
   - [ ] $\mathcal{Y_{0}}$ integrals
      - [x] (MAPLE) check $uWu$ sommation vs matrix moltiplication 
      - [x] (MAPLE) check $uZQ$ sommation vs matrix moltiplication 
      - [x] (MAPLE) check $QRu$ sommation vs matrix moltiplication 
      - [ ] function for $y_{1}^{4}$ in `ypowers_mod.f90`
      - [ ] function for $\tilde{\mathcal{P}}_{0}^{M}$ in `polynomials_mod.f90`
      - [x] (MAPLE) check $uWu$ w/ subs vs Fortran 
      - [x] (MAPLE) check $uZQ$ w/ subs vs Fortran 
      - [x] (MAPLE) check $QRu$ w/ subs vs Fortran 
      - [x] (MAPLE) check $QVQ$ w/ subs vs Fortran 
   - [ ] $\mathcal{Y_{1}}$ integrals
   - [ ] $\mathcal{Y_{2}}$ integrals
- [ ] Write the effective potential terms subroutines
- [ ] Write the time derivative subroutines
- [ ] Update the basis set module

## DONE
- [x] Copy relevant files from `Multicon_Fortran` repository
   - [x] Clean `Makefile`
- [x] Diagonalization subroutine
   - [x] check Diagonalization subroutine
- [x] Momenta subroutine
   - [x] check Momenta subroutine
   - [x] check `momcoeff` error
- [x] Matrix potential function
   - [x] check
- [x] `potential_mod.f90` at Report
