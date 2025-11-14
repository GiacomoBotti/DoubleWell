# Double Well 

## Objective
The objective of this repository is to write a code for a MacLachlan variational principle dynamics of a double well potential ($x$ coordinate) coupled with an harmonic bath ($\mathbf{y}$ coordinate):

$$V = \frac{x^{4}}{16\eta} - \frac{x^{2}}{2} + V_{c} + \sum_{i} \frac{y_{i}^{2}}{2}$$

The coupling potential can be either a **Well and Chain** coupling,

$$V_{c} = \gamma xy_{1} + \gamma \sum_{i=2} y_{i}y_{i-1}$$

or an **Hub and Spoke** coupling,

$$V_{c} = \frac{\gamma}{2}\sum_{i=1}^{F_{B}} xy_{i}^{2}$$

The idea is to represent the double well mode with two Hermite-Gauss basis functions and the bath with a single multidimensional Gaussian. The evolution will be performed with the Heller-Karplus equations of motions. More informations of the Report USC, Chapt. V.

## TO DO
- [ ] Check every polynomial and generating sum in 3D w/ MAPLE
- [ ] $V_{0}$ polynomials subroutines
   - [ ] $\mathcal{Y_{0}}$ integrals
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
