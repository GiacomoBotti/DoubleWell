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

Actually, two Hermite-Gauss basis functions are not enough. This requires a change of framework, in which the integral in $dx$ must be computed numerically. For this reason, I moved everything I did with two Hermite-Gauss basis function in a dead branch named `OLD_two_hermite`.

## NOTICE!

Since the matrix potential is initialized with the public parameter `nv`, if this parameter is changed a
```
make .clean
```
is required before recompiling with
```
make compile
```
The same goes for the Hermite polynomials basis set.

## NEW FRAMEWORK

### TO RUN

- [x] Uncoupled DW + HO
   - [x] Start at bottom
   - [x] Start at top
- [x] Numerical Stability tests
   - [x] No BOT
   - [x] BOT, only propagator on proj basis
   - [x] Full
- [ ] Coupled system:
   - [ ] QTAG
   - [ ] Increasing $N_{h}$
   - [ ] Quartic well
   - [ ] Cubic well
   - [ ] Ammonia inversion 

### TO DO
- [x] An easy way to do Gaussian-Average (Coalson) dynamics in the code
- [ ] Test it against analytical results (I need new SG, I think)
- [ ] Clean up the code
- [ ] Set `Nsq` back to normal
- [ ] `Makefile` at Report
- [x] `constants.f90` at Report
- [x] `inversion_module.f90` at Report
- [ ] `basisset_module.f90` at Report
- [ ] `potential_module.f90` at Report
- [ ] `matrix_module.f90` at Report
- [ ] `integrals_module.f90` at Report
- [ ] `kinetic_module.f90` at Report
- [ ] `normalization_module.f90` at Report
- [ ] `effectivepot_module.f90` at Report
- [ ] `eofmotion_module.f90` at Report
- [ ] `kinetic_module.f90` at Report
- [ ] `observable_module.f90` at Report
- [ ] `check_module.f90` at Report
- [ ] `main.f90` at Report
   

### DONE

- [x] Work out the integrals
- [x] Sketch the code
- [x] Hermite coefficient matrix
   - [x] Some smart way to code it
- [x] Check Hermite polynomials with MAPLE (`dw_check.mw`)
- [x] Hermite product matrix
- [x] Check Hermite product matrix with MAPLE (`dw_check.mw`)
- [x] Y0 integral
   - [x] subroutine
   - [x] check with MAPLE: determinant
   - [x] check with MAPLE: integral 
- [x] X integrals
   - [x] Gx
   - [x] check with MAPLE
   - [x] Xn 
   - [x] check with MAPLE
- [x] vector of coefficients
- [x] $[ V_{0}]_{ij}$
   - [x] $<V_{x}>$ portion
      - [x] $N^{2}$
      - [x] Coefficient Vector
      - [x] $\mathbf{c}^{T} X^{(0)} \mathbf{c}$
      - [x] compute 
      - [x] Test with MAPLE
   - [x] $<V_{xy}>$ portion 
      - [x] $\mathbb{A}^{-1}$
      - [x] Check w/ $\mathbb{A}^{-1}\mathbb{A} = \mathbb{I}$ 
      - [x] $\mathbb{A}^{-1} \mathbf{a}$
      - [x] Check $\mathbb{A}\mathbb{A}^{-1} \mathbf{a} =\mathbf{a}$
      - [x] compute 
      - [x] Check w/ MAPLE
   - [x] $<V_{y}>$ portion
      - [x] update $\mathbf{u}^{T} \mathbb{W} \mathbf{u}$
         - [x] check?
      - [x] update $\mathbf{q}^{T} \mathbb{V} \mathbf{q}$
         - [x] check?
      - [x] compute $\mathbb{V} \mathbb{A}^{-1} \mathbf{a}$
         - [x] check?
      - [x] compute $(\mathbb{A}^{-1} \mathbf{a})^{T} \mathbb{V} \mathbb{A}^{-1} \mathbf{a}$
         - [x] check?
      - [x] compute $\mathbf{q}^{T} \mathbb{V} \mathbb{A}^{-1} \mathbf{a}$
         - [x] check?
      - [x] compute 
   - [x] compute
- [x] $[ \mathbf{V}_{1}]_{ij}$
- [x] $[ \mathbb{V}_{2}]_{ij}$
- [x] Masses Matrix
- [x] Equations of motion in subroutine
- [x] $\mathbb{S}^{(00)}$
- [x] Adapt `normalization`
- [x] Adapt `total_normalization`
- [x] Adapt `bot_evo`
- [x] Minimal test: No hermite?
- [x] $\mathbb{H}^{(00)}$
   - [x] $\mathbb{T}^{(00)}$
      - [x] $d_{x} H_{i} d_{x} H_{j}$ integral
      - [x] check $d_{x} H_{i} d_{x} H_{j}$ integral
      - [x] check $d_{x} H_{i} d_{x} H_{j}$ vs MAPLE 
      - [x] $\partial_{x} \ln G d_{x} H_{i}$ integral
      - [x] check $\partial_{x} \ln G d_{x} H_{i}$ integral
      - [x] check $\partial_{x} \ln G d_{x} H_{i}$ vs MAPLE 
      - [x] $\vert \partial_{x} \ln G \vert^{2}$ integral
      - [x] check $\vert \partial_{x} \ln G \vert^{2}$ integral
      - [x] $\nabla_{y} \ln G \nabla_{y} \ln G$ integral 
      - [x] check $\nabla_{y} \ln G \nabla_{y} \ln G$ integral 
      - [x] make $\nabla_{y} \ln G \nabla_{y} \ln G$ real again 
      - [x] check $\nabla_{y} \ln G \nabla_{y} \ln G$ vs MAPLE 
      - [x] assemble everything
    - [x] assemble $\mathbb{H}^{(00)}$
    - [x] test $\mathbb{H}^{(00)}$
- [x] Debug $\mathbb{H}^{(00)}$ by checking step-by-step $\mathbb{T}^{(00)}$
- [x] Debug *dynamics* with $\tilde{\mathbb{B}}$ complex
- [x] Time shifted overlap matrix $\mathscr{S}_{mn} (t_{i},t_{j})$:
    - [x] Test in MAPLE
    - [x] Time shifted Hermite polynomials product function
    - [x] Debug w/ MAPLE
    - [x] Time shifted $N^{2} = N_{i} \times N_{j}$
    - [x] Debug w/ MAPLE using `check_module.f90`
    - [x] $S_{b}$ function
    - [x] Debug w/ MAPLE
    - [x] Adapt Cavalieri-Simpson
    - [x] Debug w/ MAPLE
- [x] Hermite normalization
   - [x] factorial
   - [x] Hermite norm
   - [x] scaled variables
   - [x] update the rest
   - [x] debug w/ MAPLE
- [x] DEBUG W/ MAPLE
- [x] Check w/ S.Garashchuk
   - [x] Hermite polynomials indeces
   - [x] Gaussian normalization
- [x] Adapt `bot`
   - [x] $\mathbb{T}(\tau,0)$ from S.G. docs
   - [x] Implement algorithm
- [x] Debug BOT
   - [x] check $dt =0$
   - [x] check with parameters' $dt = 0$
   - [x] Remove $N^{2} =1$ and $N_{i}N_{j} =1$
- [x] `\langle x \rnagle` and `\langle y \rangle` 
- [x] Create a `do` loop for print trajectory:

```
do i = 1,nstep_print
  do j = 1,10
    does evolution
  end do
  print
end do
```

- [x] **URGENT:** move the calculation of `Xn` integrals to `eofmotion_mod` and feed the integrals to `fun_V1` and `fun_V2` as input
- [x] **URGENT:** do the same for `energy`: move the integrals above and feed them to `fun_V0` and `kin_energy` (and then to `fun_Ka` and `fun_Kb`)
- [x] **URGENT:** feed them directly as $\mathbf{c}\mathcal{X}\mathbf{c}$
- [x] **URGENT:** test against a dynamics w/out this mod
- [x] Masking vector and matrix to "scale" dynamics of certain parameters
- [x] Reaction probability
- [x] Mirror correlation
- [x] Gif script

## OLD FRAMEWORK (TWO HERMITE)

### TO DO
- [ ] `check_mod.f90` at Report
- [ ] `quadratic_mod.f90` at Report
- [ ] `polynomials_mod.f90` at Report
- [ ] `ypowers_mod.f90` at Report
- [x] `tupowers_mod.f90` at Report
- [ ] $V_{0}$ polynomials subroutines
   - [x] $\mathcal{Y_{0}}$ integrals
      - [x] (MAPLE) check $uWu$ sommation vs matrix moltiplication 
      - [x] (MAPLE) check $uZQ$ sommation vs matrix moltiplication 
      - [x] (MAPLE) check $QRu$ sommation vs matrix moltiplication 
      - [x] function for $y_{1}^{4}$ in `ypowers_mod.f90`
         - [x] function for $(\mathbb{T}_{1,\cdot}\mathbf{u})$ in `tupowers_mod.f90`
         - [x] function for $(\mathbb{T}_{1,\cdot}\mathbf{u})^{2}$ in `tupowers_mod.f90`
         - [x] function for $(\mathbb{T}_{1,\cdot}\mathbf{u})^{3}$ in `tupowers_mod.f90`
         - [x] function for $(\mathbb{T}_{1,\cdot}\mathbf{u})^{4}$ in `tupowers_mod.f90`
      - [x] function for $\tilde{\mathcal{P}}_{0}^{M}$ in `polynomials_mod.f90`
      - [x] (MAPLE) check $uWu$ w/ subs vs Fortran 
      - [x] (MAPLE) check $uZQ$ w/ subs vs Fortran 
      - [x] (MAPLE) check $QRu$ w/ subs vs Fortran 
      - [x] (MAPLE) check $QVQ$ w/ subs vs Fortran 
   - [ ] $\mathcal{Y_{1}}$ integrals
      - [ ] dummy function for $\mathbb{T}_{1,\cdot}\mathbf{u} \cdot QVQ$
      - [ ] dummy function for $\mathbb{T}_{1,\cdot}\mathbf{u} \cdot uWu$
      - [ ] function for $\mathbb{T}_{1,\cdot}\mathbf{u} \cdot uZQ$
         - [ ] check with MAPLE
      - [ ] function for $\mathbb{T}_{1,\cdot}\mathbf{u} \cdot QRu$
         - [ ] check with MAPLE
      - [x] function for $y_{1}^{5}$ in `ypowers_mod.f90`
         - [x] function for $(\mathbb{T}_{1,\cdot}\mathbf{u})^5$ in `tupowers_mod.f90`
   - [ ] $\mathcal{Y_{2}}$ integrals
      - [x] function for $y_{1}^{6}$ in `ypowers_mod.f90`
         - [x] function for $(\mathbb{T}_{1,\cdot}\mathbf{u})^6$ in `tupowers_mod.f90`
- [ ] $\mathbf{V}_{1}$ polynomials subroutines
   - [ ] $\mathcal{Y_{0}}$ integrals
      - [x] function for $y_{1}^{3}$ in `ypowers_mod.f90`
   - [ ] $\mathcal{Y_{1}}$ integrals
   - [ ] $\mathcal{Y_{2}}$ integrals
- [ ] $\mathbb{V}_{2}$ polynomials subroutines
   - [ ] $\mathcal{Y_{0}}$ integrals
      - [x] function for $y_{1}^{2}$ in `ypowers_mod.f90`
   - [ ] $\mathcal{Y_{1}}$ integrals
   - [ ] $\mathcal{Y_{2}}$ integrals

### DONE
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
