# DOUBLEWELL

`DOUBLEWELL` is a Fortran code for the variational propagation of a multiconfigurational Gaussian wavefunction on a quartic double-well system coupled to a harmonic bath. The wavefunction combines a Gaussian in the nuclear coordinates with a finite Hermite-polynomial expansion, and the Gaussian parameters and expansion coefficients are evolved using equations derived from the McLachlan variational principle.

The current source is configured for one active double-well coordinate coupled to `nv` bath coordinates. The active coordinate is treated with a quartic potential, while the bath is harmonic and can include nearest-neighbour coupling. The code evaluates the required Gaussian/Hermite integrals analytically for the bath coordinates and numerically for the system one, propagates the variational parameters, and writes observables including the norm, energy, Gaussian centers and widths, expansion coefficients, autocorrelation function, reaction probability, cross correlation, and coordinate moments.

## Requirements

The code requires:

- a Fortran compiler supporting Fortran 90;
- BLAS;
- LAPACK.

The supplied Makefile uses `gfortran` and links against `-llapack -lblas`.

## Compilation

Compile with:

```bash
make compile
```

This produces:

```text
doublewell.x
```

The Makefile currently uses:

```text
-fno-automatic -O3 -ftree-vectorize -fcheck=bounds -ffast-math
```

To remove compiled objects, module files, and executables, the supplied Makefile defines:

```bash
make .clean
```
## Source structure

The program is divided into modules responsible for the different parts of the variational calculation.

| File | Purpose |
|---|---|
| `main.f90` | Main program, input handling, initialization, basis projection, and launch of the dynamics. |
| `constants.f90` | Physical and numerical constants. |
| `parameters_mod.f90` | Model dimensions, potential parameters, grid parameters, and dynamics-control flags. |
| `potential_mod.f90` | Construction of the harmonic bath potential matrix and optional 2D potential output. |
| `basisset_mod.f90` | Hermite-polynomial basis and associated coefficient matrices. |
| `matrix_mod.f90` | Matrix manipulation, diagonalization, determinants, traces, and extraction of Gaussian-width blocks. |
| `inversion_mod.f90` | Matrix inversion and linear-system routines. |
| `integrals_mod.f90` | Gaussian/Hermite overlap, coordinate, transition, and probability integrals. |
| `quadratic_mod.f90` | Quadratic forms entering the Gaussian integrals. |
| `effectivepot_mod.f90` | Effective-potential matrix elements required by the variational equations. |
| `kinetic_mod.f90` | Active-coordinate and bath kinetic-energy matrix elements. |
| `normalization_mod.f90` | Wavefunction normalization and overlap matrix construction. |
| `observable_mod.f90` | Energy, coordinate moments, and wavefunction plotting. |
| `eofmotion_mod.f90` | Equations of motion and parameter propagators. |
| `BOT_mod.f90` | Coefficient propagation and basis-overlap transformations. |
| `evolution_mod.f90` | Main time-evolution driver and observable evaluation. |
| `check_mod.f90` | Diagnostic routines for matrices, integrals, normalization, projection, and related quantities. |
| `print.f90` | Output headers and run information. |

## Model dimensions

The principal compile-time dimensions are defined in `parameters_mod.f90`:

```fortran
integer, parameter :: nv = 4
integer, parameter :: nh = 9
integer, parameter :: max_x = nh + 1
integer, parameter :: maxorder = 8
```

`nv` is the number of bath coordinates. The complete Gaussian therefore has `nv + 1` coordinates: one active coordinate followed by the bath coordinates.

`nh` is the number of Hermite-polynomial basis functions used in the multiconfigurational expansion. Changing `nv` or `nh` requires recompilation.

## Potential

For the active coordinate `x`, the code uses a quartic double-well contribution of the form

```text
V(x) = x^4 / (16 eta) + sigma x^2 / 2
```

with linear coupling between the active coordinate and the first bath coordinate controlled by `gamma`.

For a one-dimensional bath, the potential written by `write_potential2D` is:

```text
V(x,y) = x^4/(16 eta) + sigma x^2/2 + gamma x y + kappa(1) y^2/2
```

For the general bath, `potential_mod.f90` constructs an `nv x nv` harmonic matrix with `kappa_const(i)` on the diagonal and `bath_const(i)` as nearest-neighbour bath couplings.

The potential parameters are supplied through the `pot_param` namelist.

## Input

The executable reads a file named:

```text
input
```

The namelists must appear in the order in which the program reads them:

```text
setup
pot_param
inp_mass
equilibrium
initial
grid
scaling_factors
```

A template input is:

```fortran
&setup
    trj = 0, 1, 100, 1, 0
    coalson = 0
    scaling = 0
    frozen = 0
    stationary = 0
/

&pot_param
    eta_const = 1.3544d0
    sigma_const = -1.d0
    gamma_const = 0.5d0
    kappa_const = 1.d0, 1.d0, 1.d0, 1.d0
    bath_const = 0.1d0, 0.1d0, 0.1d0, 0.d0
/

&inp_mass
    masses = 1.d0, 1.d0, 1.d0, 1.d0, 1.d0
/

&equilibrium
    qeq = 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
    peq = 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
    ceq = (1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0)
    Beq = ...
/

&initial
    q0 = 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
    p0 = 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
    Bcmplx = ...
/

&grid
    lwb = -12.d0
    hgb = 12.d0
    gstep = 500
/

&scaling_factors
    scalv = 0.d0, 1.d0, 1.d0, 1.d0, 1.d0
    scalm = ...
/
```

The dimensions of the arrays depend on the compile-time values of `nv` and `nh`. With the supplied settings, `qeq`, `peq`, `q0`, `p0`, and `masses` contain five elements, `Beq` and `Bcmplx` are `5 x 5` complex matrices, and `ceq` contains nine complex coefficients.

### `setup`

The `setup` namelist controls the trajectory and optional propagation modes:

```fortran
&setup
    trj = first, last, nstep, nprint, nback
    coalson = 0
    scaling = 0
    frozen = 0
    stationary = 0
/
```

The default value of `trj` is:

```text
0  1  100  1  0
```

The first three entries define the propagation interval and number of integration steps. The time step is computed internally as:

```text
h = (last - first) / nstep
```

`trj(4)` controls the number of internal propagation steps performed between successive output points. `trj(5)` controls the optional initial backpropagation used by the evolution routine.

The remaining flags modify the variational dynamics:

| Flag | Meaning |
|---|---|
| `coalson` | Enables the Gaussian-average/Coalson treatment when nonzero. |
| `scaling` | Enables selective scaling of the propagated Gaussian parameters. |
| `frozen` | Sets the Gaussian-width scaling matrix to zero, producing frozen-Gaussian width dynamics. |
| `stationary` | Freezes both centers/momenta and widths through the scaling vectors/matrices; the source notes that this currently applies to the SCP propagator. |

All four flags default to zero.

### `pot_param`

The potential namelist is:

```fortran
&pot_param
    eta_const = ...
    sigma_const = ...
    gamma_const = ...
    kappa_const = ...
    bath_const = ...
/
```

The defaults defined in the source are:

```text
eta_const   = 1.3544
sigma_const = -1.0
gamma_const = 0.5
kappa_const(:) = 1.0
bath_const(:)  = 0.1
bath_const(nv) = 0.0
```

`kappa_const` contains the harmonic bath force constants. `bath_const` contains the nearest-neighbour bath couplings.

### `inp_mass`

The mass vector has dimension `nv + 1`:

```fortran
&inp_mass
    masses = ...
/
```

The first element corresponds to the active coordinate and the remaining elements to the bath coordinates. All masses default to `1.d0`.

The mass vector is converted internally into the inverse-mass matrix used by the equations of motion.

### `equilibrium`

The `equilibrium` namelist defines the reference wavefunction:

```fortran
&equilibrium
    qeq = ...
    peq = ...
    ceq = ...
    Beq = ...
/
```

`qeq` and `peq` are the equilibrium Gaussian centers and momenta, `ceq` is the complex Hermite-basis coefficient vector, and `Beq` is the complex Gaussian width matrix.

Before reading the namelist, the program initializes:

```text
qeq = 0
peq = 0
ceq(1) = 1
ceq(2:nh) = 0
```

and sets `Beq` to a diagonal matrix with:

```text
Beq(i,i) = sqrt(masses(i))
```

The width matrix is complex (`complex*16`), so complex values can be supplied directly in the namelist.

### `initial`

The `initial` namelist defines the wavepacket from which the dynamics starts:

```fortran
&initial
    q0 = ...
    p0 = ...
    Bcmplx = ...
/
```

If these values are omitted, they are initialized from `qeq`, `peq`, and `Beq`.

After reading the initial parameters, the program does not simply copy the equilibrium coefficient vector. Instead, it projects the equilibrium expansion onto the initial Gaussian/Hermite basis through `c_update`, producing the actual initial coefficient vector `c0`.

### `grid`

The grid used to print the wavefunction is controlled by:

```fortran
&grid
    lwb = -12.d0
    hgb = 12.d0
    gstep = 500
/
```

`lwb` and `hgb` are the lower and upper grid bounds and `gstep` is the number of grid steps.

When `nv = 1`, the code additionally writes a two-dimensional representation of the model potential to `pot.dat`.

### `scaling_factors`

Selective propagation is controlled by:

```fortran
&scaling_factors
    scalv = ...
    scalm = ...
/
```

`scalv` has dimension `nv + 1` and controls the propagation of the center/momentum variables. `scalm` is an `(nv + 1) x (nv + 1)` matrix controlling the propagation of the Gaussian width matrix.

The internal default is one for all entries except:

```text
scalv(1)   = 0
scalm(1,1) = 0
```

These user-supplied values are activated only when `scaling = 1`. The `frozen` and `stationary` flags can override the width/parameter scaling.

## Wavefunction representation

The code represents the evolving state using:

- a Gaussian center vector `q`;
- a Gaussian momentum vector `p`;
- a complex Gaussian width matrix `B`;
- a vector of `nh` complex coefficients multiplying Hermite-polynomial basis functions.

The active coordinate is the first component of `q` and `p`; components `2:nv+1` are bath coordinates.

The equilibrium and initial Gaussian/Hermite bases need not be identical. The initial coefficient vector is obtained by projecting the equilibrium wavefunction onto the requested initial basis before propagation begins.

## Dynamics

The equations of motion for the Gaussian parameters are implemented in `eofmotion_mod.f90`. `KarplusTimeDer` evaluates the time derivatives, and the module contains several propagation schemes, including Runge-Kutta, self-consistent propagation, velocity-type propagation, PECE, and self-consistent PECE routines.

The current `bot_evo` production path propagates the coefficients and Gaussian parameters through `pece_coef`. Several alternative propagators remain present in the source and are commented at the call site, making it straightforward to switch the propagation scheme during method development.

Before the forward trajectory, `bot_evo` can perform an initial backward propagation controlled by `trj(5)`. The equilibrium wavefunction and propagated state are normalized using the nonorthogonal Hermite/Gaussian overlap matrix.

At each output step, the code recomputes the normalization, energy, autocorrelation, reaction probability, cross correlation, and coordinate moments.

## Gaussian averaging and constrained dynamics

Several switches allow modified dynamics without changing the equations elsewhere in the code.

### Gaussian-average / Coalson mode

Setting `coalson` to a nonzero value activates the Gaussian-average treatment. The integer value is also stored as `coalmode`.

### Scaled dynamics

With:

```fortran
scaling = 1
```

the vectors supplied through `scalv` and `scalm` are used to selectively evolve or freeze components of the Gaussian parameters.

### Frozen Gaussian

With:

```fortran
frozen = 1
```

the width scaling matrix is set to zero, preventing propagation of the Gaussian width. The code also sets the internal frozen-Gaussian factor `ffact`.

### Stationary-width/parameter mode

With:

```fortran
stationary = 1
```

both `scalvec` and `scalmat` are set to zero. The source currently labels this option as applicable to the SCP propagator.

## Outputs

The main evolution writes a set of plain-text files.

| File | Contents |
|---|---|
| `trajectory_BOT.dat` | Time, normalization, energy, active-coordinate center and momentum, selected width-matrix elements, and diagnostic overlap/Hamiltonian elements. |
| `coefficients_BOT.dat` | Time-dependent Hermite coefficients, including magnitudes and complex components. |
| `qbath_BOT.dat` | Bath-coordinate centers. |
| `pbath_BOT.dat` | Bath momenta. |
| `energy_BOT.dat` | Reserved energy-component output file. |
| `phase_BOT.dat` | Phases of the expansion coefficients. |
| `correlation_BOT.dat` | Autocorrelation function: real part, imaginary part, squared modulus, and modulus. |
| `reaction_BOT.dat` | Reaction-probability observable. |
| `crosscorr_BOT.dat` | Cross correlation with the reflected reference state. |
| `moments_BOT.dat` | Active-coordinate and bath coordinate moments. |
| `wfx_BOT.dat` | Wavefunction representation along the active coordinate. |
| `wfy_BOT.dat` | Wavefunction representation along the bath coordinate/grid used by `plot_wfn`. |

The output files contain header information describing the propagation range, number of steps, time step, and initial normalization where applicable.

### `trajectory_BOT.dat`

The trajectory file contains columns corresponding to:

```text
time
N
E/N
q(1)
p(1)
Re[B(1,1)]
Im[B(1,1)]
Re[B(2,2)]
Im[B(2,2)]
|B(1,2)|^2
S00M(nh,nh)
Re[H00M(nh,nh)]
```

### `correlation_BOT.dat`

The correlation output contains:

```text
time
Re[C(t)]
Im[C(t)]
|C(t)|^2
|C(t)|
```

The reference state is the initial propagated wavefunction.

### `crosscorr_BOT.dat`

The cross correlation is evaluated against a reflected equilibrium reference obtained by changing the sign of the equilibrium coordinate vector.

### `moments_BOT.dat`

The moment output contains the active-coordinate first moment, its variance-like quantity, and the corresponding bath-coordinate moments evaluated by the observable routines.

## Running

The executable expects both `input` and `banner.txt` in the working directory:

```text
doublewell.x
input
banner.txt
```

Run with:

```bash
./doublewell.x
```

At startup the program prints `banner.txt`, reads the input namelists, constructs the potential and Hermite matrices, initializes the mass matrix and reference wavefunction, projects the initial coefficients, prepares the dynamics options, and launches `bot_evo`.

A typical directory after a run contains:

```text
input
banner.txt
doublewell.x
trajectory_BOT.dat
coefficients_BOT.dat
qbath_BOT.dat
pbath_BOT.dat
energy_BOT.dat
phase_BOT.dat
correlation_BOT.dat
reaction_BOT.dat
crosscorr_BOT.dat
moments_BOT.dat
wfx_BOT.dat
wfy_BOT.dat
```

## Diagnostic routines

`check_mod.f90` contains a collection of development and validation routines for checking:

- Hessian/matrix diagonalization;
- momentum matrices;
- extraction of Gaussian width matrices;
- potential matrices;
- Hermite-polynomial matrices;
- Gaussian integrals;
- effective-potential terms;
- normalization;
- shifted overlaps;
- kinetic and overlap matrices;
- basis projection.

The corresponding calls are present but commented out in `main.f90`. They can be enabled individually when validating changes to the implementation.

## Notes

- `nv` and `nh` are compile-time parameters; changing them requires recompilation.
- The code assumes one active coordinate plus `nv` bath coordinates.
- Input arrays must have dimensions consistent with the compiled values of `nv` and `nh`.
- `Beq` and `Bcmplx` are complex matrices. Fortran namelist complex values should therefore be supplied using complex-number syntax such as `(1.0d0,0.0d0)`.
- The order of the namelists in `input` matters because the program reads them sequentially from the same open unit.
- `banner.txt` is required by the startup command `cat banner.txt`.
- The current production evolution path in `evolution_mod.f90` uses `pece_coef`; alternative propagators are retained in the source as commented calls.
- BLAS and LAPACK are required for the linear-algebra routines used throughout the code.
