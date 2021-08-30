# ATPESC-codes

SUNDIALS+AMReX example codes for ATPESC 2019

## Advection-Diffusion Example

$$\frac{\partial u}{\partial t} + \vec{a} \cdot \nabla u -  \nabla \cdot ( D \nabla u ) = 0$$

where $$u$$ is the chemical concentration, $$\vec{a}$$ is the advection speed, and
$$D$$ is the diffusion coefficient.

### Problem Options

The most general set of problem options are shown below -- these apply to the general program
`Advection-Diffusion.exe` contained in `Advection-Diffusion.h` and `Advection-Diffusion.cpp`.
Each hands-on lesson will only accept an appropriate subset of these options for that lesson.


| Option               | Type   | Description                                        | Default  |
| ---------------------|--------|----------------------------------------------------|----------|
| `n_cell`             | `int`  | number of cells on each side of the square domain  | 256      |
| `max_grid_size`      | `int`  | max size of boxes in box array                     | 64       |
| `plot_int`           | `int`  | enable (1) or disable (0) plots                    | 0        |
| `stepper`            | `int`  | use CVODE (0) or ARKStep (1)                       | 0        |
| `cvode_method`       | `int`  | use BDF (0) or Adams (1) methods in CVODE          | 0        |
| `arkode_order`       | `int`  | ARKStep method order                               | 4        |
| `nls_method`         | `int`  | use Newton (0) or fixed-point (1) solver           | 0        |
| `nls_max_iter`       | `ìnt`  | maximum number of nonlinear iterations             | 3        |
| `nls_fp_acc`         | `int`  | number of fixed-point acceleration vectors         | 0        |
| `ls_max_iter`        | `int`  | maximum number of linear iterations                | 5        |
| `rhs_adv`            | `int`  | advection: disable (0), implicit (1), explicit (2) | 1        |
| `rhs_diff`           | `int`  | diffusion: disable (0), implicit (1), explicit (2) | 1        |
| `rtol`               | `Real` | relative tolerance                                 | 1e-4     |
| `atol`               | `Real` | absolute tolerance                                 | 1e-9     |
| `fixed_dt`           | `Real` | use a fixed time step size (if `fixed_dt` > 0.0)   | -1.0     |
| `tfinal`             | `Real` | final integration time                             | 1e4      |
| `dtout`              | `Real` | output frequency                                   | `tfinal` |
| `max_steps`          | `int`  | maximum number of steps between outputs            | 1000     |
| `write_diag`         | `int`  | output ARKStep diagnostics to a file               | 1        |
| `advCoeffx`          | `Real` | advection speed in the x-direction                 | 5e-4     |
| `advCoeffy`          | `Real` | advection speed in the y-direction                 | 5e-4     |
| `diffCoeffx`         | `Real` | diffusion coefficient in the x-direction           | 2e-5     |
| `diffCoeffy`         | `Real` | diffusion coefficient in the y-direction           | 2e-5     |
| `use_preconditioner` | `int`  | use preconditioning (1) or not (0)                 | 0        |

If preconditioning is enabled, then additional options may be set (see AMReX documentation of
the `MLMG` solver for descriptions):

| Option                    | Type   | Default  |
| --------------------------|--------|----------|
| mlmg.agglomeration        | `int`  | 1        |
| mlmg.consolidation        | `int`  | 1        |
| mlmg.max_coarsening_level | `int`  | 1000     |
| mlmg.linop_maxorder       | `int`  | 2        |
| mlmg.max_iter             | `int`  | 1000     |
| mlmg.max_fmg_iter         | `int`  | 1000     |
| mlmg.verbose              | `int`  | 0        |
| mlmg.bottom_verbose       | `int`  | 0        |
| mlmg.use_hypre            | `int`  | 1        |
| mlmg.hypre_interface      | `int`  | 3        |
| mlmg.use_petsc            | `int`  | 0        |
| mlmg.tol_rel              | `Real`  | 1.0e-6  |


## Building

The programs in this repository require at least SUNDIALS version 4.0, and AMReX version X.  These
can be installed locally via Spack using the commands
```bash
spack install sundials
spack install amrex@develop dimensions=2
```
To build the examples using the provided `Makefile` set the environment variables
`AMREX_INSTALL_DIR` and `SUNDIALS_INSTALL_DIR` to the directories for the AMReX
and Sundials installations respectively.  Also set the environment
variable `MPICXX` to the MPI C++ wrapper to use for compilation.  For
example, if AMReX was installed in `~/apps/amrex` and and Sundials in
`~/apps/sundials`, and the MPI C++ wrapper is just `mpicxx`, then the required
environment variables can be set with the commands
```bash
export AMREX_INSTALL_DIR=~/apps/amrex
export SUNDIALS_INSTALL_DIR=~/apps/sundials
export MPICXX=mpicxx
```
for sh/bash/ksh/zsh shells, or with
```tcsh
set AMREX_INSTALL_DIR ~/apps/amrex
set SUNDIALS_INSTALL_DIR ~/apps/sundials
set MPICXX mpicxx
```
for csh/tcsh shells.

If any of these environment variables are left unset, they default to
valid values for compilation on Cooley at the ATPESC 2019 workshop.

After the environment variables are set run `make` to build the ATPESC hands-on lesson
executables, `HandsOn1.exe`, `HandsOn2.exe`, and `HandsOn3.exe`.  The more general program,
`Advection-Diffusion.exe`, can be compiled with `make Advection-Diffusion.exe`.
