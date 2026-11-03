# SUNDIALS ★★★

Use the SUNDIALS (ARKODE) backend of `amrex::TimeIntegrator` for explicit,
implicit, implicit-explicit (IMEX) and multirate (MRI) time integration of an
advection-diffusion equation. The tutorial runs each method type in turn and
compares the number of steps, right-hand side evaluations and the error
against the exact solution.

**Key APIs:** `amrex::TimeIntegrator<MultiFab>` with `IntegratorTypes::Sundials`,
`set_rhs`, `set_imex_rhs`, `set_fast_rhs`, `set_time_step`,
`set_fast_time_step`, `set_adaptive_step`, `set_post_step_action`,
`integration.sundials.*` parameters

## The problem

On the periodic unit domain, solve

    dphi/dt = -u . grad(phi) + kappa lap(phi),   phi(t=0) = prod_d sin(2 pi x_d)

with centered differences, as in [RungeKutta](../RungeKutta). Diffusion is
stiff: explicit methods need `dt` below about `dx^2/(2 dim kappa)`. How each
method type uses the right-hand side:

| `integration.sundials.type` | Right-hand side | Step |
|-----------------------------|-----------------|------|
| `ERK`     | `set_rhs(advection + diffusion)`, explicit | adaptive, limited by stability |
| `DIRK`    | `set_rhs(advection + diffusion)`, implicit (Newton-GMRES) | `dt` |
| `IMEX-RK` | `set_imex_rhs(diffusion, advection)`: implicit, explicit | `dt` |
| `EX-MRI`  | `set_rhs(advection)` slow, `set_fast_rhs(diffusion)` fast | `dt`, fast `fast_dt` |

`methods` selects the types to run. Other SUNDIALS settings, such as
`integration.sundials.method` or `integration.sundials.fast_type = DIRK`,
are read from the inputs file (see the comments in `inputs` and the AMReX
documentation on Time Integration). `integration.rel_tol` and
`integration.abs_tol` control the adaptive step and the implicit solves.
Increase `n_cell` to make diffusion stiffer and the implicit, IMEX and MRI
methods more attractive.

## Build and run

AMReX must be built with SUNDIALS (>= 6). CMake, from the repository root
(SUNDIALS is downloaded and built unless `SUNDIALS_ROOT` is given):

```
cmake -S . -B build -DAMReXTutorials_TOPICS=TimeIntegration -DAMReX_SUNDIALS=ON
cmake --build build --target TimeIntegration_SUNDIALS -j 8
cd build/TimeIntegration/SUNDIALS && ./SUNDIALS inputs
```

GNU Make, with AMReX cloned next to this repository and SUNDIALS installed
(the `GNUmakefile` sets `USE_SUNDIALS = TRUE`; `SUNDIALS_DIR` or
`AMREX_SUNDIALS_HOME` gives its location). With `USE_MPI=TRUE`, SUNDIALS
must be built with MPI.

```
make -j8 SUNDIALS_DIR=/path/to/sundials/install
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

After the SUNDIALS statistics of each run (silence them with `amrex.verbose=0`):

```
method      steps  fast steps  RHS evals  relative error
ERK            49           0        199       2.646e-08
DIRK            5           0        236       6.983e-06
IMEX-RK         5           0        160       1.054e-06
EX-MRI          5          60        275       3.529e-04
PASSED
```
