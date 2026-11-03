# RungeKutta ★★

Advance an advection-diffusion equation with `amrex::TimeIntegrator` and its
native explicit Runge-Kutta methods, chosen at run time. The tutorial runs the
same problem with several time steps and prints the error and the observed
order of accuracy, which it checks against the order of the method.

**Key APIs:** `amrex::TimeIntegrator<MultiFab>`, `set_rhs`,
`set_post_step_action`, `set_time_step`, `evolve`, `integration.type`,
`integration.rk.type`, `ParmParse::queryAdd`

## The problem

On the periodic unit domain, solve

    dphi/dt = -u . grad(phi) + kappa lap(phi),   phi(t=0) = prod_d sin(2 pi x_d)

with second-order centered differences. The right-hand side function fills
the ghost cells of the state it is given and then evaluates the stencil.
For this initial condition the semi-discrete system has an exact solution
(each sine is advected at a slightly modified speed and decays at a modified
rate), so the measured error is the time integration error alone.

For AMR codes, AMReX also has built-in Runge-Kutta functions that fill the
coarse/fine boundary at the RK stage times; see the
[TimeIntegration docs](https://amrex-codes.github.io/amrex/tutorials_html/TimeIntegration.html)
and [AMR/WaveAmrLevel](../../AMR/WaveAmrLevel).

## Run-time parameters

| Parameter             | Default    | Meaning |
|-----------------------|------------|---------|
| `integration.type`    | RungeKutta | `ForwardEuler`, `RungeKutta`, or `SUNDIALS` (AMReX built with SUNDIALS) |
| `integration.rk.type` | 4          | 1 forward Euler, 2 trapezoid, 3 SSPRK3, 4 RK4, 0 user Butcher tableau (`integration.rk.nodes`, `weights`, `tableau`) |
| `dt`, `nruns`         | 0.01, 3    | step of the first run; each run halves it |
| `n_cell`, `stop_time` | 16, 0.25   | problem size and final time |

For a larger problem, raise `n_cell`. The explicit methods need
`dt` below about `dx^2/(2 dim kappa)`.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=TimeIntegration
cmake --build build --target TimeIntegration_RungeKutta -j 8
cd build/TimeIntegration/RungeKutta && mpiexec -n 2 ./RungeKutta inputs integration.rk.type=3
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
integration.type = RungeKutta, integration.rk.type = 4

        dt     steps  RHS evals      max error   rate
  1.000e-02      25        100  6.794e-07
  5.000e-03      50        200  4.212e-08   4.01
  2.500e-03     100        400  2.635e-09   4.00

Expected order of accuracy: 4
PASSED
```
