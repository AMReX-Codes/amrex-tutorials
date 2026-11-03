# WaveAmrLevel ★★★

An AMR code built on `Amr` and `AmrLevel` that solves the scalar wave
equation $\partial^2 u/\partial t^2 = \nabla^2 u$ as the first-order system

$$\frac{\partial u}{\partial t} = v, \qquad \frac{\partial v}{\partial t} = \nabla^2 u,$$

with a fourth-order Laplacian and Runge-Kutta time stepping (RK3 or RK4)
through `AmrLevel::RK`. The domain is $[0,1]$ in x, with $u = 0$ at
$x = 0$ and $1$, and periodic in y and z. The velocity $v$ starts as a pulse
in the middle that splits in two, reflects at the walls, and comes back;
the fine level follows the pulses. The solution does not depend on y and
z, and d'Alembert's formula gives the exact $v$ at any time, which the
program compares with the result.

**Key APIs:** `Amr`, `AmrLevel`, `LevelBld`, `StateDescriptor`,
`AmrLevel::RK`, `GpuBndryFuncFab`, `cell_quartic_interp`,
`FourthOrderInterpFromFineToCoarse`, `FillPatch`, `FillCoarsePatch`,
`TagBoxArray`

## Code structure

| File | Contents |
|------|----------|
| `AmrLevelWave.H` | The `AmrLevel` subclass, initial data and exact solution |
| `AmrLevelWave.cpp` | State setup, boundary function, time stepping, tagging, synchronization |
| `main.cpp` | `LevelBld`, the run, the error check |

- **State data.** `variableSetUp` defines one cell-centered `StateData`
  with the components `u` and `v`, two ghost cells for the stencil, and the
  fourth-order `cell_quartic_interp` interpolater for the coarse/fine
  boundary. This interpolater requires `amr.ref_ratio = 2`. The x
  boundaries use a `GpuBndryFuncFab` for the Dirichlet condition.
- **Time stepping.** `advance` calls `AmrLevel::RK` with a function that
  computes the right-hand side. `AmrLevel::RK` fills the ghost cells at
  every stage, including those at the coarse/fine boundary at the stage
  times, which on fine levels need data from the coarse level's RK stages
  (see `AMReX_RungeKutta.H` and `FillPatcher`).
- **Synchronization.** After each step, `post_timestep` replaces the coarse
  data under the fine level with `FourthOrderInterpFromFineToCoarse`. The
  code does not reflux, so the integral of $v$ is not conserved exactly
  across levels.
- **Tagging.** `errorEst` refines where $|v| > 0.25$.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=AMR
cmake --build build --target AMR_WaveAmrLevel -j 8
cd build/AMR/WaveAmrLevel && mpiexec -n 2 ./WaveAmrLevel inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Useful options: `wave.rk_order` (3 or 4), `wave.cfl`, `amr.n_cell`,
`amr.max_level`, `amr.v=1` to print every step, and
`amr.plot_files_output=1` to write plotfiles of `u` and `v`. The program
passes if the max error of $v$ is below $10^{-4}$.

## Expected output

```
Finished 320 level-0 steps at time 1 with 2 levels
Max error of v: 1.063265662e-05
PASSED
```
