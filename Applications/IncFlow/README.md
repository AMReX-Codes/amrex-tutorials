# IncFlow (Blueprint)

A small incompressible Navier-Stokes solver built on `amrex::AmrCore`, in the
spirit of [incflo](https://github.com/AMReX-Fluids/incflo) but with one
method for each step and no options. It solves

    du/dt + div(u u) = -grad(p) + nu lap(u),   div(u) = 0

with constant density, adaptive mesh refinement (all levels advance with
the same time step), and periodic or no-slip wall boundaries. It runs in 2D
and 3D, on CPUs and GPUs, and uses only AMReX: the projections are written
directly with the MLMG linear solvers. Copy this directory to start your
own `AmrCore` flow solver. It comes with three problems:

- `TaylorGreen`: decaying Taylor-Green vortices on one level, checked
  against the exact solution.
- `DoubleShearLayer`: two periodic shear layers that roll up into vortices
  (Kelvin-Helmholtz), with a fine level that follows the vorticity. The mean
  velocity must be conserved to round-off, also across levels.
- `LidDrivenCavity`: no-slip walls and a moving lid, refined near the lid.
  The velocity must stay bounded.

**Key APIs:** `AmrCore`, `FillPatchTwoLevels`, `InterpFromCoarseLevel`,
`GpuBndryFuncFab`, `MLPoisson`, `MLNodeLaplacian`, `MLABecLaplacian`,
`MLMG::solve`, `MLMG::getFluxes`, `MLMG::apply`, `computeDivergence`,
`average_down_faces`, `average_down`, `amrex_calc_xslope_extdir`,
`AMRErrorTag`, `ParmParse::queryAdd`, `ParallelFor`

## Algorithm

Each step is a predictor-corrector (Heun) method of lines:

1. `R(u) = -div(u u) + nu lap(u)` without the pressure gradient:
   - Extrapolate the cell-centered velocity to the faces with limited
     slopes, and upwind the normal component to get a face velocity `umac`.
   - **MAC projection:** solve `lap(phi) = div(umac)` with `MLPoisson` and
     subtract `grad(phi)`, so that `umac` is divergence free.
   - Upwind all velocity components with `umac` and take the divergence of
     the fluxes. Fine fluxes replace coarse ones at coarse/fine faces.
   - The viscous term is computed by `MLMG::apply` with an `MLABecLaplacian`,
     which handles walls and coarse/fine boundaries.
2. Predictor: `u* = u^n + dt R(u^n)`, then project.
3. Corrector: `u^{n+1} = u^n + dt/2 (R(u^n) + R(u*))`, then project.

**Nodal projection:** `u` is `u_div_free + dt grad(p)`, so solving
`div(grad(phi)) = div(u)` with `MLNodeLaplacian`, phi on the nodes, and
subtracting `grad(phi)` gives the new velocity and the pressure `p = phi/dt`.
Both projections solve all levels together.

## Code structure

```
Source/                 the solver, shared by all problems
  main.cpp              Initialize, InitData, Evolve, check the result
  IncFlow.H, .cpp       the AmrCore subclass: parameters, boundary conditions,
                        level creation (AmrCore callbacks), tagging, FillPatch
  Advance.cpp           time step, advection and viscous terms
  Advection_K.H         GPU kernels: slopes, face values, upwinding
  Projection.cpp        MAC and nodal projections
  IO.cpp                plotfiles, checkpoint and restart
Exec/<Problem>/         one directory per problem
  prob.H                initial velocity prob_initdata
  prob.cpp              default parameters and the check of the result
  inputs, GNUmakefile, Make.package
```

`AmrCore` owns the grid hierarchy and calls these functions of `IncFlow`:

| Function | Called |
|----------|--------|
| `MakeNewLevelFromScratch` | at start-up, to create a level and set its initial velocity |
| `MakeNewLevelFromCoarse` | when regridding creates a new level, filled from the coarser level |
| `RemakeLevel` | when regridding changes the grids of a level, filled from old data |
| `ClearLevel` | when regridding removes a level |
| `ErrorEst` | to tag cells where the vorticity exceeds `incflow.vort_tag` |

The velocity is the only state. Ghost cells are filled by `FillPatch` from
the same level, the coarser level, and the wall velocity (`ext_dir` boundary
type, set by `WallVelocityFill`). The pressure is only output.

## Making your own code

- **New problem:** copy `Exec/TaylorGreen` to `Exec/MyProblem`. Edit
  `prob_initdata` in `prob.H`, and the defaults and the check in `prob.cpp`.
  Add it to the problem list in `CMakeLists.txt`.
- **Wall velocity:** `xlo.velocity`, `yhi.velocity`, ... set the velocity
  of the walls in non-periodic directions; the normal component is zero.
- **More physics:** a scalar can be advected like a velocity component in
  `ComputeAdvection`. Forces such as gravity are added to `R` in
  `ComputeRHS`. Inflow and outflow need Dirichlet values in
  `WallVelocityFill` and Dirichlet pressure boundaries in the projections.

## Build and run

CMake, from the repository root (targets `Applications_IncFlow_<Problem>`,
binaries `IncFlow_<Problem>`):

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Applications [-DAMReX_SPACEDIM=2]
cmake --build build --target Applications_IncFlow_TaylorGreen -j 8
cd build/Applications/IncFlow && mpiexec -n 2 ./IncFlow_TaylorGreen ../../../Applications/IncFlow/Exec/TaylorGreen/inputs
```

CMake, standalone after copying this directory (AMReX is fetched, or use a
local checkout):

```
cmake -S . -B build [-DAMReX_SPACEDIM=2] [-DAMReX_GPU_BACKEND=CUDA] \
      [-DFETCHCONTENT_SOURCE_DIR_AMREX=/path/to/amrex]
cmake --build build -j 8
```

GNU Make, with AMReX cloned next to this repository:

```
cd Exec/TaylorGreen
make -j8 [DIM=2] [USE_CUDA=TRUE]
mpiexec -n 2 ./IncFlow3d.gnu.MPI.ex inputs
```

The `inputs` files stop after 10 steps; remove `max_step` for a longer
run, e.g., the shear layers roll up by `stop_time=1`, and `stop_time=20`
brings the cavity close to steady state. The nodal projection needs
`amr.ref_ratio = 2`. Restart with `amr.restart=chk00005`.

## Expected output

Taylor-Green vortices in 3D:

```
Step 10 time = 0.08138020833 dt = 0.008138020833 levels = 1 kinetic energy = 0.05493709093
Relative max error in velocity: 0.001081376713
PASSED
```

Double shear layer in 3D:

```
Step 10 time = 0.07762177993 dt = 0.007698007641 levels = 2 kinetic energy = 0.1077343874
Max |velocity| = 1.016523926, max |mean velocity| = 5.488773498e-17
PASSED
```

Lid-driven cavity in 3D:

```
Step 10 time = 0.08138020833 dt = 0.008138020833 levels = 2 kinetic energy = 0.005179961541
Max |velocity| = 0.4264155568, kinetic energy = 0.005179961541
PASSED
```
