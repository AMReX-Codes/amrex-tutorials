# AdvectionAmrCore ★★★ / Blueprint

A complete adaptive mesh refinement code built on `AmrCore`. It advects a
scalar phi with a prescribed swirling velocity field. The grids follow phi as
it moves. It covers tagging, regridding, time stepping with and without
subcycling, refluxing, plotfiles, and checkpoint/restart. Use it as a
starting point for your own `AmrCore` application.

The flow reverses at t = 1, so phi returns to its initial position at
t = 2. The program checks that the total mass of phi stays constant to
round-off (1e-12 in double precision). Refluxing and conservative
interpolation make this possible.

**Key APIs:** `AmrCore` (`InitFromScratch`, `regrid`, `MakeNewLevelFromScratch`,
`MakeNewLevelFromCoarse`, `RemakeLevel`, `ClearLevel`, `ErrorEst`),
`AMRErrorTag`, `FillPatchSingleLevel`, `FillPatchTwoLevels`,
`InterpFromCoarseLevel`, `FluxRegister`, `average_down`, `ParReduce`,
`WriteMultiLevelPlotfile`, `VisMF`, `ParmParse::queryAdd`

## Code structure

| File | Contents |
|------|----------|
| `main.cpp` | Defaults for required parameters, the run, the mass check |
| `AmrCoreAdv.H/.cpp` | Level data, `AmrCore` callbacks, FillPatch, time stepping |
| `Advance.cpp` | Single-level advance, fluxes, time step estimate |
| `IO.cpp` | Plotfiles and checkpoint/restart |
| `Kernels.H` | Initial data, velocity field, limited upwind flux |

- **`AmrCore` callbacks.** `AmrCore` creates the grids and calls back into
  the class to fill new levels with data. `MakeNewLevelFromScratch` uses the
  initial data. `MakeNewLevelFromCoarse` interpolates from the level below.
  `RemakeLevel` copies data onto new grids. `ErrorEst` marks the cells that
  need refinement. Here one `AMRErrorTag` tags cells where
  `phi > adv.phierr[lev]`.
- **Time stepping.** `TimeStepWithSubcycling` advances level `lev` and then
  calls itself `ref_ratio` times for `lev+1`. Every `amr.regrid_int` steps it
  regrids the finer levels first. `TimeStepNoSubcycling` (`adv.do_subcycle=0`)
  regrids, then advances every level with a dt that is stable on all of them.
- **Single-level advance.** `AdvanceLevel` uses Heun's method (second-order
  Runge-Kutta). It computes upwind fluxes from MC-limited linear
  reconstructions. Before each stage, `FillPatch` fills the ghost cells.
  Fine-level ghost cells are interpolated in space and time from the level
  below.
- **Synchronization.** Each level passes dt * area * flux to the flux
  registers. After the fine steps, `Reflux` corrects the coarse cells next to
  the fine level. `average_down` then overwrites the covered coarse cells.

To make your own code:

- Change `Kernels.H` for the physics.
- Add components to `phi_new`/`phi_old` (and the `BCRec`s) for more
  variables.
- Use `BCType::ext_dir` with a functor in place of `FabFillNoOp` for inflow
  boundaries.
- Add `AMRErrorTag`s, or write your own loop in `ErrorEst`, for other
  refinement criteria.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=AMR
cmake --build build --target AMR_AdvectionAmrCore -j 8
cd build/AMR/AdvectionAmrCore && mpiexec -n 2 ./AdvectionAmrCore inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8 DIM=2
mpiexec -n 2 ./main2d.gnu.MPI.ex inputs
```

The inputs file stops after 10 steps. For the full problem, use
`max_step=1000`. To restart, use `amr.restart=chk00005`. Other useful
options are `adv.do_subcycle=0`, `adv.do_reflux=0` (mass is then not
conserved), `amr.max_level` and `amr.n_cell`. View the `plt*` plotfiles with
yt, ParaView, VisIt or AMReXplorer.

## Expected output (3D)

```
Coarse step 1 time = 0.01140770953 dt = 0.01140770953 levels = 3 mass = 0.007850917675
...
Coarse step 10 time = 0.1146053625 dt = 0.01155880142 levels = 3 mass = 0.007850917675
Relative change in mass: 8.83832208e-16
PASSED
```
