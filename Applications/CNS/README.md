# CNS (Blueprint)

A compressible Euler solver built on `amrex::Amr` and `amrex::AmrLevel`, with
adaptive mesh refinement, subcycling in time, refluxing and GPU support. It
uses a second-order Godunov method with characteristic slope limiting and an
approximate Riemann solver, advanced in time with Runge-Kutta. It runs in 2D
and 3D and comes with two problems: a Sod shock tube and a Rayleigh-Taylor
instability. Copy this directory to start your own `AmrLevel` code.

At the end of a run, the code checks that total mass changed only by
round-off, if no mass can enter or leave the domain.

**Key APIs:** `Amr`, `AmrLevel`, `LevelBld`, `StateDescriptor`
(`desc_lst`), `DeriveList` (`derive_lst`), `AmrLevel::RK`, `FluxRegister`,
`average_down`, `AMRErrorTag`, `GpuBndryFuncFab`, `ParmParse::queryAdd`,
`ParallelFor`, `ParReduce`

## Code structure

```
Source/                 the solver, shared by all problems
  main.cpp              Initialize, build Amr, time loop, mass check
  CNS.H, CNS.cpp        the AmrLevel subclass: one object per level
  CNS_setup.cpp         state data, BCs, derived variables, LevelBld
  CNS_advance.cpp       advance() and the flux computation
  CNS_hydro_K.H         GPU kernels: primitive variables, slopes, Riemann solver
  CNS_K.H               GPU kernels: time step, temperature
  CNS_derive.H/.cpp     derived plot variables (pressure, velocity)
  CNS_index.H           component indices
  CNS_parm.H/.cpp       physical parameters (EOS, gravity)
Exec/<Problem>/         one directory per problem
  cns_prob.H            ProbParm and the initial condition cns_initdata
  cns_prob.cpp          default parameters and reading of prob.*
  inputs, GNUmakefile, Make.package
```

`Amr` owns the grid hierarchy and calls the virtual functions of `CNS`:

| Function | Called |
|----------|--------|
| `variableSetUp` (static) | once, to define the state data and derived variables |
| `initData` | for each level at start-up |
| `computeInitialDt`, `computeNewDt` | on level 0, to choose `dt` for all levels |
| `advance` | for each level and each (sub)step |
| `post_timestep` | after a level has caught up with the next finer level: reflux and average down |
| `errorEst` | when regridding, to tag cells for refinement |
| `init(old)`, `init()` | when regridding, to fill new grids from old data or from the coarser level |
| `post_init`, `postCoarseTimeStep` | after initialization and after each coarse step |
| `restart` | when restarting from a checkpoint |

`AmrLevel` handles checkpoint, restart and plotfile I/O of the state data,
and `FillPatch`, which fills ghost cells from the same level, the coarser
level (interpolated in space and time) and the physical boundary conditions.

The state is registered in `CNS::variableSetUp` with `desc_lst.addDescriptor`
(cell-centered, `NUM_GROW` ghost cells, conservative interpolation) and
`desc_lst.setComponent` (names, `BCRec` boundary types, and the boundary
fill function). Derived quantities such as pressure are registered with
`derive_lst.add` and computed only when needed, e.g., for
`amr.derive_plot_vars`.

`CNS::advance` calls `AmrLevel::RK`, which runs the Runge-Kutta stages and
fills ghost cells for each stage. `compute_dSdt` computes the fluxes and
adds them to the `FluxRegister`s. In `post_timestep` the coarse level is
corrected with the fine fluxes (`Reflux`) and overwritten with the average of
the fine data (`average_down`), which together keep the scheme conservative.

## Making your own code

- **New problem:** copy `Exec/Sod` to `Exec/MyProblem`. Edit `ProbParm` and
  `cns_initdata` in `cns_prob.H`, and the defaults and `prob.*` parameters in
  `cns_prob.cpp`. Add it to the problem list in `CMakeLists.txt`.
- **New physics:** add components in `CNS_index.H` and their names and BCs in
  `CNS_setup.cpp`, and change the kernels in `CNS_hydro_K.H` and
  `compute_dSdt`. Inflow boundaries need values in `CnsFillExtDir`
  and removal of the check in `CNS::read_params`.
- **New AmrLevel code:** keep the structure of `main.cpp`, `CNS_setup.cpp`
  (`LevelBld`) and the `AmrLevel` overrides in `CNS.cpp`, and replace
  `advance`, the kernels and the tagging criteria. Kernels capture small
  parameter structs such as `Parm` by value; large tables belong in GPU
  memory, e.g., `Gpu::DeviceVector`.

## Build and run

CMake, from the repository root (targets `Applications_CNS_Sod` and
`Applications_CNS_RT`, binaries `CNS_Sod` and `CNS_RT`):

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Applications [-DAMReX_SPACEDIM=2]
cmake --build build --target Applications_CNS_Sod -j 8
cd build/Applications/CNS && mpiexec -n 2 ./CNS_Sod ../../../Applications/CNS/Exec/Sod/inputs
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
cd Exec/Sod
make -j8 [DIM=2] [USE_CUDA=TRUE]
mpiexec -n 2 ./CNS3d.gnu.MPI.ex inputs
```

Without an inputs file, each problem runs to its `stop_time` with the
defaults in `cns_prob_defaults`. The `inputs` files stop after 20 steps;
remove `max_step` and increase `amr.n_cell` and `amr.max_level` for a
real run. Boundary types (`cns.lo_bc`, `cns.hi_bc`) are 0 interior
(periodic), 2 outflow, 3 symmetry, 4 slip wall and 5 no-slip wall. Type 1
(inflow) is rejected until `CnsFillExtDir` supplies the inflow values.

## Expected output

```
STEP = 20 TIME = 0.0460694... DT = 0.0021319...
...
Relative change in total mass: 1.97e-15
```
