# StaticHierarchy ★★

The core multilevel operations of AMR on a fixed two-level hierarchy, with no
regridding and no `AmrCore`. Level 0 covers the unit square/cube, and level 1
refines its middle half by a factor of 2. A Gaussian diffuses across the
coarse/fine interface, and both levels take the same time step (no
subcycling). The domain has zero-gradient walls, so the total mass must stay
constant. The program checks this and fails if the relative change exceeds
round-off (1e-12 in double precision).

Each time step:

1. `FillPatchSingleLevel` fills the ghost cells of level 0, including the
   physical boundaries through `PhysBCFunct` with `BCType::foextrap`.
2. `FillPatchTwoLevels` fills the ghost cells of level 1 from level 1 where
   possible, and otherwise by conservative interpolation (`cell_cons_interp`)
   from level 0.
3. Each level computes fluxes F = -D grad(phi) and updates phi -= dt div(F).
4. The `FluxRegister` collects dt * area * F from the coarse side (`CrseInit`)
   and the fine side (`FineAdd`) of the coarse/fine boundary. `Reflux` corrects
   the coarse cells next to the fine level, so they see the fine fluxes.
5. `average_down` replaces the covered coarse cells with fine averages.

Run with `do_reflux=0` to see mass conservation fail.

**Key APIs:** `FillPatchSingleLevel`, `FillPatchTwoLevels`, `PhysBCFunct`,
`GpuBndryFuncFab`, `BCRec`, `cell_cons_interp`, `FluxRegister`,
`average_down`, `WriteMultiLevelPlotfile`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=AMR
cmake --build build --target AMR_StaticHierarchy -j 8
cd build/AMR/StaticHierarchy && mpiexec -n 2 ./StaticHierarchy inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8 DIM=2
mpiexec -n 2 ./main2d.gnu.MPI.ex inputs
```

Parameters: `n_cell` (coarse cells per direction), `max_grid_size`, `nsteps`,
`cfl`, `diff_coef`, `do_reflux`, `plot` (write `plt_static` at the end).

## Expected output (3D)

```
Initial mass: 0.005550189972
Final mass:   0.005550189972 after 100 steps, time = 0.003662109375
Relative change in mass: 1.719036495e-15
PASSED
```
