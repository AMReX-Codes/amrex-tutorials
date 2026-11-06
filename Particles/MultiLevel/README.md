# MultiLevel ★★★

Particles on an AMR hierarchy. A particle container defined on several
levels keeps each particle on the finest level that covers it. On a fixed
two-level hierarchy (level 1 refines the middle half of the domain in each
direction) the program initializes random particles, deposits their weights
on both levels, moves the particles across the coarse/fine boundary with
`Redistribute`, and regrids level 1, which is what an `AmrCore` code does
when its grids change.

- The container is constructed with the `Geometry`, `DistributionMapping`
  and `BoxArray` of every level and the refinement ratios. A code built on
  `AmrCore` uses `AmrParticleContainer` (or `AmrTracerParticleContainer`),
  which takes the hierarchy from the `AmrCore` and follows its regridding.
- `InitRandom` and `Redistribute` place each particle on the finest level
  whose grids cover it. Loops over particles go level by level, with
  `ParIter(pc, lev)`, and `Redistribute()` after a move handles the
  periodic boundaries, the grids, the MPI ranks and the levels at once.
- The multi-level `ParticleToMesh(pc, GetVecOfPtrs(rho), lev_min, lev_max,
  f)` deposits on every level and takes care of the coarse/fine boundaries,
  including the particles near them, and the coarse density under level 1
  is the average of the fine one. With the default `vol_weight = true`, the
  result is a density, so the integral of the level-0 density is the total
  weight.
- Regridding a level means giving the container the new `BoxArray` and
  `DistributionMapping` with `SetParticleBoxArray` and
  `SetParticleDistributionMap`, then calling `Redistribute`. A container
  built from an `AmrCore` takes the new grids from it and only needs the
  `Redistribute` after `regrid`.

After the regrid, the test checks that no level-0 particle lies under
level 1, that every level-1 particle lies in its region, and that the
deposited mass equals the total weight.

**Key APIs:** `ParticleContainerPureSoA(geom, dmap, ba, rr)`,
`NumberOfParticlesAtLevel`, `ReduceSum(pc, lev, f)`, `Redistribute`,
`ParIter(pc, lev)`, multi-level `ParticleToMesh`, `SetParticleBoxArray`,
`SetParticleDistributionMap`, `WriteMultiLevelPlotfile`, `WritePlotFile`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Particles
cmake --build build --target Particles_MultiLevel -j 8
cd build/Particles/MultiLevel && mpiexec -n 2 ./MultiLevel inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

The plotfile `plt_multilevel` holds the density on both levels and the
particles; see
[Visualization](https://amrex-codes.github.io/amrex/tutorials_html/VisualizationAndDebugging.html).

## Expected output

```
Initialized 100000 particles on 2 levels
  level 0: 87584 particles, 0 under level 1; level 1: 12416 particles, 12416 in its region
  fraction on level 1: 0.12416 (expected about 0.125)
Deposited mass on level 0 (fine data averaged down): 1
After 10 steps:
  level 0: 87550 particles, 0 under level 1; level 1: 12450 particles, 12450 in its region
After regridding level 1:
  level 0: 87475 particles, 0 under level 1; level 1: 12525 particles, 12525 in its region
Deposited mass on level 0: 1
PASSED
```

The counts depend on the random numbers, which differ between runs with
different numbers of MPI ranks.
