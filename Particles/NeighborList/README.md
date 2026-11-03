# NeighborList ★★★

A small molecular dynamics code with short-range forces, built on a
neighbor particle container and neighbor lists.

- `MDParticleContainer` derives from `NeighborParticleContainerPureSoA`.
  `fillNeighbors` copies particles within one cell of each tile (including
  periodic images) into that tile, after its real particles.
  `buildNeighborList` then lists, for each real particle, the nearby real
  and neighbor particles that pass a `CheckPair` test.
- Particles are soft spheres of diameter `sigma` with the potential
  `U(r) = k/2 (sigma - r)^2` for `r < sigma`. The force kernel loops over
  `getNeighbors(i)` in a `ParallelFor`.
- Time stepping uses velocity Verlet. The list radius is one cell, larger
  than `sigma`, so the list stays valid for `num_rebuild` steps. In between,
  `updateNeighbors` refreshes the positions of the neighbor copies. At each
  rebuild, a local `Redistribute` moves particles to their new tiles.

The test checks that the total momentum is conserved to round-off (pair
forces are equal and opposite, also across tiles and MPI ranks) and that
the relative change in total energy is below 1%.

**Key APIs:** `amrex::NeighborParticleContainerPureSoA`, `fillNeighbors`,
`updateNeighbors`, `buildNeighborList`, `NeighborList::data`,
`getNeighbors`, `Redistribute`, `amrex::ParticleReduce`,
`amrex::RandomNormal`

## Build and run

CMake, from the repository root (2D or 3D):

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Particles
cmake --build build --target Particles_NeighborList -j 8
cd build/Particles/NeighborList && mpiexec -n 2 ./NeighborList inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Use `n_cell` and `nsteps` for a larger run. The time step is
`cfl*sigma/max_speed`, and `2*num_rebuild*cfl*sigma` must be smaller than
the skin `1 - sigma`.

## Expected output

```
Number of particles: 4096, dt = 0.004380251835
Step 20: kinetic 5710.24762, potential 436.2572452, total 6146.504865
...
Step 100: kinetic 5223.101324, potential 923.2618196, total 6146.363143
Relative momentum change: 2.895367347e-18
Relative energy change:   5.024214763e-05
PASSED
```
