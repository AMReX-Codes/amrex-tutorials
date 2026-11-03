# Sorting ★★

Organize particles within a tile: bin them by cell, sort them in memory,
partition them, and copy subsets or transformed copies to other containers.

1. `DenseBins` builds a counting sort of the particles of a tile by cell
   without moving them. `permutationPtr()` lists particle indices bin by
   bin and `offsetsPtr()[b]` is where bin `b` starts, so a kernel over
   cells can loop over the particles in each cell. Here it counts the
   particles per cell and compares with `ParticleContainer::Increment`.
2. `SortParticlesByCell` and `SortParticlesByBin` reorder the particle data
   so that particles in the same cell (or group of cells) are adjacent in
   memory, which helps cache and GPU memory access during deposition.
3. `partitionParticles` moves the particles that satisfy a predicate
   (here: speed above `vcut`) to the front of the tile.
4. `filterParticles` copies the particles that satisfy a predicate to
   another tile, and `transformParticles` copies particles while changing
   them (here: reversing the velocity).

The test checks each step: per-cell counts, sorted order, partition, the
number of filtered particles, and the reversed velocities.

**Key APIs:** `amrex::DenseBins`, `amrex::getParticleCell`,
`SortParticlesByCell`, `SortParticlesByBin`, `amrex::partitionParticles`,
`amrex::filterParticles`, `amrex::transformParticles`,
`DefineAndReturnParticleTile`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Particles
cmake --build build --target Particles_Sorting -j 8
cd build/Particles/Sorting && mpiexec -n 2 ./Sorting inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
Number of particles: 131072
DenseBins: max difference in particles per cell = 0
Out-of-order neighbors: initial 65496, after SortParticlesByCell 0, after SortParticlesByBin 0
partitionParticles: 68497 fast particles, 0 misplaced
filterParticles: 62575 slow particles
transformParticles: velocity sum 1123.82041 -> -1123.82041
PASSED
```

The counts depend on the random numbers, which differ between runs with
different numbers of MPI ranks or threads, and on GPUs.
