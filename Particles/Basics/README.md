# Basics ★

A first particle code: create a pure struct-of-arrays (SoA) particle
container with named components, add particles on the device, move them
through a periodic domain, and use reductions, plotfiles and checkpoints.

- Every particle stores `x y z vx vy vz w`; positions are the first
  `AMREX_SPACEDIM` real components. A component `dist` (distance traveled)
  is added at runtime.
- Each tile is filled with `nppc` random particles per cell in a
  `ParallelForRNG` kernel. Each process reserves a block of unique ids.
- Particles move in straight lines. After each step, `Redistribute` moves
  them to their new grid and MPI rank and wraps them around periodic
  boundaries.
- `ParticleReduce` sums the weight, kinetic energy and path length in one pass.
- A plotfile with the particle count per cell plus the particles, and a
  particle checkpoint that is read back with `Restart`.

The test checks that the number of particles, total weight and kinetic
energy are conserved, that the total path length is `nsteps*dt*sum|v|`, and
that the restarted particles match.

**Key APIs:** `amrex::ParticleContainerPureSoA`, `SetSoACompileTimeNames`,
`AddRealComp`, `DefineAndReturnParticleTile`, `ParticleTileData`,
`amrex::SetParticleIDandCPU`, `ParIter`, `Redistribute`,
`amrex::ParticleReduce`, `Increment`, `WritePlotFile`, `Checkpoint`, `Restart`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Particles
cmake --build build --target Particles_Basics -j 8
cd build/Particles/Basics && mpiexec -n 2 ./Basics inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Use `n_cell`, `nppc` and `nsteps` to make the problem larger. Set
`particles.do_tiling = 1` for better OpenMP performance on CPUs.

## Expected output

```
Real components: x y z vx vy vz w dist
Initial: 65536 particles, total weight 65653.57918, kinetic energy 32840.52706
Final:   65536 particles, total weight 65653.57918, kinetic energy 32840.52706
Restart: 65536 particles, total weight 65653.57918
Relative errors: weight 8.865877785e-16, energy 2.215542278e-16, path length 3.696117777e-16
PASSED
```

The particles in `plt00020/particles` can be viewed with yt or ParaView.
