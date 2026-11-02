# ParticleMesh ★★

Particle-mesh interpolation: deposit a particle quantity onto the mesh and
gather a mesh field to the particles, with cloud-in-cell (linear) weights.

- Particles are placed randomly with `InitRandom`, each with the same mass.
- `ParticleToMesh` deposits the mass density onto a cell-centered
  `MultiFab`. The callable is run once per particle; it uses
  `ParticleInterpolator::Linear`, which adds the weighted value to the
  2^dim cells around the particle with atomic adds. Contributions to ghost
  cells are summed into the valid cells that own them, including across
  periodic boundaries.
- `MeshToParticle` interpolates a vector field from the mesh to the
  particles. The mesh field needs filled ghost cells.
- The particles and the mesh use different `BoxArray`s
  (`particle_max_grid_size` and `max_grid_size`). The two functions copy
  data between the layouts as needed.

The test checks that the deposited mass equals the total particle mass, and
that the gather reproduces a linear field exactly (to round-off).

**Key APIs:** `amrex::ParticleToMesh`, `amrex::MeshToParticle`,
`amrex::ParticleInterpolator::Linear`, `InitRandom`, `amrex::ReduceSum`,
`amrex::ReduceMax`, `WritePlotFile`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Particles
cmake --build build --target Particles_ParticleMesh -j 8
cd build/Particles/ParticleMesh && mpiexec -n 2 ./ParticleMesh inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
Number of particles: 100000
Total particle mass: 1
Deposited mass:      1
Relative mass error:   4.369837825e-13
Relative gather error: 6.401010382e-16
PASSED
```

`plt00000` holds the density `rho`, the field `Ex Ey Ez`, and the particles.

## Cloud-in-cell weights

In 1D, a particle at `x` lies between the cell centers `x_i` and `x_{i+1}`.
With `f = (x - x_i)/dx`, it gives a fraction `1-f` of its mass to cell `i`
and `f` to cell `i+1`; gathering uses the same weights. In 2D and 3D the
weights are products of the 1D weights.
