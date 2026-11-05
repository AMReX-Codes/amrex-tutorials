# EBTracers ★★

Tracer particles in a flow around an embedded boundary. A scalar is carried
by particles: each particle has a fixed weight, and depositing the weights
on the mesh gives the scalar field. This is the Lagrangian counterpart of
[EB/Advection](../../EB/Advection), which solves for the same kind of
scalar on the mesh.

1. As in EB/Advection, a uniform velocity is made to flow around a sphere
   in a periodic domain with a MAC projection (`MLEBABecLap`). The
   velocity is kept at face centers, where `mac_interpolate` expects it,
   and set to zero on covered faces.
2. Particles are placed at random positions in every cell, and those inside
   the sphere are removed. Their weights sample a Gaussian blob upstream of
   the sphere.
3. Each step moves the particles with the midpoint method, interpolating
   the face velocity to their positions with `mac_interpolate`.
4. The interpolation does not know about the EB, so a few particles near
   the sphere end up inside it, up to about a cell deep. After each step,
   such particles are reflected back across the EB. The signed distance to
   the EB, filled on the nodes by `FillSignedDistance`, gives the depth
   (`-s`) and the normal (the gradient of `s`), and the particle moves by
   `2|s|` along the normal. This works for any EB geometry.
5. The weights are deposited on the mesh with cloud-in-cell weights by
   `ParticleToMesh`.

The particles move forward to `stop_time` and then back with the velocity
reversed. The program checks that the deposited weight equals the total
particle weight, and that the particles that were never reflected return
to their initial positions: the flow is steady, so they retrace their paths
up to the time stepping error, and the mean distance must be below 0.001
cells.

**Key APIs:** `ParticleContainerPureSoA`, `mac_interpolate`,
`FillSignedDistance`, `ParticleToMesh`, `ParticleInterpolator::Linear`,
`ReduceSum` over particles, `MLEBABecLap`, `MLMG::getFluxes`,
`EB_WriteSingleLevelPlotfile`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Particles
cmake --build build --target Particles_EBTracers -j 8
cd build/Particles/EBTracers && mpiexec -n 2 ./EBTracers inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Set `plot_int` to write plotfiles with the deposited scalar `phi`, the
velocity and the particles. Other options are `n_cell`, `nppc`, `radius`,
`velocity`, `stop_time` and `cfl`.

## Expected output (3D)

```
129246 particles, total weight 0.007691420379
32 steps forward and back with dt = 0.009375
Back at t = 0:
  particles reflected off the EB: 210
  mean distance of the other particles from their start: 3.902734996e-05 cells
  relative error of the deposited weight: 2.526048763e-14
PASSED
```

The particle positions are random, so the numbers differ slightly with the
number of MPI ranks.
