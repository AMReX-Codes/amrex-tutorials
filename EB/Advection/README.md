# Advection ★★★

Conservative finite-volume advection of a scalar around an embedded
boundary on a periodic domain, with redistribution for the small cut cells.

1. A uniform velocity is made discretely divergence free around a spherical
   obstacle with a MAC projection, solved with `MLEBABecLap`. Afterwards the
   sum of `area fraction * u` over the faces of every cut cell is zero, and
   nothing flows through the EB.
2. Fluxes are computed with second-order upwinding and limited slopes,
   dropping to first order next to the EB. The conservative update of a cut
   cell is `(sum of area fraction * flux) / (volume fraction * dx)`.
3. A small cut cell would limit the time step. `ApplyRedistribution` moves
   part of its update to its neighbors without changing the total. The time
   step comes from the full cell size, and the time integrator is two-stage
   SSP Runge-Kutta.

Choose `redistribution_type = StateRedist` (default), `FluxRedist`, or
`NoRedist`. In 3D, `NoRedist` is unstable with the default time step.

The state has two components: a Gaussian blob and a constant. The program
checks that there are no NaNs, that the total mass (sum of volume fraction
times state) is conserved to roundoff, and that the constant stays
constant. That last check holds only if the velocity is divergence free.

**Key APIs:** `EBFArrayBoxFactory`, `EBDataArrays`, `EB_computeDivergence`,
`MLEBABecLap`, `MLMG::getFluxes`, `ApplyRedistribution`, `EB_set_covered`,
`EB_average_face_to_cellcenter`, `EB_WriteSingleLevelPlotfile`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=EB
cmake --build build --target EB_Advection -j 8
cd build/EB/Advection && mpiexec -n 2 ./Advection inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

The CI run is short. To watch the blob go around the obstacle, run with
`n_cell = 64 stop_time = 1 plot_int = 10`.

## Expected output

```
max |div u| before projection: 881.403378
max |div u| after projection:  1.707217148e-11
Redistribution: StateRedist, dt = 0.006040948875, initial mass = 0.007724640317
Steps: 17, time = 0.1
...
Relative change in mass: 3.817692201e-15 (tolerance 2.220446049e-13)
Max deviation of the constant: 1.914024494e-13 (tolerance 1e-08)
PASSED
```
