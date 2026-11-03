# Poisson ★★

Solve the Poisson equation `-lap(phi) = f` on `[0,1]^D` outside a spherical
obstacle with the cell-centered EB solver `MLEBABecLap`. The boundary
conditions are Dirichlet on the domain boundary and inhomogeneous Dirichlet
on the EB.

The exact solution is `phi = sin(pi x) sin(pi y) sin(pi z)`. The right-hand
side is evaluated at the cell centroids and the EB values at the EB
centroids, both available from the `EBFArrayBoxFactory`. The program prints
the max-norm error at the cell centers and fails if it exceeds `10 dx^2`.
The scheme is second order: the error drops by about 4 when `n_cell`
doubles. The plotfile `plt` has `phi`, `exact`, `error` and `vfrac`.

`MLEBABecLap` solves `(a alpha - b div(beta grad)) phi = rhs`. Here `a = 0`,
`b = 1` and `beta = 1`. Without `setEBDirichlet` or `setEBHomogDirichlet`,
the EB boundary condition is homogeneous Neumann. Multigrid needs the EB on
coarsened grids, so `EB2::Build` is called with a large
`max_coarsening_level`.

**Key APIs:** `EB2::SphereIF`, `EB2::Build`, `EBFArrayBoxFactory`,
`EBData_t::centroid`, `EBData_t::bndrycent`, `MLEBABecLap`,
`MLEBABecLap::setEBDirichlet`, `MLMG`, `EB_WriteSingleLevelPlotfile`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=EB
cmake --build build --target EB_Poisson -j 8
cd build/EB/Poisson && mpiexec -n 2 ./Poisson inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
MLMG: Final Iter. 12 resid, resid/resid0 = 1.220027457e-08, 4.927121892e-12
...
Max error: 0.005236049738 (tolerance 0.009765625)
PASSED
```
