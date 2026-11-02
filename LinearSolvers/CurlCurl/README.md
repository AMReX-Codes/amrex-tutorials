# CurlCurl ★★★

Solve the curl-curl equation for an edge-centered vector field,

```
curl(alpha*curl(E)) + beta*E = f,
```

with `MLCurlCurl`. Equations like this arise in implicit time stepping of
Maxwell's equations and in frequency-domain electromagnetics.

- `E` is an `Array<MultiFab,3>`, always with three components. `E_d` lives
  on edges along direction `d`: cell-centered in `d`, nodal in the other
  directions. In 2D, `E_z` is nodal in both x and y.
- `LinOpBCType::Dirichlet` fixes the tangential components on the
  boundary, e.g., a perfect conductor. The values come from the initial
  guess passed to `solve`, here zero.
- `setScalars(alpha, beta)` sets constant coefficients. `setAlpha` (nodal)
  and `setBeta` (edge-centered) set variable ones.
- `prepareRHS` makes the right-hand side consistent on edges shared by
  several boxes.
- The solver is `MLMGT<Array<MultiFab,3>>`, MLMG templated on the field
  type.

The curl of a gradient is zero, so with a small `beta` the operator is
nearly singular. MLCurlCurl's smoother is designed for this. Try
`beta = 1.e-3`. The exact solution `E_d = prod_{e != d} sin(pi x_e)`
vanishes tangentially on the walls, and the error is checked to be
O(dx^2). Works in 1D, 2D and 3D.

**Key APIs:** `MLCurlCurl`, `setScalars`, `prepareRHS`,
`MLMGT<Array<MultiFab,3>>`, edge `BoxArray`s with `amrex::convert`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=LinearSolvers
cmake --build build --target LinearSolvers_CurlCurl -j 8
cd build/LinearSolvers/CurlCurl && mpiexec -n 2 ./CurlCurl inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
MLMG: Initial rhs               = 20.7392088
MLMG: Initial residual (resid0) = 20.7392088
MLMG: Final Iter. 11 resid, resid/bnorm = 3.530260528e-10, 1.702215625e-11
Ex: max error = 0.0007648012602
Ey: max error = 0.0007648012617
Ez: max error = 0.0007648012716
PASSED
```
