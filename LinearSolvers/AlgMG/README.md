# AlgMG ★★★

AlgMG is AMReX's algebraic multigrid solver. It solves a sparse linear
system $A x = b$ and builds its coarse levels from the matrix itself, so it
works where geometric multigrid struggles, e.g., for strongly varying or
anisotropic coefficients. This tutorial uses it in two ways:

1. **Standalone**: assemble the matrix as an `SpMatrix` and solve with
   `AlgMG`.
2. **In MLMG**: solve the same problem with `MLABecLaplacian` and `MLMG`,
   using AlgMG as the bottom solver or as the multigrid type.

The problem is $-\nabla\cdot(\beta\nabla\phi) = f$ on the unit domain with
$\phi = 0$ on the boundary. The coefficient $\beta$ is a checkerboard of
blocks of `block` cells alternating between 1 and `jump` (default $10^4$).
The right-hand side is the discrete operator applied to
$\phi = \prod_d \sin(\pi x_d)$, so the solvers should recover that $\phi$
to within the tolerance.

**Key APIs:** `AlgPartition`, `AlgVector`, `SpMatrix::define`, `SpMV`,
`AlgMG`, `AlgMG::setKrylovSolver`, `MLMG::setMultigridType`,
`MultigridType::algebraic/hybrid`, `BottomSolver::algmg`,
`MLMG::setAlgMGOptions`, `MLMG::usedAlgMG`, `MLMG::setThrowException`

## Standalone AlgMG

`AlgPartition` splits the rows (one per cell) evenly among the MPI ranks,
and `AlgVector` holds a distributed vector. The matrix is assembled from
CSR arrays of the local rows: each row has room for the $2d+1$ entries of
the stencil, and missing neighbors at the domain boundary get column -1,
which `CsrValid{false}` tells `SpMatrix::define` to drop. AlgMG then solves
with V-cycles accelerated by conjugate gradients (`algmg.krylov`).

## AlgMG in MLMG

`MLMG` runs four variants on the same operator:

- `geometric`: the default geometric multigrid. It diverges on this
  problem, because the coarse grids cannot represent the coefficient.
- `geometric + AlgMG bottom`: `setBottomSolver(BottomSolver::algmg)`. AlgMG
  solves the bottom level, but the geometric levels above it are the
  problem, so this diverges too.
- `algebraic`: `setMultigridType(MultigridType::algebraic)`. AlgMG solves
  the whole level in every MLMG iteration.
- `hybrid`: `setMultigridType(MultigridType::hybrid)`. Geometric multigrid
  first; when it stalls or diverges, MLMG switches to AlgMG and restarts
  from the best iterate. `usedAlgMG()` reports the switch.

With `setThrowException(true)`, a failed solve throws instead of aborting.
`setAlgMGOptions` takes a callback that receives the `AlgMG` object MLMG
builds, to change its options; the tolerances and verbosity come from the
MLMG bottom-solver settings. Use `geometric` when it converges, since it is
cheaper, and `hybrid` when the coefficients may be difficult.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=LinearSolvers
cmake --build build --target LinearSolvers_AlgMG -j 8
cd build/LinearSolvers/AlgMG && mpiexec -n 2 ./AlgMG inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Try `jump=1`, where geometric multigrid converges, `jump=100`, or
`algmg.krylov=none` for plain V-cycles. `algmg.verbose=1` prints the
AlgMG hierarchy, and `mlmg.tol_rel` sets the tolerance of both solvers.

## Expected output (3D)

```
Checkerboard of 4-cell blocks, jump 10000, 32^3 cells
...
AlgMG standalone: 7 levels, 9 iterations, relative residual 8.6e-11, max error 7.3e-10
...
MLMG geometric: did not converge after 88 iterations
MLMG geometric + AlgMG bottom: did not converge after 88 iterations
MLMG algebraic: 3 iterations, max error 6.6e-11
MLMG hybrid: 9 iterations, max error 9.8e-10, switched to AlgMG
PASSED
```
