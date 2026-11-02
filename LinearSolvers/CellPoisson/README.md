# CellPoisson ★

Solve Poisson's equation `lap(phi) = f` for cell-centered `phi` with
`MLPoisson` and the multigrid solver `MLMG`. This is the minimal workflow
for every AMReX linear solver:

1. Build the operator from `Geometry`, `BoxArray` and `DistributionMapping`.
2. Set the domain boundary types with `setDomainBC`, then the boundary
   values with `setLevelBC`. The values live in the ghost cells of a
   `MultiFab` (here the solution itself), and represent values on the
   domain faces.
3. Build an `MLMG` solver, set its options, and call `solve`.

The boundary type in each direction comes from the inputs, e.g.
`bc = dirichlet neumann periodic`. The exact solution
`phi = prod_d cos(2 pi x_d)` works with all three. Without a Dirichlet
direction, `phi` is determined only up to a constant, so the program
removes the mean before comparing.

The problem is solved on two grids. Halving the cell size reduces the
error by a factor of four, showing second-order accuracy. Works in 1D, 2D
and 3D.

**Key APIs:** `MLPoisson`, `LinOpBCType`, `setDomainBC`, `setLevelBC`,
`MLMG`, `MLMG::solve`, `setVerbose`, `setMaxIter`

`MLMG::solve(sol, rhs, tol_rel, tol_abs)` stops when
`max|rhs - L(sol)| <= max(tol_abs, tol_rel*R)`, where `R` is the larger of
`max|rhs|` and the initial residual. Set `tol_abs` to zero unless you know
the scale of the residual.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=LinearSolvers
cmake --build build --target LinearSolvers_CellPoisson -j 8
cd build/LinearSolvers/CellPoisson && mpiexec -n 2 ./CellPoisson inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
Boundary conditions: dirichlet neumann periodic
MLMG: Initial rhs               = 116.7325837
MLMG: Initial residual (resid0) = 2821.164753
MLMG: Final Iter. 7 resid, resid/resid0 = 4.106625795e-08, 1.455649051e-11
...
n_cell = 32: max error = 0.003249084214
n_cell = 64: max error = 0.0008196247914
Convergence rate = 1.98699761
PASSED
```
