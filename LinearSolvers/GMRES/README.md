# GMRES ★★

Use AMReX's GMRES with multigrid as the preconditioner, in two ways.

1. **Your own operator.** `amrex::GMRES<V,M>` works with any vector type
   `V` and operator class `M`. `ConvDiffOp` implements the
   convection-diffusion operator `lap(u) - c.grad(u)` without storing a
   matrix. It provides the matrix-vector product (`apply`), vector
   operations (`dotProduct`, `linComb`, ...), and a preconditioner
   (`precond`) that runs one `MLMG` V-cycle of `MLPoisson`, i.e., of the
   Laplacian alone. The operator is not symmetric, so CG cannot be used.
   With the preconditioner GMRES needs about 40 iterations, independent
   of the grid size. Without it (`use_precond = 0`) it needs about 300.
2. **An MLMG operator.** `GMRESMLMG` takes an `MLMG` object and solves its
   system with GMRES, using MLMG itself as the preconditioner. This can
   help when MLMG alone converges slowly.

The preconditioner must be a fixed linear operator. That is why it does a
fixed number of V-cycles and uses the smoother instead of BiCGStab as the
bottom solver. Both problems use homogeneous Dirichlet boundaries and the
exact solution `u = prod_d sin(pi x_d)`. Works in 1D, 2D and 3D.

**Key APIs:** `amrex::GMRES`, `GMRES::define`, `setRestartLength`,
`setMaxIters`, `getNumIters`, `getStatus`, `MLMG::precond`,
`MLMG::setBottomSolver`, `GMRESMLMG`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=LinearSolvers
cmake --build build --target LinearSolvers_GMRES -j 8
cd build/LinearSolvers/GMRES && mpiexec -n 2 ./GMRES inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
GMRES: Solve Time = 0.807119752
Convection-diffusion with GMRES: 40 iterations, status 0, max error = 0.001317131374
GMRES: Solve Time = 0.104936725
Poisson with GMRESMLMG: 7 iterations, max error = 0.0008006773422
PASSED
```
