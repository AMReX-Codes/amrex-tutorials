# CustomLinOp ★★★

Write your own linear operator and solve with `MLMG`. `MultiCompNodeLap`
is a nodal operator with `N` coupled components,

```
L(phi)_n = sum_m C_nm lap(phi_m),
```

where `C` is a constant, invertible `N x N` matrix from the inputs
(`N <= 4`).

`MultiCompNodeLap` derives from `MLNodeLinOp`, which provides the
multigrid hierarchy, ghost node exchange, the Dirichlet node mask
(`m_dirichlet_mask`) and dot products. The derived class implements:

- `getNComp`: the number of components.
- `Fapply`: `out = L(in)`, with `in`'s ghost nodes already filled.
- `Fsmooth`: block Jacobi, `sol += omega * D^{-1} (rhs - L(sol))`, where
  the diagonal block `D = C * diag(lap)` couples the components at a node.
- `restriction` (full weighting) and `interpolation` (linear) between
  multigrid levels, which must handle different box layouts after
  agglomeration.
- `prepareForSolve`: builds the masks.

It is single-level only. A multi-level operator would also implement
`averageDownSolutionRHS` and `reflux` for the coarse/fine boundaries (see
[Multi-Component Operators](https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html#multi-component-operators)
in the AMReX documentation).

The program checks the error against the exact solution
`phi_n = prod_d sin((n+1) pi x_d)`, which is zero on the Dirichlet
boundary. Works in 1D, 2D and 3D.

**Key APIs:** `MLNodeLinOp`, `Fapply`, `Fsmooth`, `restriction`,
`interpolation`, `getNComp`, `buildMasks`, `nodalSync`,
`setDirichletNodesToZero`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=LinearSolvers
cmake --build build --target LinearSolvers_CustomLinOp -j 8
cd build/LinearSolvers/CustomLinOp && mpiexec -n 2 ./CustomLinOp inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Try three components: `coeff = 1 0 0  0 2 0.5  0 0.3 1`.

## Expected output

```
Solving for 2 components with C = 2 1 0.5 1
MLMG: Initial rhs               = 139.3718454
MLMG: Initial residual (resid0) = 139.3718454
MLMG: Final Iter. 18 resid, resid/bnorm = 9.651316191e-09, 6.924867904e-11
Component 0: max error = 0.0008035776406
Component 1: max error = 0.00321896441
PASSED
```
