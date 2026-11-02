# ABecLaplacian ★★

Solve the variable-coefficient equation

```
a*alpha(x)*phi - b*div(beta(x)*grad(phi)) = f
```

with `MLABecLaplacian` as one composite solve over all levels of a fixed
AMR hierarchy. Each finer level refines the middle half of the level below
(`max_level = 1` by default). `beta` is face-centered and rises steeply but
smoothly from 1 to about `beta_in` (10 by default) inside a sphere covered by
the fine level.

- `setScalars(a, b)`, `setACoeffs(lev, alpha)` with cell-centered `alpha`,
  and `setBCoeffs(lev, beta)` with an array of face-centered `MultiFab`s.
- Boundary types per direction from the inputs: `dirichlet`, `neumann`
  (`LinOpBCType::inhomogNeumann`) or `periodic`. The ghost cells passed to
  `setLevelBC` hold `phi` on the face for Dirichlet, and `d(phi)/dx_d` on
  the face for Neumann (same sign on the low and high sides).
- `setMaxOrder` sets the order of the stencil at Dirichlet and coarse/fine
  boundaries.

After the solve, `MLMG` has averaged the fine solution onto the covered
coarse cells. The error is checked on every level against the exact
solution `phi = prod_d sin(2 pi x_d + pi/4)`. Without a Dirichlet
direction, the mean of `phi` is only controlled by the `a*alpha` term, and
the error is larger. Works in 1D, 2D and 3D.

**Key APIs:** `MLABecLaplacian`, `setScalars`, `setACoeffs`, `setBCoeffs`,
`LinOpBCType::inhomogNeumann`, `setMaxOrder`, `MLMG::solve` with
`GetVecOfPtrs`/`GetVecOfConstPtrs`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=LinearSolvers
cmake --build build --target LinearSolvers_ABecLaplacian -j 8
cd build/LinearSolvers/ABecLaplacian && mpiexec -n 2 ./ABecLaplacian inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Try `max_level = 2`, `bc = neumann neumann neumann`, or a larger `n_cell`.

## Expected output

```
MLMG: Initial rhs               = 987.8425329
MLMG: Initial residual (resid0) = 1525.720905
MLMG: Final Iter. 8 resid, resid/resid0 = 5.350898391e-08, 3.507127924e-11
Level 0: 8 boxes, 32768 cells, max error = 0.004386008405
Level 1: 8 boxes, 32768 cells, max error = 0.00145993704
PASSED
```
