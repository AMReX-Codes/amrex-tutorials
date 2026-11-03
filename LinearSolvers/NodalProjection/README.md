# NodalProjection ★★

Project a cell-centered velocity `u` onto a divergence-free field, the core
step of incompressible flow solvers:

```
div(sigma*grad(phi)) = div(u),    u_new = u - sigma*grad(phi)
```

with `phi` on the nodes, solved by `MLNodeLaplacian`.

1. `compDivergence` computes the nodal `div(u)`. The velocity needs one
   ghost cell.
2. `MLMG` solves for `phi`. With walls (`LinOpBCType::Neumann`, i.e.,
   `u.n = 0`) and periodic boundaries only, the system is singular, and
   MLMG makes it solvable.
3. `updateVelocity` subtracts `sigma*grad(phi)`.

Here `sigma = 1` is passed to the constructor. For variable density, use
`sigma = 1/rho` with `setSigma`.

The initial velocity is a divergence-free field plus a gradient. The
domain is periodic in x with walls elsewhere. This is an approximate
projection: the nodal Laplacian is not exactly the product of the discrete
divergence and gradient, so `div(u_new)` is O(dx^2), not zero. The program
checks that the remaining divergence and the error in the projected
velocity are O(dx^2). Works in 2D and 3D.

**Key APIs:** `MLNodeLaplacian`, `compDivergence`, `updateVelocity`,
nodal `BoxArray` (`convert(ba, IntVect::TheNodeVector())`)

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=LinearSolvers
cmake --build build --target LinearSolvers_NodalProjection -j 8
cd build/LinearSolvers/NodalProjection && mpiexec -n 2 ./NodalProjection inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
MLMG: Initial rhs               = 58.36656387
MLMG: Initial residual (resid0) = 58.36656387
MLMG: Final Iter. 6 resid, resid/bnorm = 3.849152108e-10, 6.594789641e-12
max |div(u)| before projection: 59.90092017
max |div(u)| after projection:  0.1477836432
max error in the divergence-free part: 0.02371584586
PASSED
```
