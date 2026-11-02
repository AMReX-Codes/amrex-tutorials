# BoundaryConditions ★★

Solve the heat equation on a non-periodic domain and fill ghost cells at
physical boundaries. Dirichlet values are set at $x = 0$ and $x = 1$, $y$
has zero-gradient boundaries, and $z$ is periodic. The solution
$\phi = \phi_{lo} + (\phi_{hi} - \phi_{lo})x + e^{-\pi^2 D t}\sin(\pi x)$
relaxes to the linear steady state, and the program checks it against this
exact solution.

**Key APIs:** `BCRec`, `BCType::ext_dir/foextrap/reflect_even/int_dir`,
`PhysBCFunct`, `GpuBndryFuncFab`, `FillBoundary`, `MFIter`,
`ParallelFor(box, [=] (IntVect const&) ...)`

## Filling ghost cells

1. `FillBoundary(geom.periodicity())` fills ghost cells that overlap other
   boxes, including periodic images.
2. `PhysBCFunct<GpuBndryFuncFab<F>>` fills ghost cells outside the domain
   according to the `BCRec`. AMReX handles `foextrap`, `reflect_even`,
   `reflect_odd` and `hoextrap`; the user functor `F` (here `DirichletFill`)
   is called for every such cell and sets the `ext_dir` values.

For `ext_dir`, the ghost cell holds the value on the boundary face, so the
Laplacian uses a half-cell distance next to that boundary.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_BoundaryConditions -j 8
cd build/Base/BoundaryConditions && mpiexec -n 2 ./BoundaryConditions inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Try `phi_lo=2 phi_hi=-1`, or `stop_time=0.5 plot_int=100` to watch the
approach to the steady state.

## Expected output

```
Finished 342 steps at time 0.05, max error = 2.437267784e-05
PASSED
```
