# Poisson ★★

Solve `lap(phi) = rhs` with FFT-based direct solvers and check the result
against a manufactured solution:

1. `FFT::Poisson` with periodic boundaries.
2. `FFT::Poisson` with homogeneous Dirichlet (`FFT::Boundary::odd`, sine
   transform) and Neumann (`FFT::Boundary::even`, cosine transform)
   boundaries, set per face with
   `Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>`.
3. 3D only: `FFT::PoissonHybrid`, which uses FFTs in `x` and `y` and a
   tridiagonal solve in `z`. The `z` boundaries must be Dirichlet or
   Neumann, and `z` may have a non-uniform cell size.

The solvers invert the standard second-order Laplacian exactly, so the
relative residual of the discrete equations is at round-off level, while the
error against the exact solution is `O(dx^2)`. If `phi` has ghost cells, the
solver fills one layer from the boundary conditions; we use them to compute
the residual. Works in 2D and 3D.

With all boundaries periodic or Neumann, the solution is unique only up to a
constant and `rhs` must have zero mean.

**Key APIs:** `FFT::Poisson`, `FFT::PoissonHybrid`, `FFT::Boundary`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=FFT
cmake --build build --target FFT_Poisson -j 8
cd build/FFT/Poisson && mpiexec -n 2 ./Poisson inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
Periodic (FFT::Poisson)
    relative residual 1.412167105e-14
    max error         0.003172687413
Dirichlet/Neumann (FFT::Poisson)
    relative residual 6.272380662e-14
    max error         0.0007345742346
Periodic x,y + Neumann/Dirichlet z (FFT::PoissonHybrid)
    relative residual 1.762746748e-14
    max error         0.003096260662
PASSED
```

Doubling `n_cell` reduces the error by about 4.
