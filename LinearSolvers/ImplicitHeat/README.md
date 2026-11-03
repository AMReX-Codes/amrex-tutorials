# ImplicitHeat ★★

Advance the heat equation `d(phi)/dt = kappa*lap(phi)` on a periodic
domain with backward Euler,

```
(1 - dt*kappa*lap) phi^{n+1} = phi^n,
```

which is `MLABecLaplacian` with `a = 1`, `alpha = 1`, `b = dt` and
`beta = kappa`. The time step is `dt_factor` (10 by default) times the
explicit stability limit `dx^2/(2*D*kappa)`; backward Euler is stable for
any `dt`.

The operator and the `MLMG` solver are built once, before the time loop,
and reused every step. Constant coefficients are set with the scalar
versions of `setACoeffs` and `setBCoeffs`. Each solve starts from the
previous solution, a good initial guess.

The initial condition `phi = prod_d sin(2 pi x_d)` is an eigenfunction of
the discrete Laplacian, so each step multiplies it by exactly
`1/(1 + dt*kappa*lambda)`. The program checks the numerical solution
against this to the solver tolerance. It also prints the error against the
PDE solution `exp(-kappa*D*(2 pi)^2*t)`, which is dominated by the O(dt)
error of backward Euler: try `dt_factor = 1 nsteps = 200`. Works in 1D, 2D
and 3D.

**Key APIs:** `MLABecLaplacian`, `setScalars`, `setACoeffs`/`setBCoeffs`
with scalars, `setLevelBC(0, nullptr)`, reusing `MLMG` across solves,
`MLMG::getNumIters`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=LinearSolvers
cmake --build build --target LinearSolvers_ImplicitHeat -j 8
cd build/LinearSolvers/ImplicitHeat && mpiexec -n 2 ./ImplicitHeat inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
Step 1: time = 0.001627604167, MLMG iterations = 7
...
Step 20: time = 0.03255208333, MLMG iterations = 7
Decay factor: numerical 0.02974420513, discrete exact 0.02974420513, PDE exact 0.0211669512
Relative error vs. discrete exact solution: 7.80386895e-11
Relative error vs. PDE exact solution:      0.399393567
PASSED
```
