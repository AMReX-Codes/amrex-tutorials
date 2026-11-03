# HeatEquation ★

Solve the heat equation $\partial\phi/\partial t = D\nabla^2\phi$ with
explicit (forward Euler) time stepping on a periodic unit domain. The
initial condition is $\phi = 1 + \prod_d \sin(2\pi x_d)$; its sine mode
decays as $\exp(-4\pi^2 d D t)$ in $d$ dimensions, and the program compares
the numerical solution with this exact solution.

**Key APIs:** `MultiFab`, `Geometry`, `FillBoundary`, `ParallelFor(mf, ...)`,
`const_arrays`, `WriteSingleLevelPlotfile`

## Algorithm

Each step fills the ghost cells with `FillBoundary` and then computes
`phi_new = phi_old + dt * D * Laplacian(phi_old)` with the standard
second-order stencil. The time step is `dt = cfl * dx^2 / (2 * dim * D)`.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_HeatEquation -j 8
cd build/Base/HeatEquation && mpiexec -n 2 ./HeatEquation inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Add `plot_int=10` to write a plotfile every 10 steps, and e.g.
`n_cell=128 stop_time=0.05` for a larger run.

## Profiling

AMReX has a built-in profiler, TinyProfiler, which is off by default in
the tutorials. Build with `-DAMReX_TINY_PROFILE=ON` (CMake) or
`make TINY_PROFILE=TRUE` (GNU Make), and at the end of the run it prints
the time spent in AMReX functions such as `FillBoundary`. To time your own
code, put `BL_PROFILE("name");` at the start of a scope; it times until the
end of the scope and does nothing when TinyProfiler is off. For GPU runs,
add `tiny_profiler.device_synchronize_around_region=1` to get accurate
times for asynchronous kernels. The extra synchronization slows the code
down, so it is not recommended for production runs.

## Expected output

```
Step 10, time = 0.00146484375
...
Finished 69 steps at time 0.01, mode amplitude = 0.3059442057
Max error = 0.001952850866, tolerance = 0.008614567486
PASSED
```
