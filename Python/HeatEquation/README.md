# HeatEquation (Python) ★

The Python version of [Base/HeatEquation](../../Base/HeatEquation): solve
the heat equation `d(phi)/dt = D Laplacian(phi)` with explicit time steps on
a periodic unit cube, entirely in Python with pyAMReX. The stencil is
written with NumPy (CPU), CuPy (CUDA/HIP) or dpnp (SYCL) array slicing on
zero-copy views of each box. The initial condition is
`1 + sin(2 pi x) sin(2 pi y) sin(2 pi z)`; its sine mode decays as
`exp(-12 pi^2 D t)`, and the final error against this exact solution is
checked.

**Key APIs:** `amr.MultiFab`, `amr.xp`, `mf.array(mfi).to_xp()`,
`mf.fill_boundary`, `geom.periodicity()`, `mf.max`, `amr.ParmParse.to_dict`,
`amr.write_single_level_plotfile`

## Run

Requires pyAMReX 26.09 or newer (`import amrex`), plus CuPy for CUDA/HIP
or dpnp for SYCL GPU builds of pyAMReX.

```
python3 main.py inputs
mpiexec -n 2 python3 main.py inputs n_cell=64 plot_int=20
```

With CMake, `-DAMReXTutorials_PYTHON=ON` registers this script as the CTest
test `Python_HeatEquation`; set `-DPython3_EXECUTABLE` to a Python that has
pyAMReX installed.

Each time step loops over the boxes in Python, but the work on each box is
done by NumPy/CuPy/dpnp, so large boxes (`max_grid_size`) keep the Python
overhead small.

## Expected output

```
Step 10, time = 0.0014648437499999998
...
Finished 69 steps at time 0.01, mode amplitude = 0.3059442057
Max error = 0.001952850866, tolerance = 0.008614567486
PASSED
```
