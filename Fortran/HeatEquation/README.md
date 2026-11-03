# HeatEquation ★

The heat equation of [Base/HeatEquation](../../Base/HeatEquation) written
entirely in Fortran with the AMReX Fortran interfaces. It solves
$\partial\phi/\partial t = D\nabla^2\phi$ with forward Euler on a periodic
unit domain, starting from $\phi = 1 + \prod_d \sin(2\pi x_d)$, and
compares the result with the exact solution, which decays as
$\exp(-4\pi^2 d D t)$ in $d$ dimensions. It runs on CPUs, with MPI and
OpenMP.

**Key APIs:** `amrex_base_module`, `amrex_init`, `amrex_parmparse`,
`amrex_box`, `amrex_boxarray`, `amrex_distromap`, `amrex_geometry`,
`amrex_multifab`, `amrex_mfiter`, `dataptr`, `fill_boundary`,
`amrex_multifab_swap`, `amrex_parallel_reduce_max`, `amrex_write_plotfile`

## Code

`main.f90` sets up the domain, the grids and two MultiFabs, and runs the
time loop. `heat.f90` has the kernels. Each kernel loops over tiles with
an `amrex_mfiter` inside an OpenMP parallel region, gets a pointer
`p(i,j,k,n)` to the whole box of data, including ghost cells, with
`dataptr`, and passes it to an ordinary Fortran loop over the tile box
`mfi%tilebox()`. The indices are the cell indices of AMReX,
zero-based, and the code is the same in 1D, 2D and 3D: indices beyond
`amrex_spacedim` are always 1.

Fortran objects such as `amrex_multifab` own C++ objects. Build them with
`amrex_*_build` and free them with `amrex_*_destroy` before
`amrex_finalize`. Assignment between them makes a shallow copy.

## Build and run

CMake, from the repository root, with a Fortran compiler:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Fortran -DAMReX_FORTRAN_INTERFACES=ON
cmake --build build --target Fortran_HeatEquation -j 8
cd build/Fortran/HeatEquation && mpiexec -n 2 ./HeatEquation inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Add `plot_int=10` to write a plotfile every 10 steps, and e.g.
`n_cell=128 stop_time=0.05` for a larger run.

## Expected output

```
Step 10, time = 1.4648437500E-03
...
Finished 69 steps at time 1.0000000000E-02, mode amplitude = 3.0594420565E-01
Max error = 1.9528508660E-03, tolerance = 8.6145674863E-03
PASSED
```
