# OpenBC ★★

Solve `lap(phi) = rho` in free space (open boundaries) for a Gaussian charge
and compare with the analytic potential. In 3D, it is
`phi = -erf(r/(sqrt(2) sigma)) / (4 pi r)`. In 2D, it is
`phi = (ln(r) + E1(r^2/(2 sigma^2))/2) / (2 pi)`, where `E1` is the
exponential integral.

The free-space solution is the convolution of `rho` with the Green's
function `G`, which is `-1/(4 pi r)` in 3D and `ln(r)/(2 pi)` in 2D.
`FFT::OpenBCSolver` computes it with FFTs on a domain doubled in each
direction (Hockney's method), so the periodic images do not interact. You
provide the Green's function with `setGreensFunction`, integrated over a
cell; here we use a `2^D`-point Gauss quadrature, which also avoids the
singularity at `r = 0`. In 3D, the tutorial also uses `FFT::PoissonOpenBC`,
which provides the same Green's function for you.

The error is second order in `dx/sigma` and is checked against a tolerance.
Works in 2D and 3D (`FFT::PoissonOpenBC` is 3D only).

**Key APIs:** `FFT::OpenBCSolver`, `setGreensFunction`,
`FFT::PoissonOpenBC`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=FFT
cmake --build build --target FFT_OpenBC -j 8
cd build/FFT/OpenBC && mpiexec -n 2 ./OpenBC inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
FFT::OpenBCSolver   relative error: 0.007061771676
FFT::PoissonOpenBC  relative error: 0.007061771676
PASSED
```

With `n_cell = 64`, the error drops to about 0.0018.
