# Basic ★

Parallel real-to-complex FFTs on a `MultiFab` with `FFT::R2C`:

1. Forward and backward transforms with the spectral data in a `cMultiFab`.
   FFTs are unnormalized: a round trip multiplies the data by the number of
   points, so we scale by `scalingFactor()`.
2. A spectral derivative `d/dx`, computed by `forwardThenBackward` with a
   lambda that multiplies each spectral point by `i*kx`.
3. A low-pass filter that zeroes all modes with `|k| > kcut`.

The test field is a smooth low-wavenumber part plus a high-wavenumber part
that the filter removes. All three results are checked against the exact
answer to round-off. Works in 2D and 3D.

The spectral data has `(nx/2+1, ny, nz)` points: only `kx >= 0` is stored,
because the transform of real data is Hermitian. In `y` and `z`, indices
above `n/2` are negative wavenumbers. A CPU build uses FFTW; GPU builds use
cuFFT, rocFFT or oneMKL.

**Key APIs:** `FFT::R2C`, `getSpectralDataLayout`, `forward`, `backward`,
`forwardThenBackward`, `scalingFactor`, `cMultiFab`, `GpuComplex`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=FFT
cmake --build build --target FFT_Basic -j 8
cd build/FFT/Basic && mpiexec -n 2 ./Basic inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
Real domain     ((0,0,0) (31,31,31) (0,0,0))
Spectral domain ((0,0,0) (16,31,31) (0,0,0))
Round trip error:        6.661338148e-16
Spectral d/dx rel. error: 3.537172464e-15
Low-pass filter error:   7.21644966e-16
PASSED
```
