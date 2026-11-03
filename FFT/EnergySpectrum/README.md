# EnergySpectrum ★★

Compute the shell-averaged kinetic energy spectrum `E(k)` of a periodic
velocity field. This is a standard diagnostic for turbulence.

The velocity field is synthetic: three Fourier modes with `|k| = 2, 5, 8`
(in units of `2 pi / L`) and known energies. All velocity components are
transformed in one batched `FFT::R2C` call. Each spectral point then adds
`|u_hat|^2 / (2 N^2)` to the shell `round(|k|)` with an atomic add on the
device, giving a histogram in a `Gpu::DeviceVector`. The local histograms
are summed across MPI ranks and written to `spectrum.txt`.

Only `kx >= 0` is stored by the real-to-complex transform. Because the
transform of real data is Hermitian, `u_hat(-k) = conj(u_hat(k))`, so points
with `0 < kx < nx/2` count twice.

The tutorial checks that the energy is in the expected shells and that
Parseval's theorem holds: `sum_k E(k)` equals the kinetic energy
`<|u|^2>/2` computed in physical space. Works in 2D and 3D.

**Key APIs:** `FFT::R2C`, `FFT::Info::setBatchSize`, `cMultiFab`,
`HostDevice::Atomic::Add`, `Gpu::DeviceVector`, `Gpu::copy`,
`ParallelDescriptor::ReduceRealSum`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=FFT
cmake --build build --target FFT_EnergySpectrum -j 8
cd build/FFT/EnergySpectrum && mpiexec -n 2 ./EnergySpectrum inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
Nonzero shells:
  E(2) = 0.25  (exact 0.25)
  E(5) = 0.0625  (exact 0.0625)
  E(8) = 0.015625  (exact 0.015625)
Kinetic energy: physical space 0.328125, spectral space 0.328125
Parseval relative error: 0
Max shell error:         8.458842092e-17
Spectrum written to spectrum.txt
PASSED
```
