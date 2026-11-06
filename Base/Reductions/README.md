# Reductions ★★

Compute sums, minima, maxima and counts over a `MultiFab` and a device
vector, on CPU and GPU, and combine the results across MPI ranks. The data,
$\phi = i + j + k$, has a simple exact sum and maximum; the program prints
them next to the computed values and checks the two sums.

**Key APIs:** `ParReduce`, `ReduceOps`, `ReduceData`, `Reducer`,
`ReduceOpSum/Min/Max`, `MultiFab::sum/min/max/norm0/norm1/norm2`,
`Reduce::Sum`, `ParallelAllReduce`, `ParallelDescriptor::ReduceRealSum`

## Which one to use

- `MultiFab` member functions for common norms; they already include the
  MPI reduction.
- `ParReduce` for several custom reductions over a `MultiFab` in one pass.
- `ReduceOps`/`ReduceData`, or the more compact `Reducer`, inside `MFIter`
  loops or over arbitrary boxes and index ranges.
- `Reduce::Sum/Min/Max` for a plain array, e.g. a `Gpu::DeviceVector`.

All except the `MultiFab` member functions are local to the MPI rank, so
finish with `ParallelAllReduce` or `ParallelDescriptor::Reduce*`.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_Reductions -j 8
cd build/Base/Reductions && mpiexec -n 2 ./Reductions inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
ParReduce:
  sum = 1523712 (exact 1523712)
  min = 0 (exact 0)
  max = 93 (exact 93)
ReduceOps/ReduceData:
  sum of squares = 79233024
  even cells = 16384
Reducer:
  sum of squares = 79233024
  even cells = 16384
MultiFab member functions:
  sum() = 1523712 (exact 1523712)
  min() = 0
  max() = 93 (exact 93)
  norm0() = 93
  norm1() = 1523712
  norm2() = 8901.293389
Reduce::Sum on a Gpu::DeviceVector:
  sum = 1001000 (exact 1001000)
PASSED
```
