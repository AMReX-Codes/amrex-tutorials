# GpuContainers ★★

Work with GPU-friendly containers and memory arenas: fill a
`Gpu::DeviceVector` with random numbers, keep the values above a threshold
with a prefix-sum based stream compaction, copy results to pinned host
memory, and check them against the same operations done on the host. The
code is the same for CPU and GPU builds.

**Key APIs:** `Gpu::DeviceVector`, `Gpu::HostVector`, `Gpu::PinnedVector`,
`Gpu::copyAsync`, `Gpu::streamSynchronize`, `ParallelForRNG`,
`amrex::Random`, `Scan::PrefixSum`, `Scan::ExclusiveSum`, `The_Arena`,
`The_Async_Arena`, `The_Pinned_Arena`

## Notes

- `DeviceVector` memory comes from `The_Arena()` (device memory on GPU), and
  `PinnedVector` from `The_Pinned_Arena()` (page-locked host memory, needed
  for truly asynchronous copies). `HostVector` is ordinary host memory.
- `copyAsync` returns before the copy is done on GPU; call
  `Gpu::streamSynchronize()` before using the data.
- `The_Async_Arena()` is for temporaries: memory freed there is reused only
  after the work already queued on the stream is done.
- `Scan::PrefixSum` takes a function that produces each input and a function
  that consumes each partial sum, so compaction needs no temporary arrays.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_GpuContainers -j 8
cd build/Base/GpuContainers && ./GpuContainers inputs
```

GNU Make, with AMReX cloned next to this repository (add `USE_CUDA=TRUE`
for a CUDA build):

```
make -j8
./main3d.gnu.MPI.ex inputs
```

## Expected output

```
Kept 24739 of 100000 values above 0.75
Compaction matches the host result: yes
ExclusiveSum total = 150000, matches the host result: yes
PASSED
```

The number of kept values depends on the random numbers.
