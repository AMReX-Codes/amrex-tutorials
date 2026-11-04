# SharedMemory ★★★

Write a GPU kernel by hand, with shared memory, and a separate CPU
implementation. The example bins the cells by the value of one field and
sums a signed integer field over each bin. Integer sums are exact, so the
two methods must agree to the bit, and they need a signed 64-bit
accumulator. With `ParallelFor` and one atomic add per cell, all threads
compete for the same few addresses in global memory, and the hardware
serializes the adds to each one. The hand-written kernel gives each thread
block private sums in shared memory, where contention is cheap, and adds
them to the global `Long` sums once per block; on CPU, each OpenMP thread
does the same with a private array. The program checks that both results
agree and times them with `BL_PROFILE` (see
[HeatEquation](../HeatEquation#profiling)).

**Key APIs:** `amrex::launch<nthreads>(nblocks, stream, f)`,
`amrex::launch<nthreads>(nblocks, shared_mem_bytes, stream, f)`,
`__shared__` (CUDA, HIP), `Gpu::Handler` (SYCL),
`Gpu::Atomic::AddNoRet`, `HostDevice::Atomic::Add`, `BoxIndexer`,
`Gpu::Device::maxBlocksPerLaunch`, `AMREX_USE_GPU`, `AMREX_USE_SYCL`

## Notes

- `ParallelFor` is portable and is the right choice almost always. A raw
  launch is for the rare kernel that needs threads of a block to cooperate:
  shared memory, block-wide synchronization, or custom reductions. The
  portable version stays as the reference to check against.
- `amrex::launch<nthreads>(nblocks, stream, f)` runs `f` on `nblocks`
  blocks of `nthreads` threads. The kernel body differs by backend: with
  CUDA and HIP it uses `threadIdx`, `blockIdx`, `__syncthreads()` and a
  static `__shared__` array. With SYCL, the shared memory is dynamic and
  its size is passed to `amrex::launch`, and that overload hands the
  kernel a `Gpu::Handler` whose `threadIdx()`, `blockDim()`, `gridDim()`,
  `globalIdx()`, `syncThreads()` and `sharedMemory()` play these roles;
  atomics on the shared memory take the `local_space` address space.
- Atomic adds to the same global address are serialized, so the time of
  the `ParallelFor` version grows as the bins get fewer. The private sums
  in shared memory are `int`, with a check that a block cannot overflow
  them, and leave one `Long` atomic per bin and block. Older AMReX
  versions implemented `Gpu::Atomic::Add` on `Long` as a compare-and-swap
  loop on CUDA and HIP, which under contention made the global version
  slower by another factor of a hundred or more.
- There is no CPU version of this launch. The CPU code is a separate path
  selected with `#ifdef AMREX_USE_GPU`; here it is an `MFIter` loop with
  tiling and a vector of private sums per OpenMP thread, merged after the
  parallel region. The cell loop is `amrex::Loop`, not `ParallelFor`,
  because it scatters into the sums.
- `BoxIndexer` converts a linear cell index into `(i,j,k)` within a box, so
  the threads of a grid-stride loop can walk a box.
- Related: `ParallelFor(Gpu::KernelInfo().setReduction(true), box, f)` passes
  a `Gpu::Handler` to the kernel for `Gpu::deviceReduceSum` and friends,
  which AMReX uses for its block reductions.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_SharedMemory -j 8
cd build/Base/SharedMemory && mpiexec -n 2 ./SharedMemory inputs
```

GNU Make, with AMReX cloned next to this repository (add `USE_CUDA=TRUE`
for a CUDA build, which makes `main3d.gnu.MPI.CUDA.ex`):

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

For timings, build with TinyProfiler on (`-DAMReX_TINY_PROFILE=ON` or
`make TINY_PROFILE=TRUE`) and use e.g. `n_cell=128 max_grid_size=128
nbins=64 nrepeat=5 tiny_profiler.print_threshold=0` on a GPU. The fewer
the bins, the more the global atomics collide and the larger the speedup
of the shared-memory version.

## Expected output

```
bin  sum
0  -652800
8  2955200
16  13470720
24  26177280
32  58570240
40  14887680
48  9210880
56  2266560
The two results are identical; total 873843200 (expected 873843200)
PASSED
```

With TinyProfiler on, the profiling tables list the timers
`binned_sum_atomic` and `binned_sum_private`. On a GPU, the shared-memory
version is faster by several times with a thousand bins and by ten to
twenty times with a few dozen bins or fewer. The table hides timers below
1% of the total unless run with `tiny_profiler.print_threshold=0`.
