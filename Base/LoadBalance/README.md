# LoadBalance ★★

Rebalance the work over MPI ranks when boxes cost different amounts. The
default `DistributionMapping` gives every rank about the same number of
cells,
which is a good choice when every cell costs the same. Here the cells
inside a ball do `work_ratio` times more work than the others, so the ranks
holding the ball are late at every synchronization point. The program
measures the time each box takes, builds a new mapping from these costs
with the knapsack or the space-filling-curve algorithm, moves the data to
it with `ParallelCopy`, and measures again.

**Key APIs:** `DistributionMapping::makeKnapSack`, `DistributionMapping::makeSFC`,
`LayoutData`, `ParallelCopy`, `amrex::second`, `ParallelAllReduce`

## Notes

- Costs are stored per box in a `LayoutData<Real>`, which is indexed like a
  `MultiFab` in an `MFIter` loop. The timers need
  `Gpu::streamSynchronize()` after each kernel, because kernels run
  asynchronously on GPUs. Any other measure of work can be used instead,
  e.g. the number of particles in each box.
- `makeKnapSack(cost, eff_cur, eff_new)` and `makeSFC(...)` return the new
  mapping and, on the I/O rank, the efficiency (average cost per rank over
  the maximum) of the current and of the proposed mapping. With
  `broadcastToAll = false`, only the I/O rank gets the proposal, so a code
  can decide whether the gain is worth moving the data. The knapsack
  algorithm balances best; the space-filling curve keeps neighboring boxes
  on the same rank, which reduces communication. With no more boxes than
  ranks there is nothing to balance, and the proposed efficiency is
  reported as 1.
- A new `MultiFab` on the new mapping receives the data with
  `ParallelCopy`. The kernel's result does not depend on the mapping, which
  the program checks.
- In an application, the costs are accumulated over some steps and the
  mapping redone every `N` steps or after regridding. A code built on
  `AmrCore` calls `SetDistributionMap(lev, dm)` and `RemakeLevel`
  ([AdvectionAmrCore](../../AMR/AdvectionAmrCore) shows `RemakeLevel`);
  `Amr`/`AmrLevel` codes
  override `AmrLevel::estimateWork` and set
  `amr.loadbalance_with_workestimates = 1` and `amr.loadbalance_level0_int`.
- With OpenMP, the boxes of a rank are timed by the threads working on
  them, so the cost of a rank is the sum of its box times, not the wall
  time.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_LoadBalance -j 8
cd build/Base/LoadBalance && mpiexec -n 2 ./LoadBalance inputs
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

Try `strategy=sfc`, more ranks (`mpiexec -n 4`), and `work_ratio=200` for a
stronger imbalance. On one rank, both mappings are the same. On a GPU, the
time of a kernel on a small box is mostly launch overhead, whatever the
work, so the costs only differ for large boxes and heavy kernels, e.g.
`n_cell=256 max_grid_size=64 work_ratio=5000`.

## Expected output

```
512 boxes on 2 ranks
Default mapping: max/avg cost per rank = 1.20104519
Efficiency of the default mapping 0.8326081387, of the knapsack mapping 0.9987378103
New mapping: max/avg cost per rank = 1.00074445
Max difference between the results: 0
PASSED
```

The numbers vary from run to run, since they are measured times.
