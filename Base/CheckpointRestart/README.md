# CheckpointRestart ★★

Write checkpoints of a heat equation run and restart from them. A
checkpoint is a directory with a `Header` (step, time and `BoxArray`) and the
`MultiFab` data written by `VisMF`. With `test_restart=1`, the program
restarts from the first checkpoint, reruns the remaining steps, and checks
that the result is bitwise identical to the uninterrupted run. It also reads
the final plotfile back with `PlotFileData` and compares it to the data in
memory.

**Key APIs:** `PreBuildDirectorHierarchy`, `MultiFabFileFullPrefix`,
`VisMF::Write`, `VisMF::Read`, `BoxArray::writeOn/readFrom`,
`ParallelDescriptor::ReadAndBcastFile`, `PlotFileData`, `ParallelCopy`

## Notes

- Only the I/O rank writes the `Header`; every rank writes its own data.
- On restart, the `BoxArray` comes from the checkpoint, and a new
  `DistributionMapping` is made, so a run may restart on a different number
  of MPI ranks.
- `PlotFileData::get` returns data with its own `DistributionMapping`;
  `ParallelCopy` moves it into the layout of another `MultiFab`.

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_CheckpointRestart -j 8
cd build/Base/CheckpointRestart
mpiexec -n 2 ./CheckpointRestart inputs
mpiexec -n 1 ./CheckpointRestart inputs restart=chk00010
mpiexec -n 2 ./CheckpointRestart inputs test_restart=1
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs test_restart=1
```

## Expected output (`test_restart=1`)

```
Writing checkpoint chk00010 at step 10
Writing checkpoint chk00020 at step 20
Finished step 20, time 0.0029296875, wrote plt00020
Restarted run vs. uninterrupted run: max diff = 0, time diff = 0
Plotfile plt00020 (step 20, time 0.0029296875) vs. memory: max diff = 0
PASSED
```
