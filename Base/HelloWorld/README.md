# HelloWorld ★

The smallest AMReX program: initialize and finalize AMReX, and print from
one or all MPI ranks.

**Key APIs:** `amrex::Initialize`, `amrex::Finalize`, `amrex::Print`,
`amrex::AllPrint`

## Build and run

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Base
cmake --build build --target Base_HelloWorld -j 8
cd build/Base/HelloWorld && mpiexec -n 2 ./HelloWorld
```

GNU Make, with AMReX cloned next to this repository:

```
make -j8
mpiexec -n 2 ./main3d.gnu.MPI.ex
```

## Expected output

```
Hello world from AMReX version 26.10
  MPI ranks:      2
  ...
```
