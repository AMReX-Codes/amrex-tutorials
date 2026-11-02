# MPMD ★★★

Couple two separately built programs that run in the same MPI job
(multiple program, multiple data) and exchange MultiFab data with
`amrex::MPMD`. The second program can be C++ (`agent.cpp`) or Python with
pyAMReX (`agent.py`).

- `main.cpp` (`MPMD_Main`) owns the domain `[0, n_cell)^D`. It fills
  component 0 with a Gaussian, sends it, and receives component 1.
- `agent.cpp` (`MPMD_Agent`) and `agent.py` own the upper half
  `[n_cell/2, n_cell)^D` with their own boxes. They receive component 0,
  set component 1 to 10 times component 0, and send it back.
- The main program checks that component 1 is 10 times component 0 where
  the agent covers the domain, and unchanged elsewhere.

`MPMD::Initialize` splits `MPI_COMM_WORLD` into one communicator per
program, and AMReX is initialized on it. `MPMD::Copier` exchanges the
BoxArrays and DistributionMappings of both programs once, and then `send`
and `recv` copy components on the intersection of the two BoxArrays. The
programs may use different boxes and numbers of ranks, but must use the
same dimension and precision. With `Copier(ba, dm, true)` on one side and
`Copier(true)` on the other, the second program adopts the BoxArray of the
first instead.

**Key APIs:** `MPMD::Initialize`, `MPMD::Finalize`, `MPMD::Copier`,
`Copier::send`, `Copier::recv`, `amrex::Initialize(argc, argv, true, comm)`;
Python: `MPMD_Initialize_without_split`, `initialize_when_MPMD`, `MPMD_Copier`

## Build and run

Requires MPI. CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=Coupling
cmake --build build --target Coupling_MPMD_Main Coupling_MPMD_Agent -j 8
cd build/Coupling/MPMD
mpiexec -n 2 ./MPMD_Main inputs : -n 2 ./MPMD_Agent inputs
```

GNU Make builds one program at a time, selected with `PROG`:

```
make -j8 PROG=main
make -j8 PROG=agent
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs : -n 2 ./agent3d.gnu.MPI.ex inputs
```

With the Python agent (needs pyAMReX built with MPI, and mpi4py; add
`dim=2` for a 2D build):

```
mpiexec -n 2 ./MPMD_Main inputs : -n 2 python3 agent.py inputs
```

The CTest tests run with `-DAMReXTutorials_PYTHON=ON` include the Python agent.

If `mpiexec` itself crashes on the two-program command line (seen with Open
MPI 5 and PMIx 5.0.9 on Ubuntu), set `PMIX_MCA_pmdl=^ompi` in the environment.
The CTest tests do this automatically when needed.

## Expected output

With the C++ agent (the Python agent prints `Agent (Python):` instead):

```
Main: 2 of 4 ranks in MPI_COMM_WORLD
Agent: 2 of 4 ranks in MPI_COMM_WORLD
Agent: received component 0 with min = 1, max = 1.92937572
Main: max error of the round trip = 0
```
