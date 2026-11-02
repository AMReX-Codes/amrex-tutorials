#!/usr/bin/env python3
"""The agent program of the MPMD pair, written in Python with pyAMReX.

It does the same as agent.cpp. Run it together with the C++ main program:

    mpiexec -n 2 ./MPMD_Main inputs : -n 2 python3 agent.py inputs

The dimension must match the C++ program; add dim=2 for a 2D build.
"""

import importlib
import sys

from mpi4py import MPI

dim = 3
for arg in sys.argv[1:]:
    if arg.startswith("dim="):
        dim = int(arg.split("=", 1)[1])
amr = importlib.import_module(f"amrex.space{dim}d")


def main():
    pp = amr.ParmParse("")
    n_cell = 32
    agent_max_grid_size = 8
    exists, value = pp.query_int("n_cell")
    if exists:
        n_cell = value
    exists, value = pp.query_int("agent_max_grid_size")
    if exists:
        agent_max_grid_size = value

    amr.Print(
        f"Agent (Python): {amr.ParallelDescriptor.NProcs()} of "
        f"{amr.MPMD_NProcs()} ranks in MPI_COMM_WORLD"
    )

    # The upper half of the main program's domain, with our own boxes
    domain = amr.Box(amr.IntVect(n_cell // 2), amr.IntVect(n_cell - 1))
    ba = amr.BoxArray(domain)
    ba.max_size(agent_max_grid_size)
    dm = amr.DistributionMapping(ba)

    mf = amr.MultiFab(ba, dm, 2, 0)

    copier = amr.MPMD_Copier(ba, dm, False)

    copier.recv(mf, 0, 1)  # receive component 0

    # Zero-copy NumPy (CPU) or CuPy (GPU) view of each box, indexed (x, y, z, comp)
    for mfi in mf:
        arr = mf.array(mfi).to_xp(copy=False)
        arr[..., 1] = 10.0 * arr[..., 0]

    amr.Print(
        f"Agent (Python): received component 0 with min = {mf.min(0)}, "
        f"max = {mf.max(0)}"
    )

    copier.send(mf, 1, 1)  # send component 1 back


# Split MPI_COMM_WORLD with mpi4py, and run AMReX on this program's ranks only.
amr.MPMD_Initialize_without_split(sys.argv[1:])
app_comm = MPI.COMM_WORLD.Split(amr.MPMD_AppNum(), amr.MPMD_MyProc())
amr.initialize_when_MPMD(sys.argv[1:], app_comm)

main()  # AMReX objects must be gone before finalize

amr.finalize()
amr.MPMD_Finalize()
