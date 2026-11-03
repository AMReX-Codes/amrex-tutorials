.. _coupling:

Coupling
========

This tutorial runs separate programs in one MPI job that exchange
``MultiFab`` data (``MPMD``). It requires MPI. Build only this topic with
``-DAMReXTutorials_TOPICS=Coupling``; see :doc:`GettingStarted`.

.. _coupling-mpmd:

MPMD
----

:Level: ★★★
:Source: `Coupling/MPMD <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Coupling/MPMD>`__
:Key APIs: ``MPMD::Initialize``, ``MPMD::Finalize``, ``MPMD::Copier``,
           ``Copier::send``, ``Copier::recv``,
           ``amrex::Initialize(argc, argv, true, comm)``;
           Python: ``MPMD_Initialize_without_split``, ``initialize_when_MPMD``,
           ``MPMD_Copier``

Two separately built programs run in the same MPI job (multiple program,
multiple data) and exchange ``MultiFab`` data with ``amrex::MPMD``. The
main program (``main.cpp``) owns the domain :math:`[0, n)^D`, fills
component 0 with a Gaussian, sends it, and receives component 1. The agent,
in C++ (``agent.cpp``) or Python with pyAMReX (``agent.py``), owns the
upper half :math:`[n/2, n)^D` with its own boxes, receives component 0,
sets component 1 to 10 times component 0, and sends it back. The main
program checks the round trip.

``MPMD::Initialize`` splits ``MPI_COMM_WORLD`` into one communicator per
program, and AMReX is initialized on it:

.. literalinclude:: ../../Coupling/MPMD/main.cpp
   :language: cpp
   :start-at: // Split MPI_COMM_WORLD into one communicator per program
   :end-at: amrex::Initialize(argc, argv, true, comm);
   :dedent:

``MPMD::Copier`` exchanges the BoxArrays and DistributionMappings of the two
programs once; ``send`` and ``recv`` then copy components on the
intersection of the two BoxArrays. The programs may use different boxes and
numbers of ranks, but must use the same dimension and precision. With
``Copier(ba, dm, true)`` on one side and ``Copier(true)`` on the other, the
second program adopts the BoxArray of the first.

.. literalinclude:: ../../Coupling/MPMD/main.cpp
   :language: cpp
   :start-at: // The Copier exchanges BoxArray and DistributionMapping
   :end-at: copier.recv(mf, 1, 1);
   :dedent:

The agent matches each ``send`` with a ``recv`` and vice versa:

.. literalinclude:: ../../Coupling/MPMD/agent.cpp
   :language: cpp
   :start-at: MPMD::Copier copier(ba, dm);
   :end-at: copier.send(mf, 1, 1);
   :dedent:

The Python agent splits the communicator with mpi4py:

.. literalinclude:: ../../Coupling/MPMD/agent.py
   :language: python
   :start-at: # Split MPI_COMM_WORLD with mpi4py
   :end-at: amr.MPMD_Finalize()

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target Coupling_MPMD_Main Coupling_MPMD_Agent -j 8
   cd build/Coupling/MPMD
   mpiexec -n 2 ./MPMD_Main inputs : -n 2 ./MPMD_Agent inputs
   mpiexec -n 2 ./MPMD_Main inputs : -n 2 python3 agent.py inputs

The Python agent needs pyAMReX built with MPI, and mpi4py; add ``dim=2``
for a 2D build. With ``-DAMReXTutorials_PYTHON=ON``, CTest also runs the
Python agent. With GNU Make, build the two programs with ``make PROG=main``
and ``make PROG=agent``.

.. code-block:: text

   Main: 2 of 4 ranks in MPI_COMM_WORLD
   Agent: 2 of 4 ranks in MPI_COMM_WORLD
   Agent: received component 0 with min = 1, max = 1.92937572
   Main: max error of the round trip = 0
   PASSED
