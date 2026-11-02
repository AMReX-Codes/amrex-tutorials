.. _particles:

Particles
=========

These tutorials use AMReX particle containers: pure struct-of-arrays (SoA)
particle data distributed over the same grids as mesh data, particle-mesh
deposition and gather, sorting and filtering within tiles, neighbor lists,
and a complete particle-in-cell code. All of them require AMReX built with
particles (``AMReX_PARTICLES``, on by default); ElectrostaticPIC also needs
the linear solvers. NeighborList and ElectrostaticPIC run in 2D and 3D only.
Build only this topic with ``-DAMReXTutorials_TOPICS=Particles``; see
:doc:`GettingStarted`.

.. _particles-basics:

Basics
------

:Level: ★
:Source: `Particles/Basics <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Particles/Basics>`__
:Key APIs: ``ParticleContainerPureSoA``, ``SetSoACompileTimeNames``,
   ``AddRealComp``, ``DefineAndReturnParticleTile``, ``ParticleTileData``,
   ``SetParticleIDandCPU``, ``ParIter``, ``Redistribute``,
   ``ParticleReduce``, ``Increment``, ``WritePlotFile``, ``Checkpoint``,
   ``Restart``

A first particle code. Each particle stores its position, velocity and
weight as separate arrays; the first ``AMREX_SPACEDIM`` real components are
the positions. A component ``dist`` (distance traveled) is added at run
time. Particles are created on the device with ``nppc`` random particles
per cell, move in straight lines through a periodic domain, and are moved to
their new grid and MPI rank by ``Redistribute``. ``ParticleReduce`` computes
the total weight, kinetic energy and path length in one pass. The program
writes a plotfile with the particle count per cell and the particles, and a
particle checkpoint that it reads back with ``Restart``.

The test checks that the number of particles, total weight and kinetic
energy are conserved, that the total path length is
:math:`n_\text{steps}\,\Delta t \sum |v|`, and that the restarted particles
match.

The components are named with an enum; the container stores all of them,
including the positions and the id/cpu word, as separate arrays.

.. literalinclude:: ../../Particles/Basics/main.cpp
   :language: cpp
   :start-at: // Real components of each particle.
   :end-before: struct Stats

Each tile is filled directly on the device. A block of unique ids is
reserved per tile, and ``ParallelForRNG`` provides a random engine.

.. literalinclude:: ../../Particles/Basics/main.cpp
   :language: cpp
   :start-at: // Reserve a block of unique particle ids
   :end-at: });
   :dedent:

The push loops over the particle tiles with ``ParIter``. Particles move
less than one cell per step, so ``Redistribute`` only communicates with
neighboring ranks.

.. literalinclude:: ../../Particles/Basics/main.cpp
   :language: cpp
   :start-at: for (PC::ParIterType pti(pc, lev); pti.isValid(); ++pti)
   :end-at: pc.Redistribute(0, -1, 0, 1);
   :dedent:

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target Particles_Basics -j 8
   cd build/Particles/Basics
   mpiexec -n 2 ./Basics inputs

Use ``n_cell``, ``nppc`` and ``nsteps`` for a larger problem, and
``particles.do_tiling=1`` for better OpenMP performance on CPUs. The
particles in ``plt00020/particles`` can be viewed with yt or ParaView (see
:doc:`VisualizationAndDebugging`).

.. _particles-particlemesh:

ParticleMesh
------------

:Level: ★★
:Source: `Particles/ParticleMesh <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Particles/ParticleMesh>`__
:Key APIs: ``ParticleToMesh``, ``MeshToParticle``,
   ``ParticleInterpolator::Linear``, ``InitRandom``, ``ReduceSum``,
   ``ReduceMax``, ``WritePlotFile``

Deposit a particle quantity onto the mesh and gather a mesh field to the
particles with cloud-in-cell (linear) weights. In 1D, a particle at
:math:`x` between the cell centers :math:`x_i` and :math:`x_{i+1}`, with
:math:`f = (x - x_i)/\Delta x`, gives a fraction :math:`1-f` of its mass to
cell :math:`i` and :math:`f` to cell :math:`i+1`; the gather uses the same
weights. In 2D and 3D the weights are products of the 1D weights.

The particles and the mesh use different ``BoxArray``\ s
(``particle_max_grid_size`` and ``max_grid_size``); ``ParticleToMesh`` and
``MeshToParticle`` copy data between the two layouts. The test checks that
the deposited mass equals the total particle mass and that the gather
reproduces a linear field to round-off.

The deposition callable runs once per particle and adds its weighted
contributions with atomic adds. Contributions to ghost cells are summed
into the valid cells that own them, including across periodic boundaries.

.. literalinclude:: ../../Particles/ParticleMesh/main.cpp
   :language: cpp
   :start-at: // Deposit the mass density with cloud-in-cell
   :end-before: const amrex::Real deposited_mass
   :dedent:

The gather needs a mesh field with filled ghost cells. The first callable
reads the mesh, the second accumulates each weighted value on the particle.

.. literalinclude:: ../../Particles/ParticleMesh/main.cpp
   :language: cpp
   :start-at: // Gather the field to the particles
   :end-before: // Largest error of the gathered field.
   :dedent:

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target Particles_ParticleMesh -j 8
   cd build/Particles/ParticleMesh
   mpiexec -n 2 ./ParticleMesh inputs

The plotfile ``plt00000`` holds the density ``rho``, the field, and the
particles.

.. _particles-sorting:

Sorting
-------

:Level: ★★
:Source: `Particles/Sorting <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Particles/Sorting>`__
:Key APIs: ``DenseBins``, ``getParticleCell``, ``SortParticlesByCell``,
   ``SortParticlesByBin``, ``partitionParticles``, ``filterParticles``,
   ``transformParticles``, ``DefineAndReturnParticleTile``

Organize the particles of a tile:

#. ``DenseBins`` builds a counting sort of the particles by cell without
   moving them. ``permutationPtr()`` lists particle indices bin by bin and
   ``offsetsPtr()[b]`` is where bin ``b`` starts, so a kernel over cells can
   loop over the particles of each cell.
#. ``SortParticlesByCell`` and ``SortParticlesByBin`` reorder the particle
   data so that particles in the same cell (or group of cells) are adjacent
   in memory, which improves memory access in deposition kernels.
#. ``partitionParticles`` moves the particles that satisfy a predicate (here:
   speed above ``vcut``) to the front of the tile.
#. ``filterParticles`` copies the particles that satisfy a predicate to
   another tile, and ``transformParticles`` copies particles while changing
   them (here: reversing the velocity).

The test checks each step: per-cell counts against
``ParticleContainer::Increment``, sorted order, the partition, the number of
filtered particles, and the reversed velocities.

Binning by cell and looping over the particles of each cell:

.. literalinclude:: ../../Particles/Sorting/main.cpp
   :language: cpp
   :start-at: amrex::DenseBins<ConstPTD> bins;
   :end-before: amrex::MultiFab::Subtract(count_bins
   :dedent:

Copying a subset of the particles, and transformed copies, into other
containers. The destination tiles must be large enough.

.. literalinclude:: ../../Particles/Sorting/main.cpp
   :language: cpp
   :start-at: PC slow_pc(geom, dm, ba);
   :end-before: const amrex::Long nslow = slow_pc
   :dedent:

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target Particles_Sorting -j 8
   cd build/Particles/Sorting
   mpiexec -n 2 ./Sorting inputs

.. _particles-neighborlist:

NeighborList
------------

:Level: ★★★
:Source: `Particles/NeighborList <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Particles/NeighborList>`__
:Key APIs: ``NeighborParticleContainerPureSoA``, ``fillNeighbors``,
   ``updateNeighbors``, ``buildNeighborList``, ``NeighborList::data``,
   ``getNeighbors``, ``Redistribute``, ``ParticleReduce``, ``RandomNormal``

A small molecular dynamics code with short-range forces. Particles are soft
spheres of diameter :math:`\sigma` with the potential
:math:`U(r) = \frac{k}{2}(\sigma - r)^2` for :math:`r < \sigma`.

``MDParticleContainer`` derives from ``NeighborParticleContainerPureSoA``.
``fillNeighbors`` copies particles within one cell of each tile (including
periodic images) into that tile, after its real particles.
``buildNeighborList`` then lists, for each real particle, the nearby real
and neighbor particles that pass a ``CheckPair`` test. The list radius is
one cell, larger than :math:`\sigma`, so the list stays valid for
``num_rebuild`` steps; in between, ``updateNeighbors`` only refreshes the
positions of the neighbor copies.

The test checks that the total momentum is conserved to round-off (pair
forces are equal and opposite, also across tiles and MPI ranks) and that the
relative change in total energy is below 1%.

The time loop uses velocity Verlet and rebuilds the neighbor list every
``num_rebuild`` steps:

.. literalinclude:: ../../Particles/NeighborList/main.cpp
   :language: cpp
   :start-at: // Velocity Verlet time stepping.
   :end-at: pc.kickDrift(0.5_rt*dt, 0.0_rt);
   :dedent:

The force kernel loops over the neighbors of each real particle:

.. literalinclude:: ../../Particles/NeighborList/MDParticleContainer.cpp
   :language: cpp
   :start-at: // For each real particle i, the list holds
   :end-at: });
   :dedent:

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target Particles_NeighborList -j 8
   cd build/Particles/NeighborList
   mpiexec -n 2 ./NeighborList inputs

Use ``n_cell`` and ``nsteps`` for a larger run. The time step is
``cfl*sigma/max_speed``, and ``2*num_rebuild*cfl*sigma`` must be smaller
than the skin ``1 - sigma``.

.. _particles-electrostaticpic:

ElectrostaticPIC
----------------

:Level: ★★★ / Blueprint
:Source: `Particles/ElectrostaticPIC <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Particles/ElectrostaticPIC>`__
:Key APIs: ``ParticleContainerPureSoA``, ``ParticleToMesh``,
   ``MeshToParticle``, ``ParticleInterpolator::Linear``, ``MLPoisson``,
   ``MLMG``, ``ParReduce``, ``Redistribute``, ``WriteSingleLevelPlotfile``,
   ``WritePlotFile``

An electrostatic particle-in-cell (PIC) code for electrons moving through a
fixed, neutralizing ion background in a periodic box, organized to be copied
as the starting point of a PIC code. Each step:

#. Deposit the electron charge with cloud-in-cell weights
   (``ParticleToMesh``) and add the ion charge.
#. Solve :math:`\nabla^2\phi = -\rho/\epsilon_0` with ``MLPoisson`` and
   ``MLMG``, and compute :math:`E = -\nabla\phi` at cell centers with
   centered differences.
#. Interpolate :math:`E` to the particles (``MeshToParticle``).
#. Leapfrog push: :math:`v \mathrel{+}= (q/m) E\,\Delta t`, then
   :math:`x \mathrel{+}= v\,\Delta t`, then ``Redistribute``.

The problem is a Langmuir (plasma) oscillation. In normalized units
(:math:`\epsilon_0 = m = n_0 = 1`, :math:`q = -1`) the plasma frequency is
:math:`\omega_p = \sqrt{n_0 q^2/(\epsilon_0 m)} = 1`. The electrons start
at rest, displaced by :math:`\xi_0 \sin(kx)`, and linear theory gives
:math:`E_x = -(q n_0 \xi_0/\epsilon_0) \sin(kx) \cos(\omega_p t)`. The code
measures the frequency from the zero crossings of the :math:`\sin(kx)` mode
of :math:`E_x`, and checks it, the initial amplitude, and the total (field
plus kinetic) energy against theory to 2%.

Code structure
^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - File
     - Contents
   * - ``PICParticleContainer.H/.cpp``
     - Pure SoA particle species: initialization, charge deposition, field
       gather, push
   * - ``FieldSolver.H/.cpp``
     - Periodic Poisson solve with MLMG and :math:`E = -\nabla\phi`
   * - ``main.cpp``
     - Setup, time loop, diagnostics and plotfiles

The time loop in ``main.cpp``:

.. literalinclude:: ../../Particles/ElectrostaticPIC/main.cpp
   :language: cpp
   :start-at: // Charge deposition, field solve and gather at time t^n.
   :end-before: const amrex::Real field_energy
   :dedent:

Charge deposition in ``PICParticleContainer.cpp``:

.. literalinclude:: ../../Particles/ElectrostaticPIC/PICParticleContainer.cpp
   :language: cpp
   :start-at: // ParticleToMesh zeros rho
   :end-before: void
   :lines: 1-14
   :dedent:

The field solver sets up ``MLPoisson`` and ``MLMG`` once, in its
constructor, and reuses them every step:

.. literalinclude:: ../../Particles/ElectrostaticPIC/FieldSolver.cpp
   :language: cpp
   :start-at: FieldSolver::FieldSolver
   :end-before: void
   :lines: 1-18

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target Particles_ElectrostaticPIC -j 8
   cd build/Particles/ElectrostaticPIC
   mpiexec -n 2 ./ElectrostaticPIC inputs

``n_cell`` and ``nppc`` take three values; 2D uses the first two. The
plotfiles ``plt00000``, ``plt00050`` and ``plt00100`` contain ``rho``,
``phi``, the electric field and the electrons.

.. code-block:: text

   Initial E mode amplitude: 0.009936866131 (linear theory 0.01)
   Measured frequency: 0.9956005796 (plasma frequency 1)
   Relative errors: frequency 0.004399420394, amplitude 0.006313386911, max energy change 0.001909800398
   PASSED

To add species, create one container per species and sum their charge
densities before the solve. For other boundary conditions, change
``setDomainBC`` in ``FieldSolver`` and the particle boundary handling; for
periodic problems, :ref:`fft-poisson` can replace MLMG when AMReX is built
with FFT support.
