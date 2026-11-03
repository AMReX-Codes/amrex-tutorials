.. _base:

Base
====

The Base tutorials cover the core of AMReX: initialization, runtime
parameters, the mesh data structures (``Box``, ``BoxArray``,
``DistributionMapping``, ``Geometry``, ``MultiFab``), GPU-portable kernels
with ``ParallelFor``, ghost cells and boundary conditions, reductions, GPU
containers, and I/O. They need no optional AMReX components; ExternalMPI
requires MPI. Build only this topic with ``-DAMReXTutorials_TOPICS=Base``
(add ``-DAMReX_FFT=OFF`` if FFTW is not installed); see :doc:`GettingStarted`.

Start with HelloWorld, MultiFab and HeatEquation, in that order.

.. _base-helloworld:

HelloWorld
----------

:Level: ★
:Source: `Base/HelloWorld <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/HelloWorld>`__
:Key APIs: ``amrex::Initialize``, ``amrex::Finalize``, ``amrex::Print``,
           ``amrex::AllPrint``

The smallest AMReX program. ``amrex::Initialize`` sets up MPI and the GPU
device; ``amrex::Finalize`` releases them. AMReX objects must be destroyed
before ``Finalize``, hence the extra scope. ``amrex::Print`` writes from
the I/O rank only, ``amrex::AllPrint`` from every rank.

.. literalinclude:: ../../Base/HelloWorld/main.cpp
   :language: cpp
   :start-at: int main

Run:

.. code-block:: console

   cmake --build build --target Base_HelloWorld -j 8
   cd build/Base/HelloWorld
   mpiexec -n 2 ./HelloWorld

.. _base-parmparse:

ParmParse
---------

:Level: ★
:Source: `Base/ParmParse <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/ParmParse>`__
:Key APIs: ``ParmParse::query``, ``get``, ``queryarr``, ``queryAdd``,
           ``queryWithParser``, ``querytable``,
           ``query_enum_case_insensitive``, ``AMREX_ENUM``,
           ``ParmParse::prettyPrintUnusedInputs``

Reads runtime parameters from an inputs file and the command line:
scalars, strings, arrays, math expressions, enums, 2D tables, TOML-like
``[section]`` headers, and unused (often misspelled) inputs. Arguments
after the inputs file, e.g. ``nsteps=200``, override the file.

- ``query`` returns whether the name was found and leaves the variable
  unchanged otherwise; ``get`` aborts if the name is missing; ``queryAdd``
  also stores the default in the database.
- Expressions without spaces (``dt = 1.0/64``) are evaluated by ``query``;
  expressions with spaces need ``queryWithParser``. Unknown symbols are
  looked up in the database, here with the ``my_constants.`` prefix.
- ``[physics]`` starts a table: ``diffusivity = 0.1`` below it is
  ``physics.diffusivity``.

The inputs file:

.. literalinclude:: ../../Base/ParmParse/inputs
   :language: none
   :start-at: nsteps = 100

Reading an enum, a table, and the entries of a section:

.. literalinclude:: ../../Base/ParmParse/main.cpp
   :language: cpp
   :start-at: // Enum
   :end-before: amrex::Print() << "physics.diffusivity
   :dedent:

Run:

.. code-block:: console

   cmake --build build --target Base_ParmParse -j 8
   cd build/Base/ParmParse
   ./ParmParse inputs nsteps=200

The output ends with the misspelled ``physics.difusivity = 0.2`` listed as
an unused input.

.. _base-multifab:

MultiFab
--------

:Level: ★
:Source: `Base/MultiFab <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/MultiFab>`__
:Key APIs: ``Box``, ``BoxArray::maxSize``, ``DistributionMapping``,
           ``Geometry``, ``RealBox``, ``MultiFab``, ``ParallelFor(mf, ...)``,
           ``MFIter``, ``TilingIfNotGPU``, ``Array4``, ``FillBoundary``,
           ``MultiFab::sum/norm0/norm1/norm2``, ``WriteSingleLevelPlotfile``

Builds the basic data structures: a ``Box`` for the index space of the
domain, a ``BoxArray`` that chops it into grids, a ``DistributionMapping``
that assigns grids to MPI ranks, and a ``Geometry`` for the physical,
periodic domain. A ``MultiFab`` holds data on the grids owned by each rank,
plus ghost cells.

.. literalinclude:: ../../Base/MultiFab/main.cpp
   :language: cpp
   :start-at: // Box: the index space
   :end-before: amrex::Print() << "Domain "
   :dedent:

The same function is written into two components in two ways.
``ParallelFor(mf, ...)`` with ``mf.arrays()`` runs one kernel over all
local boxes; it is the simplest, and on GPUs it is efficient for many small
boxes. An ``MFIter`` loop gives one ``Array4`` per box (per tile on CPU,
with ``TilingIfNotGPU()`` and OpenMP); use it when work is organized per
box, e.g. with temporary ``FArrayBox``\ es or several kernels on the same
box.

.. literalinclude:: ../../Base/MultiFab/main.cpp
   :language: cpp
   :start-at: auto const& ma = mf.arrays();
   :end-at: });
   :dedent:

.. literalinclude:: ../../Base/MultiFab/main.cpp
   :language: cpp
   :start-at: for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU());
   :end-at: mf.FillBoundary(geom.periodicity());
   :dedent:

The ``MFIter`` loop is preceded by ``#pragma omp parallel if
(amrex::Gpu::notInLaunchRegion())`` so that tiles are distributed over
OpenMP threads in CPU builds. ``FillBoundary`` fills ghost cells from
neighboring boxes and periodic images; :ref:`base-communication` covers it
and the other communication functions in more detail. The program then
checks valid and ghost cells against the exact function and writes the
plotfile ``plt_multifab``.

Run:

.. code-block:: console

   cmake --build build --target Base_MultiFab -j 8
   cd build/Base/MultiFab
   mpiexec -n 2 ./MultiFab inputs

.. _base-heatequation:

HeatEquation
------------

:Level: ★
:Source: `Base/HeatEquation <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/HeatEquation>`__
:Key APIs: ``MultiFab``, ``Geometry``, ``FillBoundary``,
           ``ParallelFor(mf, ...)``, ``const_arrays``,
           ``WriteSingleLevelPlotfile``

Solves :math:`\partial\phi/\partial t = D\nabla^2\phi` with forward Euler
on the periodic unit domain. The initial condition is
:math:`\phi = 1 + \prod_d \sin(2\pi x_d)`; its sine mode decays as
:math:`\exp(-4\pi^2 d D t)` in :math:`d` dimensions, and the program
compares the result with this exact solution. The time step is
:math:`\Delta t = \mathrm{cfl}\,\Delta x^2/(2 d D)`.

Each step fills the ghost cells and applies the standard second-order
Laplacian:

.. literalinclude:: ../../Base/HeatEquation/main.cpp
   :language: cpp
   :start-at: // Fill ghost cells, here all
   :end-before: time += dt_step;
   :dedent:

Run:

.. code-block:: console

   cmake --build build --target Base_HeatEquation -j 8
   cd build/Base/HeatEquation
   mpiexec -n 2 ./HeatEquation inputs

Add ``plot_int=10`` to write a plotfile every 10 steps, and e.g.
``n_cell=128 stop_time=0.05`` for a larger run. To see where the time goes,
build with TinyProfiler; see :ref:`profiling`.

.. _base-communication:

Communication
-------------

:Level: ★★
:Source: `Base/Communication <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/Communication>`__
:Key APIs: ``FillBoundary``, ``FillBoundary_nowait``,
           ``FillBoundary_finish``, ``ParallelCopy``, ``ParallelAdd``,
           ``SumBoundary``, ``Geometry::periodicity``

The three basic ways to exchange ``MultiFab`` data between boxes and MPI
ranks. Each valid cell holds a unique number, so every result is checked
exactly.

``FillBoundary`` copies valid data into the ghost cells of neighboring
boxes, and with ``geom.periodicity()`` also across periodic boundaries. It
can fill fewer ghost cells than the ``MultiFab`` has. The non-blocking
version returns after starting the communication, so that work that does
not need the ghost cells can be done before ``FillBoundary_finish``:

.. literalinclude:: ../../Base/Communication/main.cpp
   :language: cpp
   :start-at: a.FillBoundary_nowait(geom.periodicity());
   :end-at: a.FillBoundary_finish();
   :dedent:

``ParallelCopy`` copies data between ``MultiFab``\ s with different
``BoxArray``\ s and ``DistributionMapping``\ s, e.g. after regridding.
The number of ghost cells of the source and the destination can be given,
here 0 and 1. The destination ``BoxArray`` above is shifted and partly
outside the domain; with ``geom.periodicity()``, those cells are filled from
their periodic images. ``ParallelAdd`` adds instead of copying:

.. literalinclude:: ../../Base/Communication/main.cpp
   :language: cpp
   :start-at: dst.ParallelCopy(a,
   :end-at: dst.ParallelAdd(a,
   :dedent:

``SumBoundary`` is the reverse of ``FillBoundary``: values in ghost cells
are added to the valid cells they overlap, as in deposition, e.g. of
particles near a box boundary. Here every valid cell deposits 1 into the
:math:`3^d` cells around it, itself included, some of them ghost cells,
and after ``SumBoundary`` every valid cell holds :math:`3^d`. The ghost
cells are not updated; call ``FillBoundary`` afterwards if they are needed.

Run:

.. code-block:: console

   cmake --build build --target Base_Communication -j 8
   cd build/Base/Communication
   mpiexec -n 2 ./Communication inputs

.. _base-boundaryconditions:

BoundaryConditions
------------------

:Level: ★★
:Source: `Base/BoundaryConditions <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/BoundaryConditions>`__
:Key APIs: ``BCRec``, ``BCType::ext_dir/foextrap/reflect_even/int_dir``,
           ``PhysBCFunct``, ``GpuBndryFuncFab``, ``FillBoundary``,
           ``MFIter``, ``ParallelFor(box, [=] (IntVect const&) ...)``

Solves the heat equation on a non-periodic domain: Dirichlet values at
:math:`x = 0` and :math:`x = 1`, zero gradient in :math:`y`, periodic in
:math:`z`. The solution
:math:`\phi = \phi_{lo} + (\phi_{hi} - \phi_{lo})x + e^{-\pi^2 D t}\sin(\pi x)`
relaxes to the linear steady state and is checked against this exact
solution.

Boundary conditions are given per component as a ``BCRec``. Periodic
directions use ``int_dir`` and must match the ``Geometry``.

.. literalinclude:: ../../Base/BoundaryConditions/main.cpp
   :language: cpp
   :start-at: amrex::Vector<amrex::BCRec> bcs(1);
   :end-before: amrex::Box domain(

Ghost cells are filled in two steps. ``FillBoundary`` fills those that
overlap other boxes, including periodic images. ``PhysBCFunct`` fills those
outside the domain: AMReX handles ``foextrap``, ``reflect_even``,
``reflect_odd`` and ``hoextrap``, and calls the user functor for every such
cell to set the ``ext_dir`` values.

.. literalinclude:: ../../Base/BoundaryConditions/main.cpp
   :language: cpp
   :start-at: struct DirichletFill
   :end-at: };

.. literalinclude:: ../../Base/BoundaryConditions/main.cpp
   :language: cpp
   :start-at: DirichletFill dirichlet_fill{};
   :end-at: (geom, bcs, amrex::GpuBndryFuncFab<DirichletFill>(dirichlet_fill));
   :dedent:

.. literalinclude:: ../../Base/BoundaryConditions/main.cpp
   :language: cpp
   :start-at: phi_old.FillBoundary(geom.periodicity());
   :end-at: physbc(phi_old
   :dedent:

For ``ext_dir``, the ghost cell holds the value on the boundary face, so the
Laplacian uses a half-cell distance next to that boundary.

Run:

.. code-block:: console

   cmake --build build --target Base_BoundaryConditions -j 8
   cd build/Base/BoundaryConditions
   mpiexec -n 2 ./BoundaryConditions inputs

Try ``phi_lo=2 phi_hi=-1``, or ``stop_time=0.5 plot_int=100`` to follow
the approach to the steady state.

.. _base-checkpointrestart:

CheckpointRestart
-----------------

:Level: ★★
:Source: `Base/CheckpointRestart <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/CheckpointRestart>`__
:Key APIs: ``PreBuildDirectorHierarchy``, ``MultiFabFileFullPrefix``,
           ``VisMF::Write``, ``VisMF::Read``, ``BoxArray::writeOn/readFrom``,
           ``ParallelDescriptor::ReadAndBcastFile``, ``PlotFileData``,
           ``ParallelCopy``

Writes checkpoints of a heat equation run and restarts from them. A
checkpoint is a directory with a ``Header`` (step, time and ``BoxArray``),
written by the I/O rank, and the ``MultiFab`` data, written by every rank
with ``VisMF``.

.. literalinclude:: ../../Base/CheckpointRestart/main.cpp
   :language: cpp
   :start-at: void WriteCheckpoint
   :end-before: // Read a checkpoint

On restart, the ``BoxArray`` comes from the checkpoint and a new
``DistributionMapping`` is made, so a run may restart on a different number
of MPI ranks.

.. literalinclude:: ../../Base/CheckpointRestart/main.cpp
   :language: cpp
   :start-at: void ReadCheckpoint
   :end-before: // Advance the periodic

With ``test_restart=1``, the program restarts from the first checkpoint,
reruns the remaining steps, and checks that the result is bitwise identical
to the uninterrupted run. It also reads the final plotfile back with
``PlotFileData``; the data come with their own ``DistributionMapping``, so
``ParallelCopy`` moves them into the layout of another ``MultiFab``.

.. literalinclude:: ../../Base/CheckpointRestart/main.cpp
   :language: cpp
   :start-at: amrex::PlotFileData pf(plotfile);
   :end-at: phi_pf.ParallelCopy
   :dedent:

Run:

.. code-block:: console

   cmake --build build --target Base_CheckpointRestart -j 8
   cd build/Base/CheckpointRestart
   mpiexec -n 2 ./CheckpointRestart inputs
   mpiexec -n 1 ./CheckpointRestart inputs restart=chk00010
   mpiexec -n 2 ./CheckpointRestart inputs test_restart=1

.. _base-reductions:

Reductions
----------

:Level: ★★
:Source: `Base/Reductions <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/Reductions>`__
:Key APIs: ``ParReduce``, ``ReduceOps``, ``ReduceData``, ``Reducer``,
           ``ReduceOpSum/Min/Max``, ``MultiFab::sum/min/max/norm0/norm1/norm2``,
           ``Reduce::Sum``, ``ParallelAllReduce``,
           ``ParallelDescriptor::ReduceRealSum``

Computes sums, minima, maxima and counts over a ``MultiFab`` and a device
vector, on CPU and GPU, and combines them across MPI ranks. The data,
:math:`\phi = i + j + k`, has simple exact answers that the program checks.

- ``MultiFab`` member functions (``sum``, ``norm0``, ...) compute common
  norms and already include the MPI reduction.
- ``ParReduce`` computes several custom reductions over a ``MultiFab`` in
  one pass.
- ``ReduceOps``/``ReduceData``, or the more compact ``Reducer``, work inside
  ``MFIter`` loops or over arbitrary boxes and index ranges.
- ``Reduce::Sum/Min/Max`` reduce a plain array, e.g. a
  ``Gpu::DeviceVector``.

All except the ``MultiFab`` member functions are local to the MPI rank and
are followed by ``ParallelAllReduce`` or ``ParallelDescriptor::Reduce*``.

``ParReduce`` with a sum, a minimum and a maximum, followed by the MPI
reductions:

.. literalinclude:: ../../Base/Reductions/main.cpp
   :language: cpp
   :start-at: auto const& a = mf.const_arrays();
   :end-at: amrex::ParallelAllReduce::Max(vmax
   :dedent:

``Reducer`` with two sums of different types:

.. literalinclude:: ../../Base/Reductions/main.cpp
   :language: cpp
   :start-at: amrex::Reducer<amrex::TypeList
   :end-at: Result_t r = reducer.getResult();
   :dedent:

Run:

.. code-block:: console

   cmake --build build --target Base_Reductions -j 8
   cd build/Base/Reductions
   mpiexec -n 2 ./Reductions inputs

.. _base-gpucontainers:

GpuContainers
-------------

:Level: ★★
:Source: `Base/GpuContainers <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/GpuContainers>`__
:Key APIs: ``Gpu::DeviceVector``, ``Gpu::PinnedVector``, ``Gpu::copyAsync``,
           ``Gpu::streamSynchronize``, ``ParallelForRNG``, ``amrex::Random``,
           ``Scan::PrefixSum``, ``Scan::ExclusiveSum``, ``The_Arena``,
           ``The_Async_Arena``, ``The_Pinned_Arena``

Fills a ``Gpu::DeviceVector`` with random numbers, keeps the values above a
threshold with a prefix-sum based stream compaction, copies the results to
pinned host memory, and checks them against the same operations on the
host. The code is the same for CPU and GPU builds.

``Scan::PrefixSum`` takes one function that produces each input and one
that consumes each partial sum, so the compaction needs no temporary
arrays:

.. literalinclude:: ../../Base/GpuContainers/main.cpp
   :language: cpp
   :start-at: // Uniform random numbers
   :end-at: kept.resize(nkept);
   :dedent:

``DeviceVector`` memory comes from ``The_Arena()`` (device memory on GPU),
and ``PinnedVector`` from ``The_Pinned_Arena()`` (page-locked host memory,
needed for truly asynchronous copies).
``copyAsync`` may return before the copy is done, so synchronize the stream
before using the data:

.. literalinclude:: ../../Base/GpuContainers/main.cpp
   :language: cpp
   :start-at: // Copy to the host.
   :end-at: amrex::Gpu::streamSynchronize();
   :dedent:

``The_Async_Arena()`` is for temporaries: on GPU, memory freed there is
reused only after the work already queued on the stream is done.

Run:

.. code-block:: console

   cmake --build build --target Base_GpuContainers -j 8
   cd build/Base/GpuContainers
   ./GpuContainers inputs

.. _base-parser:

Parser
------

:Level: ★★
:Source: `Base/Parser <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/Parser>`__
:Key APIs: ``Parser``, ``IParser``, ``setConstant``, ``registerVariables``,
           ``compile<N>``, ``ParmParse::makeParser``

Evaluates functions given as strings at runtime, here an initial condition
and a tagging criterion from the inputs file, inside GPU kernels.
``amrex::Parser`` handles floating-point expressions and ``amrex::IParser``
integer ones. Expressions support ``+ - * / ^ **``, comparisons, ``and``,
``or``, ``if(cond, a, b)``, ``min``, ``max`` and the usual math functions;
``IParser`` also has floor division ``//``.

A parser is built on the host: set the constants, register the variables in
the order they will be passed, and ``compile<N>()`` with ``N`` variables.
``ParmParse::makeParser`` instead resolves unknown symbols from the
``ParmParse`` database.

.. literalinclude:: ../../Base/Parser/main.cpp
   :language: cpp
   :start-at: // 1. Parser built by hand
   :end-at: auto const ftag = iparser.compile<3>();
   :dedent:

The returned executor is a small object captured by value in device
lambdas; the ``Parser`` must stay alive while it is used.

.. literalinclude:: ../../Base/Parser/main.cpp
   :language: cpp
   :start-at: amrex::ParallelFor(phi, [=] AMREX_GPU_DEVICE
   :end-at: ta[b](i,j,k) = int(ftag(i, j, k));
   :dedent:

Run:

.. code-block:: console

   cmake --build build --target Base_Parser -j 8
   cd build/Base/Parser
   ./Parser inputs
   ./Parser inputs ic.function="x*y*exp(-t)"

.. _base-compiletimeoptions:

CompileTimeOptions
------------------

:Level: ★★★
:Source: `Base/CompileTimeOptions <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/CompileTimeOptions>`__
:Key APIs: ``ParallelFor(TypeList<CompileTimeOptions<...>...>, {...}, box, f)``,
           ``CompileTimeOptions``, ``if constexpr``,
           ``AMREX_ENUM``, ``BL_PROFILE``

Turns runtime options into compile-time constants inside kernels. A kernel
with two runtime options, the order of a Laplacian stencil and whether to
add an absorption term, is written once with runtime branches and once with
``CompileTimeOptions``. AMReX compiles a specialized kernel for every
combination and launches the one that matches the runtime values. The
program checks that both versions agree and times them with ``BL_PROFILE``
(see :ref:`profiling`).

.. literalinclude:: ../../Base/CompileTimeOptions/main.cpp
   :language: cpp
   :start-at: amrex::TypeList<amrex::CompileTimeOptions
   :end-at: });
   :dedent:
   :prepend: amrex::ParallelFor(

The kernel receives one extra argument per option, of type
``std::integral_constant<int, value>``, so ``opt.value`` can be used in
``if constexpr`` and as a template argument. Specialized kernels can drop
unused code and use fewer registers, which matters most for large GPU
kernels; the cost is compile time and binary size, which grow with the
number of combinations. ``AnyCTO`` does the same for other launch functions
such as ``ParallelForRNG``.

Run:

.. code-block:: console

   cmake --build build --target Base_CompileTimeOptions -j 8
   cd build/Base/CompileTimeOptions
   ./CompileTimeOptions inputs

For timings, build with TinyProfiler on and use e.g. ``n_cell=256
nrepeat=20 max_grid_size=128`` and an optimized GPU build.

.. _base-tagparallelfor:

TagParallelFor
--------------

:Level: ★★★
:Source: `Base/TagParallelFor <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/TagParallelFor>`__
:Key APIs: ``ParallelFor(Vector<Tag>, ncomp, f)``,
           ``ParallelFor(TagVector<Tag>, ncomp, f)``, ``TagVector``,
           ``Array4BoxOrientationTag``, ``Orientation``,
           ``OrientationIter``, ``adjCell``, ``BL_PROFILE``

Runs many small kernels as one launch. Filling the ghost cells outside the
domain boundary touches a thin box per domain face of every box at the
boundary: many small boxes of different sizes, each too small to fill a
GPU. The program fills them with a reflecting boundary condition in three
ways, checks that the results agree, and times them with ``BL_PROFILE``
(see :ref:`profiling`): one ``ParallelFor`` per box and face inside an
``MFIter`` loop, the baseline, not recommended on GPU for many small boxes
since no launch fills the device; one ``ParallelFor`` over a ``Vector`` of
tags, which packs all the boxes into a single launch; and the same with a
``TagVector`` built once and reused.

A tag is a small struct with a ``box()`` member (``size()`` for 1D work)
that carries what the kernel needs, here the ``Array4``, the thin box and
the face. The tags are built on the host:

.. literalinclude:: ../../Base/TagParallelFor/main.cpp
   :language: cpp
   :start-at: // One tag per (local box, domain face it touches)
   :end-before: // Largest difference

The kernel receives ``(i, j, k, n, tag)``. ``ParallelFor`` over a
``Vector`` of tags copies the tags to the device and synchronizes the
stream on every call, so for repeated use the ``TagVector`` is built once:

.. literalinclude:: ../../Base/TagParallelFor/main.cpp
   :language: cpp
   :start-at: // 3. The TagVector is built once
   :end-before: // Expected result
   :dedent:

On CPU, the launch is a loop over the tags that opens its own OpenMP
parallel region, so it must not be called inside one. AMReX
provides ``Array4BoxTag``, ``Array4BoxValTag``,
``Array4BoxOrientationTag``, ``Array4PairTag``, ``Array4CopyTag`` and
``VectorTag``, and any struct with ``box()`` or ``size()`` works. Use
``ParallelFor(mf, ...)`` (:ref:`base-multifab`) when the work is the same
on every box, and tags for irregular sets of boxes: physical boundary
conditions, corrections at coarse/fine faces, operations on some boxes
only.

Run:

.. code-block:: console

   cmake --build build --target Base_TagParallelFor -j 8
   cd build/Base/TagParallelFor
   mpiexec -n 2 ./TagParallelFor inputs

For timings, build with TinyProfiler on and an optimized GPU build, and
use e.g. ``n_cell=256 nrepeat=20 tiny_profiler.print_threshold=0``; the
difference between the methods grows with the number of boxes.

.. _base-sharedmemory:

SharedMemory
------------

:Level: ★★★
:Source: `Base/SharedMemory <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/SharedMemory>`__
:Key APIs: ``amrex::launch<nthreads>(nblocks, stream, f)``,
           ``amrex::launch<nthreads>(nblocks, shared_mem_bytes, stream, f)``,
           ``__shared__`` (CUDA, HIP), ``Gpu::Handler`` (SYCL),
           ``Gpu::Atomic::AddNoRet``, ``HostDevice::Atomic::Add``,
           ``BoxIndexer``, ``Gpu::Device::maxBlocksPerLaunch``,
           ``AMREX_USE_GPU``, ``AMREX_USE_SYCL``

A GPU kernel written by hand, with shared memory, and a separate CPU
implementation. ``ParallelFor`` is portable and the right choice almost
always; a raw launch is for the rare kernel whose threads must cooperate
within a block. The example bins the cells by the value of one field and
sums a signed integer field over each bin. Integer sums are exact, so the
two methods must agree to the bit, and they need a signed 64-bit
accumulator. With ``ParallelFor`` and one atomic add per cell, all threads
compete for the same few addresses in global memory:

.. literalinclude:: ../../Base/SharedMemory/main.cpp
   :language: cpp
   :start-at: void binned_sum_atomic
   :end-before: // 2. Each thread block

Atomic adds to the same address are serialized, so the fewer the bins, the
slower this version. The hand-written kernel gives each thread block
private ``int`` sums in shared memory, where contention is cheap, and adds
them to the global ``Long`` sums once per block, so the global atomics
number one per bin and block instead of one per cell. ``amrex::launch``
runs the kernel on ``nblocks`` blocks of ``nthreads`` threads; the threads
of a block walk the cells of the box in a grid-stride loop, with
``BoxIndexer`` converting the linear index to ``(i,j,k)``. With CUDA and
HIP the body uses ``threadIdx``, ``blockIdx``, ``__syncthreads()`` and a
static ``__shared__`` array:

.. literalinclude:: ../../Base/SharedMemory/main.cpp
   :language: cpp
   :start-at: amrex::launch<nthreads>(nblocks, amrex::Gpu::gpuStream(),
   :end-at: });
   :dedent:

With SYCL, the shared memory is dynamic with its size passed to
``amrex::launch``, and that overload hands the kernel a ``Gpu::Handler``
whose ``threadIdx()``, ``blockDim()``, ``gridDim()``, ``globalIdx()``,
``syncThreads()`` and ``sharedMemory()`` play these roles; atomics on the
shared memory take the ``local_space`` address space; see the
``AMREX_USE_SYCL`` branch in the source. There is no CPU version of
this launch. The CPU path, selected with ``#ifdef AMREX_USE_GPU``, is an
``MFIter`` loop with a vector of private sums per OpenMP thread, merged
after the parallel region; its cell loop is ``amrex::Loop``, not
``ParallelFor``, because it scatters into the sums:

.. literalinclude:: ../../Base/SharedMemory/main.cpp
   :language: cpp
   :start-at: // Private sums per OpenMP thread
   :end-at: } // end of merge

The program checks that both results agree and times them with
``BL_PROFILE`` (see :ref:`profiling`).

Run:

.. code-block:: console

   cmake --build build --target Base_SharedMemory -j 8
   cd build/Base/SharedMemory
   mpiexec -n 2 ./SharedMemory inputs

For timings, build with TinyProfiler on and use a GPU build with e.g.
``n_cell=128 max_grid_size=128 nbins=64 nrepeat=5
tiny_profiler.print_threshold=0``. The fewer the bins, the more the global
atomics collide, and the larger the speedup of the shared-memory version.

.. _base-externalmpi:

ExternalMPI
-----------

:Level: ★★
:Source: `Base/ExternalMPI <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Base/ExternalMPI>`__
:Key APIs: ``amrex::Initialize(argc, argv, build_parm_parse, comm, func)``,
           ``ParmParse::add``, ``ParmParse::addfile``,
           ``ParallelDescriptor::NProcs``

Uses AMReX inside an application that manages MPI itself. The application
calls ``MPI_Init``, splits ``MPI_COMM_WORLD`` into even and odd ranks, and
initializes AMReX on each group's communicator. AMReX uses a duplicate of
the given communicator and does not finalize MPI; the application calls
``MPI_Finalize`` after ``amrex::Finalize``.

With ``build_parm_parse = false``, AMReX ignores ``argc``/``argv``. The
callback runs before AMReX reads its own parameters, so it is the place to
set them in code.

.. literalinclude:: ../../Base/ExternalMPI/main.cpp
   :language: cpp
   :start-at: auto add_parameters = [&] ()
   :end-at: amrex::Initialize(argc, argv, false, group_comm, add_parameters);
   :dedent:

Run (requires ``AMReX_MPI=ON``):

.. code-block:: console

   cmake --build build --target Base_ExternalMPI -j 8
   cd build/Base/ExternalMPI
   mpiexec -n 4 ./ExternalMPI inputs
