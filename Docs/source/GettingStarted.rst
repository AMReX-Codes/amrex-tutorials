.. _getting_started:

Getting Started
===============

The tutorials need a C++20 compiler and, for MPI builds (the default), an
MPI library. GPU builds need CUDA, ROCm (HIP) or oneAPI (SYCL). CPU builds
with CMake also need FFTW (or ``-DAMReX_FFT=OFF``, which skips the FFT
tutorials). Every tutorial except the Python scripts builds with CMake or
with GNU Make.

.. code-block:: console

   git clone https://github.com/AMReX-Codes/amrex-tutorials.git
   cd amrex-tutorials

Building with CMake
-------------------

From the top of the repository, CMake configures all tutorials at once.
By default, AMReX is downloaded and built as part of the tutorials
(a superbuild with ``FetchContent``):

.. code-block:: console

   cmake -S . -B build
   cmake --build build -j 8
   ctest --test-dir build

Each tutorial ``<Topic>/<Name>`` is the target ``<Topic>_<Name>``. Its
executable, named ``<Name>``, and its ``inputs`` file are placed in
``build/<Topic>/<Name>``. The Applications have one target per problem,
e.g. ``Applications_CNS_Sod``; see their READMEs.

.. code-block:: console

   cmake --build build --target Base_HeatEquation -j 8
   cd build/Base/HeatEquation
   mpiexec -n 2 ./HeatEquation inputs

Runtime parameters can be overridden on the command line after the inputs
file, e.g. ``./HeatEquation inputs n_cell=64``.

Options
^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Option
     - Meaning
   * - ``-DAMReXTutorials_TOPICS="Base;FFT"``
     - Build only some topics
   * - ``-DAMReX_SPACEDIM=2``
     - Build for 1, 2 or 3 dimensions (default 3)
   * - ``-DAMReXTutorials_DIM=2``
     - Dimension to use when AMReX is built for several (default the last)
   * - ``-DAMReX_GPU_BACKEND=CUDA -DCMAKE_CUDA_ARCHITECTURES=native``
     - GPU build (also ``HIP``, ``SYCL``)
   * - ``-DAMReX_OMP=ON``, ``-DAMReX_MPI=OFF``
     - OpenMP on, MPI off
   * - ``-DAMReX_PRECISION=SINGLE``
     - Single precision ``amrex::Real``
   * - ``-DCMAKE_BUILD_TYPE=Debug``
     - Debug build (default ``Release``); see :doc:`VisualizationAndDebugging`
   * - ``-DFETCHCONTENT_SOURCE_DIR_AMREX=/path/to/amrex``
     - Use a local AMReX source tree instead of downloading it
   * - ``-DAMReX_ROOT=/path/to/install``
     - Use an installed AMReX
   * - ``-DAMReX_TESTING=ON``
     - Trap floating-point exceptions in the tests (CPU, GCC or Clang)
   * - ``-DAMReX_TINY_PROFILE=ON``
     - Build with TinyProfiler; see :ref:`profiling`
   * - ``-DAMReX_SUNDIALS=ON``
     - Build the :ref:`SUNDIALS <timeintegration-sundials>` tutorial
       (SUNDIALS is downloaded unless ``SUNDIALS_ROOT`` is given)
   * - ``-DAMReX_FORTRAN_INTERFACES=ON``
     - Build the :doc:`Fortran` tutorials (needs a Fortran compiler)
   * - ``-DAMReXTutorials_PYTHON=ON``
     - Test the :doc:`Python` tutorials (needs the ``amrex`` Python module)
   * - ``-DAMReXTutorials_TORCH=ON``
     - Build the :doc:`ML` tutorial (needs LibTorch)

In a superbuild, AMReX is configured with EB, FFT, particles and linear
solvers enabled. Any AMReX CMake option, e.g. ``-DAMReX_FFT=OFF``,
can be passed on the command line. With an installed AMReX, the tutorials
use its GPU backend. Tutorials that need a component the AMReX build does
not have (e.g. EB, FFT, MPI or SUNDIALS), or that do not support the chosen
dimension, are skipped with a message at configure time.

Tests
^^^^^

Each tutorial is registered as a CTest test that runs it with its
``inputs`` file, on 2 MPI ranks in MPI builds (the MPMD test needs 4;
with Open MPI, add ``-DMPIEXEC_PREFLAGS=--oversubscribe`` if needed).
Tutorials with an exact or reference answer check it and return a non-zero
exit code if the error is too large. With GCC or Clang on CPUs, the tests
trap floating-point exceptions (``amrex.fpe_trap_invalid=1
amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1``) when the tutorials are
configured with ``-DAMReX_TESTING=ON`` or AMReX is built with
``-DAMReX_FPE=ON``.

.. code-block:: console

   ctest --test-dir build --output-on-failure
   ctest --test-dir build -R LinearSolvers

On a GPU, AMReX preallocates 3/4 of the device memory at startup. If the
test ranks share one GPU, they may run out of memory; set a smaller initial
size through the ``AMREX_DEFAULT_INIT`` environment variable:

.. code-block:: console

   AMREX_DEFAULT_INIT="amrex.the_arena_init_size=268435456" ctest --test-dir build

Building with GNU Make
----------------------

Clone AMReX next to this repository, or point ``AMREX_HOME`` to an AMReX
source tree, and run ``make`` in a tutorial directory:

.. code-block:: console

   cd ..
   git clone https://github.com/AMReX-Codes/amrex.git
   cd amrex-tutorials/Base/HeatEquation
   make -j8
   mpiexec -n 2 ./main3d.gnu.MPI.ex inputs

The ``GNUmakefile`` sets the build options, which can also be given on the
command line, e.g. ``make -j8 DIM=2 USE_OMP=TRUE``:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Variable
     - Meaning
   * - ``AMREX_HOME``
     - Path to AMReX (default: ``amrex`` next to this repository)
   * - ``DIM``
     - 1, 2 or 3
   * - ``COMP``
     - Compiler family: ``gnu``, ``llvm``, ``intel-llvm``, ...
   * - ``USE_MPI``, ``USE_OMP``
     - MPI (default on) and OpenMP
   * - ``USE_CUDA``, ``USE_HIP``, ``USE_SYCL``
     - GPU backend
   * - ``DEBUG``
     - Debug build with bounds checking
   * - ``TINY_PROFILE``
     - TinyProfiler (default off); see :ref:`profiling`

The name of the executable reflects the options, e.g.
``main3d.gnu.DEBUG.MPI.OMP.ex``. ``make clean`` removes the build.

Running the tutorials
---------------------

Runtime parameters come from the ``inputs`` file and the command line;
arguments after the inputs file override it, e.g.
``mpiexec -n 2 ./HeatEquation inputs n_cell=128``. See
:ref:`base-parmparse`.

The ``inputs`` files are set up so that each tutorial runs in a few seconds,
also as a test in CI: small domains, few steps, and small grids. In
production runs, ``max_grid_size`` and ``blocking_factor`` are usually much
larger, because larger boxes mean less overhead per box and better
performance, especially on GPUs.

Many tutorials write plotfiles; :ref:`viewing_plotfiles` lists tools that
read them.

Your own project
----------------

The :doc:`BuildSystems` directory contains three standalone projects to
copy when starting a new code: a CMake superbuild that fetches AMReX
(:ref:`buildsystems-cmakesuperbuild`), a CMake project that links to an
installed AMReX (:ref:`buildsystems-cmakeinstalledamrex`), and a GNU Make
project (:ref:`buildsystems-gnumake`).

For all AMReX build options, see `Building AMReX
<https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html>`__ in the
AMReX documentation.
