.. _buildsystems:

BuildSystems
============

These are standalone projects that show how to build your own code against
AMReX. Each one is a small program meant to be copied as the start of a new
project. They are not part of the tutorials' CMake build; build them from
their own directories. See :doc:`GettingStarted` for building the tutorials
themselves.

.. _buildsystems-cmakesuperbuild:

CMakeSuperbuild
---------------

:Source: `BuildSystems/CMakeSuperbuild <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/BuildSystems/CMakeSuperbuild>`__
:Key APIs: ``FetchContent_Declare``, ``FetchContent_MakeAvailable``,
           ``AMReX::amrex``, ``setup_target_for_cuda_compilation``

AMReX is downloaded with ``FetchContent`` and built together with the
application. Most AMReX applications are built this way. AMReX options are
ordinary CMake cache variables; the project sets defaults before AMReX is
added, and they can be overridden on the command line.

.. literalinclude:: ../../BuildSystems/CMakeSuperbuild/CMakeLists.txt
   :language: cmake
   :start-at: cmake_minimum_required

Build and run:

.. code-block:: console

   cd BuildSystems/CMakeSuperbuild
   cmake -S . -B build
   cmake --build build -j 8
   mpiexec -n 2 ./build/MyApp

- AMReX options, e.g., ``-DAMReX_SPACEDIM=2 -DAMReX_OMP=ON -DAMReX_EB=ON``.
- GPU: ``-DAMReX_GPU_BACKEND=CUDA -DCMAKE_CUDA_ARCHITECTURES=native``, or
  ``HIP`` with ``-DAMReX_AMD_ARCH=gfx90a``, or ``SYCL`` with Intel's
  ``icpx`` compiler.
- A local AMReX source tree instead of a download:
  ``-DFETCHCONTENT_SOURCE_DIR_AMREX=/path/to/amrex``.
- A fixed AMReX version: set ``GIT_TAG`` to a release such as ``26.10``.

.. _buildsystems-cmakeinstalledamrex:

CMakeInstalledAMReX
-------------------

:Source: `BuildSystems/CMakeInstalledAMReX <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/BuildSystems/CMakeInstalledAMReX>`__
:Key APIs: ``find_package(AMReX)``, ``AMReX::amrex_3d``,
           ``setup_target_for_cuda_compilation``

The project links to an AMReX installed with CMake.
``find_package(AMReX REQUIRED 3D EB ...)`` fails if the installed AMReX
lacks a requested component. AMReX provides the targets
``AMReX::amrex_1d``, ``AMReX::amrex_2d`` and ``AMReX::amrex_3d`` for each
dimension it was built for. The project inherits AMReX's compiler, MPI,
OpenMP and GPU settings.

.. literalinclude:: ../../BuildSystems/CMakeInstalledAMReX/CMakeLists.txt
   :language: cmake
   :start-at: cmake_minimum_required

Install AMReX once, then build the project:

.. code-block:: console

   cmake -S /path/to/amrex -B amrex-build -DAMReX_SPACEDIM=3 -DCMAKE_INSTALL_PREFIX=$HOME/amrex-install
   cmake --build amrex-build -j 8 --target install

   cd BuildSystems/CMakeInstalledAMReX
   cmake -S . -B build -DAMReX_ROOT=$HOME/amrex-install
   cmake --build build -j 8
   mpiexec -n 2 ./build/MyApp

.. _buildsystems-gnumake:

GNUMake
-------

:Source: `BuildSystems/GNUMake <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/BuildSystems/GNUMake>`__
:Key APIs: ``AMREX_HOME``, ``Make.defs``, ``Make.package``, ``Make.rules``

A project built with AMReX's GNU Make system. The application sources are
in ``Source/`` and the build happens in ``Exec/``, a layout that scales to
codes with several problem setups (``Exec/ProblemA``, ...).
``Source/Make.package`` lists the application sources, and ``AMREX_DIRS``
lists the AMReX source directories to compile.

.. literalinclude:: ../../BuildSystems/GNUMake/Exec/GNUmakefile
   :language: make

.. literalinclude:: ../../BuildSystems/GNUMake/Source/Make.package
   :language: make

Build and run:

.. code-block:: console

   cd BuildSystems/GNUMake/Exec
   make -j8 AMREX_HOME=/path/to/amrex
   mpiexec -n 2 ./main3d.gnu.MPI.ex

Options are set in ``Exec/GNUmakefile`` or on the command line:

.. list-table::
   :header-rows: 1

   * - Variable
     - Values
   * - ``DIM``
     - 1, 2, 3
   * - ``COMP``
     - gnu, llvm, intel-llvm
   * - ``DEBUG``
     - TRUE, FALSE
   * - ``USE_MPI``, ``USE_OMP``
     - TRUE, FALSE
   * - ``USE_CUDA``, ``USE_HIP``, ``USE_SYCL``
     - TRUE, FALSE (e.g., ``USE_CUDA=TRUE CUDA_ARCH=80``)

The name of the executable records the build options, e.g.,
``main3d.gnu.MPI.ex`` for a 3D GCC build with MPI.
