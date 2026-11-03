AMReX Tutorials
===============

These tutorials show how to write codes with `AMReX
<https://github.com/AMReX-Codes/amrex>`__, a software framework for
block-structured adaptive mesh refinement. Each tutorial is a small,
self-contained program that runs on CPUs and GPUs (the Fortran ones on CPUs
only, PyTorchInference not with HIP or SYCL), with MPI and OpenMP, and
builds with CMake or GNU Make (the Python ones are scripts). Most of them
check their result against an exact or reference answer, and all of them
run as tests in CI (GPU builds are compiled only).

The source is in the `amrex-tutorials repository
<https://github.com/AMReX-Codes/amrex-tutorials>`__. For the concepts and
the full API, see the `AMReX documentation
<https://amrex-codes.github.io/amrex/docs_html/>`__.

Organization
------------

Tutorials are organized by topic. Most live in ``<Topic>/<Name>/``
with a ``README.md``, an ``inputs`` file, a ``CMakeLists.txt`` and a
``GNUmakefile``; the Applications, BuildSystems and Python ones differ.
Each has a level:

- ★ introductory
- ★★ intermediate
- ★★★ advanced
- **Blueprint**: a complete code meant to be copied as the starting point
  of an application

New users should start with :doc:`GettingStarted`, then
:ref:`base-helloworld`, :ref:`base-multifab` and :ref:`base-heatequation`
in :doc:`Base`.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Topic
     - Contents
   * - :doc:`Base`
     - Core data structures, kernels, reductions, I/O, ``ParmParse``
   * - :doc:`AMR`
     - Adaptive mesh refinement with ``AmrCore`` and ``AmrLevel``
   * - :doc:`LinearSolvers`
     - MLMG, GMRES and related solvers
   * - :doc:`EB`
     - Embedded boundaries
   * - :doc:`Particles`
     - Particle containers and particle-mesh operations
   * - :doc:`FFT`
     - ``amrex::FFT`` transforms and FFT-based solvers
   * - :doc:`TimeIntegration`
     - ``amrex::TimeIntegrator`` and SUNDIALS
   * - :doc:`Coupling`
     - Multiple programs in one MPI job (MPMD)
   * - :doc:`Fortran`
     - Applications written in Fortran with the AMReX Fortran interfaces
   * - :doc:`Python`
     - pyAMReX
   * - :doc:`ML`
     - Machine-learning inference
   * - :doc:`Applications`
     - Complete solvers combining several topics
   * - :doc:`BuildSystems`
     - Standalone projects that build against AMReX

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   GettingStarted

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   Base
   AMR
   LinearSolvers
   EB
   Particles
   FFT
   TimeIntegration
   Coupling
   Fortran
   Python
   ML
   Applications
   BuildSystems

.. toctree::
   :maxdepth: 1
   :caption: Tools

   VisualizationAndDebugging

Questions and help
------------------

Questions are welcome on `AMReX's GitHub Discussions
<https://github.com/AMReX-Codes/amrex/discussions>`__. To add a tutorial,
see `CONTRIBUTING.md
<https://github.com/AMReX-Codes/amrex-tutorials/blob/main/CONTRIBUTING.md>`__.

AMReX and these tutorials are distributed under the 3-clause BSD license;
see `LICENSE <https://github.com/AMReX-Codes/amrex-tutorials/blob/main/LICENSE>`__.
