.. _fortran:

Fortran
=======

AMReX has Fortran interfaces to its core data structures, so that a whole
application can be written in Fortran. The Fortran tutorials run on CPUs,
with MPI and OpenMP, and need a Fortran compiler. They are built only with
``-DAMReX_FORTRAN_INTERFACES=ON`` and are skipped in GPU builds. With GNU
Make, ``USE_F_INTERFACES = TRUE`` and the ``Src/F_Interfaces/Base`` package
are set in the ``GNUmakefile``. Build only this topic with
``-DAMReXTutorials_TOPICS=Fortran``; see :doc:`GettingStarted`.

.. _fortran-heatequation:

HeatEquation
------------

:Level: ★
:Source: `Fortran/HeatEquation <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Fortran/HeatEquation>`__
:Key APIs: ``amrex_base_module``, ``amrex_init``, ``amrex_parmparse``,
           ``amrex_boxarray``, ``amrex_distromap``, ``amrex_geometry``,
           ``amrex_multifab``, ``amrex_mfiter``, ``dataptr``,
           ``fill_boundary``, ``amrex_parallel_reduce_max``,
           ``amrex_write_plotfile``

The problem of :ref:`base-heatequation` written in Fortran: forward Euler
for :math:`\partial\phi/\partial t = D\nabla^2\phi` on the periodic unit
domain, checked against the exact solution.

The main program sets up the grids and the data like its C++ counterpart.
The physical domain and the periodicity are given to
``amrex_geometry`` through the ``geometry.*`` runtime parameters:

.. literalinclude:: ../../Fortran/HeatEquation/main.f90
   :language: fortran
   :start-at: ! amrex_geometry takes the physical domain
   :end-at: call amrex_multifab_build(phi_new, ba, dm, 1, 1)
   :dedent:

Kernels loop over tiles with an ``amrex_mfiter`` in an OpenMP parallel
region. ``dataptr`` returns a pointer ``p(i,j,k,n)`` to the data of a box,
indexed by AMReX's cell indices, which is passed to a plain Fortran loop.
Indices beyond ``amrex_spacedim`` are 1, so the same code works in 1D, 2D
and 3D:

.. literalinclude:: ../../Fortran/HeatEquation/heat.f90
   :language: fortran
   :start-at: ! Fill ghost cells, here all
   :end-at: !$omp end parallel
   :dedent:

Fortran objects such as ``amrex_multifab`` own C++ objects. They are built
with ``amrex_*_build`` and freed with ``amrex_*_destroy``, and those in the
main program must be destroyed before ``amrex_finalize``.

Configure from the repository root, build and run:

.. code-block:: console

   cmake -S . -B build -DAMReXTutorials_TOPICS=Fortran -DAMReX_FORTRAN_INTERFACES=ON
   cmake --build build --target Fortran_HeatEquation -j 8
   cd build/Fortran/HeatEquation
   mpiexec -n 2 ./HeatEquation inputs

.. code-block:: text

   Finished 69 steps at time 1.0000000000E-02, mode amplitude = 3.0594420565E-01
   Max error = 1.9528508660E-03, tolerance = 8.6145674863E-03
   PASSED
