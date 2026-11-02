.. _eb:

EB
==

These tutorials cover embedded boundaries (EB): building cut-cell
geometries from implicit functions and STL files, solving elliptic problems
with boundary conditions on the EB, and conservative cut-cell advection.
All of them need AMReX built with EB support (``Poisson`` and
``Advection`` also need the linear solvers) and run in 2D and 3D only. The
CMake superbuild enables both by default. Build only
this topic with ``-DAMReXTutorials_TOPICS=EB``; see :doc:`GettingStarted`.

.. _eb-geometry:

Geometry
--------

:Level: ★
:Source: `EB/Geometry <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/EB/Geometry>`__
:Key APIs: ``EB2::SphereIF``, ``EB2::BoxIF``, ``EB2::CylinderIF``,
   ``EB2::PlaneIF``, ``EB2::ParserIF``, ``EB2::makeUnion``,
   ``EB2::makeIntersection``, ``EB2::makeDifference``, ``EB2::rotate``,
   ``EB2::translate``, ``EB2::makeShop``, ``EB2::Build``,
   ``makeEBFabFactory``, ``EBFArrayBoxFactory::getEBDataArrays``,
   ``EBData_t``, ``EB_WriteSingleLevelPlotfile``, ``WriteEBSurface``

Build EB geometries on :math:`[-1,1]^D` and query the cut-cell data. The
program builds each geometry listed in ``geometries``:

- ``sphere``: a sphere (a circle in 2D) with the fluid inside.
- ``csg``: constructive solid geometry, a box with a cylindrical hole,
  rotated and translated, plus a hemisphere made by intersecting a ball
  with a plane.
- ``parser``: an implicit function given as a runtime expression,
  ``parser.function``.
- ``stl``: a triangulated surface from ``octahedron.stl``, read through the
  ``eb2.*`` runtime parameters.

An implicit function is negative in the fluid and positive in the body, and
the CSG operations act on the bodies:

.. literalinclude:: ../../EB/Geometry/main.cpp
   :language: cpp
   :start-at: // Body 1: a box with a cylindrical hole
   :end-at: max_coarsening_level);
   :dedent:

A ``Parser`` expression in :math:`x, y, z` becomes an implicit function
with ``EB2::ParserIF``. Without a ``GeometryShop``, ``EB2::Build`` reads the
geometry from the ``eb2.*`` parameters (e.g., ``eb2.geom_type = stl``):

.. literalinclude:: ../../EB/Geometry/main.cpp
   :language: cpp
   :start-at: Parser parser(function);
   :end-at: max_coarsening_level);
   :dedent:

For each geometry, the program counts regular, cut and covered cells, and
sums the fluid volume, the EB surface area, and the fluid area of the plane
:math:`x = 0`. The EB data (cell flags, volume and area fractions,
centroids, ...) are accessed in a kernel through the factory's
``EBDataArrays``:

.. literalinclude:: ../../EB/Geometry/main.cpp
   :language: cpp
   :start-at: // The factory gives access to the EB data
   :end-at: });
   :dedent:

The results are compared with the exact values; the program fails if a
relative error exceeds 2% at ``n_cell = 32`` (the tolerance scales with
:math:`\Delta x^2`). It writes a plotfile ``plt_<geometry>`` with the cell
type and the volume fraction, and in 3D the EB surface of
``surface_geometry`` to ``eb.pvtp`` for ParaView.

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target EB_Geometry -j 8
   cd build/EB/Geometry
   mpiexec -n 2 ./Geometry inputs

Try ``geometries=parser parser.function="x^2+y^2+z^2-0.25"``, or a finer
grid, ``n_cell=128``.

.. _eb-poisson:

Poisson
-------

:Level: ★★
:Source: `EB/Poisson <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/EB/Poisson>`__
:Key APIs: ``EB2::SphereIF``, ``EB2::Build``, ``EBFArrayBoxFactory``,
   ``EBData_t::centroid``, ``EBData_t::bndrycent``, ``MLEBABecLap``,
   ``MLEBABecLap::setEBDirichlet``, ``MLMG``, ``EB_WriteSingleLevelPlotfile``

Solve :math:`-\nabla^2\phi = f` on :math:`[0,1]^D` outside a spherical
obstacle with the cell-centered EB solver ``MLEBABecLap``, with Dirichlet
conditions on the domain boundary and inhomogeneous Dirichlet conditions on
the EB. The exact solution is :math:`\phi = \prod_d \sin(\pi x_d)`. The
program prints the max-norm error at cell centers and fails if it exceeds
:math:`10\,\Delta x^2`; the scheme is second order. The plotfile ``plt``
has ``phi``, ``exact``, ``error`` and ``vfrac``.

Multigrid needs the EB on coarsened grids, so ``EB2::Build`` is called with
a large ``max_coarsening_level``. MultiFabs built with the factory carry the
EB data:

.. literalinclude:: ../../EB/Poisson/main.cpp
   :language: cpp
   :start-at: // A spherical obstacle with the fluid outside.
   :end-at: grids, dmap, {2,2,2}, EBSupport::full);
   :dedent:

In cut cells, the right-hand side is evaluated at the cell centroid and the
EB value at the EB centroid. Both are stored relative to the cell center in
units of :math:`\Delta x`:

.. literalinclude:: ../../EB/Poisson/main.cpp
   :language: cpp
   :start-at: GpuArray<Real,AMREX_SPACEDIM> xc, xb;
   :end-at: phieb_a[b](i,j,k) = exact_solution(AMREX_D_DECL(xb[0],xb[1],xb[2]));
   :dedent:

``MLEBABecLap`` solves :math:`(a\alpha - b\nabla\cdot\beta\nabla)\phi =
\mathrm{rhs}`, here with :math:`a = 0` and :math:`b = \beta = 1`. Without
``setEBDirichlet`` or ``setEBHomogDirichlet``, the EB condition is
homogeneous Neumann (see :ref:`linearsolvers-abeclaplacian` for the
non-EB operator):

.. literalinclude:: ../../EB/Poisson/main.cpp
   :language: cpp
   :start-at: // Operator: (a alpha - b div beta grad) phi = rhs
   :end-at: mlmg.solve({&phi}, {&rhs}, tol_rel, tol_abs);
   :dedent:

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target EB_Poisson -j 8
   cd build/EB/Poisson
   mpiexec -n 2 ./Poisson inputs

The error drops by about 4 when ``n_cell`` doubles:

.. code-block:: text

   MLMG: Final Iter. 12 resid, resid/resid0 = 1.220027457e-08, 4.927121892e-12
   ...
   Max error: 0.005236049738 (tolerance 0.009765625)
   PASSED

.. _eb-advection:

Advection
---------

:Level: ★★★
:Source: `EB/Advection <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/EB/Advection>`__
:Key APIs: ``EBFArrayBoxFactory``, ``EBDataArrays``, ``EB_computeDivergence``,
   ``MLEBABecLap``, ``MLMG::getFluxes``, ``ApplyRedistribution``,
   ``EB_set_covered``, ``EB_average_face_to_cellcenter``,
   ``EB_WriteSingleLevelPlotfile``

Conservative finite-volume advection of a scalar around a spherical
obstacle on a periodic domain, with redistribution for small cut cells.

1. A uniform velocity is made discretely divergence free around the
   obstacle with a MAC projection solved by ``MLEBABecLap``. Afterwards,
   the sum of area fraction times velocity over the faces of every cell is
   zero, and nothing flows through the EB.
2. Fluxes use second-order upwinding with limited slopes, dropping to first
   order next to the EB. A cut cell with volume fraction :math:`\kappa` and
   face area fractions :math:`\alpha` is updated with

   .. math::

      \frac{dS}{dt} = -\frac{1}{\kappa}\sum_d
      \frac{(\alpha u S)_{d,\mathrm{hi}} - (\alpha u S)_{d,\mathrm{lo}}}{\Delta x_d}.

3. A small cut cell would limit the time step. ``ApplyRedistribution``
   moves part of its update to its neighbors without changing the total, so
   the time step is based on the full cell size. Time integration is
   two-stage SSP Runge-Kutta.

The projection solves :math:`\nabla^2\phi = \nabla\cdot\mathbf{u}` with the
default homogeneous Neumann condition on the EB, and ``getFluxes`` returns
:math:`-\nabla\phi` on face centroids:

.. literalinclude:: ../../EB/Advection/main.cpp
   :language: cpp
   :start-at: MLEBABecLap op({geom}, {ba}, {dm}, LPInfo(), {&factory});
   :end-before: EB_computeDivergence(divu, GetArrOfConstPtrs(u), geom, true);
   :dedent:

The conservative divergence divides by the volume fraction:

.. literalinclude:: ../../EB/Advection/main.cpp
   :language: cpp
   :start-at: // Conservative divergence, (sum of a u S)
   :end-at: m_dSdt.FillBoundary(m_geom.periodicity());
   :dedent:

Redistribution is applied box by box where the box has cut cells:

.. literalinclude:: ../../EB/Advection/main.cpp
   :language: cpp
   :start-at: const FabType type = flags[mfi].getType(amrex::grow(bx,nghost));
   :end-before: } else {
   :dedent:

The state has two components, a Gaussian blob and a constant. The program
checks that there are no NaNs, that the total mass
(:math:`\sum \kappa S\,\Delta V`) is conserved to round-off, and that the
constant stays constant, which holds only if the velocity is divergence
free.

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target EB_Advection -j 8
   cd build/EB/Advection
   mpiexec -n 2 ./Advection inputs

Choose ``redistribution_type`` = ``StateRedist`` (default), ``FluxRedist``
or ``NoRedist``; in 3D, ``NoRedist`` is unstable with the default time
step. To watch the blob go around the obstacle, run with ``n_cell=64
stop_time=1 plot_int=10``.

.. code-block:: text

   max |div u| before projection: 881.403378
   max |div u| after projection:  1.707217148e-11
   Redistribution: StateRedist, dt = 0.006040948875, initial mass = 0.007724640317
   Steps: 17, time = 0.1
   ...
   Relative change in mass: 3.817692201e-15 (tolerance 2.220446049e-13)
   Max deviation of the constant: 1.914024494e-13 (tolerance 1e-08)
   PASSED
