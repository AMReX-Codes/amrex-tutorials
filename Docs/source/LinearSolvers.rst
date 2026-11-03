.. _linearsolvers:

LinearSolvers
=============

AMReX's geometric multigrid solver ``MLMG`` with cell-centered, nodal and
edge-centered operators, GMRES, the algebraic multigrid solver AlgMG, and a
user-defined operator. These
tutorials require AMReX built with ``AMReX_LINEAR_SOLVERS=ON`` (the
default). They run in 1D, 2D and 3D, except NodalProjection (2D and 3D).
Build only this topic with ``-DAMReXTutorials_TOPICS=LinearSolvers``; see
:doc:`GettingStarted`.

All tutorials here read the ``MLMG`` options ``mlmg.verbose`` and
``mlmg.tol_rel`` (and some ``mlmg.max_iter``) from the inputs.
``MLMG::solve(sol, rhs, tol_rel, tol_abs)`` stops when the max norm of the
residual is at most ``max(tol_abs, tol_rel*norm0)``, where ``norm0`` is by
default the larger of the max norms of the right-hand side and the initial
residual. Set ``tol_abs`` to zero unless you know the scale of the residual.

.. _linearsolvers-cellpoisson:

CellPoisson
-----------

:Level: ★
:Source: `LinearSolvers/CellPoisson <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/LinearSolvers/CellPoisson>`__
:Key APIs: ``MLPoisson``, ``LinOpBCType``, ``setDomainBC``, ``setLevelBC``,
           ``MLMG``, ``MLMG::solve``, ``setVerbose``, ``setMaxIter``

Solve :math:`\nabla^2\phi = f` for a cell-centered :math:`\phi` with
``MLPoisson``. This is the workflow common to all AMReX linear solvers:

#. Build the operator from ``Geometry``, ``BoxArray`` and
   ``DistributionMapping`` (one of each per AMR level).
#. Set the domain boundary types with ``setDomainBC``, then the boundary
   values with ``setLevelBC``. The values are taken from the ghost cells of
   a ``MultiFab`` (here the solution itself) and represent values on the
   domain faces.
#. Build an ``MLMG`` solver, set its options, and call ``solve``.

The boundary type in each direction comes from the inputs, e.g.,
``bc = dirichlet neumann periodic``. The exact solution
:math:`\phi = \prod_d \cos(2\pi x_d)` works with all three. Without a
Dirichlet direction, :math:`\phi` is determined only up to a constant, so
the program removes the mean before comparing. The problem is solved on two
grids, and the error drops by a factor of four when the cell size is halved.

.. literalinclude:: ../../LinearSolvers/CellPoisson/main.cpp
   :language: cpp
   :start-at: MLPoisson mlpoisson({geom}, {ba}, {dm});
   :end-at: mlmg.solve({&phi}, {&rhs}, p.tol_rel, p.tol_abs);
   :dedent:

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target LinearSolvers_CellPoisson -j 8
   cd build/LinearSolvers/CellPoisson
   mpiexec -n 2 ./CellPoisson inputs

With ``mlmg.verbose = 1``, MLMG prints the norms and the number of
iterations:

.. code-block:: text

   MLMG: Initial rhs               = 116.7325837
   MLMG: Initial residual (resid0) = 2821.164753
   MLMG: Final Iter. 7 resid, resid/resid0 = 4.106625795e-08, 1.455649051e-11
   ...
   n_cell = 32: max error = 0.003249084214
   n_cell = 64: max error = 0.0008196247914
   Convergence rate = 1.98699761

.. _linearsolvers-abeclaplacian:

ABecLaplacian
-------------

:Level: ★★
:Source: `LinearSolvers/ABecLaplacian <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/LinearSolvers/ABecLaplacian>`__
:Key APIs: ``MLABecLaplacian``, ``setScalars``, ``setACoeffs``, ``setBCoeffs``,
           ``LinOpBCType::inhomogNeumann``, ``setMaxOrder``, ``MLMG::solve``
           with ``GetVecOfPtrs``/``GetVecOfConstPtrs``

Solve the variable-coefficient equation

.. math::

   a\,\alpha(\mathbf{x})\,\phi - b\,\nabla\cdot\left(\beta(\mathbf{x})\nabla\phi\right) = f

with ``MLABecLaplacian`` as one composite solve over all levels of a fixed
AMR hierarchy. Each finer level refines the middle half of the level below
(``max_level = 1`` by default):

.. literalinclude:: ../../LinearSolvers/ABecLaplacian/main.cpp
   :language: cpp
   :start-at: // Fixed AMR hierarchy
   :end-before: // Cell-centered phi (with one ghost cell
   :dedent:

:math:`\alpha` is cell-centered and :math:`\beta` is face-centered.
:math:`\beta` changes steeply by a factor ``beta_in`` (10 by default) across
a sphere covered by the fine level. The boundary type per direction is
``dirichlet``, ``neumann`` (``LinOpBCType::inhomogNeumann``) or ``periodic``.
The ghost cells passed to ``setLevelBC`` hold :math:`\phi` on the face for
Dirichlet, and :math:`\partial\phi/\partial x_d` on the face for Neumann.
``setMaxOrder`` sets the order of the stencil at Dirichlet and coarse/fine
boundaries.

.. literalinclude:: ../../LinearSolvers/ABecLaplacian/main.cpp
   :language: cpp
   :start-at: // The operator for all levels of the composite solve.
   :end-at: mlmg.solve(GetVecOfPtrs(phi)
   :dedent:

After the solve, ``MLMG`` has averaged the fine solution onto the covered
coarse cells, and the error is checked on every level against
:math:`\phi = \prod_d \sin(2\pi x_d + \pi/4)`.

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target LinearSolvers_ABecLaplacian -j 8
   cd build/LinearSolvers/ABecLaplacian
   mpiexec -n 2 ./ABecLaplacian inputs

Try ``max_level=2``, ``bc="neumann neumann neumann"``, or a larger
``n_cell``.

.. _linearsolvers-implicitheat:

ImplicitHeat
------------

:Level: ★★
:Source: `LinearSolvers/ImplicitHeat <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/LinearSolvers/ImplicitHeat>`__
:Key APIs: ``MLABecLaplacian``, ``setScalars``, ``setACoeffs``/``setBCoeffs``
           with scalars, ``setLevelBC(0, nullptr)``, reusing ``MLMG`` across
           solves, ``MLMG::getNumIters``

Advance the heat equation :math:`\partial\phi/\partial t = \kappa\nabla^2\phi`
on a periodic domain with backward Euler,

.. math::

   (1 - \Delta t\,\kappa\nabla^2)\,\phi^{n+1} = \phi^n,

which is ``MLABecLaplacian`` with :math:`a = 1`, :math:`\alpha = 1`,
:math:`b = \Delta t` and :math:`\beta = \kappa`. The time step is
``dt_factor`` (10 by default) times the explicit stability limit
:math:`\Delta x^2/(2D\kappa)`; backward Euler is stable for any
:math:`\Delta t`. Compare with the explicit :ref:`base-heatequation`.

The operator and the solver are built once and reused every step. Each
solve starts from the previous solution:

.. literalinclude:: ../../LinearSolvers/ImplicitHeat/main.cpp
   :language: cpp
   :start-at: // Build the operator and the solver once
   :end-before: // The mode is an eigenfunction
   :dedent:

The initial condition :math:`\phi = \prod_d \sin(2\pi x_d)` is an
eigenfunction of the discrete Laplacian, so each step multiplies it by a
known factor, and the program checks the solution against this to the
solver tolerance. It also prints the error against the PDE solution, which
is dominated by the :math:`O(\Delta t)` error of backward Euler.

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target LinearSolvers_ImplicitHeat -j 8
   cd build/LinearSolvers/ImplicitHeat
   mpiexec -n 2 ./ImplicitHeat inputs

Try ``dt_factor=1 nsteps=200`` to reduce the time error.

.. _linearsolvers-nodalprojection:

NodalProjection
---------------

:Level: ★★
:Source: `LinearSolvers/NodalProjection <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/LinearSolvers/NodalProjection>`__
:Key APIs: ``MLNodeLaplacian``, ``compDivergence``, ``updateVelocity``,
           ``setSigma``, nodal ``BoxArray`` (``convert(ba, IntVect::TheNodeVector())``)

Project a cell-centered velocity :math:`\mathbf{u}` onto a divergence-free
field, the core step of incompressible flow solvers:

.. math::

   \nabla\cdot(\sigma\nabla\phi) = \nabla\cdot\mathbf{u}, \qquad
   \mathbf{u}^\mathrm{new} = \mathbf{u} - \sigma\nabla\phi,

with :math:`\phi` on the nodes, solved by ``MLNodeLaplacian``. Here
:math:`\sigma = 1` is passed to the constructor; for variable density, use
:math:`\sigma = 1/\rho` with ``setSigma``. The domain is periodic in
:math:`x` with walls (``LinOpBCType::Neumann``, i.e.,
:math:`\mathbf{u}\cdot\mathbf{n} = 0`) elsewhere. The system is singular,
and MLMG handles it.

``compDivergence`` computes the nodal :math:`\nabla\cdot\mathbf{u}` (the
velocity needs one ghost cell), and ``updateVelocity`` subtracts
:math:`\sigma\nabla\phi`:

.. literalinclude:: ../../LinearSolvers/NodalProjection/main.cpp
   :language: cpp
   :start-at: // Nodal Poisson operator with constant sigma
   :end-at: linop.updateVelocity({&vel}, {&phi});
   :dedent:

This is an approximate projection: the nodal Laplacian is not exactly the
product of the discrete divergence and gradient, so
:math:`\nabla\cdot\mathbf{u}^\mathrm{new}` is :math:`O(\Delta x^2)`, not zero.
The program checks that the remaining divergence and the error in the
projected velocity are :math:`O(\Delta x^2)`.

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target LinearSolvers_NodalProjection -j 8
   cd build/LinearSolvers/NodalProjection
   mpiexec -n 2 ./NodalProjection inputs

.. _linearsolvers-gmres:

GMRES
-----

:Level: ★★
:Source: `LinearSolvers/GMRES <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/LinearSolvers/GMRES>`__
:Key APIs: ``amrex::GMRES``, ``GMRES::define``, ``setRestartLength``,
           ``setMaxIters``, ``getNumIters``, ``getStatus``, ``MLMG::precond``,
           ``MLMG::setBottomSolver``, ``GMRESMLMG``

Use AMReX's GMRES with multigrid as the preconditioner, in two ways.

**Your own operator.** ``amrex::GMRES<V,M>`` works with any vector type
``V`` and operator class ``M``. ``ConvDiffOp`` implements the
convection-diffusion operator :math:`\nabla^2 u - \mathbf{c}\cdot\nabla u`
without storing a matrix. It provides the matrix-vector product (``apply``),
vector operations (``dotProduct``, ``linComb``, ...), and a preconditioner
(``precond``) that runs one ``MLMG`` V-cycle of ``MLPoisson``, i.e., of the
Laplacian alone. The operator is not symmetric, so CG cannot be used. With
the preconditioner GMRES needs about 40 iterations, independent of the grid
size; without it (``use_precond = 0``) about 300.

.. literalinclude:: ../../LinearSolvers/GMRES/ConvDiffOp.cpp
   :language: cpp
   :start-at: void ConvDiffOp::precond

The preconditioner must be a fixed linear operator. That is why it does a
fixed number of V-cycles and uses the smoother instead of BiCGStab as the
bottom solver:

.. literalinclude:: ../../LinearSolvers/GMRES/main.cpp
   :language: cpp
   :start-at: MLMG mlmg(mlpoisson);
   :end-before: MultiFab::Subtract(u, exact, 0, 0, 1, 0);
   :dedent:

**An MLMG operator.** ``GMRESMLMG`` takes an ``MLMG`` object and solves its
system with GMRES, using MLMG itself as the preconditioner. This can help
when MLMG alone converges slowly.

.. literalinclude:: ../../LinearSolvers/GMRES/main.cpp
   :language: cpp
   :start-at: GMRESMLMG gmsolver(mlmg);
   :end-at: gmsolver.solve(u, lap_rhs, tol_rel, Real(0));
   :dedent:

Both problems use homogeneous Dirichlet boundaries and the exact solution
:math:`u = \prod_d \sin(\pi x_d)`.

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target LinearSolvers_GMRES -j 8
   cd build/LinearSolvers/GMRES
   mpiexec -n 2 ./GMRES inputs

Runtime parameters: ``n_cell``, ``velocity`` (the vector :math:`\mathbf{c}`),
``use_precond``, ``verbose`` and ``tol_rel``.

.. _linearsolvers-algmg:

AlgMG
-----

:Level: ★★★
:Source: `LinearSolvers/AlgMG <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/LinearSolvers/AlgMG>`__
:Key APIs: ``AlgPartition``, ``AlgVector``, ``SpMatrix::define``, ``SpMV``,
           ``AlgMG``, ``AlgMG::setKrylovSolver``, ``MLMG::setMultigridType``,
           ``BottomSolver::algmg``, ``MLMG::setAlgMGOptions``,
           ``MLMG::usedAlgMG``

AlgMG is AMReX's algebraic multigrid solver. It solves a sparse linear
system :math:`A x = b` and builds its coarse levels from the matrix itself,
so it works where geometric multigrid struggles, e.g. for strongly varying
or anisotropic coefficients. The problem here is
:math:`-\nabla\cdot(\beta\nabla\phi) = f` on the unit domain with
:math:`\phi = 0` on the boundary, where :math:`\beta` is a checkerboard of
blocks of ``block`` cells alternating between 1 and ``jump`` (default
:math:`10^4`). The right-hand side is the discrete operator applied to
:math:`\phi = \prod_d \sin(\pi x_d)`.

**Standalone.** ``AlgPartition`` splits the rows, one per cell, evenly
among the MPI ranks, and ``AlgVector`` is a distributed vector. The matrix
is assembled from CSR arrays of the local rows. Each row has room for the
:math:`2d+1` entries of the stencil; a missing neighbor at the domain
boundary gets column -1, which ``CsrValid{false}`` tells
``SpMatrix::define`` to drop. AlgMG then solves with V-cycles accelerated
by conjugate gradients:

.. literalinclude:: ../../LinearSolvers/AlgMG/main.cpp
   :language: cpp
   :start-at: SpMatrix<Real> A;
   :end-at: amg.solve(x, b);
   :dedent:

**In MLMG.** The same problem with ``MLABecLaplacian`` is solved with four
``MLMG`` variants:

.. literalinclude:: ../../LinearSolvers/AlgMG/main.cpp
   :language: cpp
   :start-at: Variant const variants[] = {
   :end-at: };
   :dedent:

Geometric multigrid diverges on this problem, because its coarse grids
cannot represent the coefficient. AlgMG as the bottom solver does not help,
since the geometric levels above the bottom are the problem. The algebraic
type solves the whole level with AlgMG in every MLMG iteration. The hybrid
type starts with geometric multigrid and, when it stalls or diverges,
switches to AlgMG and restarts from the best iterate; ``usedAlgMG()``
reports the switch. ``setAlgMGOptions`` takes a callback that receives the
``AlgMG`` object MLMG builds, to change its options:

.. literalinclude:: ../../LinearSolvers/AlgMG/main.cpp
   :language: cpp
   :start-at: mlmg.setAlgMGOptions(
   :end-at: });
   :dedent:

Use the geometric type when it converges, since it is cheaper, and the
hybrid type when the coefficients may be difficult. With
``setThrowException(true)``, a failed solve throws an exception instead of
aborting.

Run:

.. code-block:: console

   cmake --build build --target LinearSolvers_AlgMG -j 8
   cd build/LinearSolvers/AlgMG
   mpiexec -n 2 ./AlgMG inputs

.. code-block:: text

   AlgMG standalone: 7 levels, 9 iterations, relative residual 8.6e-11, max error 7.3e-10
   MLMG geometric: did not converge after 88 iterations
   MLMG geometric + AlgMG bottom: did not converge after 88 iterations
   MLMG algebraic: 3 iterations, max error 6.6e-11
   MLMG hybrid: 9 iterations, max error 9.8e-10, switched to AlgMG
   PASSED

Runtime parameters: ``n_cell``, ``max_grid_size``, ``block``, ``jump``,
``max_iter``, ``algmg.krylov``, ``algmg.smoother``, ``algmg.verbose``,
``mlmg.verbose`` and ``mlmg.tol_rel`` (the tolerance of both solvers). Try ``jump=1``, where geometric multigrid converges.

.. _linearsolvers-curlcurl:

CurlCurl
--------

:Level: ★★★
:Source: `LinearSolvers/CurlCurl <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/LinearSolvers/CurlCurl>`__
:Key APIs: ``MLCurlCurl``, ``setScalars``, ``prepareRHS``,
           ``MLMGT<Array<MultiFab,3>>``, edge ``BoxArray``\ s with ``amrex::convert``

Solve the curl-curl equation for an edge-centered vector field,

.. math::

   \nabla\times(\alpha\nabla\times\mathbf{E}) + \beta\mathbf{E} = \mathbf{f},

with ``MLCurlCurl``. Equations like this arise in implicit time stepping of
Maxwell's equations and in frequency-domain electromagnetics.

:math:`\mathbf{E}` is an ``Array<MultiFab,3>``, always with three
components. :math:`E_d` lives on edges along direction :math:`d`:
cell-centered in :math:`d` and nodal in the other directions. In 2D,
:math:`E_z` is nodal in both :math:`x` and :math:`y`. The initial guess also
supplies the tangential components on Dirichlet boundaries (here zero, a
perfect conductor):

.. literalinclude:: ../../LinearSolvers/CurlCurl/main.cpp
   :language: cpp
   :start-at: // Edge-centered fields: an Array of three MultiFabs.
   :end-at: E[dir].setVal(0);
   :dedent:

``setScalars(alpha, beta)`` sets constant coefficients; ``setAlpha`` (nodal)
and ``setBeta`` (edge-centered) set variable ones. ``prepareRHS`` makes the
right-hand side consistent on edges shared by several boxes. The solver is
``MLMGT<Array<MultiFab,3>>``, MLMG templated on the field type:

.. literalinclude:: ../../LinearSolvers/CurlCurl/main.cpp
   :language: cpp
   :start-at: MLCurlCurl mlcc({geom}, {ba}, {dm});
   :end-at: mlmg.solve({&E}, {&rhs}, tol_rel, Real(0));
   :dedent:

The curl of a gradient is zero, so with a small :math:`\beta` the operator is
nearly singular; MLCurlCurl's smoother is designed for this. The exact
solution :math:`E_d = \prod_{e \ne d} \sin(\pi x_e)` vanishes tangentially
on the walls, and the error is checked to be :math:`O(\Delta x^2)`.

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target LinearSolvers_CurlCurl -j 8
   cd build/LinearSolvers/CurlCurl
   mpiexec -n 2 ./CurlCurl inputs

Try ``beta=1.e-3``.

.. _linearsolvers-customlinop:

CustomLinOp
-----------

:Level: ★★★
:Source: `LinearSolvers/CustomLinOp <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/LinearSolvers/CustomLinOp>`__
:Key APIs: ``MLNodeLinOp``, ``Fapply``, ``Fsmooth``, ``restriction``,
           ``interpolation``, ``getNComp``, ``buildMasks``, ``nodalSync``,
           ``setDirichletNodesToZero``

Write your own linear operator and solve with ``MLMG``.
``MultiCompNodeLap`` is a nodal operator with :math:`N` coupled components,

.. math::

   L(\phi)_n = \sum_m C_{nm} \nabla^2\phi_m,

where :math:`C` is a constant, invertible :math:`N \times N` matrix from the
inputs (:math:`N \le 4`). It derives from ``MLNodeLinOp``, which provides
the multigrid hierarchy, ghost node exchange, the Dirichlet node mask
(``m_dirichlet_mask``) and dot products. The derived class implements:

.. literalinclude:: ../../LinearSolvers/CustomLinOp/MultiCompNodeLap.H
   :language: cpp
   :start-at: [[nodiscard]] int getNComp
   :end-before: private:
   :dedent:

``restriction`` (full weighting) and ``interpolation`` (linear) must handle
different box layouts on the two multigrid levels, which occur after
agglomeration. The smoother is block Jacobi, where the diagonal block
:math:`D = C\,\mathrm{diag}(\nabla^2)` couples the components at a node:

.. literalinclude:: ../../LinearSolvers/CustomLinOp/MultiCompNodeLap.cpp
   :language: cpp
   :start-at: // Jacobi with the diagonal blocks
   :end-at: nodalSync(amrlev, mglev, sol);
   :dedent:

The operator is single-level only. A multi-level operator would also
implement ``averageDownSolutionRHS`` and ``reflux`` for the coarse/fine
boundaries (see `Multi-Component Operators
<https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html#multi-component-operators>`__
in the AMReX documentation). The program checks the error against the exact solution
:math:`\phi_n = \prod_d \sin((n+1)\pi x_d)`, which is zero on the Dirichlet
boundary.

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target LinearSolvers_CustomLinOp -j 8
   cd build/LinearSolvers/CustomLinOp
   mpiexec -n 2 ./CustomLinOp inputs

Try three components: ``coeff="1 0 0  0 2 0.5  0 0.3 1"``.
