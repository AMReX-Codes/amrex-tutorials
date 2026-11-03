.. _amr:

AMR
===

These tutorials cover block-structured adaptive mesh refinement: filling
ghost cells across levels, synchronizing levels with refluxing and
averaging down, and tagging, regridding and subcycling in time, with
``AmrCore`` and with ``Amr``/``AmrLevel``. StaticHierarchy and
AdvectionAmrCore run in 2D and 3D, and WaveAmrLevel in 1D, 2D and 3D.
WaveAmrLevel needs AMReX built with ``AMReX_AMRLEVEL=ON`` (the default). Build only this topic with ``-DAMReXTutorials_TOPICS=AMR``; see
:doc:`GettingStarted`. For complete codes, see :ref:`applications-cns` on
``Amr``/``AmrLevel`` and :ref:`applications-incflow` on ``AmrCore``.

.. _amr-statichierarchy:

StaticHierarchy
---------------

:Level: ★★
:Source: `AMR/StaticHierarchy <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/AMR/StaticHierarchy>`__
:Key APIs: ``FillPatchSingleLevel``, ``FillPatchTwoLevels``, ``PhysBCFunct``,
   ``GpuBndryFuncFab``, ``BCRec``, ``cell_cons_interp``, ``FluxRegister``,
   ``average_down``, ``WriteMultiLevelPlotfile``

The core multilevel operations of AMR on a fixed two-level hierarchy, with
no regridding and no ``AmrCore``. Level 0 covers the unit square/cube, and
level 1 refines its middle half by a factor of 2. A Gaussian diffuses,
:math:`\partial\phi/\partial t = \nabla\cdot(D\nabla\phi)`, across the
coarse/fine interface. Both levels take the same time step (no subcycling),
limited by the fine cell size.

The walls are zero-gradient (``BCType::foextrap``), so the total mass must
stay constant. The program checks this and fails if the relative change
exceeds round-off (:math:`10^{-12}` in double precision). Run with
``do_reflux=0`` to see conservation fail.

The physical boundary conditions are described by a ``BCRec`` per component
and applied by ``PhysBCFunct`` (see also :ref:`base-boundaryconditions`).
The ``FluxRegister`` is defined on the fine level's grids:

.. literalinclude:: ../../AMR/StaticHierarchy/main.cpp
   :language: cpp
   :start-at: // Homogeneous Neumann walls
   :end-before: // No subcycling
   :dedent:

Each time step fills ghost cells, updates both levels, corrects the coarse
cells next to the fine level with the fine fluxes (``Reflux``), and replaces
the covered coarse cells with the average of the fine cells
(``average_down``). ``FillPatchTwoLevels`` copies from level 1 where
possible and otherwise interpolates conservatively from level 0:

.. literalinclude:: ../../AMR/StaticHierarchy/main.cpp
   :language: cpp
   :start-at: // Level 0 ghost cells
   :end-at: amrex::average_down(phi[1], phi[0], 0, 1, ratio);
   :dedent:

Run (from the top of the repository, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target AMR_StaticHierarchy -j 8
   cd build/AMR/StaticHierarchy
   mpiexec -n 2 ./StaticHierarchy inputs

Parameters: ``n_cell`` (coarse cells per direction), ``max_grid_size``,
``nsteps``, ``cfl``, ``diff_coef``, ``do_reflux``, ``plot`` (write
``plt_static`` at the end).

.. code-block:: text

   Initial mass: 0.005550189972
   Final mass:   0.005550189972 after 100 steps, time = 0.003662109375
   Relative change in mass: 1.719036495e-15
   PASSED

.. _amr-advectionamrcore:

AdvectionAmrCore
----------------

:Level: ★★★ / Blueprint
:Source: `AMR/AdvectionAmrCore <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/AMR/AdvectionAmrCore>`__
:Key APIs: ``AmrCore`` (``InitFromScratch``, ``regrid``,
   ``MakeNewLevelFromScratch``, ``MakeNewLevelFromCoarse``, ``RemakeLevel``,
   ``ClearLevel``, ``ErrorEst``), ``AMRErrorTag``, ``FillPatchSingleLevel``,
   ``FillPatchTwoLevels``, ``InterpFromCoarseLevel``, ``FluxRegister``,
   ``average_down``, ``ParReduce``, ``WriteMultiLevelPlotfile``, ``VisMF``,
   ``ParmParse::queryAdd``

A complete adaptive code built on ``AmrCore``: a scalar :math:`\phi` is
advected, :math:`\partial\phi/\partial t + \nabla\cdot(\mathbf{u}\phi) = 0`,
by a swirling velocity derived from the stream function
:math:`\psi = \sin^2(\pi x)\sin^2(\pi y)\cos(\pi t/2)/\pi`. The grids
follow :math:`\phi` as it moves. The flow reverses at :math:`t = 1`, so
:math:`\phi` returns to its initial position at :math:`t = 2`. The code
covers tagging, regridding, time stepping with and without subcycling,
refluxing, plotfiles, and checkpoint/restart. It checks that the total mass
stays constant to round-off, which refluxing and conservative interpolation
guarantee. Use it as a starting point for your own ``AmrCore`` application.

Code structure
^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - File
     - Contents
   * - ``main.cpp``
     - Defaults for parameters ``AmrCore`` requires (``queryAdd``), the run,
       the mass check
   * - ``AmrCoreAdv.H``, ``AmrCoreAdv.cpp``
     - Level data, ``AmrCore`` callbacks, ``FillPatch``, time stepping
   * - ``Advance.cpp``
     - Single-level advance, fluxes, time step estimate
   * - ``IO.cpp``
     - Plotfiles and checkpoint/restart
   * - ``Kernels.H``
     - Initial data, velocity field, limited upwind flux

To adapt it: change ``Kernels.H`` for the physics; add components to
``phi_new``/``phi_old`` (and the ``BCRec``\ s) for more variables; use
``BCType::ext_dir`` with a functor in place of ``FabFillNoOp`` for inflow
boundaries; add ``AMRErrorTag``\ s, or a custom loop in ``ErrorEst``, for
other refinement criteria.

AmrCore callbacks
^^^^^^^^^^^^^^^^^

``AmrCore`` owns the grids and geometry of all levels and calls back into
the class to create level data. ``MakeNewLevelFromScratch`` uses the
initial data, ``MakeNewLevelFromCoarse`` interpolates from the level below,
``RemakeLevel`` moves existing data onto new grids, and ``ErrorEst`` tags
cells for refinement, here with one ``AMRErrorTag`` that tags
:math:`\phi > \mathtt{adv.phierr[lev]}`:

.. literalinclude:: ../../AMR/AdvectionAmrCore/AmrCoreAdv.cpp
   :language: cpp
   :start-after: phi[b](i,j,k) = adv::initial_phi(x, y);
   :end-before: AmrCoreAdv::GetData
   :lines: 4-39

Time stepping
^^^^^^^^^^^^^

With subcycling, ``TimeStepWithSubcycling`` advances level ``lev`` and then
calls itself ``ref_ratio`` times for ``lev+1``. Every ``amr.regrid_int``
steps it regrids the finer levels first. After the fine steps, ``Reflux``
and ``average_down`` synchronize the levels. ``TimeStepNoSubcycling``
(``adv.do_subcycle=0``) regrids, then advances every level with a ``dt``
that is stable on all of them.

.. literalinclude:: ../../AMR/AdvectionAmrCore/AmrCoreAdv.cpp
   :language: cpp
   :start-at: // Advance level lev by one step, then its finer levels
   :end-before: // Regrid, then advance all levels with the same dt.

``AdvanceLevel`` uses Heun's method (two-stage SSP Runge-Kutta) with upwind
fluxes from MC-limited linear reconstructions. Before each stage,
``FillPatch`` fills the ghost cells; fine-level ghost cells are interpolated
in space and time from the level below. Each level passes
:math:`\Delta t \cdot \text{area} \cdot F` to the flux registers it
shares with its coarser and finer neighbors:

.. literalinclude:: ../../AMR/AdvectionAmrCore/Advance.cpp
   :language: cpp
   :start-at: // Stage 0: phi_new = phi_old
   :end-before: // Sborder and flux are freed
   :dedent:

Running
^^^^^^^

Run (from the top of the repository, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target AMR_AdvectionAmrCore -j 8
   cd build/AMR/AdvectionAmrCore
   mpiexec -n 2 ./AdvectionAmrCore inputs

The ``inputs`` file stops after 10 coarse steps; use ``max_step=1000`` for
the full problem and ``amr.restart=chk00005`` to restart from a checkpoint.
Other useful options are ``adv.do_subcycle=0``, ``adv.do_reflux=0`` (mass
is then not conserved), ``amr.max_level`` and ``amr.n_cell``.

.. code-block:: text

   Coarse step 1 time = 0.01140770953 dt = 0.01140770953 levels = 3 mass = 0.007850917675
   ...
   Coarse step 10 time = 0.1146053625 dt = 0.01155880142 levels = 3 mass = 0.007850917675
   Relative change in mass: 8.83832208e-16
   PASSED

.. _amr-waveamrlevel:

WaveAmrLevel
------------

:Level: ★★★
:Source: `AMR/WaveAmrLevel <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/AMR/WaveAmrLevel>`__
:Key APIs: ``Amr``, ``AmrLevel``, ``LevelBld``, ``StateDescriptor``,
           ``AmrLevel::RK``, ``GpuBndryFuncFab``, ``cell_quartic_interp``,
           ``FourthOrderInterpFromFineToCoarse``

An AMR code built on ``Amr`` and ``AmrLevel`` that solves the scalar wave
equation :math:`\partial^2 u/\partial t^2 = \nabla^2 u` as the first-order
system :math:`\partial u/\partial t = v`, :math:`\partial v/\partial t =
\nabla^2 u`, with a fourth-order Laplacian and Runge-Kutta time stepping.
The domain is :math:`[0,1]` in x with :math:`u = 0` at :math:`x = 0` and
:math:`1`, and periodic in y and z. The velocity :math:`v` starts as a pulse
in the middle, which splits in two, reflects at the walls and comes back,
and the fine level follows the pulses. The solution does not depend on y
and z, so d'Alembert's formula gives the exact :math:`v`, which the program
compares with the result.

``Amr`` drives the run and creates one ``AmrLevelWave`` object per level
through a ``LevelBld``. ``variableSetUp`` defines the state data: two
components, two ghost cells for the stencil, a boundary function for the x
boundaries, and the fourth-order ``cell_quartic_interp`` interpolater for
the coarse/fine boundary. This interpolater requires ``amr.ref_ratio = 2``.

The time step uses AMReX's built-in Runge-Kutta support. ``AmrLevel::RK``
advances the state with RK3 or RK4 (``wave.rk_order``) and fills the
ghost cells of every stage, including those at the coarse/fine boundary at
the stage times, from the coarse level's RK stage data. The code only
provides the right-hand side:

.. literalinclude:: ../../AMR/WaveAmrLevel/AmrLevelWave.cpp
   :language: cpp
   :start-at: AmrLevel::RK(rk_order, State_Type
   :end-at: });
   :dedent:

After each step, ``post_timestep`` replaces the coarse data under the fine
level with a fourth-order average of the fine data. The code does not
reflux, so the integral of :math:`v` is not conserved exactly across
levels.

Run:

.. code-block:: console

   cmake --build build --target AMR_WaveAmrLevel -j 8
   cd build/AMR/WaveAmrLevel
   mpiexec -n 2 ./WaveAmrLevel inputs

.. code-block:: text

   Finished 320 level-0 steps at time 1 with 2 levels
   Max error of v: 1.063265662e-05
   PASSED

Useful options are ``wave.rk_order``, ``wave.cfl``, ``amr.n_cell``,
``amr.max_level``, ``amr.v=1`` and ``amr.plot_files_output=1``.
