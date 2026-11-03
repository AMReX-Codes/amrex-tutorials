.. _timeintegration:

TimeIntegration
===============

These tutorials advance a ``MultiFab`` in time with
``amrex::TimeIntegrator``, which takes a right-hand side function and
provides native explicit Runge-Kutta methods and, when AMReX is built with
SUNDIALS, the explicit, implicit, IMEX and multirate methods of ARKODE. The
method is selected with ``integration.*`` runtime parameters. RungeKutta has
no extra requirements; SUNDIALS requires ``-DAMReX_SUNDIALS=ON`` (SUNDIALS is
downloaded and built unless ``SUNDIALS_ROOT`` is given). Build only this
topic with ``-DAMReXTutorials_TOPICS=TimeIntegration``; see
:doc:`GettingStarted`.

AMReX also has built-in Runge-Kutta functions, ``RungeKutta::RK2``,
``RK3`` and ``RK4`` in ``AMReX_RungeKutta.H``, which advance a
``MultiFab`` by one step given a right-hand side function and a function
that fills ghost cells. They are designed for AMR: together with
``FillPatcher``, they fill the ghost cells of a fine level at the
coarse/fine boundary at the intermediate times of the RK stages, using the
coarse level's stage data. For codes
built on ``AmrLevel``, ``AmrLevel::RK`` does all of this; see
:ref:`amr-waveamrlevel`.

.. note::

   ``TimeIntegrator`` advances a single ``MultiFab``. In an AMR code, the
   right-hand side function of a fine level must fill the ghost cells at the
   coarse/fine boundary at each stage time, and ``TimeIntegrator`` leaves
   that to the user. Interpolating the coarse data in time between the old
   and new states is in general not consistent with the RK stages, so for
   methods of order higher than two it can reduce the accuracy.

Both tutorials solve the same advection-diffusion problem on the periodic
unit domain,

.. math::

   \frac{\partial\phi}{\partial t} = -\mathbf{u}\cdot\nabla\phi + \kappa\nabla^2\phi,
   \qquad \phi(t=0) = \prod_d \sin(2\pi x_d),

with second-order centered differences. For this initial condition the
semi-discrete system has an exact solution (each sine is advected at a
slightly modified speed and decays at a modified rate), so the measured
error is the time integration error alone.

.. _timeintegration-rungekutta:

RungeKutta
----------

:Level: ★★
:Source: `TimeIntegration/RungeKutta <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/TimeIntegration/RungeKutta>`__
:Key APIs: ``TimeIntegrator<MultiFab>``, ``set_rhs``,
   ``set_post_step_action``, ``set_time_step``, ``evolve``,
   ``integration.type``, ``integration.rk.type``, ``ParmParse::queryAdd``

Advance the problem with the native explicit Runge-Kutta methods of
``TimeIntegrator``, chosen at run time. The tutorial runs the problem
``nruns`` times, halving the time step each time, and prints the error and
the observed order of accuracy, which it checks against the order of the
method.

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Parameter
     - Default
     - Meaning
   * - ``integration.type``
     - ``RungeKutta``
     - ``ForwardEuler``, ``RungeKutta``, or ``SUNDIALS`` (AMReX built with SUNDIALS)
   * - ``integration.rk.type``
     - 4
     - 1 forward Euler, 2 trapezoid, 3 SSPRK3, 4 RK4, 0 user Butcher tableau
       (``integration.rk.nodes``, ``weights``, ``tableau``)
   * - ``dt``, ``nruns``
     - 0.01, 3
     - step of the first run; each run halves it
   * - ``n_cell``, ``stop_time``
     - 16, 0.25
     - problem size and final time

The right-hand side function receives the integrator's stage data, so it
fills the ghost cells before applying the stencil in a ``ParallelFor``:

.. literalinclude:: ../../TimeIntegration/RungeKutta/main.cpp
   :language: cpp
   :start-at: // rhs = -vel.grad(phi) + kappa*lap(phi)
   :end-at: phi.FillBoundary(geom.periodicity());

The integrator is constructed from the state, given the right-hand side and
a time step, and advanced to the final time:

.. literalinclude:: ../../TimeIntegration/RungeKutta/main.cpp
   :language: cpp
   :start-at: // The integrator allocates its stage data like phi.
   :end-at: integrator.evolve(phi, stop_time);
   :dedent:

``queryAdd`` stores the defaults in ``ParmParse`` when they are not given in
the inputs, where ``TimeIntegrator`` reads them:

.. literalinclude:: ../../TimeIntegration/RungeKutta/main.cpp
   :language: cpp
   :start-at: std::string int_type = "RungeKutta";
   :end-at: ParmParse("integration.rk").queryAdd("type", rk_type);
   :dedent:

Run (from the top of the repository, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target TimeIntegration_RungeKutta -j 8
   cd build/TimeIntegration/RungeKutta
   mpiexec -n 2 ./RungeKutta inputs integration.rk.type=3

With the default RK4:

.. code-block:: text

           dt     steps  RHS evals      max error   rate
     1.000e-02      25        100  6.794e-07
     5.000e-03      50        200  4.212e-08   4.01
     2.500e-03     100        400  2.635e-09   4.00

   Expected order of accuracy: 4
   PASSED

The explicit methods need ``dt`` below about :math:`\Delta x^2/(2 d \kappa)`
in :math:`d` dimensions.

.. _timeintegration-sundials:

SUNDIALS
--------

:Level: ★★★
:Source: `TimeIntegration/SUNDIALS <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/TimeIntegration/SUNDIALS>`__
:Key APIs: ``TimeIntegrator<MultiFab>`` with ``IntegratorTypes::Sundials``,
   ``set_rhs``, ``set_imex_rhs``, ``set_fast_rhs``, ``set_time_step``,
   ``set_fast_time_step``, ``set_adaptive_step``, ``set_post_step_action``,
   ``set_post_fast_step_action``, ``integration.sundials.*``

Use the SUNDIALS (ARKODE) backend of ``TimeIntegrator`` for explicit,
implicit, implicit-explicit (IMEX) and multirate (MRI) time integration of
the problem above, here with a larger diffusion coefficient. Diffusion is
stiff: explicit methods need ``dt`` below about
:math:`\Delta x^2/(2 d \kappa)`. The tutorial runs each method type listed in
``methods`` in turn and compares the number of steps, right-hand side
evaluations and the error against the exact solution.

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - ``integration.sundials.type``
     - Right-hand side
     - Step
   * - ``ERK``
     - ``set_rhs(advection + diffusion)``, explicit
     - adaptive, limited by stability
   * - ``DIRK``
     - ``set_rhs(advection + diffusion)``, implicit (Newton-GMRES)
     - ``dt``
   * - ``IMEX-RK``
     - ``set_imex_rhs(diffusion, advection)``: implicit, explicit
     - ``dt``
   * - ``EX-MRI``
     - ``set_rhs(advection)`` slow, ``set_fast_rhs(diffusion)`` fast
     - ``dt``, fast ``fast_dt``

The integrator is created with
``TimeIntegrator<MultiFab> integrator(IntegratorTypes::Sundials, phi, time)``,
which selects the SUNDIALS backend in the code instead of with
``integration.type``. The method type is set with
``ParmParse("integration.sundials").add("type", method)``, and the
right-hand side functions and step sizes depend on it:

.. literalinclude:: ../../TimeIntegration/SUNDIALS/main.cpp
   :language: cpp
   :start-at: if (method == "ERK" || method == "DIRK") {
   :end-before: int steps = 0;
   :dedent:

Other SUNDIALS settings, such as ``integration.sundials.method`` or
``integration.sundials.fast_type = DIRK``, are read from the inputs file
(see the comments in ``inputs`` and the `Time Integration
<https://amrex-codes.github.io/amrex/docs_html/TimeIntegration_Chapter.html>`__
chapter of the AMReX documentation). ``integration.rel_tol`` and
``integration.abs_tol`` control the adaptive step and the implicit solves.
Increase ``n_cell`` to make diffusion stiffer and the implicit, IMEX and MRI
methods more attractive.

Run (from the top of the repository, with the build configured with
``-DAMReX_SUNDIALS=ON``, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target TimeIntegration_SUNDIALS -j 8
   cd build/TimeIntegration/SUNDIALS
   mpiexec -n 2 ./SUNDIALS inputs

After the SUNDIALS statistics of each run (silenced with
``amrex.verbose=0``):

.. code-block:: text

   method      steps  fast steps  RHS evals  relative error
   ERK            49           0        199       2.646e-08
   DIRK            5           0        236       6.983e-06
   IMEX-RK         5           0        160       1.054e-06
   EX-MRI          5          60        275       3.529e-04
   PASSED

With GNU Make, the ``GNUmakefile`` sets ``USE_SUNDIALS = TRUE``; pass the
SUNDIALS installation with ``SUNDIALS_DIR=/path/to/sundials/install``.
