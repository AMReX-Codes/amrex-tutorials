.. _visualization_and_debugging:

Visualization, Debugging and Profiling
======================================

.. _viewing_plotfiles:

Viewing plotfiles
-----------------

Tutorials write AMReX plotfiles with ``WriteSingleLevelPlotfile`` or
``WriteMultiLevelPlotfile``. A plotfile is a directory, e.g. ``plt00010``,
with a ``Header`` (variable names, time, geometry and the ``BoxArray`` of
each level) and one subdirectory per level. Several tools read it directly.
For details, see `Visualization
<https://amrex-codes.github.io/amrex/docs_html/Visualization.html>`__ in the
AMReX documentation.

yt
^^

`yt <https://yt-project.org>`__ reads plotfiles, including multi-level and
particle data, from Python:

.. code-block:: python

   import yt

   ds = yt.load("plt00010")
   ds.print_stats()
   slc = yt.SlicePlot(ds, "z", ("boxlib", "phi"))
   slc.annotate_grids()
   slc.save("phi.png")

   # Uniform-resolution array of level 0
   cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge,
                         dims=ds.domain_dimensions)
   phi = cg["boxlib", "phi"].to_ndarray()

pyAMReX
^^^^^^^

`pyAMReX <https://pyamrex.readthedocs.io>`__ reads plotfiles with
``PlotFileData``, which returns each level as a ``MultiFab``. Each box can
be viewed as a NumPy (CPU), CuPy (CUDA/HIP) or dpnp (SYCL) array without
copying:

.. code-block:: python

   import amrex.space3d as amr

   amr.initialize([])
   pf = amr.PlotFileData("plt00010")
   print(pf.varNames(), pf.time(), pf.finestLevel())
   mf = pf.get(0, "phi")          # level 0, one variable
   for mfi in mf:
       a = mf.array(mfi).to_xp()  # indexed as a[i, j, k, comp]
       print(mfi.validbox(), a.max())
   del mf, pf
   amr.finalize()

In C++, ``amrex::PlotFileData`` does the same; see
:ref:`base-checkpointrestart`.

ParaView and VisIt
^^^^^^^^^^^^^^^^^^

`ParaView <https://www.paraview.org>`__ opens a plotfile directory with
the "AMReX/BoxLib Grid Reader"; to open a time series, select the plotfile
directories as a group. Particles in the plotfile are read with the
"AMReX/BoxLib Particles Reader". `VisIt <https://visit-dav.github.io>`__
opens the ``Header`` file inside a plotfile directory. To open a series in
VisIt, list the headers in a file with the ``.visit`` extension:

.. code-block:: console

   ls -1v plt*/Header | tee movie.visit

AMReXplorer
^^^^^^^^^^^

`AMReXplorer <https://github.com/AMReX-Codes/amrexplorer>`__ is a desktop
viewer for 2D and 3D plotfiles, as well as FAB and MultiFab files. It shows
the data values of individual cells, which is convenient when debugging,
and plays plotfile sequences as animations, e.g. ``amrexplorer plt00010``
or ``amrexplorer plt*``.

Debugging
---------

AMReX has several aids for debugging that are controlled at compile time
and with runtime parameters of the ``amrex.`` prefix. See `Debugging
<https://amrex-codes.github.io/amrex/docs_html/Debugging.html>`__ and
`Runtime Parameters
<https://amrex-codes.github.io/amrex/docs_html/RuntimeParameters.html>`__ in
the AMReX documentation for more.

Debug builds
^^^^^^^^^^^^

Build with ``-DCMAKE_BUILD_TYPE=Debug`` (CMake) or ``DEBUG=TRUE`` (GNU Make).
This defines ``AMREX_DEBUG``, which

- turns on ``AMREX_ASSERT`` and bounds checking in ``Array4``: an access
  outside the box of an ``Array4``, e.g. reading a ghost cell that the
  ``MultiFab`` does not have, aborts with the offending index and the box;
- initializes floating-point ``MultiFab``, ``FArrayBox`` and
  ``Gpu::DeviceVector`` data to signaling NaNs (``amrex.init_snan``), so
  that the use of uninitialized data can be trapped;
- turns off optimization and adds debug information.

``Array4`` bounds checking alone, in an optimized build, is enabled with
``-DAMReX_BOUND_CHECK=ON`` or ``BOUND_CHECK=TRUE``, and assertions with
``-DAMReX_ASSERTIONS=ON`` or ``USE_ASSERTION=TRUE``.

Floating-point exceptions
^^^^^^^^^^^^^^^^^^^^^^^^^

These runtime parameters make the program stop at the first
floating-point exception, with a backtrace, instead of propagating NaNs and
infinities:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Parameter
     - Traps
   * - ``amrex.fpe_trap_invalid=1``
     - Invalid operations that produce NaN, e.g. ``0/0``, ``sqrt(-1)``,
       or arithmetic on a signaling NaN
   * - ``amrex.fpe_trap_zero=1``
     - Division by zero
   * - ``amrex.fpe_trap_overflow=1``
     - Overflow
   * - ``amrex.init_snan=1``
     - Not a trap itself: initializes new floating-point data to
       signaling NaNs (default on in debug builds), so that
       ``fpe_trap_invalid`` catches reads of uninitialized data

For example:

.. code-block:: console

   mpiexec -n 2 ./HeatEquation inputs amrex.fpe_trap_invalid=1 amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1

Trapping works on CPUs only; GPU kernels do not trap. Vectorized code
raises floating-point exceptions in inactive SIMD lanes, where the source
code does not (see :ref:`base-simd`), so the traps are a debugging aid for
scalar runs and the CTest tests do not enable them. With clang,
``-DAMReX_TESTING=ON`` or ``-DAMReX_FPE=ON`` (CMake superbuild) or
``TEST=TRUE`` (GNU Make) compiles AMReX and the tutorials so that the
compiler avoids such exceptions where it can, at a cost in performance;
GCC gets no such flag.

Signals and backtraces
^^^^^^^^^^^^^^^^^^^^^^

By default (``amrex.signal_handling=1``), AMReX handles segmentation
faults, floating-point exceptions, ``amrex::Abort``, failed assertions and
``Ctrl-C``. It prints an error message and writes a backtrace of the call
stack to a file ``Backtrace.<rank>`` (``Backtrace.<rank>.<thread>`` with
OpenMP) for each rank that failed. Line numbers in backtraces need debug
information (``-g``); ``amrex/Tools/Backtrace/parse_bt.py`` converts the
addresses in a backtrace file if they were not resolved.

To catch the failure in a debugger such as ``gdb`` instead, turn off
AMReX's signal handling:

.. code-block:: console

   gdb --args ./HeatEquation inputs amrex.signal_handling=0 amrex.throw_exception=1

With ``amrex.throw_exception=1`` (and signal handling off), ``amrex::Abort``
throws a ``std::runtime_error``, so a debugger can stop where it was
called.

Another useful parameter is ``amrex.verbose=1``, which reports, e.g., the
number of MPI ranks and GPU devices.

GPU kernels run asynchronously, so an error may be reported long after the
kernel that caused it. Adding ``amrex::Gpu::streamSynchronize()`` after
suspect kernels, or running with ``CUDA_LAUNCH_BLOCKING=1`` on NVIDIA GPUs,
helps locate it. Both slow the code down, so use them for debugging only
and remove the extra synchronizations afterwards. ``compute-sanitizer``
(CUDA) finds out-of-bounds and race errors in device code.

.. _profiling:

Profiling with TinyProfiler
^^^^^^^^^^^^^^^^^^^^^^^^^^^

TinyProfiler is off by default in the tutorials. Build with
``-DAMReX_TINY_PROFILE=ON`` (CMake) or ``TINY_PROFILE=TRUE`` (GNU Make) to
turn it on. AMReX functions are already instrumented, and code regions are
timed with ``BL_PROFILE("name")``, which measures until the end of the
enclosing scope. Without TinyProfiler, ``BL_PROFILE`` does nothing. At
``amrex::Finalize``, TinyProfiler prints two tables, sorted by exclusive
and by inclusive time, with the number of calls and the minimum, average
and maximum time over MPI ranks:

.. code-block:: text

   TinyProfiler total time across processes [min...avg...max]: 1.234 ... 1.235 ... 1.236

   --------------------------------------------------------------------------------
   Name                          NCalls  Excl. Min  Excl. Avg  Excl. Max   Max %
   --------------------------------------------------------------------------------
   compile_time_options               5     0.4101     0.4102     0.4103  33.20%
   ...

A large difference between the minimum and maximum time indicates load
imbalance. On GPUs, kernels are asynchronous, so a region's time may be
attributed to a later region that synchronizes;
``tiny_profiler.device_synchronize_around_region=1`` adds synchronization
for accurate per-region timing. The extra synchronization slows the code
down, so it is not recommended for production runs.
``tiny_profiler.output_file=<name>`` writes the tables to a file, and
``tiny_profiler.enabled=0`` turns the profiler off. See
`AMReX Profiling Tools <https://amrex-codes.github.io/amrex/docs_html/AMReX_Profiling_Tools.html>`__
for more.
