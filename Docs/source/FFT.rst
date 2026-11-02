.. _fft:

FFT
===

Parallel fast Fourier transforms with ``amrex::FFT`` and the FFT-based
Poisson solvers built on them. The transforms work directly on
``MultiFab`` data distributed over MPI ranks. A CPU build uses FFTW; GPU
builds use cuFFT, rocFFT or oneMKL. These tutorials require AMReX built
with ``AMReX_FFT=ON`` (the default in the superbuild) and run in 2D and 3D.
Build only this topic with ``-DAMReXTutorials_TOPICS=FFT``; see
:doc:`GettingStarted`.

.. _fft-basic:

Basic
-----

:Level: ★
:Source: `FFT/Basic <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/FFT/Basic>`__
:Key APIs: ``FFT::R2C``, ``getSpectralDataLayout``, ``forward``, ``backward``,
           ``forwardThenBackward``, ``scalingFactor``, ``cMultiFab``, ``GpuComplex``

Real-to-complex transforms of a ``MultiFab`` on a periodic unit domain: a
forward and backward round trip, a spectral derivative :math:`d/dx`, and a
low-pass filter that zeroes all modes with :math:`|k| > k_\mathrm{cut}`.
The test field is a smooth low-wavenumber part plus a high-wavenumber part
that the filter removes, and all three results are checked against the exact
answer to round-off.

The spectral data has :math:`(n_x/2+1, n_y, n_z)` points: only
:math:`k_x \ge 0` is stored, because the transform of real data is
Hermitian. In :math:`y` and :math:`z`, indices above :math:`n/2` are
negative wavenumbers. FFTs are unnormalized, so a round trip multiplies the
data by the number of points; ``scalingFactor()`` returns its inverse.

An ``R2C`` object does the expensive setup once and is reused. The spectral
data can be stored in a ``cMultiFab`` with the layout given by
``getSpectralDataLayout``:

.. literalinclude:: ../../FFT/Basic/main.cpp
   :language: cpp
   :start-at: // R2C does the expensive setup
   :end-before: int const nx = n_cell;
   :dedent:

``forwardThenBackward`` calls a lambda at every spectral point between the
two transforms, without copying the spectral data into a ``cMultiFab``.
The spectral derivative multiplies by :math:`i k_x`:

.. literalinclude:: ../../FFT/Basic/main.cpp
   :language: cpp
   :start-at: // 2. Spectral derivative d/dx
   :end-before: MultiFab dphi_exact
   :dedent:

The low-pass filter maps the indices :math:`j` and :math:`k` to signed
wavenumbers:

.. literalinclude:: ../../FFT/Basic/main.cpp
   :language: cpp
   :start-at: // 3. Low-pass filter
   :end-before: MultiFab phi_low_exact
   :dedent:

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target FFT_Basic -j 8
   cd build/FFT/Basic
   mpiexec -n 2 ./Basic inputs

Runtime parameters: ``n_cell``, ``max_grid_size`` and ``kcut``.

.. _fft-poisson:

Poisson
-------

:Level: ★★
:Source: `FFT/Poisson <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/FFT/Poisson>`__
:Key APIs: ``FFT::Poisson``, ``FFT::PoissonHybrid``, ``FFT::Boundary``

Solve :math:`\nabla^2\phi = f` with FFT-based direct solvers and check the
result against a manufactured solution:

#. ``FFT::Poisson`` with periodic boundaries.
#. ``FFT::Poisson`` with homogeneous Dirichlet (``FFT::Boundary::odd``, sine
   transform) and Neumann (``FFT::Boundary::even``, cosine transform)
   boundaries, set per face.
#. 3D only: ``FFT::PoissonHybrid``, which uses FFTs in :math:`x` and
   :math:`y` and a tridiagonal solve in :math:`z`. The :math:`z` boundaries
   must be Dirichlet or Neumann, and :math:`z` may have a non-uniform cell
   size.

The solvers invert the standard second-order Laplacian exactly, so the
residual of the discrete equations is at round-off level, while the error
against the exact solution is :math:`O(\Delta x^2)`. If ``phi`` has ghost
cells, the solver fills one layer from the boundary conditions. With all
boundaries periodic or Neumann, the solution is unique only up to a constant
and the right-hand side must have zero mean.

The boundary types are given as a pair of ``FFT::Boundary`` values (low and
high face) per direction, with a matching ``Geometry`` periodicity:

.. literalinclude:: ../../FFT/Poisson/main.cpp
   :language: cpp
   :start-at: // 2. Homogeneous Dirichlet (odd)
   :end-before: #if (AMREX_SPACEDIM == 3)
   :dedent:

``FFT::PoissonHybrid`` is set up in the same way:

.. literalinclude:: ../../FFT/Poisson/main.cpp
   :language: cpp
   :start-at: // 3. FFT in x and y
   :end-before: #endif
   :dedent:

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target FFT_Poisson -j 8
   cd build/FFT/Poisson
   mpiexec -n 2 ./Poisson inputs

Doubling ``n_cell`` reduces the error by about a factor of four.

.. _fft-openbc:

OpenBC
------

:Level: ★★
:Source: `FFT/OpenBC <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/FFT/OpenBC>`__
:Key APIs: ``FFT::OpenBCSolver``, ``setGreensFunction``, ``FFT::PoissonOpenBC``

Solve :math:`\nabla^2\phi = \rho` in free space (open boundaries) for a
Gaussian charge and compare with the analytic potential. In 3D,

.. math::

   \phi = -\frac{\mathrm{erf}\left(r/(\sqrt{2}\sigma)\right)}{4\pi r},

and in 2D, :math:`\phi = \left(\ln r + E_1\left(r^2/(2\sigma^2)\right)/2\right)/(2\pi)`,
where :math:`E_1` is the exponential integral.

The free-space solution is the convolution of :math:`\rho` with the Green's
function :math:`G`, which is :math:`-1/(4\pi r)` in 3D and
:math:`\ln(r)/(2\pi)` in 2D. ``FFT::OpenBCSolver`` computes it with FFTs on
a domain doubled in each direction (Hockney's method), so that periodic
images do not interact. The user provides :math:`G` integrated over a cell
with ``setGreensFunction``; here with a :math:`2^D`-point Gauss quadrature,
which also avoids the singularity at :math:`r = 0`:

.. literalinclude:: ../../FFT/OpenBC/main.cpp
   :language: cpp
   :start-at: // The open-BC solution is the convolution
   :end-at: solver.solve(phi, rho);

In 3D, ``FFT::PoissonOpenBC`` provides the same Green's function:

.. literalinclude:: ../../FFT/OpenBC/main.cpp
   :language: cpp
   :start-at: FFT::PoissonOpenBC<MultiFab> poisson(geom);
   :end-at: poisson.solve(phi2, rho);
   :dedent:

The error is second order in :math:`\Delta x/\sigma` and is checked against
a tolerance.

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target FFT_OpenBC -j 8
   cd build/FFT/OpenBC
   mpiexec -n 2 ./OpenBC inputs

Runtime parameters: ``n_cell``, ``max_grid_size`` and the charge width
``sigma``. With ``n_cell = 64``, the relative error drops from about 0.007
to 0.0018.

.. _fft-energyspectrum:

EnergySpectrum
--------------

:Level: ★★
:Source: `FFT/EnergySpectrum <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/FFT/EnergySpectrum>`__
:Key APIs: ``FFT::R2C``, ``FFT::Info::setBatchSize``, ``cMultiFab``,
           ``HostDevice::Atomic::Add``, ``Gpu::DeviceVector``, ``Gpu::copy``,
           ``ParallelDescriptor::ReduceRealSum``

Compute the shell-averaged kinetic energy spectrum :math:`E(k)` of a periodic
velocity field, a standard diagnostic for turbulence. The velocity field is
synthetic: three Fourier modes with :math:`|k| = 2, 5, 8` (in units of
:math:`2\pi/L`) and known energies. The program checks that the energy is in
the expected shells and that Parseval's theorem holds:
:math:`\sum_k E(k) = \langle|\mathbf{u}|^2\rangle/2`.

All velocity components are transformed in one batched call:

.. literalinclude:: ../../FFT/EnergySpectrum/main.cpp
   :language: cpp
   :start-at: // Batched forward transform
   :end-at: r2c.forward(vel, vel_hat);
   :dedent:

Each spectral point adds :math:`|\hat{u}|^2/(2N^2)` to the shell
:math:`\mathrm{round}(|k|)` with an atomic add, giving a histogram on the
device. Because only :math:`k_x \ge 0` is stored, points with
:math:`0 < k_x < n_x/2` count twice. The local histograms are then summed
across MPI ranks:

.. literalinclude:: ../../FFT/EnergySpectrum/main.cpp
   :language: cpp
   :start-at: // Histogram on the device
   :end-at: ParallelDescriptor::ReduceRealSum
   :dedent:

Run (from a CMake build directory, see :doc:`GettingStarted`):

.. code-block:: console

   cmake --build build --target FFT_EnergySpectrum -j 8
   cd build/FFT/EnergySpectrum
   mpiexec -n 2 ./EnergySpectrum inputs

The spectrum is written to ``spectrum.txt`` (set with ``outfile``):

.. code-block:: text

   Nonzero shells:
     E(2) = 0.25  (exact 0.25)
     E(5) = 0.0625  (exact 0.0625)
     E(8) = 0.015625  (exact 0.015625)
   Kinetic energy: physical space 0.328125, spectral space 0.328125
