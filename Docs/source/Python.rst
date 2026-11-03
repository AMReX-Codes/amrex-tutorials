.. _python:

Python
======

These tutorials use `pyAMReX <https://pyamrex.readthedocs.io>`__, the
Python bindings of AMReX, with zero-copy NumPy (CPU), CuPy (CUDA/HIP) or
dpnp (SYCL) views of AMReX data. They are Python scripts written for 3D
(``amrex.space3d``) and need pyAMReX 26.09 or newer (``import amrex``),
plus CuPy for CUDA/HIP or dpnp for SYCL GPU builds of pyAMReX; see the
pyAMReX documentation for installation. As with the C++ programs,
arguments after the script name are passed to AMReX, so
``python3 main.py inputs n_cell=64`` reads the inputs file and overrides
``n_cell``. The scripts are not compiled; with ``-DAMReXTutorials_PYTHON=ON``
(and ``-DPython3_EXECUTABLE`` pointing to a Python that has pyAMReX), CTest runs
them as the tests ``Python_<Name>``. See :doc:`GettingStarted`. The
:ref:`coupling-mpmd` tutorial also has a pyAMReX program.

.. _python-multifab:

MultiFab
--------

:Level: ★
:Source: `Python/MultiFab <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Python/MultiFab>`__
:Key APIs: ``amr.initialize``, ``amr.Box``, ``amr.BoxArray``,
           ``amr.DistributionMapping``, ``amr.Geometry``, ``amr.MultiFab``,
           ``amr.xp``, ``mf.array(mfi).to_xp()``, ``mf.to_xp()``,
           ``mf.fill_boundary``, ``mf.min``, ``mf.max``, ``mf.sum``,
           ``amr.write_single_level_plotfile``, ``amr.ParmParse``

The Python version of :ref:`base-multifab`: build a ``MultiFab``, fill it
through zero-copy views, fill ghost cells, compute reductions, and write a
plotfile.

``for mfi in mf:`` loops over the boxes owned by this MPI rank, and
``mf.array(mfi).to_xp(copy=False)`` returns a NumPy, CuPy or dpnp array
that shares memory with the box's data. The view includes the ghost cells, is
indexed ``[x, y, z, component]`` with indices starting at 0, so the valid
cells are ``[ng:ng+nx, ng:ng+ny, ng:ng+nz]`` for a box of ``nx*ny*nz``
cells. ``mf.to_xp()`` returns the views of all local boxes as a list. Pass
``order="C"`` for ``[component, z, y, x]`` indexing, which is usually faster
on GPUs.

.. literalinclude:: ../../Python/MultiFab/main.py
   :language: python
   :start-at: # Loop over the boxes owned by this rank.
   :end-at: phi[..., 1] = 2.0 * phi[..., 0]
   :dedent:

Command-line arguments, such as an inputs file, are passed to AMReX, and
all AMReX objects must be destroyed before ``amr.finalize``:

.. literalinclude:: ../../Python/MultiFab/main.py
   :language: python
   :start-at: if __name__ == "__main__":
   :end-at: sys.exit(status)

Run:

.. code-block:: console

   cd Python/MultiFab
   python3 main.py inputs
   mpiexec -n 2 python3 main.py inputs n_cell=64

.. code-block:: text

   min including ghost cells before FillBoundary: 0.0
   min including ghost cells after  FillBoundary: 1.0
   max(phi)      = 1.9293757195853387
   integral(phi) = 1.0055683279968117 (exact: 1.0055683279968317)
   max(2 phi)    = 3.8587514391706774
   Wrote plt00000
   PASSED

.. _python-heatequation:

HeatEquation
------------

:Level: ★
:Source: `Python/HeatEquation <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/Python/HeatEquation>`__
:Key APIs: ``amr.MultiFab``, ``amr.xp``, ``mf.array(mfi).to_xp()``,
           ``mf.fill_boundary``, ``geom.periodicity()``, ``mf.max``,
           ``amr.ParmParse.to_dict``,
           ``amr.write_single_level_plotfile``

The Python version of :ref:`base-heatequation`: the heat equation
:math:`\partial\phi/\partial t = D\nabla^2\phi` with explicit time steps on
a periodic unit cube. The initial condition is
:math:`1 + \sin(2\pi x)\sin(2\pi y)\sin(2\pi z)`; its sine mode decays
as :math:`\exp(-12\pi^2 D t)`, and the final error against this exact
solution is checked.

The stencil is written with array slicing on the zero-copy view of each
box. The loop over boxes is in Python, but the work on each box is done by
NumPy, CuPy or dpnp, so large boxes (``max_grid_size``) keep the Python
overhead small.

.. literalinclude:: ../../Python/HeatEquation/main.py
   :language: python
   :start-at: # Fill ghost cells, here all from periodic images and neighbors.
   :end-before: time += dt_step
   :dedent:

Run:

.. code-block:: console

   cd Python/HeatEquation
   python3 main.py inputs
   mpiexec -n 2 python3 main.py inputs n_cell=64 plot_int=20

.. code-block:: text

   Step 10, time = 0.0014648437499999998
   ...
   Finished 69 steps at time 0.01, mode amplitude = 0.3059442057
   Max error = 0.001952850866, tolerance = 0.008614567486
   PASSED
