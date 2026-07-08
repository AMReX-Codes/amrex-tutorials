.. _tutorials_python:

======
Python
======

These examples show how to use AMReX from Python, via `pyAMReX <https://github.com/AMReX-Codes/pyamrex/>`__.
AMReX applications can also be interfaced to Python with the same logic.

Installation
============

In order to run the Python tutorials, you need to have pyAMReX installed.
Please see the `pyAMReX documentation <https://pyamrex.readthedocs.io>`__ for installation details.

Alternatively, you can build the ExampleCodes in this repository with ``-DTUTORIAL_PYTHON=ON`` added to the CMake configuration options,
then install with ``cmake --build build --target pyamrex_pip_install``, and pyAMReX will be installed for you.

Running
=======

Python tutorials are written so they run the same way as their C++ counterparts:
command line arguments after the script name are forwarded to AMReX, so an
:ref:`inputs file <amrex_docs:sec:basics:parmparse>` and ``key=value`` overrides can be passed as usual:

.. code-block:: sh

   python3 main.py inputs

   # with runtime parameter override(s)
   python3 main.py inputs nsteps=20

   # MPI-parallel
   mpiexec -n 2 python3 main.py inputs

Tutorials
=========

Guided tutorials:

- :download:`MultiFab <../../GuidedTutorials/MultiFab/main.py>`: define, fill and plot a MultiFab
- :download:`Heat Equation <../../GuidedTutorials/HeatEquation/Source/main.py>`: explicit heat equation solve with ghost cell exchanges and runtime parameters
  (run with :download:`inputs <../../GuidedTutorials/HeatEquation/Exec/inputs>`)

Example codes:

- :download:`MPMD Case-2 <../../ExampleCodes/MPMD/Case-2/main.py>`: a Python app coupled to a C++ app via AMReX MPMD (see :ref:`tutorials_mpmd`)
