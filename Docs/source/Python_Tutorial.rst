.. _tutorials_python:

======
Python
======

These examples show how to use AMReX from Python.
AMReX applications can also be interfaced to Python with the same logic.

In order to run the Python tutorials, you need to have pyAMReX installed.
Please see `pyAMReX <https://github.com/AMReX-Codes/pyamrex/>`__ for more details.

Alternatively, you can build the ExampleCodes in this repository with ``-DTUTORIAL_PYTHON=ON`` added to the CMake configuration options,
then install with ``cmake --build build --target pyamrex_pip_install``, and pyamrex will be installed for you.

Once pyAMReX is installed, you can run the following Guided Tutorial Examples:

- :download:`MultiFab <../../GuidedTutorials/MultiFab/main.py>`
- :download:`Heat Equation <../../GuidedTutorials/HeatEquation/Source/main.py>`

