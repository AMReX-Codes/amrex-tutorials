.. _tutorials_sundials:

=============================
SUNDIALS and Time Integrators
=============================

This example code shows how to use the AMReX TimeIntegrator class
with SUNDIALS backend for integration. It also include two example
of AMReX native integrators, Forward Euler, and Explicit Runge Kutta.
Each integration type can be chosen by selecting the corresponding
inputs file:

  - ``inputs_forward_euler`` -- Native AMReX Forward Euler integrator

  - ``inputs_rk3`` -- Native AMReX Explicit Runge Kutta

  - ``inputs_sundials_erk`` -- SUNDIALS backend

Both Runge Kutta and SUNDIALS have additional options which can
be set by modifying the inputs file. Please see each respective inputs
file or `AMReX User Guide:Time Integration`_ for more details.

.. _`AMReX User Guide:Time Integration`: https://amrex-codes.github.io/amrex/docs_html/TimeIntegration_Chapter.html#
