.. _tutorials_sundials:

=============================
SUNDIALS and Time Integrators
=============================

These example codes demonstrate how to use the AMReX TimeIntegrator class
with SUNDIALS backend for integration.

The first example code at ``amrex-tutorials/ExampleCodes/SUNDIALS/Single-Rate``
solves the heat equation:

.. math:: \frac{\partial\phi}{\partial t} = \nabla^2\phi.

The inputs file contains a template for single process time integration strategies.

The second example code at ``amrex-tutorials/ExampleCodes/SUNDIALS/Reaction-Diffusion``
solves the reaction-diffusion equation, where :math:`R` and :math:`D` are
user-supplied reaction and diffusion coefficients:

.. math:: \frac{\partial\phi}{\partial t} = D \nabla^2\phi - R \phi.

The inputs file contains a template for MRI approaches, where the diffusion process
can be treated as a "fast" partition relative to the reaction process.

The third example code at ``amrex-tutorials/ExampleCodes/SUNDIALS/AdvDiff-MLMGPrecon``
solves the advection-diffusion equation, where :math:`A` and :math:`D` are
user-supplied advection and diffusion coefficients:

.. math:: \frac{\partial\phi}{\partial t} = D \nabla^2\phi - A \nabla\phi.

The inputs file contains a template for IMEX approaches, where diffusion is
treated implicitly with a SUNDIALS preconditioner that uses AMReX MLMG. The
example can also use MLMG's HYPRE bottom solver interface when AMReX is built
with HYPRE support.

Please see the inputs file or
`AMReX User Guide:Time Integration`_ for more details.

.. _`AMReX User Guide:Time Integration`: https://amrex-codes.github.io/amrex/docs_html/TimeIntegration_Chapter.html#
