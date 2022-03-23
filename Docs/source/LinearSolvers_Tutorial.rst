.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _tutorials_linearsolvers:

LinearSolvers
==========================

There are five examples in the ``ExampleCodes/LinearSolvers`` directory.

``ABecLaplacian_C`` demonstrates how to solve with cell-centered data in a C++ framework.
This example shows how to use either hypre or PETSc as a bottom-solver (or to solve
the equation at the finest level if you set the "max coarsening level" to 0.

Step-by-step instructions for configuring AMReX to use HYPRE for this example are availble
:ref:`here <HYPRE>`.

``ABecLaplacian_F`` demonstrates how to solve with cell-centered data using the Fortran interfaces.

``NodalPoisson`` demonstrates how to set up and solve a variable coefficient Poisson equation
with the rhs and solution data on nodes.

``MultiComponent`` demonstrates how to solve using a multi-component operator with nodal data.
The operator is of the form :math:`D(\mathbf{\phi})_i = \alpha_{ij}\nabla^2\phi_j`, where :math:`\mathbf{\phi}` is the multi-component solution variable and :math:`\alpha_{ij}` is a matrix of coefficients.
The operator (``MCNodalLinOp``) is implemented using the "reflux-free" method.

``NodeTensorLap`` solves the equation :math:`\nabla \cdot ( \sigma \nabla \phi ) = \text{RHS}`, where
:math:`\phi` and :math:`\text{RHS}` are nodal MultiFabs and :math:`\sigma` is a tensor constant. It is
assumed that the tensor is symmetric with rank equal to ``AMREX_SPACEDIM``.

