.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _tutorials_eb:

EB
==========================

``amrex-tutorials/ExampleCodes/EB/CNS`` is an AMR code for solving compressible
Navier-Stokes equations with the embedded boundary approach.

``amrex-tutorials/ExampleCodes/EB/Poisson`` is a single-level code that is a proxy for
solving the electrostatic Poisson equation for a grounded sphere with a point
charge inside.

``amrex-tutorials/ExampleCodes/EB/MacProj`` is a single-level code that computes a divergence-free
flow field around a sphere.  A MAC projection is performed on an initial velocity
field of :math:`(1,0,0)`.
