.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _tutorials_fftw:

FFTW
==========================

These tutorials demonstrate how to call fftw3 (CPU) or cuFFT (GPU) to solve for and manipulate Fourier transform data using a single MPI rank.

There are two FFTW tutorials, ``Basic`` and ``Poisson``:

- ``Basic`` tutorial: The tutorial found in ``amrex-tutorials/ExampleCodes/FFTW/Basic`` is
  useful if the objective is to simply take a forward FFT of data, and the DFT's ordering in k-space
  matters to the user.  This tutorial initializes a 3D or 2D :cpp:`MultiFab`, takes a forward FFT,
  and then redistributes the data in k-space where the center cell in the domain corresponds to the k=0 mode.
  The results are written to a plot file.

- ``Poisson`` tutorial: This tutorial: ``amrex-tutorials/ExampleCodes/FFTW/Poisson``
  solves a Poisson equation with periodic boundary conditions.  In it, both a forward FFT and reverse FFT
  are called to solve the equation, however, no reordering of the DFT data in k-space is performed.

We note that both fftw and cufft assume a row-major ordering of data; since a :cpp:`MultiFab` is column major,
the output to the spectral array is spatially-transposed.

.. toctree::
   :maxdepth: 1

.. _section:fftw_tutorial:fftw_basic:

Basic
--------------------------

This tutorial initializes a 3D or 2D :cpp:`MultiFab`, takes a forward FFT,
and then redistributes the data in k-space where the center cell in the domain corresponds to the k=0 mode.
The results are written to a plot file.

.. _section:fftw_tutorial:fftw_pois:

Poisson
--------------------------

In this test case we set up a right hand side (rhs), call the forward transform,
modify the coefficients, then call the backward solver and output the solution
to the discrete Poisson equation.
