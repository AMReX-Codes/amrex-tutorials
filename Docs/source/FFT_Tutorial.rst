.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _tutorials_fft:

FFT
==========================

These tutorials demonstrate how to use the amrex::FFT classes to solve for and manipulate Fourier transform data.
For more information on the amrex::FFT class refer to the AMReX documentation 
`here <https://amrex-codes.github.io/amrex/docs_html/FFT_Chapter.html>`_.

There are two FFT tutorials, ``Basic`` and ``Poisson``.

.. toctree::
   :maxdepth: 1

.. _section:fft_tutorial:fft_basic:

Basic
--------------------------

The tutorial found in ``amrex-tutorials/ExampleCodes/FFT/Basic`` demonstrates how
to take a forward FFT, manipulate or access the spectral data (in this example by copying the spectral
data into a :cpp:`MultiFab` for plotfile visualziation) and then taking the inverse FFT.

.. _section:fft_tutorial:fft_pois:

Poisson
--------------------------

This tutorial: ``amrex-tutorials/ExampleCodes/FFT/Poisson``
solves a Poisson equation with periodic boundary conditions.  It relies on the :cpp:`AMReX_FFT_Poisson.H`
routines to solve the equation by using the ``forwardThenBackward`` function which takes the forward and inverse
transform of the data with spectral data scaling in between as required by the Poisson equation.
