


AMReX Guided Tutorials and Example Codes
========================================

Welcome to AMReX Tutorials and Example Codes. Here you will find
a progression of hands on demonstrations and useful stand-alone examples
designed to assist users in learning and writing their own AMReX code.
For additional details on the topics presented here, please see the
`AMReX Source Documentation`_. The resources are divided into two categories:

.. _`AMReX Source Documentation`: https://amrex-codes.github.io/amrex/docs_html/


Guided Tutorials
----------------

Guided tutorials provide a gentle introduction to AMReX
features by focusing on key concepts in a progressive way. They are designed
to be followed from start to finish. They will focus on a few specific goals
of larger importance to using AMReX rather than comment on every line of code.

.. toctree::
   :maxdepth: 1

   Guided Tutorials  <GuidedTutorials>

Example Codes
-------------

Example codes are stand-alone examples that demonstrate how to use
different parts of the AMReX functionality. They are aimed at
users who are comfortable with the basics of AMReX.
They present a straightforward application of AMReX features, and
provide a starting place for user's to develop their unique applications.

The example codes listed below are provided in the ``AMReX-Codes/amrex-tutorials`` repo,
under the directory ``amrex-tutorials/ExampleCodes``. The examples are
sorted by the following categories:

- :ref:`AMR<tutorials_amr>` -- Examples of adaptive mesh refinement.
- :ref:`Basic<tutorials_basic>` -- Fundamental operations supported by AMReX.
- :ref:`Blueprint<tutorials_blueprint>` -- Convert AMReX mesh data into an in-memory
  Conduit Mesh Blueprint for use with the ALPINE Ascent in situ visualization
  and analysis tool.
- :ref:`EB<tutorials_eb>`  -- Examples of embedded boundaries.
- :ref:`ForkJoin<tutorials_forkjoin>` -- Parallel execution and subgrouping of MPI ranks.
- :ref:`GPU<tutorials_gpu>`  -- Offload work to the GPUs using AMReX tools.
- :ref:`Linear Solvers<tutorials_linearsolvers>`  -- Examples of several linear solvers.
- :ref:`MUI<tutorials_mui>`  -- Incorporates the MxUI/MUI (Multiscale Universal interface) frame into AMReX.
- :ref:`Particles<tutorials_particles>`  -- Basic usage of AMReX's particle data structures.
- :ref:`SDC<tutorials_sdc>`  -- Example usage of a "Multi-Implicit" Spectral Deferred Corrections (MISDC) integrator
  to solve a scalar advection-diffusion-reaction equation.
- :ref:`SENSEI<tutorials_sensei>`  -- In situ data analysis and visualization through a unified interface.
- :ref:`SUNDIALS<tutorials_sundials>`  -- Time integration with SUNDIALS backend and native AMReX types.
- :ref:`SWFFT<tutorials_swfft>`  -- Demonstrates how to call the SWFFT wrapper to the FFTW3 (A distributed memory
  implementation of the discrete Fourier transform) solver.


.. toctree::
   :hidden:

   AMR_Tutorial
   Basic_Tutorial
   Blueprint_Tutorial
   EB_Tutorial
   ForkJoin_Tutorial
   GPU_Tutorial
   LinearSolvers_Tutorial
   ML_Tutorial
   MUI_Tutorial
   Particles_Tutorial
   SDC_Tutorial
   SENSEI_Tutorial
   SUNDIALS_Tutorial
   SWFFT_Tutorial


.. _`AMR`:  AMR_Tutorial.html

.. _`Basic`:  Basic_Tutorial.html

.. _`Blueprint`:  Blueprint_Tutorial.html

.. _`EB`:  EB_Tutorial.html

.. _`ForkJoin`:  ForkJoin_Tutorial.html

.. _`GPU`:  GPU_Tutorial.html

.. _`Linear Solvers`:  LinearSolvers_Tutorial.html

.. _`MUI`: MUI_Tutorial.html

.. _`Particles`: Particles_Tutorial.html

.. _`SDC`: SDC_Tutorial.html

.. _`SENSEI`: SENSEI_Tutorial.html

.. _`SWFFT`: SWFFT_Tutorial.html



|
|

Additional Questions and Help
-----------------------------

Didn't find what you were looking for? Have questions we didn't answer?
Please let us know how we can improve by posting on `AMReX's GitHub Discussions`_.

.. _`AMReX's GitHub Discussions`: https://github.com/AMReX-Codes/amrex/discussions

The copyright notice of AMReX is included in the AMReX home directory as README.md.

Your use of this software is under the 3-clause BSD license -- the license agreement is included in the
AMReX home directory as license.txt.

