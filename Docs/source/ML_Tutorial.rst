.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

Tutorials/ML
==========================

In ML/PYTORCH there is an example of using a pre-trained Pytorch model in the
AMReX framework.  In this example we initialize data on a MultiFab, copy
the data into a pytorch tensor, call the model, and load the result back
into a MultiFab.  The example works on both the host and GPU.
A model is included as an example (`/Exec/model.pt`).
See ML/PYTORCH/README.md for libtorch download instructions.
