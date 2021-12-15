.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

ML/PYTORCH
==========================

In ``amrex-tutorials/ExampleCodes/ML/PYTORCH`` there is an example demonstrating the usage of a pre-trained PyTorch model in the AMReX framework.  In this example, we are using an ML model to solve a beta-decay problem. We start by initializing data on a MultiFab, then copy the data into a PyTorch tensor, call the model, and finally load the result back into a MultiFab.  The program works on both the host (CPU) and GPU.

**Running AMReX application with PyTorch model**
------------------------------------------------

Below is a step-by-step guide to succesfully run an AMReX program that uses PyTorch model. It is based on ``ML/PYTORCH/README.md`` and will require the model to have been saved as a TorchScript. For more information on TorchScript, please visit `here <https://pytorch.org/tutorials/beginner/Intro_to_TorchScript_tutorial.html>`.

1. Before compiling, either a CPU or GPU version of LibTorch (PyTorch C++ library) must be downloaded into ``ML/PYTORCH/``. An example of downloading the CUDA 11.1 version of ``libtorch`` and renaming it to ``libtorch_cuda`` is shown here:

   .. highlight:: console

   ::
      wget https://download.pytorch.org/libtorch/cu111/libtorch-cxx11-abi-shared-with-deps-1.9.0%2Bcu111.zip
      unzip libtorch-cxx11-abi-shared-with-deps-1.9.0+cu111.zip
      mv libtorch libtorch_cuda

   You can also check PyTorch `website <https://pytorch.org/get-started/locally/>` to download the latest version of LibTorch.

2. Go to ``ML/PYTORCH/Exec`` to compile the executable. If using GPU, compile with ``USE_CUDA=TRUE``. Run ``make`` and it should result in an executable named ``main2d.gnu.MPI.CUDA.ex``

3. Then you can run the example: ``./main2d.gnu.MPI.CUDA.ex inputs``.

**Beta Decay**
----------------------

In this example, the ML model is a regression model pre-trained to solve the two-component ODE system describing beta decay. The input is a timestep ``dt`` and output is the two-component solution of the ODE system at time ``t = dt``.

**Pre-trained Model**
---------------------
The TorchScript model that is included in this example is located at (``ML/PYTORCH/Exec/model.pt``). If you wish to change the model, edit the ``model_file`` parameter in ``ML/PYTORCH/Exec/inputs`` to your desired PyTorch model file location.

