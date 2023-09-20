.. _tutorials_ml:

ML/PYTORCH
==========

The overall goal of this machine learning tutorial is to accelerate computationally expensive point-wise kernels/routines within an AMReX simulation.
This tutorial demonstrates how to interface a pre-trained PyTorch machine learning model to an AMReX simulation by querying inputs from and supplying outputs to an AMReX MultiFab.
PyTorch is a commonly used machine learning package with a C++ API library called LibTorch.
Located in the directory ``amrex-tutorials/ExampleCodes/ML/PYTORCH``, this example uses a machine learning model to solve a radioactive beta decay problem.
Here we use a 1-input, 2-output model to illustrate the interface between the PyTorch model and a MultiFab.

**Beta Decay Reaction**
--------------

In this example, the machine learning model is a regression model pre-trained to solve a two-component ODE system describing beta decay.

.. math:: \frac{\partial X_0}{\partial t} = -X_0
.. math:: \frac{\partial X_1}{\partial t} = X_0
.. math:: X_0(0) = 1; ~~~ X_1(0) = 0

In the context of the pytorch model, the input is a time step ``dt`` and output is the two-component solution of the ODE system at time ``t = dt``.
          
**Pre-trained Model**
---------------------
The TorchScript model that is included in this example is located at ``ML/PYTORCH/Exec/model.pt``.
If you wish to change the model, edit the ``model_file`` parameter in ``inputs``.

**Running an AMReX application with a PyTorch model**
-----------------------------------------------------
To begin, we initialize a MultiFab full of data representing different ``dt`` values, then copy this data into a PyTorch tensor, then call the pre-trained model to compute the outputs, and finally load the result back into a MultiFab.
The model can be evaluated on the CPU or GPU.

Below is a step-by-step guide to successfully run an AMReX program that uses a PyTorch model. It will require the model to have been saved as a TorchScript. In this example the TorchScript file is ``model.pt``. For more information on TorchScript, please see their `intro tutorial <https://pytorch.org/tutorials/beginner/Intro_to_TorchScript_tutorial.html>`_.

   1. Before compiling, either a CPU or CUDA version of LibTorch (PyTorch C++ library) must be downloaded into ``ML/PYTORCH/``. To download the CPU-only version of ``libtorch`` and rename it to ``libtorch_cpu``:

      .. code-block:: console

         wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcpu.zip
         unzip libtorch-cxx11-abi-shared-with-deps-2.0.1+cpu.zip
         mv libtorch libtorch_cpu

      Similarly, the CUDA 11.8 version of ``libtorch`` can be downloaded and renamed to ``libtorch_cuda``:

      .. code-block:: console

         wget https://download.pytorch.org/libtorch/cu118/libtorch-cxx11-abi-shared-with-deps-2.0.1%2Bcu118.zip
         unzip libtorch-cxx11-abi-shared-with-deps-2.0.1+cu118.zip
         mv libtorch libtorch_cuda

      You can also check the website, `PyTorch <https://pytorch.org/get-started/locally/>`_ to download the latest version of LibTorch.

   2. Go to ``ML/PYTORCH/Exec`` to compile the executable.
      Run ``make`` and optionally ``USE_CUDA=TRUE`` and it should result in an executable named, e.g., ``main2d.gnu.MPI.CUDA.ex``

   3. Then you can run the example, e.g., ``./main2d.gnu.MPI.CUDA.ex inputs`` or ``mpiexec -n 4 ./main2d.gnu.MPI.ex inputs``.
      There will be two plotfiles, ``plt_inputs`` (containing ``dt``) and ``plt_outputs`` (containing ``X_0`` and ``X_1`` at the final time).
