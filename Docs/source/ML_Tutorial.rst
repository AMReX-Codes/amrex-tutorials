ML/PYTORCH
==========

The overall goal of machine learning models in this context is to accelerate computationally expensive kernels/routines as part of an AMReX simulation.
This tutorial demonstrates how to interface a pre-trained PyTorch machine learning model to an AMReX simulation by querying inputs from and supplying outputs to an AMReX MultiFab.
Here we use a 1-input, 2-output model to illustrate the interface between the PyTorch model and a MultiFab.

PyTorch is a commonly used machine learning package with a C++ API library called LibTorch.
Located in the directory ``amrex-tutorials/ExampleCodes/ML/PYTORCH``, this example uses a machine learning model to solve a radioactive beta decay problem.
To begin, we initialize data on a MultiFab, then copy the data into a PyTorch tensor, then we call the pre-trained model to compute the outputs, and finally we load the result back into a MultiFab.
The program runs on either only the CPU or both the CPU and GPU.

**Running an AMReX application with a PyTorch model**
-----------------------------------------------------

Below is a step-by-step guide to successfully run an AMReX program that uses a PyTorch model. It is based on ``ML/PYTORCH/README.md`` and will require the model to have been saved as a TorchScript. In this example the TorchScript file is ``model.pt``. For more information on TorchScript, please see their `intro tutorial <https://pytorch.org/tutorials/beginner/Intro_to_TorchScript_tutorial.html>`_.

   1. Before compiling, either a CPU or GPU version of LibTorch (PyTorch C++ library) must be downloaded into ``ML/PYTORCH/``. An example of downloading the CUDA 11.1 version of ``libtorch`` and renaming it to ``libtorch_cuda`` is shown here:

      .. code-block:: console

         wget https://download.pytorch.org/libtorch/cu111/libtorch-cxx11-abi-shared-with-deps-1.9.0%2Bcu111.zip
         unzip libtorch-cxx11-abi-shared-with-deps-1.9.0+cu111.zip
         mv libtorch libtorch_cuda


      You can also check the website, `PyTorch <https://pytorch.org/get-started/locally/>`_ to download the latest version of LibTorch.

   2. Go to ``ML/PYTORCH/Exec`` to compile the executable. If using GPU, compile with ``USE_CUDA=TRUE``. Run ``make`` and it should result in an executable named ``main2d.gnu.MPI.CUDA.ex``

   3. Then you can run the example: ``./main2d.gnu.MPI.CUDA.ex inputs``.

**Beta Decay**
--------------

In this example, the machine learning model is a regression model pre-trained to solve a two-component ODE system describing beta decay. The input is a time step ``dt`` and output is the two-component solution of the ODE system at time ``t = dt``.

**Pre-trained Model**
---------------------
The TorchScript model that is included in this example is located at ``ML/PYTORCH/Exec/model.pt``. If you wish to change the model, edit the ``model_file`` parameter in ``ML/PYTORCH/Exec/inputs`` to your desired PyTorch model file location.

