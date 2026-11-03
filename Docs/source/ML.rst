.. _ml:

ML
==

This topic shows how to evaluate machine-learning models on AMReX data. The
tutorial needs LibTorch, the C++ API of PyTorch, and is built only with
``-DAMReXTutorials_TORCH=ON``, with ``CMAKE_PREFIX_PATH`` or ``Torch_DIR``
pointing to LibTorch (an unpacked LibTorch download, or the one in the
``torch`` Python package: ``python3 -c "import torch;
print(torch.utils.cmake_prefix_path)"``). It supports CPU and CUDA builds;
use a CUDA LibTorch with ``-DAMReX_GPU_BACKEND=CUDA``. Build only this topic
with ``-DAMReXTutorials_TOPICS=ML``; see :doc:`GettingStarted`.

.. _ml-pytorchinference:

PyTorchInference
----------------

:Level: ★★★
:Source: `ML/PyTorchInference <https://github.com/AMReX-Codes/amrex-tutorials/tree/main/ML/PyTorchInference>`__
:Key APIs: ``torch::jit::load``, ``torch::from_blob``, ``c10::InferenceMode``,
           ``c10::cuda::CUDAStreamGuard``, ``Gpu::DeviceVector``, ``Box::index``,
           ``MFIter``

A PyTorch model replaces a pointwise kernel: for every cell it maps two
input fields to two output fields. The cells of each box are gathered into
a batch, wrapped in a tensor without a copy, passed through the model, and
the result is scattered back into a ``MultiFab``. In CUDA builds, the data
stay on the GPU and torch runs on AMReX's CUDA stream.

``make_model.py`` writes ``model.pt``, a TorchScript model with one dense
layer, :math:`y = \tanh(W x + b)`, with fixed weights, so that the program
can compare the model output with this formula. A trained network is saved
and used the same way, as long as it maps a ``(cells, inputs)`` tensor to a
``(cells, outputs)`` tensor.

.. literalinclude:: ../../ML/PyTorchInference/make_model.py
   :language: python
   :pyobject: Model

The model is loaded on the device where the MultiFab data live (the CPU,
or the current CUDA device), and converted to the tensor type of
``amrex::Real``:

.. literalinclude:: ../../ML/PyTorchInference/main.cpp
   :language: cpp
   :start-at: // Load the TorchScript model written by make_model.py.
   :end-at: model.eval();
   :dedent:

For each box, the cells are gathered into an ``(ncells, nin)`` batch, the
model is evaluated, and the output is scattered back:

.. literalinclude:: ../../ML/PyTorchInference/main.cpp
   :language: cpp
   :start-at: // Gather the cells of the box into an (ncells, nin) batch
   :end-before: // batch and output are freed at the end of this iteration.
   :dedent:

Configure from the repository root, build and run (see :doc:`GettingStarted`):

.. code-block:: console

   cmake -S . -B build -DAMReXTutorials_TOPICS=ML -DAMReXTutorials_TORCH=ON \
         -DCMAKE_PREFIX_PATH=/path/to/libtorch
   cmake --build build --target ML_PyTorchInference -j 8
   cd build/ML/PyTorchInference
   mpiexec -n 2 ./PyTorchInference inputs

``model.pt`` is in the repository; regenerate it with
``python3 make_model.py`` (needs the ``torch`` Python package). With GNU
Make, pass ``PYTORCH_ROOT=/path/to/libtorch``.

.. code-block:: text

   Loaded model.pt on cpu
   Output range: [-0.7585439934, 0.9635137063] and [-0.9863230868, 0.9632540745]
   Max error vs. tanh(W x + b): 1.110223025e-16
   PASSED
