# PyTorchInference ★★★

Evaluate a PyTorch model on `MultiFab` data with LibTorch, the C++ API of
PyTorch. The model replaces a pointwise kernel: for every cell it maps two
input fields to two output fields. The cells of each box are gathered into a
batch, wrapped in a tensor without a copy, passed through the model, and the
result is scattered back into a `MultiFab`. In CUDA builds, the data stay on
the GPU and torch runs on AMReX's CUDA stream.

**Key APIs:** `torch::jit::load`, `torch::from_blob`, `c10::InferenceMode`,
`c10::cuda::CUDAStreamGuard`, `Gpu::DeviceVector`, `Box::index`, `MFIter`

## The model

`make_model.py` writes `model.pt`, a TorchScript model with one dense layer,
`y = tanh(W x + b)`, with fixed weights. The tutorial compares the model
output with this formula computed in an AMReX kernel. A trained network is
saved and used the same way, as long as it maps a `(cells, inputs)` tensor
to a `(cells, outputs)` tensor. The model is converted to `amrex::Real` when
it is loaded.

```
python3 make_model.py    # needs the torch Python package
```

## Build and run

LibTorch is needed: an unpacked LibTorch download from pytorch.org, or the
one in the `torch` Python package (`python3 -c "import torch;
print(torch.utils.cmake_prefix_path)"` gives its CMake directory). Use a
CUDA LibTorch with `AMReX_GPU_BACKEND=CUDA`.

CMake, from the repository root:

```
cmake -S . -B build -DAMReXTutorials_TOPICS=ML -DAMReXTutorials_TORCH=ON \
      -DCMAKE_PREFIX_PATH=/path/to/libtorch
cmake --build build --target ML_PyTorchInference -j 8
cd build/ML/PyTorchInference && mpiexec -n 2 ./PyTorchInference inputs
```

GNU Make, with AMReX cloned next to this repository (LibTorch built with the
pre-C++11 ABI also needs `-D_GLIBCXX_USE_CXX11_ABI=0`):

```
make -j8 PYTORCH_ROOT=/path/to/libtorch
mpiexec -n 2 ./main3d.gnu.MPI.ex inputs
```

## Expected output

```
Loaded model.pt on cpu
Output range: [-0.7585439934, 0.9635137063] and [-0.9863230868, 0.9632540745]
Max error vs. tanh(W x + b): 1.110223025e-16
PASSED
```
