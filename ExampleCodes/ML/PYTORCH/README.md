This tutorial presents an example of using a pre-trained Pytorch model in the AMReX framework.
A model is included as an example (`/Exec/model.pt`).

You will first need to download the appropriate Pytorch library (libtorch) in the current directory. Note that we have renamed the folder to indicate whether it uses CPU or CUDA.

```shell
# CPU only
wget https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.9.0%2Bcpu.zip
unzip libtorch-cxx11-abi-shared-with-deps-1.9.0+cpu.zip
mv libtorch libtorch_cpu

# CUDA 11.1
wget https://download.pytorch.org/libtorch/cu111/libtorch-cxx11-abi-shared-with-deps-1.9.0%2Bcu111.zip
unzip libtorch-cxx11-abi-shared-with-deps-1.9.0+cu111.zip
mv libtorch libtorch_cuda

# CUDA 10.2
wget https://download.pytorch.org/libtorch/cu102/libtorch-cxx11-abi-shared-with-deps-1.9.0%2Bcu102.zip
unzip libtorch-cxx11-abi-shared-with-deps-1.9.0+cu102.zip
mv libtorch libtorch_cuda
```

The source code is located in the directory `/Source`. To compile the example, go to the `/Exec` directory and enter
`make -j`.

Then you can run the example.

`mpiexec -n 1 ./main2d.gnu.MPI.ex inputs`


<!---
Torchscript code based on https://pytorch.org/tutorials/advanced/cpp_export.html 
-->
