# required dependencies
module load gpu
module load PrgEnv-gnu
module load craype
module load craype-x86-milan
module load craype-accel-nvidia80
module load cudatoolkit
module load cmake/3.24.3

# necessary to use CUDA-Aware MPI and run a job
export CRAY_ACCEL_TARGET=nvidia80

# optimize CUDA compilation for A100
export AMREX_CUDA_ARCH=8.0

# optimize CPU microarchitecture for AMD EPYC 3rd Gen (Milan/Zen3)
# note: the cc/CC/ftn wrappers below add those
export CXXFLAGS="-march=znver3"
export CFLAGS="-march=znver3"

# compiler environment hints
export CC=cc
export CXX=CC
export FC=ftn
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=CC
