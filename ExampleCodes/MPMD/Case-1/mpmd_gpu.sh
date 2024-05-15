#!/bin/bash
#SBATCH -N 3 # Total number of nodes
#SBATCH -n 12 # Total number of tasks
#SBATCH -c 4 # number of processors per MPI task
#SBATCH -C gpu
#SBATCH -G 12 # Total number of GPUs
#SBATCH -q debug
#SBATCH -J mpmd_test
#SBATCH -t 00:05:00
#SBATCH -A mpxxx

source ./perlmutter_gpu.profile

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
# Taken from WarpX
export MPICH_OFI_NIC_POLICY=GPU
GPU_AWARE_MPI="amrex.use_gpu_aware_mpi=1"

#run the application:
srun --multi-prog --cpu_bind=cores --gpu-bind=single:1 ./mpmd_gpu.conf
