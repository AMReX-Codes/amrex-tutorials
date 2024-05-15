#!/bin/bash
#SBATCH -N 1 # Total number of nodes
#SBATCH -n 12 # Total number of tasks
#SBATCH -c 4 # number of processors per MPI task
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -J mpmd_test
#SBATCH -t 00:05:00
#SBATCH -A mpxxx

# Activate the virtual environment
source /path/to/pyamrex-gpu/bin/activate

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application:
srun --multi-prog --cpu_bind=cores ./mpmd_cpu.conf

deactivate
