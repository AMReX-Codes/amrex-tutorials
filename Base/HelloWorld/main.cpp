#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_OpenMP.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

int main (int argc, char* argv[])
{
    // Initialize MPI and the GPU device, if enabled.
    amrex::Initialize(argc, argv);
    {
        // amrex::Print only prints on the I/O process (rank 0).
        amrex::Print() << "Hello world from AMReX version "
                       << amrex::Version() << "\n";

        amrex::Print() << "  MPI ranks:      " << amrex::ParallelDescriptor::NProcs() << "\n"
                       << "  OpenMP threads: " << amrex::OpenMP::get_max_threads() << "\n"
                       << "  GPU backend:    "
#if defined(AMREX_USE_CUDA)
                       << "CUDA"
#elif defined(AMREX_USE_HIP)
                       << "HIP"
#elif defined(AMREX_USE_SYCL)
                       << "SYCL"
#else
                       << "none"
#endif
                       << "\n";

        // amrex::AllPrint prints on every process.
        amrex::AllPrint() << "  Rank " << amrex::ParallelDescriptor::MyProc()
                          << " says hello\n";
    }
    amrex::Finalize();
}
