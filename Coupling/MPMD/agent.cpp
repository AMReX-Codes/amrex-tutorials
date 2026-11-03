#include <AMReX.H>
#include <AMReX_MPMD.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using namespace amrex;

// The second program of the MPMD pair. It works on the upper half of the
// main program's domain, with its own BoxArray and DistributionMapping.
int main (int argc, char* argv[])
{
    MPI_Comm comm = MPMD::Initialize(argc, argv);
    amrex::Initialize(argc, argv, true, comm);
    {
        int n_cell = 32;
        int agent_max_grid_size = 8;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("agent_max_grid_size", agent_max_grid_size);
        }

        amrex::Print() << "Agent: " << ParallelDescriptor::NProcs() << " of "
                       << MPMD::NProcs() << " ranks in MPI_COMM_WORLD\n";

        Box domain(IntVect(n_cell/2), IntVect(n_cell-1));
        BoxArray ba(domain);
        ba.maxSize(agent_max_grid_size);
        DistributionMapping dm(ba);

        MultiFab mf(ba, dm, 2, 0);

        MPMD::Copier copier(ba, dm);

        copier.recv(mf, 0, 1); // receive component 0

        auto const& ma = mf.arrays();
        ParallelFor(mf, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k)
        {
            ma[bi](i,j,k,1) = 10._rt*ma[bi](i,j,k,0);
        });
        Gpu::streamSynchronize();

        amrex::Print() << "Agent: received component 0 with min = " << mf.min(0)
                       << ", max = " << mf.max(0) << "\n";

        copier.send(mf, 1, 1); // send component 1 back
    }
    amrex::Finalize();
    MPMD::Finalize();
}
