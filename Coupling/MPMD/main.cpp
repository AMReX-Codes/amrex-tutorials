#include <AMReX.H>
#include <AMReX_MPMD.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParReduce.H>
#include <AMReX_Print.H>

#include <cmath>
#include <limits>

using namespace amrex;

// The first program of the MPMD pair. It owns the whole domain, sends
// component 0 to the agent program, and receives component 1 back on the
// part of the domain that the agent covers.
int main (int argc, char* argv[])
{
    // Split MPI_COMM_WORLD into one communicator per program, and run AMReX
    // on this program's communicator only.
    MPI_Comm comm = MPMD::Initialize(argc, argv);
    amrex::Initialize(argc, argv, true, comm);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
        }

        amrex::Print() << "Main: " << ParallelDescriptor::NProcs() << " of "
                       << MPMD::NProcs() << " ranks in MPI_COMM_WORLD\n";

        Box domain(IntVect(0), IntVect(n_cell-1));
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        // Component 0: a Gaussian. Component 1: -1, to be overwritten by the agent.
        MultiFab mf(ba, dm, 2, 0);
        auto const& ma = mf.arrays();
        Real const dx = 1._rt/static_cast<Real>(n_cell);
        ParallelFor(mf, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k)
        {
            IntVect iv(AMREX_D_DECL(i,j,k));
            Real r2 = 0._rt;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                Real x = (static_cast<Real>(iv[d]) + 0.5_rt)*dx - 0.5_rt;
                r2 += x*x;
            }
            ma[bi](i,j,k,0) = 1._rt + std::exp(-r2/0.01_rt);
            ma[bi](i,j,k,1) = -1._rt;
        });

        // The Copier exchanges BoxArray and DistributionMapping with the other
        // program and builds the communication pattern for the intersection
        // of the two BoxArrays. Both programs must construct it.
        MPMD::Copier copier(ba, dm);

        copier.send(mf, 0, 1); // send component 0
        copier.recv(mf, 1, 1); // receive into component 1

        // The agent covers the upper half of the domain in every direction
        // and returns 10 times what it received.
        Box const agent_box(IntVect(n_cell/2), IntVect(n_cell-1));
        auto const& cma = mf.const_arrays();
        Real err = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, mf, IntVect(0),
            [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) -> GpuTuple<Real>
            {
                Real expected = agent_box.contains(IntVect(AMREX_D_DECL(i,j,k)))
                    ? 10._rt*cma[bi](i,j,k,0) : -1._rt;
                return { std::abs(cma[bi](i,j,k,1) - expected) };
            });
        ParallelDescriptor::ReduceRealMax(err);

        amrex::Print() << "Main: max error of the round trip = " << err << "\n";
        if (err > 100._rt*std::numeric_limits<Real>::epsilon()) {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    MPMD::Finalize();
    return status;
}
