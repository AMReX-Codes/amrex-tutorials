#include "MyKernels.H"

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(63));
        amrex::BoxArray ba(domain);
        ba.maxSize(32);
        amrex::DistributionMapping dm(ba);

        amrex::MultiFab mf(ba, dm, 1, 0);
        fill(mf);

        amrex::Print() << "Sum = " << mf.sum(0) << "\n";
    }
    amrex::Finalize();
}
