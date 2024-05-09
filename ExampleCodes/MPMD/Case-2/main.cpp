
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <mpi.h>
#include <AMReX_MPMD.H>

int main(int argc, char* argv[])
{
    // Initialize amrex::MPMD to establish communication across the two apps
    MPI_Comm comm = amrex::MPMD::Initialize(argc, argv);
    amrex::Initialize(argc,argv,true,comm);
    {
        amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";
        // Number of data components at each grid point in the MultiFab
        int ncomp = 2;
        // how many grid cells in each direction over the problem domain
        int n_cell = 32;
        // how many grid cells are allowed in each direction over each box
        int max_grid_size = 16;
        //BoxArray -- Abstract Domain Setup
        // integer vector indicating the lower coordindate bounds
        amrex::IntVect dom_lo(0,0,0);
        // integer vector indicating the upper coordindate bounds
        amrex::IntVect dom_hi(n_cell-1, n_cell-1, n_cell-1);
        // box containing the coordinates of this domain
        amrex::Box domain(dom_lo, dom_hi);
        // will contain a list of boxes describing the problem domain
        amrex::BoxArray ba(domain);
        // chop the single grid into many small boxes
        ba.maxSize(max_grid_size);
        // Distribution Mapping
        amrex::DistributionMapping dm(ba);
        // Create an MPMD Copier that
        // sends the BoxArray information to the other (python) application
        auto copr = amrex::MPMD::Copier(ba,dm,true);
        //Define MuliFab
        amrex::MultiFab mf(ba, dm, ncomp, 0);
        //Geometry -- Physical Properties for data on our domain
        amrex::RealBox real_box ({0., 0., 0.}, {1. , 1., 1.});
        amrex::Geometry geom(domain, &real_box);
        //Calculate Cell Sizes
        amrex::GpuArray<amrex::Real,3> dx = geom.CellSizeArray();  //dx[0] = dx dx[1] = dy dx[2] = dz
        //Fill only the first component of the MultiFab
        for(amrex::MFIter mfi(mf); mfi.isValid(); ++mfi){
            const amrex::Box& bx = mfi.validbox();
            const amrex::Array4<amrex::Real>& mf_array = mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){

                amrex::Real x = (i+0.5) * dx[0];
                amrex::Real y = (j+0.5) * dx[1];
                amrex::Real z = (k+0.5) * dx[2];
                amrex::Real r_squared = ((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))/0.01;

                mf_array(i,j,k,0) = 1.0 + std::exp(-r_squared);

            });
        }
        // Send ONLY the first populated MultiFab component to the other app
        copr.send(mf,0,1);
        // Receive ONLY the second MultiFab component from the other app
        copr.recv(mf,1,1);
        //Plot MultiFab Data
        WriteSingleLevelPlotfile("plt_cpp_001", mf, {"comp0","comp1"}, geom, 0., 0);
    }
    amrex::Finalize();
    amrex::MPMD::Finalize();
}
