#include <iostream>
#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void initialize_domain(const Box& bx, Array4<Real> fabarray)
{
    auto lo = bx.loVect3d();
    auto hi = bx.hiVect3d();
    for (int k = lo[2]; k <= hi[2]; k++) {
        for (int j = lo[1]; j <= hi[1]; j++) {
            for (int i = lo[0]; i <= hi[0]; i++) {
                fabarray(i, j, k) = static_cast<Real>(i + 100.0*j + 10000.0*k);
            }
        }
    }
}

int main(int argc,
         char* argv[])
{
    Initialize(argc, argv);

    {
        int number_cells = 32;
        std::vector<int> ncells {AMREX_D_DECL(number_cells, number_cells, number_cells)};
        int max_grid_size = 16;

        // Read input parameters
        {
            ParmParse pp;
            if (!pp.queryarr("n_cells", ncells, 0, AMREX_SPACEDIM))
                Print() << "n_cells not specified, so using 32 cells in each dimension.\n";

            if (!pp.query("max_grid_size", max_grid_size))
                Print() << "max_grid_size not specified, so using 16.\n";
        }

        BoxArray ba;
        Geometry geo;

        // Define BoxArray and Geometry for our domain
        {
            // Define index space
            Box bx(IntVect(AMREX_D_DECL(0, 0, 0)),
                   IntVect(AMREX_D_DECL(ncells[0]-1, ncells[1]-1, ncells[2]-1)),
                   IntVect(AMREX_D_DECL(0, 0, 0)));
            ba.define(bx);
            ba.maxSize(max_grid_size);

            // Define physical space
            RealBox rbox(AMREX_D_DECL(0.0, 0.0, 0.0),
                         AMREX_D_DECL(1.0, 1.0, 1.0));

            // Cartesian coordinate system
            int coordinates = 0;

            // Fully periodic domain
            std::array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(1,1,1)};

            // Define Geometry
            geo.define(bx, &rbox, coordinates, is_periodic.data());
        }

        // Construct DistributionMapping
        DistributionMapping dm {ba};

        // 1 component, no ghost cells
        int num_components = 1;
        int num_ghost_cell = 0;

        // Build MultiFab
        MultiFab state(ba, dm, num_components, num_ghost_cell);

        // Initialize MultiFab
        for (MFIter mfi(state); mfi.isValid(); ++mfi) {
            const Box &bx = mfi.validbox();
            Array4<Real> fabarray = state.array(mfi);
            initialize_domain(bx, fabarray);
        }

        // Write MultiFab to plotfile
        std::string pfname = "plt_" + std::to_string(ncells[0]);
#if (AMREX_SPACEDIM >= 2)
        pfname += "_" + std::to_string(ncells[1]);
#endif
#if (AMREX_SPACEDIM == 3)
        pfname += "_" + std::to_string(ncells[2]);
#endif
        pfname += "_" + std::to_string(max_grid_size);

        const Vector<std::string> varnames {"phi"};

        WriteSingleLevelPlotfile(pfname, state, varnames, geo, 0.0, 0);

    }

    amrex::Finalize();
}
