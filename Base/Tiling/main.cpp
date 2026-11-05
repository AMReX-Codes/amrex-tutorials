#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <limits>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        using amrex::Real;
        using namespace amrex::literals;

        int n_cell = 64;
        int max_grid_size = 32;
        int nrepeat = 10;
        amrex::IntVect tile_size(AMREX_D_DECL(1024000, 8, 8));
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("nrepeat", nrepeat);
            amrex::Vector<int> ts;
            if (pp.queryarr("tile_size", ts, 0, AMREX_SPACEDIM)) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { tile_size[d] = ts[d]; }
            }
        }

        // Part 1: the boxes an MFIter returns, on one box of 16 cells with
        // one ghost cell, owned by rank 0 (the I/O rank, so Print shows them).
        {
            amrex::Box bx(amrex::IntVect(0), amrex::IntVect(15));
            amrex::BoxArray ba(bx);
            amrex::DistributionMapping dm(amrex::Vector<int>{0});
            amrex::MultiFab mf(ba, dm, 1, 1);

            int niter = 0;
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) { // no tiling
                ++niter;
                amrex::Print() << "Without tiling, tilebox " << mfi.tilebox()
                               << " = validbox " << mfi.validbox() << "\n";
            }
            amrex::Print() << "  (" << niter << " iteration)\n";

            amrex::Print() << "With tile size " << tile_size << ":\n";
            amrex::Box validbox, fabbox;
            for (amrex::MFIter mfi(mf, tile_size); mfi.isValid(); ++mfi) {
                amrex::Print() << "  tile " << mfi.LocalTileIndex() << " of box " << mfi.index()
                               << ": tilebox " << mfi.tilebox() << "\n"
                               << "      growntilebox()  " << mfi.growntilebox() << "\n"
                               << "      nodaltilebox(" << AMREX_SPACEDIM-1 << ") "
                               << mfi.nodaltilebox(AMREX_SPACEDIM-1) << "\n";
                validbox = mfi.validbox(); // the same for every tile
                fabbox = mfi.fabbox();     // the valid box with its ghost cells
            }
            amrex::Print() << "  validbox " << validbox << " and fabbox " << fabbox
                           << " are the same for every tile\n";
        }

        // The domain for the rest: n_cell cells in each direction, periodic.
        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        amrex::DistributionMapping dm(ba);
        amrex::RealBox real_box({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
        amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian,
                             {AMREX_D_DECL(1, 1, 1)});
        amrex::Print() << "\nDomain " << domain << " in " << ba.size()
                       << " boxes, tile size " << tile_size << "\n";

        // Part 2: count how many times each cell is written for different
        // choices of the box of an iteration.  The adds are atomic so that
        // the count is exact even where the boxes of different tiles
        // overlap; without the atomics, such overlaps would be a data race
        // with OpenMP.
        amrex::iMultiFab visits(ba, dm, 1, 1);
        int vmin, vmax, gmin, gmax; // min and max of the counts

        // Correct: each tile does its own cells.
        visits.setVal(0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(visits, tile_size); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.tilebox();
            auto const& a = visits.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::HostDevice::Atomic::Add(&a(i,j,k), 1);
            });
        }
        vmin = visits.min(0);   vmax = visits.max(0);   // valid cells
        gmin = visits.min(0,1); gmax = visits.max(0,1); // ghost cells included
        amrex::Print() << "  tilebox          : visits per valid cell " << vmin << " to " << vmax
                       << ", ghost cells included " << gmin << " to " << gmax << "\n";
        bool tilebox_ok = (vmin == 1 && vmax == 1 && gmin == 0);

        // Wrong: every tile does the whole box, so each cell is done once
        // per tile, and with OpenMP several threads write the same cells.
        visits.setVal(0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(visits, tile_size); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.validbox();
            auto const& a = visits.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::HostDevice::Atomic::Add(&a(i,j,k), 1);
            });
        }
        vmin = visits.min(0);   vmax = visits.max(0);   // valid cells
        gmin = visits.min(0,1); gmax = visits.max(0,1); // ghost cells included
        amrex::Print() << "  validbox         : visits per valid cell " << vmin << " to " << vmax
                       << ", ghost cells included " << gmin << " to " << gmax << "\n";

        // Wrong: a tile box grown by hand overlaps the neighboring tiles.
        visits.setVal(0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(visits, tile_size); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = amrex::grow(mfi.tilebox(), 1);
            auto const& a = visits.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::HostDevice::Atomic::Add(&a(i,j,k), 1);
            });
        }
        vmin = visits.min(0);   vmax = visits.max(0);   // valid cells
        gmin = visits.min(0,1); gmax = visits.max(0,1); // ghost cells included
        amrex::Print() << "  grow(tilebox, 1) : visits per valid cell " << vmin << " to " << vmax
                       << ", ghost cells included " << gmin << " to " << gmax << "\n";

        // Correct: growntilebox grows only on the sides of the tile that are
        // on the boundary of the box, so the tiles still do not overlap and
        // together they cover the ghost cells too.
        visits.setVal(0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(visits, tile_size); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.growntilebox();
            auto const& a = visits.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::HostDevice::Atomic::Add(&a(i,j,k), 1);
            });
        }
        vmin = visits.min(0);   vmax = visits.max(0);   // valid cells
        gmin = visits.min(0,1); gmax = visits.max(0,1); // ghost cells included
        amrex::Print() << "  growntilebox()   : visits per valid cell " << vmin << " to " << vmax
                       << ", ghost cells included " << gmin << " to " << gmax << "\n";
        bool growntilebox_ok = (gmin == 1 && gmax == 1);

        // Faces in direction fdir: the face-centered boxes of the tiles of
        // the cell-centered data do not overlap either; each face belongs
        // to one tile.  The MFIter runs over the cell-centered visits, so
        // that the tiles are cell-centered; fvisits has the same boxes up
        // to the index type.
        int const fdir = AMREX_SPACEDIM-1;
        amrex::iMultiFab fvisits(amrex::convert(ba, amrex::IntVect::TheDimensionVector(fdir)),
                                 dm, 1, 0);
        fvisits.setVal(0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(visits, tile_size); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.nodaltilebox(fdir);
            auto const& a = fvisits.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::HostDevice::Atomic::Add(&a(i,j,k), 1);
            });
        }
        int fmin = fvisits.min(0), fmax = fvisits.max(0);
        amrex::Print() << "  nodaltilebox(" << fdir << ")  : visits per face "
                       << fmin << " to " << fmax << "\n";
        bool nodaltilebox_ok = (fmin == 1 && fmax == 1);

        // Wrong: converting the tile box by hand gives the face between two
        // neighboring tiles to both of them.
        fvisits.setVal(0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(visits, tile_size); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = amrex::surroundingNodes(mfi.tilebox(), fdir);
            auto const& a = fvisits.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::HostDevice::Atomic::Add(&a(i,j,k), 1);
            });
        }
        fmin = fvisits.min(0); fmax = fvisits.max(0);
        amrex::Print() << "  surroundingNodes(tilebox, " << fdir << "): visits per face "
                       << fmin << " to " << fmax << "\n";

        // Part 3: a stencil with and without tiling.  TilingIfNotGPU() is
        // true on CPU and false on GPU, where tiling would only make the
        // kernels smaller.  The tile size is fabarray.mfiter_tile_size.
        // Both loops share the work among the OpenMP threads: the tiled
        // one by tiles, the untiled one by boxes.
        amrex::MultiFab phi(ba, dm, 1, 1);
        amrex::MultiFab lap_tiled(ba, dm, 1, 0);
        amrex::MultiFab lap_untiled(ba, dm, 1, 0);
        auto const prob_lo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();
        auto const dxinv = geom.InvCellSizeArray();
        {
            // phi = prod_d sin(2 pi x_d) at the cell centers
            auto const& pa = phi.arrays();
            amrex::ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                amrex::IntVect iv(AMREX_D_DECL(i,j,k));
                Real r = 1.0_rt;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    Real x = prob_lo[d] + (iv[d] + 0.5_rt) * dx[d];
                    r *= std::sin(2.0_rt * amrex::Math::pi<Real>() * x);
                }
                pa[b](i,j,k) = r;
            });
            phi.FillBoundary(geom.periodicity()); // ghost cells, periodic
        }

        for (int irep = 0; irep < nrepeat; ++irep)
        {
            {
                BL_PROFILE("laplacian_tiled");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (amrex::MFIter mfi(phi, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    amrex::Box const& bx = mfi.tilebox();
                    auto const& p = phi.const_array(mfi);
                    auto const& l = lap_tiled.array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        l(i,j,k) = AMREX_D_TERM(
                              (p(i+1,j,k) - 2.0_rt*p(i,j,k) + p(i-1,j,k)) * dxinv[0]*dxinv[0],
                            + (p(i,j+1,k) - 2.0_rt*p(i,j,k) + p(i,j-1,k)) * dxinv[1]*dxinv[1],
                            + (p(i,j,k+1) - 2.0_rt*p(i,j,k) + p(i,j,k-1)) * dxinv[2]*dxinv[2]);
                    });
                }
                amrex::Gpu::streamSynchronize(); // so that the timer includes the GPU kernels
            }
            {
                BL_PROFILE("laplacian_untiled");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (amrex::MFIter mfi(phi); mfi.isValid(); ++mfi) { // one box per iteration
                    amrex::Box const& bx = mfi.validbox();
                    auto const& p = phi.const_array(mfi);
                    auto const& l = lap_untiled.array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        l(i,j,k) = AMREX_D_TERM(
                              (p(i+1,j,k) - 2.0_rt*p(i,j,k) + p(i-1,j,k)) * dxinv[0]*dxinv[0],
                            + (p(i,j+1,k) - 2.0_rt*p(i,j,k) + p(i,j-1,k)) * dxinv[1]*dxinv[1],
                            + (p(i,j,k+1) - 2.0_rt*p(i,j,k) + p(i,j,k-1)) * dxinv[2]*dxinv[2]);
                    });
                }
                amrex::Gpu::streamSynchronize();
            }
        }
        // The two results are the same up to rounding.
        amrex::MultiFab::Subtract(lap_untiled, lap_tiled, 0, 0, 1, 0);
        Real const diff = lap_untiled.norm0();
        Real const eps = std::numeric_limits<Real>::epsilon();
        amrex::Print() << "  Laplacian with and without tiling: max difference " << diff << "\n";
        bool laplacian_ok = (diff <= 100.0_rt * eps * dxinv[0]*dxinv[0]);

        if (tilebox_ok && growntilebox_ok && nodaltilebox_ok && laplacian_ok) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
